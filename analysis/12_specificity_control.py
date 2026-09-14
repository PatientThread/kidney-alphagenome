#!/usr/bin/env python3
"""
12_specificity_control.py

The control the review asked for: does the KIDNEY output track carry
kidney-specific predictive information, or would any tissue's track do?

THE QUESTION. Our benchmark scored each kidney eQTL using the kidney cortex
track. That shows the model predicts kidney eQTL effects. It does not show the
kidney track is doing the work. If the same variants were predicted just as well
using liver, brain or whole blood tracks, then what the model captures is shared
regulatory effect, not kidney biology, and the entire premise that kidney output
coverage matters would be undermined.

THE DESIGN. Re-score every kidney benchmark variant, but extract the predicted
effect on the same target gene from EVERY GTEx tissue track the model exposes.
Then compute, for each of those tissue tracks, the correlation with the OBSERVED
KIDNEY eQTL effect. Kidney's own track should rank at or near the top if it
carries specific information.

WHAT EACH OUTCOME MEANS, decided before running.

  kidney track clearly top        the kidney output carries specific information
                                  and its scarcity matters as the paper argues
  kidney track mid-pack           predictions are driven by shared cross-tissue
                                  regulatory signal; the paper's premise weakens
                                  substantially and must be rewritten
  all tracks similar              the model is not resolving tissue at this locus
                                  set at all, which is itself a reportable result

This is a control with a real chance of undermining the paper. That is why it is
worth running.

Output: results/specificity_control.json, results/specificity_by_tissue.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(HERE))
from ag_auth import get_key, masked  # noqa: E402

SEQ_1MB = 1_048_576
KIDNEY_TRACK = "Kidney_Cortex"


def main() -> None:
    from alphagenome.models import dna_client, variant_scorers
    from alphagenome.data import genome

    bench = pd.read_csv(RESULTS / "benchmark_variants" / "kidney_cortex.tsv.gz",
                        sep="\t")
    print(f"kidney benchmark variants: {len(bench)}")
    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]

    rows, fails = [], 0
    t0 = time.time()
    for i, r in enumerate(bench.itertuples()):
        chrom, pos, ref, alt = str(r.variant).split("_")
        pos = int(pos)
        gene = str(r.molecular_trait_id).split(".")[0]
        start = max(0, pos - SEQ_1MB // 2)
        try:
            res = client.score_variant(
                interval=genome.Interval(chromosome=chrom, start=start,
                                         end=start + SEQ_1MB),
                variant=genome.Variant(chromosome=chrom, position=pos,
                                       reference_bases=ref, alternate_bases=alt),
                variant_scorers=[scorer])
        except Exception as exc:                                # noqa: BLE001
            fails += 1
            print(f"  FAIL {r.variant}: {type(exc).__name__}")
            continue

        ad = res[0]
        gcol = next((c for c in ("gene_id", "gene_name", "gene")
                     if c in ad.obs.columns), None)
        gm = (ad.obs[gcol].astype(str).str.split(".").str[0] == gene).values
        if not gm.any():
            fails += 1
            continue
        X = np.asarray(ad.X)
        gi = np.where(gm)[0]

        if "gtex_tissue" not in ad.var.columns:
            raise SystemExit("no gtex_tissue column in result metadata")
        rec = {"variant": r.variant, "gene": gene, "observed_beta": r.beta}
        for t in ad.var["gtex_tissue"].dropna().unique():
            cols = np.where(ad.var["gtex_tissue"].values == t)[0]
            rec[t] = float(np.nanmean(X[gi[:, None], cols]))
        rows.append(rec)
        time.sleep(0.15)
        if (i + 1) % 20 == 0:
            print(f"  {i+1}/{len(bench)}  ({time.time()-t0:.0f}s)")

    d = pd.DataFrame(rows)
    print(f"\nscored {len(d)}/{len(bench)}, {fails} failures\n")

    tissues = [c for c in d.columns
               if c not in ("variant", "gene", "observed_beta")]
    res = []
    for t in tissues:
        s = d[[t, "observed_beta"]].dropna()
        if len(s) < 10:
            continue
        rho, p = stats.spearmanr(s[t], s["observed_beta"])
        agree = float((np.sign(s[t]) == np.sign(s["observed_beta"])).mean())
        res.append({"tissue_track": t, "n": len(s), "spearman": float(rho),
                    "p": float(p), "sign_agreement": agree})
    m = pd.DataFrame(res).sort_values("spearman", ascending=False).reset_index(drop=True)
    m["rank"] = m.index + 1
    m.to_csv(RESULTS / "specificity_by_tissue.csv", index=False)

    kr = m[m["tissue_track"] == KIDNEY_TRACK]
    krank = int(kr["rank"].iloc[0]) if len(kr) else None
    krho = float(kr["spearman"].iloc[0]) if len(kr) else float("nan")

    print("=" * 74)
    print("PREDICTING KIDNEY eQTL EFFECTS USING EACH TISSUE'S OWN TRACK")
    print("=" * 74)
    print(f"  {'rank':>4}  {'tissue track':<38}{'rho':>8}{'sign':>8}")
    for _, r in m.head(10).iterrows():
        star = "  <-- KIDNEY" if r["tissue_track"] == KIDNEY_TRACK else ""
        print(f"  {int(r['rank']):>4}  {r['tissue_track']:<38}"
              f"{r['spearman']:>8.3f}{r['sign_agreement']:>8.2f}{star}")
    if krank and krank > 10:
        print("   ...")
        print(f"  {krank:>4}  {KIDNEY_TRACK:<38}{krho:>8.3f}"
              f"{float(kr['sign_agreement'].iloc[0]):>8.2f}  <-- KIDNEY")
    print()
    print(f"  kidney cortex track ranks {krank} of {len(m)} tissue tracks")
    print(f"  best track: {m.iloc[0]['tissue_track']} (rho {m.iloc[0]['spearman']:.3f})")
    print(f"  spread across tracks: {m['spearman'].min():.3f} to "
          f"{m['spearman'].max():.3f}, sd {m['spearman'].std():.3f}")

    pct = 100 * (1 - (krank - 1) / len(m))
    out = {"n_variants": int(len(d)), "n_tissue_tracks": int(len(m)),
           "kidney_rank": krank, "kidney_rho": krho,
           "kidney_percentile": round(pct, 1),
           "best_track": m.iloc[0]["tissue_track"],
           "best_rho": float(m.iloc[0]["spearman"]),
           "spread_sd": float(m["spearman"].std()),
           "range": [float(m["spearman"].min()), float(m["spearman"].max())]}
    (RESULTS / "specificity_control.json").write_text(json.dumps(out, indent=2))

    print()
    print("  READING:")
    if krank and krank <= 5:
        print("  The kidney track is at or near the top. It carries kidney-specific")
        print("  predictive information, and its scarcity is consequential.")
    elif krank and krank <= len(m) / 2:
        print("  The kidney track is above the middle but not clearly best. Some")
        print("  specificity, much shared signal. State this plainly.")
    else:
        print("  The kidney track does NOT outperform other tissues' tracks on")
        print("  kidney eQTLs. Predictions are driven by shared cross-tissue")
        print("  regulatory signal. The paper's premise that kidney OUTPUT")
        print("  coverage drives kidney performance is NOT supported and the")
        print("  manuscript must be rewritten around that.")
    print(f"\n  Wrote {RESULTS}/specificity_control.json")


if __name__ == "__main__":
    main()
