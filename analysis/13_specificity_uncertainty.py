#!/usr/bin/env python3
"""
13_specificity_uncertainty.py

Second-round review asked what the cross-track comparison can actually support.

THE OBJECTION, and it is correct. Our benchmark comprises eQTLs DETECTED in
kidney. It does not establish that those effects are kidney-SPECIFIC. Many
regulatory effects are shared across tissues, so a model that predicts shared
effects correctly would perform similarly across tissue tracks without failing at
tissue specificity. Ranking 6th of 55 is also near the top, not mid-pack.

So the earlier reading, that the model "does not resolve tissue identity", is not
established by that ranking. Two things are needed before any statement about
similarity of the OUTPUTS.

  A. PAIRED UNCERTAINTY. A difference in correlation between two tracks on the
     same 58 observations needs an interval. We bootstrap the 58 pairs and report
     the paired difference between the kidney track and each comparator, so the
     question becomes whether kidney is distinguishable from the best track
     rather than whether it ranked first.

  B. DO THE TRACKS ACTUALLY MAKE THE SAME PREDICTIONS? Similar correlation with
     a third quantity does not mean similar predictions. We therefore correlate
     the tissue tracks' PREDICTION VECTORS against each other directly. If those
     pairwise correlations are near 1, the outputs really are near-duplicates on
     these variants. If they are not, the tracks differ and the similar
     performance has some other explanation.

Analysis B is the one that can distinguish the two explanations, and it is
reported whichever way it comes out.

Output: results/specificity_uncertainty.json,
        results/specificity_pairwise.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
KID = "Kidney_Cortex"
RNG = np.random.default_rng(20260914)
B = 10000


def main() -> None:
    # rebuild the per-variant matrix produced by script 12
    src = RESULTS / "specificity_predictions.tsv.gz"
    if not src.exists():
        raise SystemExit(
            "results/specificity_predictions.tsv.gz not found. Re-run script 12; "
            "it must persist the per-variant prediction matrix, not only the "
            "per-tissue summary.")
    d = pd.read_csv(src, sep="\t")
    tissues = [c for c in d.columns
               if c not in ("variant", "gene", "observed_beta")]
    obs = d["observed_beta"].values
    n = len(d)
    print(f"{n} kidney variant-gene pairs, {len(tissues)} tissue tracks\n")

    # ------------------------------------------------------------------ A
    print("=" * 76)
    print("A. PAIRED DIFFERENCE FROM THE KIDNEY TRACK (bootstrap over the 58 pairs)")
    print("=" * 76)
    idx = RNG.integers(0, n, size=(B, n))
    def rho_boot(v):
        return np.array([stats.spearmanr(v[i], obs[i]).statistic for i in idx])

    kb = rho_boot(d[KID].values)
    rows = []
    for t in tissues:
        if t == KID:
            continue
        tb = rho_boot(d[t].values)
        diff = tb - kb
        lo, hi = np.nanpercentile(diff, [2.5, 97.5])
        rows.append({"tissue_track": t,
                     "rho": float(stats.spearmanr(d[t], obs).statistic),
                     "diff_vs_kidney": float(np.nanmean(diff)),
                     "lo": float(lo), "hi": float(hi),
                     "excludes_zero": bool(lo > 0 or hi < 0)})
    pr = pd.DataFrame(rows).sort_values("diff_vs_kidney", ascending=False)
    n_better = int(pr["excludes_zero"].sum())
    print(f"  kidney rho {stats.spearmanr(d[KID], obs).statistic:.3f}")
    print(f"  {'track':<34}{'rho':>7}{'diff':>8}{'95% CI of difference':>26}")
    for _, r in pr.head(6).iterrows():
        print(f"  {r['tissue_track']:<34}{r['rho']:>7.3f}{r['diff_vs_kidney']:>+8.3f}"
              f"   [{r['lo']:>+6.3f}, {r['hi']:>+6.3f}]"
              + ("  *" if r["excludes_zero"] else ""))
    print(f"\n  tracks whose difference from kidney excludes zero: "
          f"{n_better} of {len(pr)}")
    if n_better == 0:
        print("  No track is distinguishable from kidney on these 58 pairs.")
        print("  The ranking is therefore not evidence that any track is better,")
        print("  and equally not evidence that kidney carries no specific signal.")

    # ------------------------------------------------------------------ B
    print()
    print("=" * 76)
    print("B. DO DIFFERENT TISSUE TRACKS MAKE THE SAME PREDICTIONS?")
    print("=" * 76)
    P = d[tissues].values
    C = np.corrcoef(np.apply_along_axis(stats.rankdata, 0, P), rowvar=False)
    iu = np.triu_indices(len(tissues), k=1)
    pair = C[iu]
    ki = tissues.index(KID)
    kid_pairs = np.delete(C[ki], ki)

    pw = pd.DataFrame({"pair_rho": pair})
    pw.to_csv(RESULTS / "specificity_pairwise.csv", index=False)

    print(f"  pairwise rank correlation between tissue tracks' prediction vectors")
    print(f"    all {len(pair):,} pairs : median {np.median(pair):.3f}, "
          f"IQR {np.percentile(pair,25):.3f}-{np.percentile(pair,75):.3f}, "
          f"min {pair.min():.3f}")
    print(f"    kidney vs others  : median {np.median(kid_pairs):.3f}, "
          f"min {kid_pairs.min():.3f}, max {kid_pairs.max():.3f}")
    print()
    med = float(np.median(pair))
    if med > 0.95:
        verdict = ("Tracks make near-identical predictions on these variants. "
                   "Similar performance reflects near-duplicate outputs.")
    elif med > 0.8:
        verdict = ("Tracks are highly but not perfectly correlated. Outputs are "
                   "largely shared with some tissue-level variation.")
    else:
        verdict = ("Tracks make substantially DIFFERENT predictions despite "
                   "similar correlation with observed kidney effects. Similar "
                   "performance does NOT imply duplicate outputs, and the "
                   "'model ignores tissue' reading is NOT supported.")
    print(f"  VERDICT: {verdict}")

    out = {
        "n_pairs": int(n), "n_tracks": int(len(tissues)),
        "bootstrap_iterations": B,
        "resampling": "with replacement over the 58 variant-gene pairs, paired across tracks",
        "kidney_rho": float(stats.spearmanr(d[KID], obs).statistic),
        "tracks_distinguishable_from_kidney": n_better,
        "best_track": pr.iloc[0]["tissue_track"],
        "best_diff": float(pr.iloc[0]["diff_vs_kidney"]),
        "best_diff_ci": [float(pr.iloc[0]["lo"]), float(pr.iloc[0]["hi"])],
        "pairwise_prediction_rho": {
            "median": med,
            "iqr": [float(np.percentile(pair, 25)), float(np.percentile(pair, 75))],
            "min": float(pair.min()), "max": float(pair.max())},
        "kidney_vs_others_rho": {
            "median": float(np.median(kid_pairs)),
            "min": float(kid_pairs.min()), "max": float(kid_pairs.max())},
        "verdict": verdict,
    }
    (RESULTS / "specificity_uncertainty.json").write_text(json.dumps(out, indent=2))
    print(f"\n  Wrote {RESULTS}/specificity_uncertainty.json")


if __name__ == "__main__":
    main()
