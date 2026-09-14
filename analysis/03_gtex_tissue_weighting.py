#!/usr/bin/env python3
"""
03_gtex_tissue_weighting.py

LEG 2, part 1: quantify the weighting bias in the published benchmark.

THE ARGUMENT. AlphaGenome's headline eQTL performance figure is a Spearman
correlation "averaged across all tissues, weighted by the number of eQTLs in
each tissue". In GTEx v8, the release the model was benchmarked on, kidney
cortex has the fewest eQTLs of any tissue. So kidney contributes the least of
any tissue to every published number, by construction, and no reader can tell
from the paper how the model does in kidney.

NAME THE RELEASE, ALWAYS. The unqualified claim is false as of v10, where
Bladder is smaller and kidney cortex is second of 50. This script computes both
releases side by side precisely so the paper cannot make that error.

This script does not need the model. It needs only the GTEx tissue table, and it
establishes:
  1. kidney cortex's rank by eQTL sample size and by eGene count, among all
     tissues that have eQTLs at all;
  2. the weight kidney cortex carries in an eGene-weighted average, expressed
     both as a raw fraction and relative to an unweighted average;
  3. the same for kidney medulla, which has never had eQTLs called in any
     release and therefore carries weight zero;
  4. how far the weighted and unweighted averages could diverge, given a
     hypothetical kidney-specific performance drop.

Point 4 is deliberately framed as a sensitivity analysis, not a result. We do
not yet know kidney's true performance. What we can show now is how much the
published averaging would hide it if it were poor.

Input : GTEx API v2 (live), cached to data/gtex_v10_tissues.json
Output: results/gtex_tissue_weighting.{json,csv}

Author: Christopher Lawrence
"""

from __future__ import annotations
import json
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

API = ("https://gtexportal.org/api/v2/dataset/tissueSiteDetail"
       "?datasetId={ds}&itemsPerPage=100")

# WHICH RELEASE MATTERS.
#   v8  is the release AlphaGenome was BENCHMARKED on. All statements about the
#       published metric must use v8.
#   v10 is the current release, and is what a reader will see if they look the
#       tissue up today.
# The two give different answers and the difference is not cosmetic:
# in v8 kidney cortex is the smallest tissue of 49 on both measures; in v10 it
# is second smallest of 50, because Bladder was added with an even smaller n.
# Saying "kidney has the fewest eQTLs of any tissue" without naming the release
# is therefore wrong half the time. Both are computed and both are reported.
RELEASES = ("gtex_v8", "gtex_v10")


def fetch(ds: str) -> list[dict]:
    cache = DATA / f"{ds}_tissues.json"
    if cache.exists():
        return json.loads(cache.read_text())["data"]
    req = urllib.request.Request(API.format(ds=ds),
                                 headers={"accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as fh:
        payload = json.load(fh)
    cache.write_text(json.dumps(payload, indent=1))
    return payload["data"]


def summarise(ds: str) -> tuple[pd.DataFrame, dict]:
    rows = fetch(ds)
    recs = []
    for t in rows:
        eq = t.get("eqtlSampleSummary") or {}
        rna = t.get("rnaSeqSampleSummary") or {}
        recs.append({
            "tissue": t["tissueSiteDetailId"],
            "label": t.get("tissueSiteDetail", ""),
            "site": t.get("tissueSite", ""),
            "eqtl_n": eq.get("totalCount"),
            "rnaseq_n": rna.get("totalCount"),
            "has_egenes": bool(t.get("hasEGenes")),
            "egenes": t.get("eGeneCount"),
            "sgenes": t.get("sGeneCount"),
            "expressed_genes": t.get("expressedGeneCount"),
        })
    df = pd.DataFrame(recs)

    # tissues that actually enter an eQTL benchmark
    eq = df[df["has_egenes"] & df["egenes"].notna()].copy()
    eq["egenes"] = eq["egenes"].astype(int)
    eq["eqtl_n"] = eq["eqtl_n"].astype(int)

    eq = eq.sort_values("egenes").reset_index(drop=True)
    eq["rank_by_egenes"] = eq.index + 1
    eq = eq.sort_values("eqtl_n").reset_index(drop=True)
    eq["rank_by_sample_size"] = eq.index + 1

    n_tissues = len(eq)
    total_egenes = int(eq["egenes"].sum())
    eq["weight_egene"] = eq["egenes"] / total_egenes
    eq["weight_unweighted"] = 1.0 / n_tissues
    eq["weight_ratio"] = eq["weight_egene"] / eq["weight_unweighted"]

    eq.sort_values("egenes").to_csv(
        RESULTS / f"gtex_tissue_weighting_{ds}.csv", index=False)

    kc = eq[eq["tissue"] == "Kidney_Cortex"]
    km = df[df["tissue"] == "Kidney_Medulla"]

    if kc.empty:
        raise SystemExit("Kidney_Cortex not found in the eQTL tissue set")
    kc = kc.iloc[0]

    # sensitivity: if kidney's true Spearman were D lower than the mean of the
    # rest, how much of that gap survives into each averaging scheme?
    sens = {}
    for drop in (0.05, 0.10, 0.20):
        # weighted average shifts by weight * drop; unweighted by 1/n * drop
        w_shift = kc["weight_egene"] * drop
        u_shift = (1 / n_tissues) * drop
        sens[f"kidney_deficit_{drop:.2f}"] = {
            "shift_in_egene_weighted_mean": round(w_shift, 5),
            "shift_in_unweighted_mean": round(u_shift, 5),
            "unweighted_reveals_x_more": round(u_shift / w_shift, 2),
        }

    out = {
        "gtex_release": ds,
        "tissues_total": int(len(df)),
        "tissues_with_eqtls": n_tissues,
        "kidney_cortex": {
            "eqtl_n": int(kc["eqtl_n"]),
            "rnaseq_n": int(kc["rnaseq_n"]),
            "egenes": int(kc["egenes"]),
            "rank_by_sample_size": int(kc["rank_by_sample_size"]),
            "rank_by_egenes": int(kc["rank_by_egenes"]),
            "is_smallest_by_sample_size": int(kc["rank_by_sample_size"]) == 1,
            "is_smallest_by_egenes": int(kc["rank_by_egenes"]) == 1,
            "weight_in_egene_weighted_mean": round(float(kc["weight_egene"]), 5),
            "weight_if_unweighted": round(1 / n_tissues, 5),
            "down_weighted_fold": round(1 / float(kc["weight_ratio"]), 2),
        },
        "kidney_medulla": {
            "rnaseq_n": int(km["rnaseq_n"].iloc[0]) if len(km) and pd.notna(
                km["rnaseq_n"].iloc[0]) else 0,
            "eqtl_n": int(km["eqtl_n"].iloc[0]) if len(km) and pd.notna(
                km["eqtl_n"].iloc[0]) else 0,
            "has_egenes": bool(km["has_egenes"].iloc[0]) if len(km) else False,
        },
        "median_eqtl_n": int(eq["eqtl_n"].median()),
        "median_egenes": int(eq["egenes"].median()),
        "smallest_five_by_sample_size": eq.nsmallest(5, "eqtl_n")[
            ["tissue", "eqtl_n", "egenes"]].to_dict(orient="records"),
        "sensitivity": sens,
    }
    return eq, out


def report(ds: str, out: dict) -> None:
    k = out["kidney_cortex"]
    n_tissues = out["tissues_with_eqtls"]
    tag = "BENCHMARK RELEASE" if ds == "gtex_v8" else "CURRENT RELEASE"
    print("=" * 74)
    print(f"{ds.upper().replace('GTEX_','GTEx ')}  ({tag})")
    print("=" * 74)
    print(f"  tissues with eQTLs called     : {n_tissues}")
    print(f"  Kidney cortex eQTL samples    : {k['eqtl_n']}   "
          f"(median {out['median_eqtl_n']})")
    print(f"  Kidney cortex eGenes          : {k['egenes']:,}   "
          f"(median {out['median_egenes']:,})")
    print(f"  rank by sample size           : {k['rank_by_sample_size']} of "
          f"{n_tissues}" + ("   <- SMALLEST" if k['is_smallest_by_sample_size'] else ""))
    print(f"  rank by eGene count           : {k['rank_by_egenes']} of "
          f"{n_tissues}" + ("   <- SMALLEST" if k['is_smallest_by_egenes'] else ""))
    print(f"  weight in PUBLISHED metric    : {k['weight_in_egene_weighted_mean']:.5f}")
    print(f"  weight if unweighted          : {k['weight_if_unweighted']:.5f}")
    print(f"  DOWN-WEIGHTED BY              : {k['down_weighted_fold']}x")
    print()
    print("  five smallest tissues by eQTL sample size:")
    for r in out["smallest_five_by_sample_size"]:
        star = "  <-- kidney" if r["tissue"] == "Kidney_Cortex" else ""
        print(f"    {r['tissue']:<34} n={r['eqtl_n']:>4}  "
              f"eGenes={r['egenes']:,}{star}")
    km = out["kidney_medulla"]
    print()
    print(f"  Kidney MEDULLA: rnaseq n={km['rnaseq_n']}, eQTL n={km['eqtl_n']}, "
          f"eGenes called={km['has_egenes']}")
    print()


def main() -> None:
    allout = {}
    for ds in RELEASES:
        eq, out = summarise(ds)
        allout[ds] = out
        report(ds, out)

    v8, v10 = allout["gtex_v8"], allout["gtex_v10"]
    print("=" * 74)
    print("THE STATEMENT THE PAPER CAN MAKE, AND THE ONE IT CANNOT")
    print("=" * 74)
    print("  CAN say, of the release the model was benchmarked on (v8):")
    print(f"    kidney cortex is the SMALLEST of {v8['tissues_with_eqtls']} tissues on")
    print(f"    both sample size (n={v8['kidney_cortex']['eqtl_n']}) and eGene count")
    print(f"    ({v8['kidney_cortex']['egenes']:,}), and an eGene-weighted average")
    print(f"    down-weights it {v8['kidney_cortex']['down_weighted_fold']}x "
          "relative to an unweighted one.")
    print()
    print("  CANNOT say it of the current release (v10):")
    print(f"    kidney cortex is rank {v10['kidney_cortex']['rank_by_sample_size']} "
          f"of {v10['tissues_with_eqtls']}; Bladder is smaller. Down-weighting is")
    print(f"    {v10['kidney_cortex']['down_weighted_fold']}x, still substantial but not extreme.")
    print()
    print("  ALWAYS NAME THE RELEASE. An unqualified 'fewest eQTLs of any")
    print("  tissue' is false as of v10 and a reviewer will catch it.")
    print()
    print("  Kidney medulla carries weight ZERO in both: never eQTL-called.")
    print()
    (RESULTS / "gtex_tissue_weighting.json").write_text(json.dumps(allout, indent=2))
    print(f"Wrote {RESULTS}/gtex_tissue_weighting.json and per-release CSVs")


if __name__ == "__main__":
    main()
