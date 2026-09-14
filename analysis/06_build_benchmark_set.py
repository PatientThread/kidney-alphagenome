#!/usr/bin/env python3
"""
06_build_benchmark_set.py

LEG 2, part 3: build the per-tissue variant set that an UNWEIGHTED benchmark
runs on, and quantify how badly kidney is served by it.

THE DESIGN. A variant-effect benchmark needs positives with real causal support,
not merely significant associations, because in a linkage block hundreds of
variants share a p-value and only one or a few are causal. SuSiE fine-mapping
gives exactly that: per gene, a credible set with a posterior inclusion
probability on each variant.

Positive set, per tissue:
    variants in a SuSiE credible set with PIP >= --pip (default 0.9)
    and credible-set size <= --max-cs (default 5)

Both filters matter. PIP alone admits variants from 200-member credible sets
where the posterior is diffuse; the size cap requires that fine-mapping actually
resolved the signal. The observed effect (beta, z) is the regression target.

WHAT THIS SCRIPT DELIBERATELY DOES NOT DO. It does not score anything. Model
predictions are a separate, gated step (see §"NEXT" below). Building and
publishing the evaluation set first, before any model is run, is what stops the
thresholds being tuned to a flattering answer. The set is frozen here.

THE KIDNEY QUESTION THIS ANSWERS ON ITS OWN. Even before any model runs, the
benchmark set tells you how much evidence kidney can possibly contribute. If
kidney yields a few dozen usable positives against thousands for whole blood,
then no per-tissue metric computed on it will be stable, and that is itself a
reportable limit on what anyone can claim about kidney performance.

Output: results/benchmark_set_summary.csv
        results/benchmark_variants/<sample_group>.tsv.gz

Author: Christopher Lawrence
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
SUSIE = DATA / "eqtl_catalogue" / "susie"
OUTDIR = RESULTS / "benchmark_variants"
OUTDIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pip", type=float, default=0.9,
                    help="minimum posterior inclusion probability")
    ap.add_argument("--max-cs", type=int, default=5,
                    help="maximum credible-set size")
    args = ap.parse_args()

    md = pd.read_csv(DATA / "eqtl_catalogue_dataset_metadata_r7.tsv", sep="\t")
    g = md[(md["study_label"] == "GTEx") & (md["quant_method"] == "ge")]
    lookup = dict(zip(g["dataset_id"], g["sample_group"]))
    sizes = dict(zip(g["dataset_id"], g["sample_size"]))

    rows = []
    for f in sorted(SUSIE.glob("*.credible_sets.tsv.gz")):
        ds = f.name.split(".")[0]
        grp = lookup.get(ds, ds)
        df = pd.read_csv(f, sep="\t")

        n_all = len(df)
        n_cs = df["cs_id"].nunique()
        keep = df[(df["pip"] >= args.pip) & (df["cs_size"] <= args.max_cs)].copy()
        keep = keep.drop_duplicates(subset=["variant", "molecular_trait_id"])

        if len(keep):
            keep["sample_group"] = grp
            keep["dataset_id"] = ds
            cols = ["sample_group", "dataset_id", "variant", "rsid",
                    "molecular_trait_id", "gene_id", "cs_id", "cs_size",
                    "pip", "beta", "se", "z", "pvalue", "region"]
            keep[[c for c in cols if c in keep.columns]].to_csv(
                OUTDIR / f"{grp}.tsv.gz", sep="\t", index=False,
                compression="gzip")

        rows.append({
            "sample_group": grp,
            "dataset_id": ds,
            "sample_size": sizes.get(ds),
            "credible_set_rows": n_all,
            "credible_sets": n_cs,
            "benchmark_variants": len(keep),
            "distinct_genes": int(keep["molecular_trait_id"].nunique()) if len(keep) else 0,
            "median_abs_beta": round(float(keep["beta"].abs().median()), 4) if len(keep) else None,
        })

    s = pd.DataFrame(rows).sort_values("benchmark_variants")
    s["pct_of_max"] = (100 * s["benchmark_variants"] /
                       s["benchmark_variants"].max()).round(1)
    s.to_csv(RESULTS / "benchmark_set_summary.csv", index=False)

    kc = s[s["sample_group"] == "kidney_cortex"]
    n_tissues = len(s)

    print("=" * 78)
    print(f"BENCHMARK SET  (PIP >= {args.pip}, credible-set size <= {args.max_cs})")
    print("=" * 78)
    print(f"  tissues built            : {n_tissues}")
    print(f"  total benchmark variants : {int(s['benchmark_variants'].sum()):,}")
    print(f"  median per tissue        : {int(s['benchmark_variants'].median()):,}")
    print()
    print("  TEN SMALLEST TISSUES BY USABLE POSITIVES")
    print(f"  {'tissue':<34}{'n':>6}{'variants':>10}{'genes':>8}{'% of max':>10}")
    for _, r in s.head(10).iterrows():
        star = "  <-- KIDNEY" if r["sample_group"] == "kidney_cortex" else ""
        print(f"  {r['sample_group']:<34}{int(r['sample_size']):>6}"
              f"{int(r['benchmark_variants']):>10}{int(r['distinct_genes']):>8}"
              f"{r['pct_of_max']:>9.1f}%{star}")
    print()
    print("  THREE LARGEST")
    for _, r in s.tail(3).iterrows():
        print(f"  {r['sample_group']:<34}{int(r['sample_size']):>6}"
              f"{int(r['benchmark_variants']):>10}{int(r['distinct_genes']):>8}"
              f"{r['pct_of_max']:>9.1f}%")
    print()

    if len(kc):
        k = kc.iloc[0]
        rank = int((s["benchmark_variants"] < k["benchmark_variants"]).sum()) + 1
        biggest = s.iloc[-1]
        print("-" * 78)
        print("  KIDNEY CORTEX")
        print(f"    usable benchmark variants : {int(k['benchmark_variants']):,}")
        print(f"    distinct genes            : {int(k['distinct_genes']):,}")
        print(f"    rank among {n_tissues} tissues     : {rank}")
        print(f"    versus largest ({biggest['sample_group']}): "
              f"{int(biggest['benchmark_variants']):,} "
              f"= {biggest['benchmark_variants'] / max(k['benchmark_variants'],1):.0f}x more")
        print()
        print("  CONSEQUENCE, and it is reportable without running any model:")
        print("  a per-tissue metric estimated on this many variants carries a")
        print("  much wider confidence interval than one estimated on thousands.")
        print("  Kidney is not merely down-weighted in the published average, it")
        print("  is the tissue where any per-tissue estimate is least stable.")
    print()
    print(f"  Wrote {RESULTS}/benchmark_set_summary.csv")
    print(f"  Wrote per-tissue variant files to {OUTDIR}/")
    print()
    print("  NEXT: these variants need MODEL PREDICTIONS. That step is gated on")
    print("  an AlphaGenome API key (free, non-commercial research use only) or")
    print("  on Atlas precomputed scores. Terms must be accepted personally, not")
    print("  on behalf of a company. Nothing here has run the model.")


if __name__ == "__main__":
    main()
