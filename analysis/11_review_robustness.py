#!/usr/bin/env python3
"""
11_review_robustness.py

Responses to the statistical objections raised in review.

THREE ANALYSES.

  A. INDEPENDENCE (Reviewer 2, point 5). Observations may share a gene, a locus
     or a credible set. Ordinary correlation and binomial intervals then
     overstate precision. We report the number of unique genes and of
     independent credible sets per tissue, and recompute kidney's interval by
     cluster bootstrap resampling whole credible sets rather than variants.

  B. SMALL-SAMPLE EXPECTATION (Reviewer 2, point 6). Kidney's interval is 3.2
     times wider than the median. The reviewer's objection is that this is the
     arithmetic consequence of n=58, not a kidney-specific defect. We test that
     directly: repeatedly downsample every other tissue to kidney's n and ask
     where kidney's observed correlation and interval width fall in the
     resulting distribution. If kidney is unremarkable once size is matched,
     the paper must say so.

  C. SELECTION OF THE MAXIMUM (Reviewer 2, point 3). Kidney has the highest
     directional agreement of 49 tissues. Selecting the maximum of 49 noisy
     estimates is a biased procedure. We ask how often a tissue with kidney's
     sample size attains the observed maximum by chance under the null that all
     tissues share a common underlying rate.

The purpose of all three is to establish what the data cannot support, not to
rescue the original claims.

Output: results/review_robustness.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCORES = RESULTS / "scores"
BENCH = RESULTS / "benchmark_variants"

RNG = np.random.default_rng(20260914)
B = 5000


def load() -> dict[str, pd.DataFrame]:
    out = {}
    for f in sorted(glob.glob(str(SCORES / "*.tsv.gz"))):
        g = os.path.basename(f).replace(".tsv.gz", "")
        d = pd.read_csv(f, sep="\t").dropna(subset=["predicted_lfc", "observed_beta"])
        b = BENCH / f"{g}.tsv.gz"
        if b.exists():
            cs = pd.read_csv(b, sep="\t")[["variant", "molecular_trait_id", "cs_id"]]
            d = d.merge(cs, left_on=["variant", "gene"],
                        right_on=["variant", "molecular_trait_id"], how="left")
        out[g] = d
    return out


def spearman_ci(rho: float, n: int) -> tuple[float, float]:
    if n < 6:
        return float("nan"), float("nan")
    r = max(min(rho, 0.999999), -0.999999)
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3))
    z = np.arctanh(r)
    return float(np.tanh(z - 1.959964 * se)), float(np.tanh(z + 1.959964 * se))


def main() -> None:
    T = load()
    out: dict = {}

    # ------------------------------------------------------------------ A
    print("=" * 76)
    print("A. INDEPENDENCE OF OBSERVATIONS")
    print("=" * 76)
    rows = []
    for g, d in T.items():
        rows.append({"tissue": g, "n": len(d),
                     "unique_genes": d["gene"].nunique(),
                     "unique_cs": d["cs_id"].nunique() if "cs_id" in d else np.nan})
    idp = pd.DataFrame(rows)
    idp["obs_per_gene"] = (idp["n"] / idp["unique_genes"]).round(3)
    k = idp[idp["tissue"] == "kidney_cortex"].iloc[0]
    print(f"  kidney: {int(k['n'])} observations, {int(k['unique_genes'])} unique genes, "
          f"{int(k['unique_cs'])} unique credible sets")
    print(f"  across all tissues, observations per gene: median "
          f"{idp['obs_per_gene'].median():.3f}, max {idp['obs_per_gene'].max():.3f}")
    near_one = (idp["obs_per_gene"] < 1.05).mean()
    print(f"  fraction of tissues with <1.05 observations per gene: {near_one:.0%}")

    kd = T["kidney_cortex"]
    if "cs_id" in kd and kd["cs_id"].notna().any():
        groups = [x for _, x in kd.groupby("cs_id")]
        boots = []
        for _ in range(B):
            pick = RNG.integers(0, len(groups), len(groups))
            s = pd.concat([groups[i] for i in pick])
            if len(s) > 5 and s["predicted_lfc"].nunique() > 2:
                boots.append(stats.spearmanr(s["predicted_lfc"],
                                             s["observed_beta"]).statistic)
        lo, hi = np.nanpercentile(boots, [2.5, 97.5])
        rho_k = stats.spearmanr(kd["predicted_lfc"], kd["observed_beta"]).statistic
        plo, phi = spearman_ci(rho_k, len(kd))
        print()
        print(f"  kidney rho {rho_k:.3f}")
        print(f"    parametric 95% CI      [{plo:.3f}, {phi:.3f}]  width {phi-plo:.3f}")
        print(f"    credible-set bootstrap [{lo:.3f}, {hi:.3f}]  width {hi-lo:.3f}")
        out["independence"] = {
            "kidney_n": int(len(kd)), "kidney_genes": int(kd["gene"].nunique()),
            "kidney_cs": int(kd["cs_id"].nunique()),
            "rho": float(rho_k),
            "parametric_ci": [round(plo, 3), round(phi, 3)],
            "cluster_bootstrap_ci": [round(float(lo), 3), round(float(hi), 3)],
            "median_obs_per_gene_all_tissues": float(idp["obs_per_gene"].median()),
        }

    # ------------------------------------------------------------------ B
    print()
    print("=" * 76)
    print("B. IS KIDNEY'S WIDE INTERVAL SIMPLY ITS SAMPLE SIZE?")
    print("=" * 76)
    n_k = len(kd)
    rho_k = stats.spearmanr(kd["predicted_lfc"], kd["observed_beta"]).statistic
    w_k = spearman_ci(rho_k, n_k)
    w_k = w_k[1] - w_k[0]

    rhos, widths = [], []
    for g, d in T.items():
        if g == "kidney_cortex" or len(d) < n_k:
            continue
        for _ in range(200):
            s = d.sample(n_k, random_state=int(RNG.integers(0, 1 << 31)))
            r = stats.spearmanr(s["predicted_lfc"], s["observed_beta"]).statistic
            if np.isfinite(r):
                rhos.append(r)
                a, b2 = spearman_ci(r, n_k)
                widths.append(b2 - a)
    rhos, widths = np.array(rhos), np.array(widths)
    pct_rho = float((rhos < rho_k).mean() * 100)
    pct_w = float((widths < w_k).mean() * 100)
    print(f"  kidney: n={n_k}, rho={rho_k:.3f}, CI width={w_k:.3f}")
    print(f"  other tissues downsampled to n={n_k} ({len(rhos):,} replicates):")
    print(f"    rho        median {np.median(rhos):.3f}, "
          f"2.5-97.5th [{np.percentile(rhos,2.5):.3f}, {np.percentile(rhos,97.5):.3f}]")
    print(f"    CI width   median {np.median(widths):.3f}, "
          f"2.5-97.5th [{np.percentile(widths,2.5):.3f}, {np.percentile(widths,97.5):.3f}]")
    print(f"  kidney's rho sits at the {pct_rho:.0f}th percentile of size-matched tissues")
    print(f"  kidney's CI width sits at the {pct_w:.0f}th percentile")
    print()
    if 5 < pct_w < 95:
        print("  VERDICT: once sample size is matched, kidney's interval width is")
        print("  UNREMARKABLE. The 3.2-fold difference is the arithmetic of n=58,")
        print("  not a kidney-specific property. The reviewer is right and the")
        print("  manuscript must say so.")
    out["downsampling"] = {
        "kidney_n": int(n_k), "kidney_rho": float(rho_k), "kidney_ci_width": float(w_k),
        "matched_rho_median": float(np.median(rhos)),
        "matched_width_median": float(np.median(widths)),
        "kidney_rho_percentile": round(pct_rho, 1),
        "kidney_width_percentile": round(pct_w, 1),
        "replicates": int(len(rhos)),
    }

    # ------------------------------------------------------------------ C
    print()
    print("=" * 76)
    print("C. IS KIDNEY'S TOP DIRECTIONAL AGREEMENT A SELECTION ARTEFACT?")
    print("=" * 76)
    agree = {g: float((np.sign(d["predicted_lfc"]) == np.sign(d["observed_beta"])).mean())
             for g, d in T.items()}
    ns = {g: len(d) for g, d in T.items()}
    k_ag = agree["kidney_cortex"]
    pooled = float(np.average(list(agree.values()), weights=list(ns.values())))
    print(f"  kidney directional agreement {k_ag:.3f} (n={ns['kidney_cortex']})")
    print(f"  pooled rate across tissues   {pooled:.3f}")

    hits = 0
    for _ in range(B):
        sim = {g: RNG.binomial(n, pooled) / n for g, n in ns.items()}
        if max(sim, key=sim.get) == "kidney_cortex":
            hits += 1
    p_top = hits / B
    print(f"  under a common underlying rate, kidney is the top tissue in "
          f"{p_top:.1%} of {B:,} simulations")
    print()
    if p_top > 0.05:
        print("  VERDICT: kidney topping the table is WELL WITHIN what chance")
        print("  produces for the smallest tissue, because small samples have the")
        print("  most variable estimates. The claim of highest directional")
        print("  agreement must be withdrawn as evidence of superiority.")
    out["max_selection"] = {"kidney_agreement": k_ag, "pooled_rate": pooled,
                            "p_kidney_is_max_by_chance": p_top}

    (RESULTS / "review_robustness.json").write_text(json.dumps(out, indent=2))
    print()
    print(f"  Wrote {RESULTS}/review_robustness.json")


if __name__ == "__main__":
    main()
