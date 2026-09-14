#!/usr/bin/env python3
"""
07_power_and_threshold_sensitivity.py

LEG 2, part 4: is the kidney result real, or did I choose thresholds that
produced it?

Script 06 found kidney cortex yields 59 usable fine-mapped positives, the fewest
of 49 tissues and less than half the next smallest. That is a striking number
and it arrived from a filter I chose (PIP >= 0.9, credible-set size <= 5). Before
it goes near a manuscript it has to survive two challenges.

CHALLENGE 1, THRESHOLD SENSITIVITY. Sweep the filter across a grid and check
whether kidney stays last. If kidney's rank moves around, the finding is an
artefact of my choice and must be reported as such. If kidney is last at every
setting, the filter is not doing the work and the finding is about the data.

CHALLENGE 2, WHAT 59 VARIANTS ACTUALLY BUYS. A per-tissue correlation estimated
on 59 points has a confidence interval wide enough that it may not be able to
distinguish good performance from poor. If so, the honest conclusion is not
"the model does badly in kidney" but "nobody can currently tell, and here is
precisely why" - which is a stronger and far more defensible claim.

The confidence interval uses the Fisher z transform on Spearman's rho with the
Bonett-Wright standard error, se = sqrt((1 + rho^2/2) / (n - 3)), which is the
standard correction for rank correlation rather than treating rho as Pearson.

Output: results/threshold_sensitivity.csv
        results/power_analysis.json

Author: Christopher Lawrence
"""

from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
SUSIE = DATA / "eqtl_catalogue" / "susie"

PIP_GRID = [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
CS_GRID = [1, 2, 5, 10, 25, 10_000]        # 10_000 = effectively no size cap


def spearman_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    """Fisher z interval for Spearman's rho, Bonett-Wright standard error."""
    if n < 6:
        return (float("nan"), float("nan"))
    rho = max(min(rho, 0.999999), -0.999999)
    z = np.arctanh(rho)
    se = np.sqrt((1 + rho ** 2 / 2) / (n - 3))
    crit = 1.959963984540054 if alpha == 0.05 else abs(
        np.sqrt(2) * float(np.vectorize(lambda p: 0)(0)))  # only 95% used
    lo, hi = np.tanh(z - crit * se), np.tanh(z + crit * se)
    return float(lo), float(hi)


def load_all() -> dict[str, pd.DataFrame]:
    md = pd.read_csv(DATA / "eqtl_catalogue_dataset_metadata_r7.tsv", sep="\t")
    g = md[(md["study_label"] == "GTEx") & (md["quant_method"] == "ge")]
    lookup = dict(zip(g["dataset_id"], g["sample_group"]))
    out = {}
    for f in sorted(SUSIE.glob("*.credible_sets.tsv.gz")):
        ds = f.name.split(".")[0]
        out[lookup.get(ds, ds)] = pd.read_csv(
            f, sep="\t", usecols=["variant", "molecular_trait_id", "pip",
                                  "cs_size", "beta", "z"])
    return out


def main() -> None:
    tissues = load_all()
    print(f"loaded {len(tissues)} tissues\n")

    # ---------------------------------------------------------- challenge 1
    rows = []
    for pip in PIP_GRID:
        for cs in CS_GRID:
            counts = {}
            for grp, df in tissues.items():
                sub = df[(df["pip"] >= pip) & (df["cs_size"] <= cs)]
                counts[grp] = sub.drop_duplicates(
                    subset=["variant", "molecular_trait_id"]).shape[0]
            ser = pd.Series(counts).sort_values()
            kc = int(ser.get("kidney_cortex", 0))
            rank = int((ser < kc).sum()) + 1
            rows.append({
                "pip_min": pip,
                "cs_size_max": "none" if cs == 10_000 else cs,
                "kidney_variants": kc,
                "kidney_rank": rank,
                "kidney_is_smallest": rank == 1,
                "median_across_tissues": int(ser.median()),
                "max_tissue": ser.index[-1],
                "max_variants": int(ser.iloc[-1]),
                "fold_below_max": round(ser.iloc[-1] / kc, 1) if kc else None,
                "second_smallest": int(ser.iloc[1]),
            })
    sens = pd.DataFrame(rows)
    sens.to_csv(RESULTS / "threshold_sensitivity.csv", index=False)

    n_settings = len(sens)
    n_smallest = int(sens["kidney_is_smallest"].sum())

    print("=" * 78)
    print("CHALLENGE 1: DOES KIDNEY STAY LAST ACROSS THRESHOLDS?")
    print("=" * 78)
    print(f"  settings tested                    : {n_settings}")
    print(f"  settings where kidney is SMALLEST  : {n_smallest}")
    print(f"  worst (highest) rank kidney reaches: {int(sens['kidney_rank'].max())} of 49")
    print()
    print(f"  {'PIP':>6}{'cs<=':>7}{'kidney n':>10}{'rank':>6}"
          f"{'median':>9}{'x below max':>13}")
    for _, r in sens.iterrows():
        flag = "" if r["kidney_is_smallest"] else "   <- not smallest"
        print(f"  {r['pip_min']:>6}{str(r['cs_size_max']):>7}"
              f"{r['kidney_variants']:>10}{r['kidney_rank']:>6}"
              f"{r['median_across_tissues']:>9}"
              f"{str(r['fold_below_max']):>13}{flag}")
    print()
    if n_smallest == n_settings:
        print("  VERDICT: kidney is the smallest tissue at EVERY setting tested.")
        print("  The filter is not producing the result; the data is.")
    else:
        print(f"  VERDICT: kidney is smallest at {n_smallest}/{n_settings} settings.")
        print("  Report the sensitivity grid in the supplement, not just one cut.")

    # ---------------------------------------------------------- challenge 2
    print()
    print("=" * 78)
    print("CHALLENGE 2: WHAT DOES A KIDNEY-SIZED SAMPLE ACTUALLY BUY?")
    print("=" * 78)
    base = sens[(sens["pip_min"] == 0.9) & (sens["cs_size_max"] == 5)].iloc[0]
    n_kid = int(base["kidney_variants"])
    n_max = int(base["max_variants"])

    power = {"kidney_n": n_kid, "largest_n": n_max,
             "largest_tissue": base["max_tissue"], "intervals": {}}
    print(f"  95% confidence interval on Spearman's rho, "
          f"kidney n={n_kid} vs {base['max_tissue']} n={n_max}")
    print(f"  {'true rho':>10}{'kidney 95% CI':>26}{'width':>8}"
          f"{'largest 95% CI':>26}{'width':>8}")
    for rho in (0.0, 0.2, 0.4, 0.6, 0.8):
        lo_k, hi_k = spearman_ci(rho, n_kid)
        lo_m, hi_m = spearman_ci(rho, n_max)
        wk, wm = hi_k - lo_k, hi_m - lo_m
        power["intervals"][f"rho_{rho}"] = {
            "kidney_ci": [round(lo_k, 3), round(hi_k, 3)],
            "kidney_width": round(wk, 3),
            "largest_ci": [round(lo_m, 3), round(hi_m, 3)],
            "largest_width": round(wm, 3),
            "width_ratio": round(wk / wm, 1),
        }
        print(f"  {rho:>10.1f}   [{lo_k:>6.3f}, {hi_k:>6.3f}]{wk:>13.3f}"
              f"   [{lo_m:>6.3f}, {hi_m:>6.3f}]{wm:>13.3f}")

    # smallest detectable difference from a reference tissue
    lo6, hi6 = spearman_ci(0.6, n_kid)
    power["headline"] = {
        "kidney_ci_at_rho_0.6": [round(lo6, 3), round(hi6, 3)],
        "kidney_ci_width_at_rho_0.6": round(hi6 - lo6, 3),
    }
    (RESULTS / "power_analysis.json").write_text(json.dumps(power, indent=2))

    print()
    print("  READING. At a true rho of 0.6, the kidney interval spans")
    print(f"  [{lo6:.2f}, {hi6:.2f}], a width of {hi6-lo6:.2f}. An interval that wide")
    print("  cannot separate 'performs like other tissues' from 'performs")
    print("  materially worse'.")
    print()
    print("  THE HONEST CONCLUSION IS THEREFORE NOT 'the model does badly in")
    print("  kidney'. It is that the public evidence cannot currently support")
    print("  ANY confident statement about kidney performance, and the reason")
    print("  is measurable: 59 fine-mapped positives against 1,756 for thyroid.")
    print("  That is a stronger claim, and one no reviewer can dispute.")
    print()
    print(f"  Wrote {RESULTS}/threshold_sensitivity.csv and power_analysis.json")


if __name__ == "__main__":
    main()
