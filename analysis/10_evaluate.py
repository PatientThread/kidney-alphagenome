#!/usr/bin/env python3
"""
10_evaluate.py

LEG 2, part 6: turn scores into the paper's central comparison.

For each tissue with scored variants, compute agreement between the model's
predicted log fold change and the observed eQTL effect:

    spearman        rank correlation, the metric the AlphaGenome paper reports
    sign_agreement  fraction of variants where predicted and observed effects
                    point the same way. This is the more interpretable quantity
                    for a clinician and is robust to the scale mismatch between
                    a predicted LFC (order 0.01) and an eQTL beta (order 1).

Then average across tissues TWO ways and report the difference:

    eGene-weighted  what the published headline does
    unweighted      one tissue, one vote

WHY BOTH. The published average weights each tissue by its eQTL count. Kidney
cortex has the fewest in the release the model was benchmarked on, so it is
down-weighted 7.71-fold. If kidney performs differently from the rest, the
published metric is the one least able to show it.

CONFIDENCE INTERVALS ARE MANDATORY HERE, NOT OPTIONAL. Kidney has 58 usable
variants against a median of 616. A point estimate without an interval invites
exactly the overreach this project is trying to avoid. Intervals use the Fisher
z transform with the Bonett-Wright standard error for rank correlation, and
sign agreement uses a Wilson interval.

Output: results/evaluation.json, results/per_tissue_metrics.csv

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
SCORES = RESULTS / "scores"

Z95 = 1.959963984540054


def spearman_ci(rho: float, n: int) -> tuple[float, float]:
    if n < 6 or not np.isfinite(rho):
        return (float("nan"), float("nan"))
    r = max(min(rho, 0.999999), -0.999999)
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3))
    z = np.arctanh(r)
    return float(np.tanh(z - Z95 * se)), float(np.tanh(z + Z95 * se))


def wilson(k: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + Z95 ** 2 / n
    c = (p + Z95 ** 2 / (2 * n)) / d
    h = Z95 * np.sqrt(p * (1 - p) / n + Z95 ** 2 / (4 * n ** 2)) / d
    return float(c - h), float(c + h)


def main() -> None:
    files = sorted(SCORES.glob("*.tsv.gz"))
    if not files:
        raise SystemExit("no scored tissues yet; run 09_score_variants.py first")

    rows = []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        d = d.dropna(subset=["predicted_lfc", "observed_beta"])
        n = len(d)
        if n < 6:
            continue
        rho, p = stats.spearmanr(d["predicted_lfc"], d["observed_beta"])
        lo, hi = spearman_ci(rho, n)
        agree = int((np.sign(d["predicted_lfc"]) ==
                     np.sign(d["observed_beta"])).sum())
        slo, shi = wilson(agree, n)
        rows.append({
            "sample_group": f.name.replace(".tsv.gz", ""),
            "n": n,
            "spearman": round(float(rho), 4),
            "spearman_lo": round(lo, 4),
            "spearman_hi": round(hi, 4),
            "spearman_ci_width": round(hi - lo, 4),
            "spearman_p": float(p),
            "sign_agreement": round(agree / n, 4),
            "sign_lo": round(slo, 4),
            "sign_hi": round(shi, 4),
            "sign_ci_width": round(shi - slo, 4),
        })

    m = pd.DataFrame(rows).sort_values("n")
    m.to_csv(RESULTS / "per_tissue_metrics.csv", index=False)

    # weights: eGene counts from the v8 tissue table (benchmark release)
    # Weights must use the SAME curated tissue map as the scoring step. A bare
    # normalisation merge silently drops the 16 hand-mapped tissues (blood, LCL,
    # the brain subregions...), which would compute the "published style" average
    # on a biased subset. Import the map rather than re-deriving it.
    import importlib.util
    _spec = importlib.util.spec_from_file_location(
        "score_mod", Path(__file__).resolve().parent / "09_score_variants.py")
    _score = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_score)

    def _key(x):
        import re
        return re.sub(r"[^a-z0-9]", "", str(x).lower())

    w_path = RESULTS / "gtex_tissue_weighting_gtex_v8.csv"
    weighted = None
    n_weighted = 0
    if w_path.exists():
        wt = pd.read_csv(w_path)
        by_key = {_key(t): e for t, e in zip(wt["tissue"], wt["egenes"])}

        def egenes_for(grp: str):
            gtex = _score.MANUAL_TISSUE_MAP.get(grp, grp)
            return by_key.get(_key(gtex))

        m["egenes"] = m["sample_group"].map(egenes_for)
        ok = m.dropna(subset=["egenes"])
        n_weighted = len(ok)
        if n_weighted >= 2:
            weighted = float(np.average(ok["spearman"], weights=ok["egenes"]))

    unweighted = float(m["spearman"].mean())

    kid = m[m["sample_group"] == "kidney_cortex"]
    others = m[m["sample_group"] != "kidney_cortex"]

    out = {
        "tissues_evaluated": int(len(m)),
        "variants_total": int(m["n"].sum()),
        "unweighted_mean_spearman": round(unweighted, 4),
        "egene_weighted_mean_spearman": (round(weighted, 4)
                                         if weighted is not None else None),
        "difference": (round(unweighted - weighted, 4)
                       if weighted is not None else None),
        "mean_sign_agreement": round(float(m["sign_agreement"].mean()), 4),
        "tissues_with_weights": int(n_weighted),
        "spearman_spread": {
            "min": round(float(m["spearman"].min()), 4),
            "max": round(float(m["spearman"].max()), 4),
            "sd": round(float(m["spearman"].std()), 4),
        },
    }
    if len(kid):
        k = kid.iloc[0]
        out["kidney"] = {
            "n": int(k["n"]),
            "spearman": float(k["spearman"]),
            "ci": [float(k["spearman_lo"]), float(k["spearman_hi"])],
            "ci_width": float(k["spearman_ci_width"]),
            "sign_agreement": float(k["sign_agreement"]),
            "sign_ci": [float(k["sign_lo"]), float(k["sign_hi"])],
            "rank_by_spearman": int((m["spearman"] < k["spearman"]).sum()) + 1,
            "median_ci_width_other_tissues": (
                round(float(others["spearman_ci_width"].median()), 4)
                if len(others) else None),
        }
    (RESULTS / "evaluation.json").write_text(json.dumps(out, indent=2))

    print("=" * 86)
    print(f"PER-TISSUE EVALUATION  ({len(m)} tissues, "
          f"{int(m['n'].sum()):,} variants)")
    print("=" * 86)
    print(f"  {'tissue':<34}{'n':>6}{'rho':>8}{'95% CI':>20}"
          f"{'width':>8}{'sign':>8}")
    for _, r in m.iterrows():
        star = "  <-- KIDNEY" if r["sample_group"] == "kidney_cortex" else ""
        print(f"  {r['sample_group']:<34}{int(r['n']):>6}{r['spearman']:>8.3f}"
              f"   [{r['spearman_lo']:>6.3f},{r['spearman_hi']:>6.3f}]"
              f"{r['spearman_ci_width']:>8.3f}{r['sign_agreement']:>8.2f}{star}")
    print()
    print("-" * 86)
    print(f"  unweighted mean rho      : {unweighted:.4f}")
    if weighted is not None:
        print(f"  eGene-weighted mean rho  : {weighted:.4f}   (the published style,"
              f" {n_weighted}/{len(m)} tissues weighted)")
        print(f"  difference               : {unweighted - weighted:+.4f}")
        sd = float(m['spearman'].std())
        print(f"  spread across tissues    : sd {sd:.4f} "
              f"(min {m['spearman'].min():.3f}, max {m['spearman'].max():.3f})")
        if abs(unweighted - weighted) < 0.01:
            print()
            print("  *** THE WEIGHTING ARGUMENT DOES NOT SURVIVE. ***")
            print("  Performance is near-uniform across tissues, so ANY weighting")
            print("  scheme gives nearly the same mean. The down-weighting of")
            print("  kidney is real but numerically inconsequential. Say so, and")
            print("  drop it as a headline claim.")
    if "kidney" in out:
        k = out["kidney"]
        print()
        print(f"  KIDNEY rho = {k['spearman']:.3f}, 95% CI "
              f"[{k['ci'][0]:.3f}, {k['ci'][1]:.3f}], width {k['ci_width']:.3f}")
        print(f"    rank {k['rank_by_spearman']} of {len(m)} tissues")
        print(f"    median CI width in other tissues: "
              f"{k['median_ci_width_other_tissues']}")
        print(f"    sign agreement {k['sign_agreement']:.2f} "
              f"[{k['sign_ci'][0]:.2f}, {k['sign_ci'][1]:.2f}]")
        if k["ci"][0] <= 0 <= k["ci"][1]:
            print("    NOTE: the kidney interval INCLUDES ZERO. On this evidence")
            print("    kidney performance is not distinguishable from chance,")
            print("    and equally not distinguishable from the other tissues.")
            print("    Report it that way. Do NOT claim the model fails in kidney.")
    print()
    print(f"  Wrote {RESULTS}/evaluation.json and per_tissue_metrics.csv")


if __name__ == "__main__":
    main()
