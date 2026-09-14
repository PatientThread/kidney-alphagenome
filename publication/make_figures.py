#!/usr/bin/env python3
"""
make_figures.py  -  publication figures for the kidney representation paper.

THREE FIGURES, EACH CARRYING ONE CLAIM.

  Figure 1  what the model contains for kidney, by modality
  Figure 2  per-tissue performance with confidence intervals (forest plot)
  Figure 3  the evidence gap: precision against sample size, kidney as outlier

DESIGN DECISIONS, stated so they are not re-litigated.

Colour does one job here: separate ONE highlighted entity (kidney) from a
neutral background population (the other tissues). It is not a categorical
palette, so the usual chroma floor does not apply to the grey - the grey is
deliberately neutral context, not a series carrying identity.

  accent  #2a78d6   kidney
  neutral #85857f   all other tissues
  surface #fcfcfb, ink #0b0b0b / #52514e

Validated: CVD separation dE 15.8 (protan), normal-vision dE 17.5, both clear of
the floors; both colours >= 3:1 against the surface.

**Kidney is never identified by colour alone.** In every figure it also carries a
larger marker, a distinct shape, and a direct text label. That is what makes the
figures safe in greyscale print and for colour-blind readers, which matters more
for a journal than it does on screen.

Output: figures/figure{1,2,3}.{png,tif} at 600 dpi (TIFF for submission).

Author: Christopher Lawrence
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.use("Agg")

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGS = ROOT / "publication" / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

ACCENT = "#2a78d6"
NEUTRAL = "#85857f"
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#dededa"

mpl.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.edgecolor": INK2,
    "axes.linewidth": 0.6,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
    "axes.labelcolor": INK,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def save(fig, name: str) -> None:
    for ext in ("png", "tif"):
        fig.savefig(FIGS / f"{name}.{ext}", dpi=600, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"} if ext == "tif" else None)
    plt.close(fig)
    print(f"  wrote {name}.png and {name}.tif")


# --------------------------------------------------------------- figure 1
def figure1() -> None:
    """Kidney share of model output tracks, by modality."""
    t = pd.read_csv(RESULTS / "tracks_by_modality.csv", index_col=0)
    t = t[t.index.astype(str) != ""]
    t = t[t["TOTAL"] > 0].copy()
    # SPLICE_SITES carries no biosample at all (4 donor/acceptor tracks, one per
    # strand): it is an annotation-level output, not a tissue-specific one. A
    # "0/4" kidney share there is meaningless, so it is excluded from an
    # organ-coverage comparison rather than shown as missing coverage.
    t = t.drop(index="SPLICE_SITES", errors="ignore")
    t["pct"] = 100 * t["KIDNEY_PARENCHYMA"] / t["TOTAL"]
    t = t.sort_values("pct")

    labels = {
        "CHIP_TF": "TF ChIP-seq", "ATAC": "ATAC-seq",
        "CHIP_HISTONE": "Histone ChIP-seq", "CAGE": "CAGE",
        "RNA_SEQ": "RNA-seq", "DNASE": "DNase-seq",
        "SPLICE_SITE_USAGE": "Splice site usage",
        "SPLICE_JUNCTIONS": "Splice junctions",
        "CONTACT_MAPS": "Contact maps", "PROCAP": "PRO-cap",
        "SPLICE_SITES": "Splice sites",
    }
    names = [labels.get(i, i) for i in t.index]

    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    y = np.arange(len(t))
    ax.barh(y, t["pct"], color=NEUTRAL, height=0.62, zorder=3)
    # highlight the two that carry the claim
    for i, idx in enumerate(t.index):
        if idx in ("CHIP_TF", "ATAC"):
            ax.barh(i, t["pct"].iloc[i], color=ACCENT, height=0.62, zorder=4)

    ax.set_yticks(y, names)
    ax.set_xlabel("Kidney share of tracks in that modality (%)")
    ax.set_xlim(0, max(5.0, t["pct"].max() * 1.35))
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    for i, (pct, n, tot) in enumerate(zip(t["pct"], t["KIDNEY_PARENCHYMA"],
                                          t["TOTAL"])):
        ax.text(pct + 0.12, i, f"{int(n)}/{int(tot)}", va="center",
                fontsize=6.8, color=INK2)

    ax.set_title("Kidney representation across AlphaGenome output modalities",
                 loc="left", pad=8, fontweight="bold")

    # No leader line: it ran from the bar tip straight through the "4/1617"
    # count label. The text names the modality outright, and sits on that row,
    # so the pointer was redundant as well as untidy.
    ax.text(2.05, list(t.index).index("CHIP_TF"),
            "all 4 kidney TF tracks are CTCF;\nHEK293 alone carries 111",
            fontsize=6.8, color=ACCENT, va="center")
    save(fig, "figure1")


# --------------------------------------------------------------- figure 2
def figure2() -> None:
    """Forest plot of per-tissue Spearman with 95% CIs."""
    m = pd.read_csv(RESULTS / "per_tissue_metrics.csv").sort_values("n")
    # .capitalize() lowercases the rest of the string, so acronyms must be
    # restored AFTERWARDS or they come out as "Lcl" and "Esophagus gej".
    pretty = (m["sample_group"].str.replace("_", " ").str.capitalize()
              .str.replace("gej", "GEJ", regex=False)
              .str.replace("Lcl", "LCL", regex=False))

    h = max(4.2, 0.135 * len(m) + 1.3)
    fig, ax = plt.subplots(figsize=(5.0, h))
    y = np.arange(len(m))

    for i, r in enumerate(m.itertuples()):
        is_kid = r.sample_group == "kidney_cortex"
        c = ACCENT if is_kid else NEUTRAL
        ax.plot([r.spearman_lo, r.spearman_hi], [i, i], color=c,
                linewidth=2.2 if is_kid else 1.2, solid_capstyle="round",
                zorder=4 if is_kid else 3)
        ax.scatter([r.spearman], [i], s=60 if is_kid else 16,
                   marker="D" if is_kid else "o", color=c,
                   edgecolor=SURFACE, linewidth=0.7, zorder=5)

    ax.set_yticks(y, pretty, fontsize=6.4)
    for tick, sg in zip(ax.get_yticklabels(), m["sample_group"]):
        if sg == "kidney_cortex":
            tick.set_color(ACCENT)
            tick.set_fontweight("bold")
            tick.set_fontsize(7.4)

    mean_rho = m["spearman"].mean()
    ax.axvline(mean_rho, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
    ax.text(mean_rho + 0.004, len(m) + 0.75, f"unweighted mean {mean_rho:.2f}",
            ha="left", fontsize=6.8, color=INK2)

    ax.set_ylim(-0.9, len(m) + 1.4)
    ax.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_title("Per-tissue prediction correlations and uncertainty",
                 loc="left", pad=10, fontweight="bold")
    # y=0 is the BOTTOM in matplotlib and the sort is ascending, so the tissue
    # with fewest variants is at the bottom. Say that, do not invert the claim.
    # scored sample size printed beside each tissue, so the ordering is
    # self-explanatory without consulting Figure 3
    xr = ax.get_xlim()
    for i, nn in enumerate(m["n"]):
        ax.text(xr[1] + 0.004 * (xr[1] - xr[0]), i, f"{int(nn)}",
                va="center", ha="left", fontsize=5.8, color=INK2)
    ax.text(xr[1] + 0.004 * (xr[1] - xr[0]), len(m) + 0.2, "n",
            va="center", ha="left", fontsize=6.2, color=INK2, fontweight="bold")
    ax.set_xlim(xr)
    ax.set_xlabel("Spearman's $\\rho$, predicted vs observed eQTL effect\n"
                  "tissues ordered by scored pairs, fewest at the bottom",
                  linespacing=1.7)
    save(fig, "figure2")




# --------------------------------------------------------------- figure 3
def figure3() -> None:
    """Precision against evidence. The money figure."""
    m = pd.read_csv(RESULTS / "per_tissue_metrics.csv")
    kid = m[m["sample_group"] == "kidney_cortex"]
    oth = m[m["sample_group"] != "kidney_cortex"]

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    ax.scatter(oth["n"], oth["spearman_ci_width"], s=26, color=NEUTRAL,
               edgecolor=SURFACE, linewidth=0.6, zorder=3, label="Other tissues")
    ax.scatter(kid["n"], kid["spearman_ci_width"], s=110, color=ACCENT,
               marker="D", edgecolor=SURFACE, linewidth=1.0, zorder=5,
               label="Kidney cortex")

    ax.set_xscale("log")
    ax.set_xlabel("Successfully scored variant-gene pairs (log scale)")
    ax.set_ylabel("Width of 95% CI on Spearman's $\\rho$")
    ax.grid(color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    if len(kid):
        k = kid.iloc[0]
        ax.annotate(f"Kidney cortex\nn = {int(k['n'])}, CI width {k['spearman_ci_width']:.2f}",
                    xy=(k["n"], k["spearman_ci_width"]),
                    xytext=(k["n"] * 1.9, k["spearman_ci_width"] - 0.02),
                    fontsize=7.2, color=ACCENT, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=ACCENT, linewidth=0.8))
        med = oth["spearman_ci_width"].median()
        ax.axhline(med, color=INK2, linewidth=0.7, linestyle=(0, (4, 3)), zorder=2)
        # Put the label at the LEFT end of the line. Both ends of the right-hand
        # side are occupied: the large-n tissues sit just below the line and the
        # mid-n tissues just above it. At the left the only mark is the kidney
        # diamond, far above, so this band is genuinely empty.
        x0 = ax.get_xlim()[0]
        ax.text(x0 * 1.06, med + 0.008,
                f"median across other tissues  {med:.2f}",
                ha="left", va="bottom", fontsize=6.8, color=INK2)

    ax.legend(loc="upper right", fontsize=7.2)
    ax.set_title("Kidney cortex has the smallest benchmark and widest interval",
                 loc="left", pad=8, fontweight="bold")
    save(fig, "figure3")


# --------------------------------------------------------------- figure 4
def figure4() -> None:
    """Sensitivity to the choice of tissue output, with paired uncertainty.

    Panel A is the ranking. Panel B is the part that matters: paired bootstrap
    differences from the kidney track, with a zero reference. Showing only the
    ranking invites the inference that the higher-ranked tracks are better, which
    the paired analysis does not support.
    """
    m = pd.read_csv(RESULTS / "specificity_by_tissue.csv")
    pr = pd.read_csv(RESULTS / "specificity_paired.csv")
    kid = m[m["tissue_track"] == "Kidney_Cortex"]
    oth = m[m["tissue_track"] != "Kidney_Cortex"]

    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(4.8, 4.3), gridspec_kw={"height_ratios": [1, 2.1]})

    # ---- A: distribution
    axA.scatter(oth["spearman"], np.zeros(len(oth)), s=22, color=NEUTRAL,
                alpha=0.75, edgecolor=SURFACE, linewidth=0.5, zorder=3,
                label="Other tissue outputs")
    axA.scatter(kid["spearman"], [0], s=85, color=ACCENT, marker="D",
                edgecolor=SURFACE, linewidth=0.9, zorder=5,
                label="Kidney cortex")
    axA.set_yticks([]); axA.set_ylim(-0.5, 0.5)
    for sp in ("left", "right", "top"):
        axA.spines[sp].set_visible(False)
    axA.set_xlabel("Spearman $\\rho$ with observed kidney eQTL effects", fontsize=7.6)
    axA.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    axA.set_axisbelow(True)
    axA.legend(loc="upper left", fontsize=6.6, ncol=2,
               bbox_to_anchor=(0.0, 1.45))
    axA.set_title("A", loc="left", fontweight="bold", fontsize=9)

    # ---- B: paired differences
    top = pr.nlargest(8, "diff_vs_kidney").iloc[::-1]
    y = np.arange(len(top))
    axB.axvline(0, color=INK2, linewidth=0.8, zorder=2)
    for i, r in enumerate(top.itertuples()):
        axB.plot([r.lo, r.hi], [i, i], color=NEUTRAL, linewidth=1.6,
                 solid_capstyle="round", zorder=3)
        axB.scatter([r.diff_vs_kidney], [i], s=26, color=NEUTRAL,
                    edgecolor=SURFACE, linewidth=0.6, zorder=4)
    axB.set_yticks(y, [t.replace("_", " ") for t in top["tissue_track"]],
                   fontsize=6.8)
    axB.set_xlabel("Difference in $\\rho$ from the kidney cortex output\n"
                   "(paired bootstrap over the same 58 pairs, 95% CI, unadjusted)",
                   fontsize=7.6, linespacing=1.6)
    axB.grid(axis="x", color=GRID, linewidth=0.5, zorder=0)
    axB.set_axisbelow(True)
    for sp in ("right", "top"):
        axB.spines[sp].set_visible(False)
    axB.text(0.0, len(top) - 0.35, "  favours kidney  |  favours other output",
             fontsize=6.2, color=INK2, ha="center")
    axB.set_title("B", loc="left", fontweight="bold", fontsize=9)

    fig.suptitle("Kidney eQTL prediction across tissue outputs",
                 x=0.02, ha="left", fontweight="bold", fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save(fig, "figure4")


def main() -> None:
    print("building figures")
    figure1()
    figure2()
    figure3()
    figure4()
    print(f"all figures in {FIGS}")


if __name__ == "__main__":
    main()
