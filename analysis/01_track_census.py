#!/usr/bin/env python3
"""
01_track_census.py

Census of kidney representation in the AlphaGenome human output-track manifest.

Input : data/alphagenome_tracks_human_raw.csv
        (parsed from the embedded dataframe in google-deepmind/alphagenome
         colabs/tissue_ontology_mapping.ipynb; cross-checks against
         data/OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto)

Output: results/track_census.json + several .csv tables, all regenerable.

INCLUSION RULE (stated explicitly because the whole paper rests on it).

  Biosamples are assigned to exactly one of four classes by curated name, NOT by
  substring matching. Substring matching on "renal" is WRONG: it captures
  "adrenal gland" (27 tracks) and "perirenal adipocyte" (2). Those are a
  different organ and adipose tissue respectively. This is the single easiest
  way to inflate a kidney count and it is why every assignment below is explicit.

  KIDNEY_PARENCHYMA  the nephron-bearing organ and its resident cell types.
                     This is the denominator that matters for kidney disease
                     heritability, and it is what "kidney" should mean.
  URINARY_TRACT      ureter, renal pelvis, urothelium. Drained by the kidney,
                     not part of the nephron. Reported separately, never pooled.
  KIDNEY_CELL_LINE   HEK293 / HEK293T. Human Embryonic Kidney 293 is a
                     transformed line of debated, probably neuronal, lineage. It
                     is not kidney tissue and must never be counted as such, but
                     it is reported because its track count is the paper's
                     sharpest contrast.
  EXCLUDED           adrenal gland, perirenal adipose. Substring false positives.

Author: Christopher Lawrence
"""

from __future__ import annotations
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

# ---------------------------------------------------------------- classification

KIDNEY_PARENCHYMA = {
    # whole organ / gross subdivisions
    "kidney",
    "left kidney",
    "right kidney",
    "cortex of kidney",
    "outer medulla of kidney",
    "renal cortex interstitium",
    "left renal cortex interstitium",
    "right renal cortex interstitium",
    # developmental
    "metanephros",
    "nephron",
    "nephron progenitor cell",
    # resident cell types
    "kidney epithelial cell",
    "kidney tubule cell",
    "nephron tubule epithelial cell",
    "epithelial cell of proximal tubule",
    "renal cortical epithelial cell",
    "mesangial cell",
    "glomerular endothelial cell",
    "glomerular visceral epithelial cell",   # = podocyte
    "kidney capillary endothelial cell",
}

URINARY_TRACT = {
    "ureter",
    "renal pelvis",
    "left renal pelvis",
    "right renal pelvis",
    "urothelial cell",
    "urothelium cell line",
}

KIDNEY_CELL_LINE = {"HEK293", "HEK293T"}

EXCLUDED_FALSE_POSITIVES = {
    "adrenal gland",          # different organ; matches /renal/ by accident
    "perirenal adipocyte cell",
    "perirenal preadipocyte",
}

# Distal nephron cell types that a kidney-disease model would need.
# Checked for ABSENCE. Names follow Cell Ontology / UBERON usage.
#
# NOTE on "renal papilla": the bare term "papilla" must NOT be used. It matches
# "fungiform papilla" (tongue, 2 tracks) and "hair follicle dermal papilla cell"
# (5 tracks), neither of which is kidney. Same class of error as "adrenal gland"
# matching /renal/. Every term below is qualified for this reason.
DISTAL_NEPHRON_EXPECTED = [
    "collecting duct",
    "thick ascending limb",
    "loop of Henle",
    "distal convoluted tubule",
    "principal cell",
    "intercalated cell",
    "macula densa",
    "juxtaglomerular",
    "parietal epithelial cell",
    "connecting tubule",
    "inner medulla of kidney",
    "renal papilla",
]


def classify(name: str) -> str:
    if not isinstance(name, str):
        return "OTHER"
    if name in KIDNEY_PARENCHYMA:
        return "KIDNEY_PARENCHYMA"
    if name in URINARY_TRACT:
        return "URINARY_TRACT"
    if name in KIDNEY_CELL_LINE:
        return "KIDNEY_CELL_LINE"
    if name in EXCLUDED_FALSE_POSITIVES:
        return "EXCLUDED_FALSE_POSITIVE"
    return "OTHER"


# ---------------------------------------------------------------- load

def load() -> pd.DataFrame:
    df = pd.read_csv(DATA / "alphagenome_tracks_human_raw.csv")
    df["biosample_name"] = df["biosample_name"].fillna("")
    df["klass"] = df["biosample_name"].map(classify)
    # modality: strip the OutputType. prefix for readability
    df["modality"] = df["output_type"].fillna("").str.replace(
        "OutputType.", "", regex=False)
    return df


# ---------------------------------------------------------------- guard

def guard_no_substring_leak(df: pd.DataFrame) -> dict:
    """Every biosample whose name contains a kidney-ish substring must have been
    explicitly classified. If a new name appears in a future release that we have
    not curated, fail loudly rather than silently miscount."""
    pat = re.compile(
        r"kidney|renal|nephr|glomer|podocyte|tubul|ureter|HEK|urothel|mesangi|"
        r"juxtaglom|henle|collecting duct|intercalated|principal cell|pelvis",
        re.I)
    hits = df[df["biosample_name"].str.contains(pat)]
    unclassified = sorted(set(hits.loc[hits["klass"] == "OTHER", "biosample_name"]))
    return {
        "substring_matching_tracks": int(len(hits)),
        "unclassified_names": unclassified,
        "ok": len(unclassified) == 0,
    }


# ---------------------------------------------------------------- census

def main() -> None:
    df = load()
    total = len(df)

    guard = guard_no_substring_leak(df)
    if not guard["ok"]:
        raise SystemExit(
            "UNCLASSIFIED kidney-like biosample names found; curate them before "
            f"reporting any count: {guard['unclassified_names']}")

    kid = df[df["klass"] == "KIDNEY_PARENCHYMA"]
    uri = df[df["klass"] == "URINARY_TRACT"]
    line = df[df["klass"] == "KIDNEY_CELL_LINE"]

    # ---- headline
    headline = {
        "total_human_tracks": total,
        "total_biosamples": int(df["biosample_name"].nunique()),
        "kidney_parenchyma_tracks": int(len(kid)),
        "kidney_parenchyma_pct": round(100 * len(kid) / total, 2),
        "kidney_parenchyma_biosamples": int(kid["biosample_name"].nunique()),
        "urinary_tract_tracks": int(len(uri)),
        "kidney_cell_line_tracks": int(len(line)),
        "excluded_false_positive_tracks": int(
            (df["klass"] == "EXCLUDED_FALSE_POSITIVE").sum()),
    }

    # ---- by modality
    by_mod = (
        pd.crosstab(df["modality"], df["klass"])
          .reindex(columns=["KIDNEY_PARENCHYMA", "URINARY_TRACT",
                            "KIDNEY_CELL_LINE", "OTHER",
                            "EXCLUDED_FALSE_POSITIVE"], fill_value=0)
    )
    by_mod["TOTAL"] = by_mod.sum(axis=1)
    by_mod["kidney_pct"] = (
        100 * by_mod["KIDNEY_PARENCHYMA"] / by_mod["TOTAL"]).round(2)
    by_mod = by_mod.sort_values("TOTAL", ascending=False)
    by_mod.to_csv(RESULTS / "tracks_by_modality.csv")

    # ---- the TF contrast: which TFs are assayed in real kidney vs in HEK293
    tf = df[df["modality"] == "CHIP_TF"]
    kid_tf = sorted(set(tf.loc[tf["klass"] == "KIDNEY_PARENCHYMA",
                               "transcription_factor"].dropna()))
    line_tf = sorted(set(tf.loc[tf["klass"] == "KIDNEY_CELL_LINE",
                                "transcription_factor"].dropna()))
    tf_contrast = {
        "kidney_parenchyma_tf_tracks": int((tf["klass"] == "KIDNEY_PARENCHYMA").sum()),
        "kidney_parenchyma_distinct_tfs": kid_tf,
        "cell_line_tf_tracks": int((tf["klass"] == "KIDNEY_CELL_LINE").sum()),
        "cell_line_distinct_tf_count": len(line_tf),
        "total_tf_tracks": int(len(tf)),
    }

    # ---- per-biosample table
    per_bio = (
        df[df["klass"].isin(["KIDNEY_PARENCHYMA", "URINARY_TRACT",
                             "KIDNEY_CELL_LINE"])]
        .groupby(["klass", "biosample_name", "biosample_type",
                  "biosample_life_stage"])
        .size().rename("n_tracks").reset_index()
        .sort_values(["klass", "n_tracks"], ascending=[True, False])
    )
    per_bio.to_csv(RESULTS / "kidney_biosamples.csv", index=False)

    # ---- per-biosample x modality, kidney parenchyma only
    pivot = pd.crosstab(kid["biosample_name"], kid["modality"])
    pivot["TOTAL"] = pivot.sum(axis=1)
    pivot.sort_values("TOTAL", ascending=False).to_csv(
        RESULTS / "kidney_biosample_by_modality.csv")

    # ---- distal nephron absence
    all_names = " || ".join(sorted(set(df["biosample_name"]))).lower()
    distal = {term: (term.lower() in all_names) for term in DISTAL_NEPHRON_EXPECTED}

    # ---- life stage skew
    def stage_mix(sub: pd.DataFrame) -> dict:
        vc = sub["biosample_life_stage"].fillna("unknown").value_counts()
        return {k: int(v) for k, v in vc.items()}

    life_stage = {
        "kidney_parenchyma": stage_mix(kid),
        "all_tracks": stage_mix(df),
    }

    # ---- GTEx kidney tissues present?
    gtex = df["gtex_tissue"].dropna()
    gtex_kidney = sorted({t for t in gtex if "kidney" in str(t).lower()})

    out = {
        "headline": headline,
        "guard": guard,
        "by_modality": by_mod.to_dict(orient="index"),
        "tf_contrast": tf_contrast,
        "distal_nephron_present": distal,
        "life_stage": life_stage,
        "gtex_kidney_tissues": gtex_kidney,
        "gtex_tissues_total": int(gtex.nunique()),
    }
    (RESULTS / "track_census.json").write_text(json.dumps(out, indent=2))

    # ---------------------------------------------------------------- report
    h = headline
    print("=" * 72)
    print("ALPHAGENOME HUMAN TRACK CENSUS  -  KIDNEY")
    print("=" * 72)
    print(f"Total human output tracks      : {h['total_human_tracks']:,}")
    print(f"Distinct biosamples            : {h['total_biosamples']:,}")
    print()
    print(f"Kidney parenchyma tracks       : {h['kidney_parenchyma_tracks']:,} "
          f"({h['kidney_parenchyma_pct']}%)")
    print(f"  across biosamples            : {h['kidney_parenchyma_biosamples']}")
    print(f"Urinary tract (not nephron)    : {h['urinary_tract_tracks']}")
    print(f"HEK293/HEK293T cell line       : {h['kidney_cell_line_tracks']}")
    print(f"Excluded substring false pos.  : {h['excluded_false_positive_tracks']}"
          "   <- 'adrenal gland', 'perirenal adipocyte'")
    print()
    print("-" * 72)
    print("BY MODALITY (kidney parenchyma / total)")
    print("-" * 72)
    for mod, row in by_mod.iterrows():
        if not mod:
            continue
        print(f"  {mod:<26} {int(row['KIDNEY_PARENCHYMA']):>5} / "
              f"{int(row['TOTAL']):>5}   {row['kidney_pct']:>5.2f}%")
    print()
    print("-" * 72)
    print("TRANSCRIPTION FACTOR CONTRAST")
    print("-" * 72)
    t = tf_contrast
    print(f"  TF ChIP tracks, kidney parenchyma : {t['kidney_parenchyma_tf_tracks']}")
    print(f"  distinct TFs assayed in kidney    : {t['kidney_parenchyma_distinct_tfs']}")
    print(f"  TF ChIP tracks, HEK293(T)         : {t['cell_line_tf_tracks']}")
    print(f"  distinct TFs in HEK293(T)         : {t['cell_line_distinct_tf_count']}")
    if t["kidney_parenchyma_tf_tracks"]:
        print(f"  ratio cell line : kidney          = "
              f"{t['cell_line_tf_tracks'] / t['kidney_parenchyma_tf_tracks']:.1f}x")
    print()
    print("-" * 72)
    print("DISTAL NEPHRON CELL TYPES PRESENT IN ANY BIOSAMPLE NAME?")
    print("-" * 72)
    for term, present in distal.items():
        print(f"  {term:<28} {'PRESENT' if present else 'ABSENT'}")
    print()
    print("-" * 72)
    print(f"GTEx tissues represented          : {out['gtex_tissues_total']}")
    print(f"GTEx kidney tissues               : {gtex_kidney or 'NONE'}")
    print()
    print(f"Wrote {RESULTS}/track_census.json and 3 CSV tables")


if __name__ == "__main__":
    main()
