#!/usr/bin/env python3
"""
04_encode_kidney_inventory.py

LEG 3, part 1: inventory the kidney data that actually exists in ENCODE, so the
paper can say what the model COULD have been trained on versus what it WAS.

WHY ENCODE. Of every resource surveyed, ENCODE has the cleanest terms. Its Data
Use Policy states that external users "may freely download, analyze and publish
results based on any ENCODE data without restrictions", with no grace period.
Nothing else in the kidney space is that unencumbered.

THE CURATION PROBLEM, AGAIN. ENCODE's `organ_slims=kidney` facet returns 592
human experiments, but that set includes HEK293 (a transformed line of debated,
probably neuronal, lineage), RCC and other renal carcinoma lines, and immortal
lines derived from kidney. Counting those as "kidney tissue" would repeat
exactly the error corrected in script 01. Biosamples are therefore classified
into tissue / primary cell / cell line, and the headline count is tissue and
primary cell only.

Output answers three questions:
  1. How many human kidney TISSUE experiments exist in ENCODE, by assay?
  2. Which transcription factors have been ChIP'd in real kidney, if any?
     (Script 01 found the model exposes only CTCF. Is that a limit of the model
     or a limit of the underlying data? This decides whether the paper's framing
     is "DeepMind under-sampled kidney" or "the field has not generated the
     data". Those are very different claims and only one is fair.)
  3. What is downloadable now for leg 3 validation.

Author: Christopher Lawrence
"""

from __future__ import annotations
import json
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

CACHE = DATA / "encode_kidney_experiments.json"

FIELDS = [
    "accession", "assay_title", "assay_term_name",
    "biosample_ontology.term_name", "biosample_ontology.classification",
    "biosample_ontology.organ_slims", "biosample_summary",
    "target.label", "status", "date_released",
    "replicates.library.biosample.life_stage",
]

URL = (
    "https://www.encodeproject.org/search/?type=Experiment"
    "&biosample_ontology.organ_slims=kidney"
    "&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
    "&status=released&limit=all&format=json"
    + "".join(f"&field={f}" for f in FIELDS)
)

# Cell lines and cancer lines that must NOT be counted as kidney tissue.
# HEK293 is the important one: it is the single biggest distorting presence in
# every "kidney" track count, including AlphaGenome's.
NON_TISSUE_LINES = {
    "HEK293", "HEK293T", "HEK293-T-REx", "293T",
    "RCC", "RCC 7860", "786-O", "Caki-2", "ACHN", "A498", "G-401",
    "HK-2", "SK-NEP-1", "G401", "RPTEC TERT1", "RPTEC/TERT1",
}


def fetch() -> list[dict]:
    if CACHE.exists():
        return json.loads(CACHE.read_text())
    req = urllib.request.Request(URL, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as fh:
        payload = json.load(fh)
    graph = payload["@graph"]
    CACHE.write_text(json.dumps(graph, indent=1))
    return graph


def main() -> None:
    graph = fetch()
    recs = []
    for e in graph:
        bo = e.get("biosample_ontology") or {}
        tgt = (e.get("target") or {}).get("label")
        recs.append({
            "accession": e.get("accession"),
            "assay": e.get("assay_title"),
            "biosample": bo.get("term_name"),
            "classification": bo.get("classification"),
            "target": tgt,
            "summary": e.get("biosample_summary", ""),
            "date_released": e.get("date_released"),
        })
    df = pd.DataFrame(recs)
    df["is_cell_line"] = (
        df["biosample"].isin(NON_TISSUE_LINES)
        | df["classification"].fillna("").str.contains("cell line", case=False)
    )

    tissue = df[~df["is_cell_line"]]
    lines = df[df["is_cell_line"]]

    df.to_csv(RESULTS / "encode_kidney_experiments.csv", index=False)

    by_assay = (pd.crosstab(df["assay"], df["is_cell_line"])
                  .rename(columns={False: "tissue_or_primary", True: "cell_line"}))
    by_assay["total"] = by_assay.sum(axis=1)
    by_assay = by_assay.sort_values("total", ascending=False)
    by_assay.to_csv(RESULTS / "encode_kidney_by_assay.csv")

    # THE DECISIVE QUESTION: TF ChIP in real kidney tissue
    tf_all = df[df["assay"].fillna("").str.contains("TF ChIP", case=False)]
    tf_tissue = tf_all[~tf_all["is_cell_line"]]
    tf_lines = tf_all[tf_all["is_cell_line"]]
    tf_tissue_targets = sorted(set(tf_tissue["target"].dropna()))
    tf_line_targets = sorted(set(tf_lines["target"].dropna()))

    out = {
        "encode_human_kidney_organ_experiments": int(len(df)),
        "tissue_or_primary_cell": int(len(tissue)),
        "cell_line": int(len(lines)),
        "distinct_biosamples_tissue": int(tissue["biosample"].nunique()),
        "distinct_biosamples_cell_line": int(lines["biosample"].nunique()),
        "top_cell_line_biosamples": lines["biosample"].value_counts().head(8).to_dict(),
        "top_tissue_biosamples": tissue["biosample"].value_counts().head(12).to_dict(),
        "assays_tissue_only": tissue["assay"].value_counts().to_dict(),
        "tf_chip": {
            "total_experiments": int(len(tf_all)),
            "in_tissue_or_primary": int(len(tf_tissue)),
            "in_cell_line": int(len(tf_lines)),
            "distinct_targets_in_tissue": tf_tissue_targets,
            "n_distinct_targets_in_tissue": len(tf_tissue_targets),
            "n_distinct_targets_in_cell_line": len(tf_line_targets),
        },
    }
    (RESULTS / "encode_kidney_inventory.json").write_text(json.dumps(out, indent=2))

    print("=" * 74)
    print("ENCODE HUMAN KIDNEY INVENTORY")
    print("=" * 74)
    print(f"  experiments under organ_slims=kidney : {out['encode_human_kidney_organ_experiments']}")
    print(f"    tissue or primary cell             : {out['tissue_or_primary_cell']}"
          f"   ({out['distinct_biosamples_tissue']} biosamples)")
    print(f"    CELL LINE (excluded from headline) : {out['cell_line']}"
          f"   ({out['distinct_biosamples_cell_line']} biosamples)")
    print()
    print("  biggest cell-line contributors (why naive counts mislead):")
    for k, v in out["top_cell_line_biosamples"].items():
        print(f"    {str(k):<34} {v:>4}")
    print()
    print("  tissue / primary-cell biosamples:")
    for k, v in out["top_tissue_biosamples"].items():
        print(f"    {str(k):<34} {v:>4}")
    print()
    print("  assays available in kidney TISSUE / primary cell:")
    for k, v in sorted(out["assays_tissue_only"].items(), key=lambda x: -x[1]):
        print(f"    {str(k):<34} {v:>4}")
    print()
    print("-" * 74)
    print("THE DECISIVE QUESTION: IS THE CTCF-ONLY FINDING THE MODEL'S FAULT?")
    print("-" * 74)
    t = out["tf_chip"]
    print(f"  TF ChIP experiments, all kidney-slim : {t['total_experiments']}")
    print(f"    in tissue / primary cell           : {t['in_tissue_or_primary']}")
    print(f"    in cell lines                      : {t['in_cell_line']}"
          f"   ({t['n_distinct_targets_in_cell_line']} distinct targets)")
    print(f"  distinct TF targets in kidney TISSUE : {t['n_distinct_targets_in_tissue']}")
    if t["distinct_targets_in_tissue"]:
        print(f"    {t['distinct_targets_in_tissue']}")
    print()
    if t["n_distinct_targets_in_tissue"] <= 2:
        print("  READING: the underlying DATA barely exists. The model's CTCF-only")
        print("  kidney TF coverage reflects the field, not a DeepMind choice.")
        print("  Frame the paper as a FIELD-LEVEL data gap with consequences for")
        print("  model users, NOT as a criticism of the model developers.")
    else:
        print("  READING: kidney TF data exists in ENCODE beyond CTCF but is not")
        print("  represented in the model. That is a selection question worth")
        print("  putting to the developers.")
    print()
    print(f"Wrote {RESULTS}/encode_kidney_inventory.json + 2 CSVs")


if __name__ == "__main__":
    main()
