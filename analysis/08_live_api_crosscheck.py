#!/usr/bin/env python3
"""
08_live_api_crosscheck.py

THIRD independent verification of the track census, this time against the LIVE
AlphaGenome API rather than a file.

The census so far rests on two static artefacts that both ship in DeepMind's
repositories: the notebook dataframe and the protobuf manifest. They reconcile
exactly once padding is removed. But both could in principle be stale relative
to what the served model actually exposes. This script asks the deployed service
directly.

If all three agree, the kidney numbers are as solid as they can be made without
DeepMind's cooperation, and the paper can say so.

Requires an API key. See docs/GETTING_MODEL_ACCESS.md. The key is read from
~/.alphagenome/api_key via ag_auth and never appears in this repository.

Run with the venv that has the alphagenome package:
    /tmp/agvenv/bin/python analysis/08_live_api_crosscheck.py

Output: data/alphagenome_live_api_tracks.csv
        results/live_api_crosscheck.json

Author: Christopher Lawrence
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

sys.path.insert(0, str(HERE))
from ag_auth import get_key, masked  # noqa: E402

# reuse the curated kidney classification, single source of truth
spec = importlib.util.spec_from_file_location("census", HERE / "01_track_census.py")
census = importlib.util.module_from_spec(spec)
spec.loader.exec_module(census)

MODALITIES = [
    "atac", "cage", "chip_histone", "chip_tf", "contact_maps", "dnase",
    "procap", "rna_seq", "splice_junctions", "splice_site_usage", "splice_sites",
]

# map API attribute names onto the OutputType.* labels used in the static CSV
API_TO_CSV = {
    "atac": "ATAC", "cage": "CAGE", "chip_histone": "CHIP_HISTONE",
    "chip_tf": "CHIP_TF", "contact_maps": "CONTACT_MAPS", "dnase": "DNASE",
    "procap": "PROCAP", "rna_seq": "RNA_SEQ",
    "splice_junctions": "SPLICE_JUNCTIONS",
    "splice_site_usage": "SPLICE_SITE_USAGE", "splice_sites": "SPLICE_SITES",
}


def main() -> None:
    from alphagenome.models import dna_client

    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    om = client.output_metadata(organism=dna_client.Organism.HOMO_SAPIENS)
    print("live output metadata retrieved\n")

    frames = []
    live_counts = {}
    for m in MODALITIES:
        df = getattr(om, m, None)
        if df is None or not isinstance(df, pd.DataFrame):
            live_counts[API_TO_CSV[m]] = 0
            continue
        d = df.copy()
        d["modality"] = API_TO_CSV[m]
        frames.append(d)
        live_counts[API_TO_CSV[m]] = len(d)

    live = pd.concat(frames, ignore_index=True)
    live.to_csv(DATA / "alphagenome_live_api_tracks.csv", index=False)

    # classify kidney the same way as everywhere else
    if "biosample_name" not in live.columns:
        raise SystemExit("live metadata has no biosample_name column; "
                         f"columns are {list(live.columns)}")
    live["biosample_name"] = live["biosample_name"].fillna("")
    live["klass"] = live["biosample_name"].map(census.classify)

    # static CSV for comparison
    static = pd.read_csv(DATA / "alphagenome_tracks_human_raw.csv")
    static["biosample_name"] = static["biosample_name"].fillna("")
    static["klass"] = static["biosample_name"].map(census.classify)
    static["modality"] = static["output_type"].fillna("").str.replace(
        "OutputType.", "", regex=False)
    static_counts = static["modality"].value_counts().to_dict()

    mods = sorted(set(live_counts) | set(static_counts))
    agree = {}
    for m in mods:
        l, s = int(live_counts.get(m, 0)), int(static_counts.get(m, 0))
        agree[m] = {"live": l, "static": s, "match": l == s}
    all_match = all(v["match"] for v in agree.values())

    kid_live = int((live["klass"] == "KIDNEY_PARENCHYMA").sum())
    kid_static = int((static["klass"] == "KIDNEY_PARENCHYMA").sum())

    # kidney TF check, live
    tf_live = live[live["modality"] == "CHIP_TF"]
    kid_tf = tf_live[tf_live["klass"] == "KIDNEY_PARENCHYMA"]
    tf_col = next((c for c in ("transcription_factor", "target", "assay")
                   if c in tf_live.columns), None)
    kid_tf_targets = (sorted(set(kid_tf[tf_col].dropna()))
                      if tf_col else ["<no tf column>"])

    out = {
        "live_total_tracks": int(len(live)),
        "static_total_tracks": int(len(static)),
        "totals_match": len(live) == len(static),
        "per_modality": agree,
        "all_modalities_match": all_match,
        "kidney_parenchyma_live": kid_live,
        "kidney_parenchyma_static": kid_static,
        "kidney_match": kid_live == kid_static,
        "kidney_pct_live": round(100 * kid_live / len(live), 2),
        "kidney_tf_tracks_live": int(len(kid_tf)),
        "kidney_tf_targets_live": kid_tf_targets,
    }
    (RESULTS / "live_api_crosscheck.json").write_text(json.dumps(out, indent=2))

    print("=" * 70)
    print("LIVE API vs STATIC MANIFEST")
    print("=" * 70)
    print(f"  live total tracks   : {len(live):,}")
    print(f"  static total tracks : {len(static):,}")
    print(f"  totals match        : {'YES' if out['totals_match'] else 'NO'}")
    print()
    print(f"  {'modality':<22}{'live':>8}{'static':>8}   match")
    for m in mods:
        v = agree[m]
        print(f"  {m:<22}{v['live']:>8}{v['static']:>8}   "
              f"{'ok' if v['match'] else 'MISMATCH'}")
    print()
    print(f"  kidney parenchyma, live   : {kid_live}  ({out['kidney_pct_live']}%)")
    print(f"  kidney parenchyma, static : {kid_static}")
    print(f"  kidney match              : {'YES' if out['kidney_match'] else 'NO'}")
    print()
    print(f"  kidney TF tracks, live    : {out['kidney_tf_tracks_live']}")
    print(f"  kidney TF targets, live   : {kid_tf_targets}")
    print()
    if all_match and out["kidney_match"]:
        print("  VERDICT: the deployed model, the shipped protobuf manifest and")
        print("  the notebook dataframe all agree. Three independent sources.")
        print("  The census is as solid as it can be made externally.")
    else:
        print("  VERDICT: sources DISAGREE. The live API is authoritative for")
        print("  what is served today. Investigate before publishing any count.")
    print()
    print(f"  Wrote {RESULTS}/live_api_crosscheck.json")


if __name__ == "__main__":
    main()
