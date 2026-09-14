#!/usr/bin/env python3
"""
09_score_variants.py

LEG 2, part 5: score the frozen benchmark variants with AlphaGenome.

WHAT IS BEING MEASURED. For each fine-mapped eQTL, the observed quantity is the
effect of the alternate allele on expression of a specific gene in a specific
tissue (the eQTL beta). The matching model quantity is the predicted log fold
change in RNA-seq signal over that same gene, in the track corresponding to that
same GTEx tissue. That is what GeneMaskLFCScorer produces, so it is the only
honest comparator. Correlating against a generic "variant effect" score would be
comparing different things.

METHOD, per variant:
  1. build a 1 Mb interval centred on the variant, the model's largest context
  2. score with the recommended RNA_SEQ scorer
  3. from the returned AnnData, take the column for the matching GTEx tissue
     track and the row for the target gene
  4. record predicted LFC alongside observed beta

The per-tissue metric is Spearman's rho between predicted LFC and observed beta,
computed on that tissue's own variants only, then averaged across tissues BOTH
ways: eQTL-count weighted (as published) and unweighted (as it should be).

TISSUE MAPPING IS CURATED, NOT FUZZY. 33 of 49 eQTL Catalogue sample groups match
an AlphaGenome gtex_tissue label after normalisation; the remaining 16 are mapped
by hand below. Fuzzy matching would silently pair, for example, "brain_spinal_cord"
with the wrong brain region. Same discipline as the kidney biosample curation.

ORDER OF WORK. Kidney is scored FIRST and alone by default. It is the smallest
tissue at 59 variants, so it is cheap, and it answers the question the paper is
actually about before any budget is spent on the other 48.

Usage:
    /tmp/agvenv/bin/python analysis/09_score_variants.py --tissues kidney_cortex
    /tmp/agvenv/bin/python analysis/09_score_variants.py --all

Output: results/scores/<sample_group>.tsv.gz
        results/scoring_log.csv

Author: Christopher Lawrence
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = ROOT / "results"
BENCH = RESULTS / "benchmark_variants"
OUT = RESULTS / "scores"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(HERE))
from ag_auth import get_key, masked  # noqa: E402

# eQTL Catalogue sample_group -> AlphaGenome gtex_tissue label.
# Only the 16 that do not match by case/punctuation normalisation are listed.
MANUAL_TISSUE_MAP = {
    "brain_spinal_cord": "Brain_Spinal_cord_cervical_c-1",
    "brain_anterior_cingulate_cortex": "Brain_Anterior_cingulate_cortex_BA24",
    "brain_putamen": "Brain_Putamen_basal_ganglia",
    "brain_frontal_cortex": "Brain_Frontal_Cortex_BA9",
    "brain_nucleus_accumbens": "Brain_Nucleus_accumbens_basal_ganglia",
    "brain_caudate": "Brain_Caudate_basal_ganglia",
    "small_intestine": "Small_Intestine_Terminal_Ileum",
    "breast": "Breast_Mammary_Tissue",
    "esophagus_gej": "Esophagus_Gastroesophageal_Junction",
    "adipose_visceral": "Adipose_Visceral_Omentum",
    "skin_not_sun_exposed": "Skin_Not_Sun_Exposed_Suprapubic",
    "skin_sun_exposed": "Skin_Sun_Exposed_Lower_leg",
    "muscle": "Muscle_Skeletal",
    "blood": "Whole_Blood",
    "LCL": "Cells_EBV-transformed_lymphocytes",
    "fibroblast": "Cells_Cultured_fibroblasts",
}

SLEEP = 0.15          # polite pacing between API calls
MAX_RETRIES = 3


def normalise(x: str) -> str:
    import re
    return re.sub(r"[^a-z0-9]", "", str(x).lower())


def build_tissue_map(live_tracks: pd.DataFrame) -> dict[str, str]:
    g = live_tracks[(live_tracks["modality"] == "RNA_SEQ")
                    & live_tracks["gtex_tissue"].notna()]
    by_norm = {normalise(t): t for t in g["gtex_tissue"].unique()}
    mapping = {}
    for f in sorted(BENCH.glob("*.tsv.gz")):
        grp = f.name.replace(".tsv.gz", "")
        if grp in MANUAL_TISSUE_MAP:
            mapping[grp] = MANUAL_TISSUE_MAP[grp]
        elif normalise(grp) in by_norm:
            mapping[grp] = by_norm[normalise(grp)]
        else:
            mapping[grp] = None
    return mapping


# The model accepts ONLY these exact context lengths. 1_000_000 is rejected;
# the "1 MB" option is 2**20 = 1_048_576. Getting this wrong fails every call
# with a clear error, which is the good kind of failure.
SEQ_1MB = 1_048_576


def score_one(client, scorer, chrom, pos, ref, alt, width=SEQ_1MB):
    from alphagenome.data import genome
    variant = genome.Variant(chromosome=chrom, position=int(pos),
                             reference_bases=ref, alternate_bases=alt)
    half = width // 2
    start = max(0, int(pos) - half)
    # keep the interval exactly `width` long even when clipped at the contig start
    interval = genome.Interval(chromosome=chrom, start=start, end=start + width)
    return client.score_variant(interval=interval, variant=variant,
                                variant_scorers=[scorer])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissues", nargs="*", default=["kidney_cortex"])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--limit", type=int, default=0,
                    help="cap variants per tissue (smoke test)")
    ap.add_argument("--force", action="store_true",
                    help="re-score tissues that already have an output file")
    args = ap.parse_args()

    from alphagenome.models import dna_client, variant_scorers

    live = pd.read_csv(ROOT / "data" / "alphagenome_live_api_tracks.csv")
    tmap = build_tissue_map(live)

    groups = ([g for g in sorted(tmap)] if args.all else args.tissues)
    unmapped = [g for g in groups if not tmap.get(g)]
    if unmapped:
        raise SystemExit(f"unmapped tissues, curate them first: {unmapped}")

    print("authenticating (key not logged)")
    client = dna_client.create(get_key())
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
    print(f"scorer: {type(scorer).__name__}\n")

    log = []
    for grp in groups:
        gtex = tmap[grp]
        done_path = OUT / f"{grp}.tsv.gz"
        # A full 49-tissue run is ~10 hours. It WILL be interrupted, so resume
        # is the default and re-scoring must be asked for explicitly.
        if done_path.exists() and not args.force and not args.limit:
            n = len(pd.read_csv(done_path, sep="\t"))
            print(f"=== {grp}: already scored ({n} variants), skipping")
            continue
        bench = pd.read_csv(BENCH / f"{grp}.tsv.gz", sep="\t")
        if args.limit:
            bench = bench.head(args.limit)
        print(f"=== {grp}  ->  {gtex}   ({len(bench)} variants)")

        rows, failures = [], 0
        t0 = time.time()
        for i, r in bench.iterrows():
            chrom, pos, ref, alt = str(r["variant"]).split("_")
            gene = str(r["molecular_trait_id"])
            got = None
            for attempt in range(1, MAX_RETRIES + 1):
                try:
                    adatas = score_one(client, scorer, chrom, pos, ref, alt)
                    got = adatas
                    break
                except Exception as exc:                      # noqa: BLE001
                    if attempt == MAX_RETRIES:
                        failures += 1
                        print(f"    FAIL {r['variant']} {gene}: "
                              f"{type(exc).__name__}: {str(exc)[:120]}")
                    else:
                        time.sleep(SLEEP * 4 * attempt)
            if got is None:
                continue

            ad = got[0]
            # rows = genes, columns = tracks
            try:
                obs = ad.obs
                var = ad.var
                gene_col = next((c for c in ("gene_id", "gene_name", "gene")
                                 if c in obs.columns), None)
                gmask = (obs[gene_col].astype(str).str.split(".").str[0]
                         == gene.split(".")[0]) if gene_col else None
                tcol = "gtex_tissue" if "gtex_tissue" in var.columns else None
                tmask = (var[tcol] == gtex) if tcol else None
                if gmask is None or tmask is None or not gmask.any() or not tmask.any():
                    failures += 1
                    continue
                import numpy as np
                sub = np.asarray(ad.X)[np.where(gmask)[0][:, None],
                                       np.where(tmask)[0]]
                pred = float(np.nanmean(sub))
            except Exception:                                  # noqa: BLE001
                failures += 1
                continue

            rows.append({
                "sample_group": grp, "gtex_tissue": gtex,
                "variant": r["variant"], "rsid": r.get("rsid"),
                "gene": gene, "pip": r.get("pip"),
                "observed_beta": r.get("beta"), "observed_z": r.get("z"),
                "predicted_lfc": pred,
            })
            time.sleep(SLEEP)
            if (i + 1) % 25 == 0:
                print(f"    {i+1}/{len(bench)}  ({time.time()-t0:.0f}s)")

        df = pd.DataFrame(rows)
        if len(df):
            df.to_csv(OUT / f"{grp}.tsv.gz", sep="\t", index=False,
                      compression="gzip")
        dt = time.time() - t0
        print(f"    scored {len(df)}/{len(bench)}, {failures} failures, {dt:.0f}s")
        log.append({"sample_group": grp, "gtex_tissue": gtex,
                    "attempted": len(bench), "scored": len(df),
                    "failures": failures, "seconds": round(dt, 1)})

    pd.DataFrame(log).to_csv(RESULTS / "scoring_log.csv", index=False)
    print(f"\nWrote {OUT}/ and {RESULTS}/scoring_log.csv")


if __name__ == "__main__":
    main()
