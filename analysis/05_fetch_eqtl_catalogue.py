#!/usr/bin/env python3
"""
05_fetch_eqtl_catalogue.py

LEG 2, part 2: pull the ground truth for an UNWEIGHTED per-tissue benchmark.

WHY THE eQTL CATALOGUE RATHER THAN GTEx DIRECTLY. Three reasons, and the third
is the one that decides it.

  1. Uniform processing. All 49 GTEx v8 tissues are re-quantified through one
     pipeline, so a per-tissue comparison is not confounded by pipeline
     differences between tissues.
  2. Licence. CC BY 4.0, commercial use permitted. GTEx's own terms are also
     unrestricted, but the Susztak atlas (the other large kidney resource) is
     licence-conflicted and cannot safely anchor a published figure.
  3. It is the SAME VINTAGE the model was benchmarked on. The catalogue's GTEx
     is v8, kidney cortex n=73, which reproduces the GTEx API v8 figures exactly
     and independently. Benchmarking against v10 would compare the model to data
     that did not exist when it was evaluated.

WHAT IS PULLED, AND WHY NOT EVERYTHING. The full nominal sumstats are ~2.4 GB
per dataset; 49 of those is ~120 GB and unnecessary. Two small files per tissue
carry what a benchmark needs:

  permuted.tsv.gz       one lead variant per gene, with the permutation p-value
                        and beta. This is the eGene set and the effect
                        directions. ~1 MB.
  credible_sets.tsv.gz  SuSiE fine-mapped credible sets. These are the variants
                        with real posterior support for being causal, which is
                        the honest positive set for a variant-effect benchmark.
                        ~300 KB.

RATE LIMITING IS NOT OPTIONAL. The catalogue's own documentation warns that
frequent requests are read as a denial-of-service attempt and the source address
is blacklisted. This script sleeps between requests, retries with backoff, and
resumes rather than re-downloading. Do not remove the sleeps.

Output: data/eqtl_catalogue/{permuted,susie}/<dataset_id>.tsv.gz
        results/eqtl_catalogue_fetch_log.csv

Author: Christopher Lawrence
"""

from __future__ import annotations
import argparse
import gzip
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
OUT = DATA / "eqtl_catalogue"
(OUT / "permuted").mkdir(parents=True, exist_ok=True)
(OUT / "susie").mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(exist_ok=True)

BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL"
SLEEP = 1.5          # seconds between requests. Do not lower.
RETRIES = 4
TIMEOUT = 300
UA = "kidney-alphagenome-research/1.0 (academic; contact via ORCID 0000-0002-8159-0879)"


def targets() -> pd.DataFrame:
    md = pd.read_csv(DATA / "eqtl_catalogue_dataset_metadata_r7.tsv", sep="\t")
    g = md[(md["study_label"] == "GTEx") & (md["quant_method"] == "ge")].copy()
    return g.sort_values("sample_size").reset_index(drop=True)


def get(url: str, dest: Path) -> tuple[bool, str]:
    if dest.exists() and dest.stat().st_size > 0:
        return True, "cached"
    for attempt in range(1, RETRIES + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as fh:
                body = fh.read()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(body)
            tmp.rename(dest)
            return True, f"ok {len(body):,}B"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False, "404"
            wait = SLEEP * (2 ** attempt)
            print(f"      HTTP {e.code}, retry {attempt}/{RETRIES} in {wait:.0f}s",
                  file=sys.stderr)
            time.sleep(wait)
        except Exception as e:                       # noqa: BLE001
            wait = SLEEP * (2 ** attempt)
            print(f"      {type(e).__name__}, retry {attempt}/{RETRIES} in {wait:.0f}s",
                  file=sys.stderr)
            time.sleep(wait)
    return False, "failed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0,
                    help="fetch only the N smallest tissues (for a dry run)")
    ap.add_argument("--kind", choices=["both", "permuted", "susie"],
                    default="both")
    args = ap.parse_args()

    g = targets()
    if args.limit:
        g = g.head(args.limit)

    print(f"GTEx gene-expression datasets to fetch: {len(g)}")
    print(f"throttle {SLEEP}s between requests\n")

    log = []
    for i, row in g.iterrows():
        ds, st = row["dataset_id"], row["study_id"]
        grp, n = row["sample_group"], int(row["sample_size"])
        tag = f"[{i+1}/{len(g)}] {grp} (n={n})"
        print(f"{tag}")

        if args.kind in ("both", "permuted"):
            url = f"{BASE}/sumstats/{st}/{ds}/{ds}.permuted.tsv.gz"
            dest = OUT / "permuted" / f"{ds}.permuted.tsv.gz"
            ok, msg = get(url, dest)
            print(f"    permuted : {msg}")
            log.append(dict(dataset_id=ds, sample_group=grp, sample_size=n,
                            kind="permuted", ok=ok, msg=msg))
            if msg != "cached":
                time.sleep(SLEEP)

        if args.kind in ("both", "susie"):
            url = f"{BASE}/susie/{st}/{ds}/{ds}.credible_sets.tsv.gz"
            dest = OUT / "susie" / f"{ds}.credible_sets.tsv.gz"
            ok, msg = get(url, dest)
            print(f"    susie    : {msg}")
            log.append(dict(dataset_id=ds, sample_group=grp, sample_size=n,
                            kind="susie", ok=ok, msg=msg))
            if msg != "cached":
                time.sleep(SLEEP)

    df = pd.DataFrame(log)
    df.to_csv(RESULTS / "eqtl_catalogue_fetch_log.csv", index=False)
    bad = df[~df["ok"]]
    total_mb = sum(p.stat().st_size for p in OUT.rglob("*.tsv.gz")) / 1e6
    print(f"\nfetched {int(df['ok'].sum())}/{len(df)} files, "
          f"{total_mb:.0f} MB on disk")
    if len(bad):
        print("FAILURES:")
        print(bad.to_string(index=False))
    else:
        print("no failures")


if __name__ == "__main__":
    main()
