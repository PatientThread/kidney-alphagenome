# Kidney representation and uncertainty in AlphaGenome regulatory variant prediction

Analysis code and data for the manuscript of the same name.

**Author:** Christopher Lawrence, Consultant Nephrologist, London, United Kingdom.
ORCID [0000-0002-8159-0879](https://orcid.org/0000-0002-8159-0879).

---

## What this repository contains

Everything needed to reproduce the paper: the frozen benchmark set, the curated
tissue mapping, the per-tissue prediction tables, the output-track manifests and
the thirteen numbered analysis scripts that produce every number and figure.

## The three findings

1. **Kidney is sparsely represented.** Kidney parenchyma accounts for 103 of
   5563 human output tracks (1.85%; 1.97% counting splice junctions per strand),
   across 20 of 715 biosample labels. All four kidney transcription-factor tracks
   assay CTCF, and cell-type-resolved accessibility is limited to seven DNase
   tracks, only one of them from an adult donor.
2. **Kidney is the least well evidenced tissue in this benchmark.** Kidney cortex
   yields 59 usable fine-mapped eQTLs against a median of 616, and is the
   smallest tissue at all 36 threshold settings tested.
3. **Kidney performance is insensitive to the choice of tissue output.** Scoring
   the same kidney variants with every tissue track in turn, the kidney track
   ranks 6th of 54, and no higher-ranked track is distinguishable from it by
   paired bootstrap over the same 58 pairs. Kidney was nonetheless better than 10
   of 53 comparators, and the tracks' prediction vectors are substantially but
   not wholly shared (median pairwise rank correlation 0.827). Whether the
   insensitivity reflects genuinely shared regulatory effects or limited tissue
   resolution cannot be separated by this design.

## What this repository does *not* claim

The manuscript explicitly withdraws three claims that earlier drafts made and
that the data do not support: that performance in kidney is *adequate*, that it
is *equivalent* to other tissues, and that missing experiments are the *decisive*
constraint. See the Discussion.

---

## Reproducing the analysis

```bash
python3 -m venv venv && ./venv/bin/pip install -r requirements.txt
```

Scripts are numbered in dependency order.

| Script | Needs network | Needs API key | Produces |
|---|---|---|---|
| `01_track_census.py` | no | no | track census, curated classification, guards |
| `02_textproto_crosscheck.py` | no | no | manifest reconciliation, padding analysis |
| `03_gtex_tissue_weighting.py` | yes (GTEx) | no | v8 and v10 tissue weighting |
| `04_encode_kidney_inventory.py` | yes (ENCODE) | no | kidney experiment inventory |
| `05_fetch_eqtl_catalogue.py` | yes (EBI) | no | downloads 49 tissues, ~150 MB |
| `06_build_benchmark_set.py` | no | no | **the frozen benchmark set** |
| `07_power_and_threshold_sensitivity.py` | no | no | 36-setting threshold sweep, power |
| `08_live_api_crosscheck.py` | yes | **yes** | live manifest verification |
| `09_score_variants.py` | yes | **yes** | per-tissue predictions (~10 h) |
| `10_evaluate.py` | no | no | per-tissue metrics, figures' input |
| `11_review_robustness.py` | no | no | independence, downsampling, max-selection |
| `12_specificity_control.py` | yes | **yes** | **the specificity control** |
| `13_specificity_uncertainty.py` | no | no | paired bootstrap, pairwise predictions |
| `publication/make_figures.py` | no | no | Figures 1-4 |

Scripts 01, 02, 06, 07, 10, 11 and 13 run offline from what is committed here,
and reproduce every headline number including the specificity result summary.

### The API key

Scripts 08, 09 and 12 require an AlphaGenome API key, obtained under the
developers' non-commercial research terms. **The key is never stored in this
repository.** It is read from `~/.alphagenome/api_key` (mode 600) or the
`ALPHAGENOME_API_KEY` environment variable, by `analysis/ag_auth.py`, which
refuses to load a key whose file permissions are looser than 600. No key
material, masked or otherwise, is written to any committed file.

### Not committed

`data/eqtl_catalogue/` (~150 MB) is regenerable with script 05 and is excluded.
Everything else needed is here.

---

## Data sources and licences

| Source | Use | Licence |
|---|---|---|
| EMBL-EBI eQTL Catalogue | benchmark ground truth, 49 GTEx tissues | CC BY 4.0 |
| GTEx | tissue sample sizes and eGene counts | no commercial restriction |
| ENCODE | kidney experiment inventory | free to download, analyse and publish |
| AlphaGenome | model predictions | non-commercial research terms; weights not redistributed |

This is a secondary analysis of public summary-level data. No individual-level
clinical data were accessed and no participants were recruited.

## Licence

Analysis code in this repository is released under the MIT licence. Data files
redistributed here remain under the licences of their original sources, listed
above. AlphaGenome model weights are not redistributed.

## Citation

Lawrence C. Kidney representation and uncertainty in AlphaGenome regulatory
variant prediction. *Manuscript in preparation.*
