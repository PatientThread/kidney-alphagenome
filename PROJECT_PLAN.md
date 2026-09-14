# Kidney representation in sequence-to-function genome models

**Paper 1 of two.** Started 13 September 2026. Christopher Lawrence.
Paper 2 is the APOL1 application and it depends on this one. See §7.

---

## 1. The claim

> Kidney is systematically under-represented in the training data of
> sequence-to-function genome models, and the under-representation is
> concentrated precisely in the cell types that carry kidney disease
> heritability.

This is a measurement paper, not a modelling paper. Its value is that it is
verifiable, reproducible from public files, and that nobody has done it.

### 1a. FRAMING, settled on day 1 after checking ENCODE. Read this before writing a word.

The obvious paper is "DeepMind under-sampled kidney". **That paper is wrong, and
writing it would be both unfair and easy to rebut.**

ENCODE holds 592 human kidney-organ experiments. Of the 228 transcription-factor
ChIP experiments among them, **221 are in cell lines, covering 218 distinct
factors. Exactly 7 are in kidney tissue or primary cells, and all 7 are CTCF**,
spanning 2011 to 2021 (ENCSR000DXD, ENCSR000DMC, ENCSR000DVH, ENCSR194WQV,
ENCSR966OUM, ENCSR799GJD, ENCSR266ZUX).

So the model's CTCF-only kidney transcription-factor coverage is not a selection
failure. **It is a faithful reflection of the fact that the field has never
ChIP-seq'd a single non-CTCF transcription factor in human kidney tissue.**

The correct claim is therefore a **field-level data gap with consequences for
anyone using these models on kidney**, not a criticism of the model developers.
That version is stronger, fairer, survives the obvious reviewer rebuttal ("that
is simply what ENCODE contains"), and points at something nephrology can
actually fix by generating data.

Two honest caveats that belong in the Discussion rather than being suppressed.
The model does not take quite everything that exists: ENCODE has one intact Hi-C
and two PRO-cap experiments in kidney, and the model exposes zero kidney contact
maps and zero kidney PRO-cap. Small, but it should be stated rather than
discovered by a reviewer.

## 2. Why this is publishable, stated honestly

Three things have to be true, and I have checked all three.

**Nobody has done it for kidney.** A search of 2,237 Europe PMC hits found no
benchmark of these models stratified by kidney tissue. The nearest work,
Kathail et al., *Genome Biol* 25:202 (2024), establishes the method of
stratified evaluation and does not do kidney.

**There is a hard external anchor.** Loeb et al., *Nat Genet* 56:2078-2092
(2024) put 56% of kidney-function heritability in tubule epithelial candidate
regulatory elements and 7% in podocyte elements. That is the number the census
gets compared against, and it is what turns a track count into a finding.

**The model's own documentation concedes the general point.** The AlphaGenome
FAQ states that "accurately capturing tissue-specific effects and long-range
genomic interactions remains challenging" and that the model "has not yet been
benchmarked for predicting individual (personal) human genomes."

## 3. Results in hand as of day 1

All produced by `analysis/01_track_census.py` and
`analysis/02_textproto_crosscheck.py`, both reproducible from `data/`.

**VERIFIED THREE WAYS, 13 Sep 2026.** The census now agrees exactly across three
independent sources: the notebook dataframe, the shipped protobuf manifest
(once its 563 padding rows are removed), and **the deployed model queried live
through the API** (`analysis/08_live_api_crosscheck.py`). All eleven modalities
match to the track, total 5,563, kidney parenchyma 103 in all three, and the
live API independently confirms four kidney transcription-factor tracks, all
CTCF. The census is as solid as it can be made without DeepMind's cooperation.

| Quantity | Value |
|---|---|
| Human output tracks, total | 5,563 |
| Kidney parenchyma tracks | 103 (1.85%) |
| Kidney parenchyma biosamples | 20 of 715 |
| Urinary tract tracks, reported separately | 25 |
| HEK293/HEK293T tracks | 116 |

By modality, kidney parenchyma against the total for that modality:

| Modality | Kidney | Total | Share |
|---|---|---|---|
| Transcription factor ChIP | 4 | 1,617 | 0.25% |
| ATAC | 1 | 167 | 0.60% |
| Histone ChIP | 8 | 1,116 | 0.72% |
| CAGE | 12 | 546 | 2.20% |
| RNA-seq | 23 | 667 | 3.45% |
| DNase | 13 | 305 | 4.26% |
| Splice site usage | 28 | 734 | 3.81% |
| Splice junctions | 14 | 367 | 3.81% |
| Contact maps | 0 | 28 | 0% |
| PRO-cap | 0 | 12 | 0% |

**The four sharpest facts.**

1. **All four kidney transcription-factor tracks are CTCF.** Not one
   kidney-relevant transcription factor (HNF1B, PAX2, WT1, PPARGC1A, ESRRB) has
   been assayed in kidney tissue. HEK293 alone carries 111 TF tracks, 27.8 times
   the kidney total, and covers 110 distinct factors.
2. **Kidney parenchyma has exactly one ATAC track**, the whole adult organ. No
   cell-type-resolved kidney chromatin accessibility exists in the model at all.
3. **The entire distal nephron is absent.** Collecting duct, thick ascending
   limb, loop of Henle, distal convoluted tubule, principal cell, intercalated
   cell, macula densa, juxtaglomerular apparatus, parietal epithelial cell,
   connecting tubule, inner medulla and renal papilla: twelve of twelve absent.
   Podocyte appears once, as a single DNase track, from a child donor.
4. **Kidney tracks are developmentally skewed.** 37% of kidney tracks are
   embryonic against 15% across all tracks, and only 25% are adult against 46%
   overall. Much of what the model knows about kidney, it learned from fetal
   tissue.

## 4. Three methodological traps, all found and fixed on day 1

These belong in the Methods section. Each one would have produced a wrong
number, and the first two would have inflated the kidney count in our favour,
which is the direction a reviewer will look hardest at.

**"adrenal gland" contains the string "renal".** Substring matching on /renal/
captures 27 adrenal tracks plus 2 perirenal adipocyte and 2 perirenal
preadipocyte tracks, 31 in total. Every biosample in this analysis is therefore
classified by an explicit curated name list, never by substring. A guard
function fails the run if a future release introduces a kidney-like name that
has not been curated.

**"papilla" matches tongue and hair follicle.** Testing for renal papilla with
the bare term returns PRESENT because of "fungiform papilla" and "hair follicle
dermal papilla cell". Qualified to "renal papilla", it is correctly ABSENT.

**The canonical manifest contains 563 padding placeholders.** The protobuf
manifest has 6,126 rows against the notebook dataframe's 5,563. The difference
is entirely rows literally named "Padding", which pad each modality's output
head to a round size: PRO-cap to 128, ATAC to 256, DNase to 384, CAGE to 640,
RNA-seq to 768, histone to 1,152, TF to 1,664. Counting them inflates every
denominator. Removing them reconciles the manifest to the dataframe exactly, in
every modality. A near-miss worth recording: identifying padding by "empty
biosample" instead of by name also discards the four real splice-site tracks,
which are annotation-derived and legitimately tissue-free, and breaks
reconciliation by exactly four.

## 5. What remains to be done

**Leg 2, the quantitative re-analysis. Part 1 is DONE**
(`analysis/03_gtex_tissue_weighting.py`).

The headline accuracy figure in the AlphaGenome paper is a Spearman correlation
averaged across tissues *weighted by the number of eQTLs in each tissue*. That
weighting is what buries kidney, and the size of the effect is now measured:

| | GTEx v8 (benchmark release) | GTEx v10 (current) |
|---|---|---|
| Tissues with eQTLs | 49 | 50 |
| Kidney cortex eQTL n | 73 (median 233) | 103 (median 304) |
| Kidney cortex eGenes | 1,260 (median 9,642) | 2,718 (median 12,574) |
| Rank by sample size | **1st, smallest** | 2nd (Bladder smaller) |
| Rank by eGenes | **1st, smallest** | 2nd |
| Weight vs unweighted | **down-weighted 7.71x** | down-weighted 4.47x |

**ALWAYS NAME THE RELEASE.** In v8, the release the model was benchmarked on,
kidney cortex is the smallest tissue of 49 on both measures. In v10 it is second
of 50, because Bladder was added. An unqualified "fewest eQTLs of any tissue" is
false today and a reviewer will catch it. The script computes both side by side
so the manuscript cannot drift.

Kidney medulla carries weight **zero** in both: 11 RNA-seq samples in v10 and
eQTLs never called, in any release.

Sensitivity, framed as a hypothetical because kidney's true performance is not
yet measured: a kidney-specific Spearman deficit of any size shows up **4.5 to
7.7 times more visibly** in an unweighted average than in the published one.

**Part 2, DONE: ground truth pulled.** Decision taken 13 Sep 2026 to re-run
against the eQTL Catalogue rather than lift per-tissue numbers from the published
figures. Three reasons, and the third decides it:

1. **Uniform processing.** All 49 GTEx tissues re-quantified through one
   pipeline, so per-tissue comparison is not confounded by pipeline differences.
2. **Licence.** CC BY 4.0, commercial use permitted. The other large kidney eQTL
   resource is licence-conflicted and cannot safely anchor a published figure.
3. **Same vintage as the benchmark.** The catalogue's GTEx is v8, kidney cortex
   n=73, which reproduces the GTEx API v8 figures exactly and independently.
   Benchmarking against v10 would compare the model to data that did not exist
   when it was evaluated.

Confirmed on disk: the catalogue carries **all 49 GTEx v8 gene-expression
datasets**, n=73 to 702, median 233, with kidney cortex the smallest. Two files
per tissue are enough and the full 2.4 GB nominal sumstats are not needed:
`permuted` (one lead variant per gene, ~0.9 MB) and `susie` (fine-mapped credible
sets, ~0.3-1 MB). Roughly 100 MB total. Fetcher is
`analysis/05_fetch_eqtl_catalogue.py`, throttled at 1.5 s with backoff and
resume, because the catalogue blacklists addresses it reads as denial-of-service.

**Part 3, benchmark set construction** (`analysis/06_build_benchmark_set.py`).
Positives are variants in a SuSiE credible set with posterior inclusion
probability at or above 0.9 **and** credible-set size at or below 5. Both filters
are needed: probability alone admits variants from 200-member sets where the
posterior is diffuse, and the size cap requires that fine-mapping actually
resolved the signal. Observed beta and z are the regression target.

**The set is frozen before any model runs.** That is deliberate. Building and
publishing the evaluation set first is what stops the thresholds being tuned to a
flattering answer, and it is the difference between a benchmark and a
demonstration.

**Part 4, DONE, and it changed the paper's conclusion**
(`analysis/07_power_and_threshold_sensitivity.py`).

Building the set produced a result bigger than the weighting argument.

| Tissue | Fine-mapped positives |
|---|---|
| Kidney cortex | **59** |
| Vagina (2nd smallest) | 146 |
| Median across 49 tissues | 616 |
| Thyroid (largest) | 1,756 |

Kidney cortex is **last of 49**, with fewer than half the next smallest, and
30-fold below the largest. Total across all tissues is 32,818 variants.

**Challenge 1, is it my thresholds?** No. Swept across 36 settings (posterior
probability 0.5 to 0.99, credible-set size 1 to uncapped), **kidney is the
smallest tissue at every single one**, always 27 to 33-fold below the largest.
The filter is not producing the result; the data is. The grid goes in the
supplement.

**Challenge 2, what do 59 variants buy?** A 95% interval on Spearman's rho at
n=59, using the Bonett-Wright standard error for rank correlation:

| True rho | Kidney (n=59) | Thyroid (n=1,756) |
|---|---|---|
| 0.0 | -0.26 to 0.26 | -0.05 to 0.05 |
| 0.4 | 0.15 to 0.60 | 0.36 to 0.44 |
| 0.6 | 0.39 to 0.75 | 0.57 to 0.63 |
| 0.8 | 0.66 to 0.89 | 0.78 to 0.82 |

Kidney's interval is roughly five to six times wider throughout. At a true rho
of 0.6 it spans 0.39 to 0.75, which cannot separate "performs like other
tissues" from "performs materially worse".

### THE CONCLUSION THIS FORCES, and it is better than the one intended

Do **not** write "the model performs poorly in kidney". We cannot show that, and
attempting it would be the same overreach as the retracted 86% in the HLA paper.

Write instead: **the public evidence cannot currently support any confident
statement about kidney performance, and the reason is measurable.** Kidney is
down-weighted 7.71-fold in the published average, and separately it is the
tissue with the fewest fine-mapped positives by a factor of thirty, so even a
correctly unweighted per-tissue benchmark yields an interval too wide to
interpret.

That claim is unassailable, it is a call to action rather than an accusation,
and it sits naturally alongside the §1a finding that the underlying assay data
does not exist either. Same lesson as the HLA paper: the negative result was the
publishable one.

## 5a. THE GATE BEFORE LEG 2 CAN FINISH

Everything to this point needed no permission. The next step does, and it is a
decision for CL personally rather than something to be arranged around.

Model predictions require either an **AlphaGenome API key** or **Atlas
precomputed scores**. Both are free for non-commercial research and both require
accepting terms **in a personal capacity, never on behalf of a company**. The
weights and outputs are restricted to non-commercial use by non-commercial
organisations, so this work is authored personally.

Scale matters for the choice. The documentation states the live API suits
analyses "requiring 1000s of predictions" and is "likely not suitable for large
scale analyses requiring more than 1 million predictions", while the precomputed
Atlas "will typically have a larger query rate". Whether 49 tissues of fine-mapped
positives fits the API or needs the Atlas is answered by the output of script 06,
which is why that script reports the exact variant count before anything is run.

**Leg 3, independent validation. Part 1 is DONE**
(`analysis/04_encode_kidney_inventory.py`).

ENCODE human kidney inventory, curated the same way as §4 so cell lines cannot
contaminate the count:

| | Experiments | Biosamples |
|---|---|---|
| Tissue or primary cell | 284 | 23 |
| Cell line (excluded from headline) | 308 | 9 |
| of which HEK293 / HEK293T | 285 | 2 |

Assays available in kidney tissue or primary cell: DNase-seq 85, long-read
RNA-seq 58, polyA RNA-seq 27, Mint-ChIP 24, histone ChIP 22, DNA methylation
array 8, total RNA-seq 8, **TF ChIP 7 (all CTCF)**, RAMPAGE 5, RRBS 4, PRO-cap 2,
snRNA-seq 2, snATAC-seq 2, ATAC-seq 2, intact Hi-C 1.

That table is the evidence base for the framing in §1a, and it is also the menu
for leg 3. DNase-seq at 85 experiments is the only kidney assay with real depth.

**Part 2, still to do:** score variants against kidney measurements the model did
not train on, separating tubule from glomerulus, which maps directly onto the
compartment asymmetry in §3.

**Data for legs 2 and 3, all open and already identified.**

| Source | Content | Licence |
|---|---|---|
| ENCODE kidney | 570 GRCh38 bigWigs | Free download, analyse and publish, no restrictions |
| GTEx v10 kidney cortex | n=103, 2,718 eGenes | No commercial restriction |
| EMBL-EBI eQTL Catalogue | GTEx kidney cortex, 5 quantifications | CC BY 4.0 |
| NephQTL2 | 240 glomerular, 311 tubulointerstitial | Not stated |
| Susztak kidney atlas | n=686 meta, tubule 356 / glomerulus 303 | **Conflicted, see §6** |

Kidney medulla has never had eQTLs called, in any release. It is 11 RNA-seq
samples and zero eQTL samples in v10. That is worth one sentence in the paper.

## 6. Risks, and what to do about each

**A competitor is close.** One group at UCSF has a nephrologist who co-authored
both the kidney heritability paper and the stratified-benchmark paper that
defines this method. The data, the clinician and the technique are already in
one lab, and they have not published. Mitigation: the descriptive leg in §3 is
complete today. Preprint it rather than waiting for legs 2 and 3.

**The Susztak licence conflict.** The download files carry CC BY 4.0 on the
repository, while the lab's own site agreement bars commercial use *and* bars
publishing any result without the principal investigator's written consent.
Mitigation: use ENCODE, GTEx and the eQTL Catalogue as primary. If the Susztak
atlas is genuinely needed, resolve it in writing first. Do not build a
validation figure on it and discover the problem at submission.

**Licence on the model itself.** The weights are non-commercial and restricted
to non-commercial organisations. This work is therefore authored personally, not
on behalf of any company. The model *code* is Apache 2.0, which is what makes a
future trained-model paper commercially clean, but that is Paper 2 territory.

## 7. Sequence, and why this order

Paper 1 is the methods and caution paper. Paper 2 is APOL1: whether predicted
regulatory variation on African haplotype backgrounds explains part of the
incomplete penetrance of the high-risk genotype.

Paper 1 is the justification for Paper 2's method. You cannot fine-map APOL1
regulation from kidney expression data, because GTEx holds 133 significant
APOL1 eQTLs and not one is in kidney, and kidney cortex is the smallest tissue
in the resource. That is exactly when a sequence model earns its place, and
saying so honestly, having measured it, is stronger than asserting it.

APOL1 is the time-critical one. The Atlas is days old and APOL1 is an obvious
target for anyone with the same idea. Paper 1 is nearly written; Paper 2 should
start as soon as leg 1 is preprinted, not after Paper 1 is accepted.

**Standing constraints that apply to both.** File before disclosing, because
Europe has no grace period. Keep commercial positioning until after 22 January
2027. Publish under the personal clinical identity, not a company banner.

## 8. Target and deadline

**UK Kidney Week 2027**, Newport, 7 to 9 June 2027. Abstracts close
**31 December 2026, hard, no extensions**, decisions early January. Categories
include Rare Disease and GN, AKI/Basic Science and Diagnostics, and
Transplantation.

Journal: the descriptive and benchmark legs together suit *Genome Biology* or
*Kidney International Reports*. Decide after leg 2.

## 9. Reproducing everything here

```
cd "kidney-alphagenome"
python3 analysis/01_track_census.py        # census, guards, tables
python3 analysis/02_textproto_crosscheck.py # independent manifest parse
```

Both write to `results/`. Neither needs network access. See `data/PROVENANCE.md`
for the origin of every input file.
