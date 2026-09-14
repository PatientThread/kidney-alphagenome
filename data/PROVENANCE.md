# Provenance of input files

Every file in `data/` with its origin, retrieval date, and licence. Nothing in
this directory was hand-edited. Retrieved 13 September 2026.

---

## alphagenome_tracks_human_raw.csv

**5,563 rows, one per human output track. 16 columns.**

Origin: the embedded dataframe in DeepMind's own notebook
`colabs/tissue_ontology_mapping.ipynb` in `github.com/google-deepmind/alphagenome`.
This is a convenience export, not the authoritative artefact, which is why it is
cross-checked against the protobuf manifest below.

Columns: `index, name, strand, assay_title, ontology_curie, biosample_name,
biosample_type, biosample_life_stage, data_source, endedness,
genetically_modified, nonzero_mean, output_type, gtex_tissue, histone_mark,
transcription_factor`.

Licence: Apache 2.0 (the repository's code licence). Note this is the CODE
licence and it does permit commercial use. It is the model *weights* that are
restricted to non-commercial use by non-commercial organisations, under a
separate custom licence. Do not conflate the two.

## OutputMetadataResponse_ORGANISM_HOMO_SAPIENS.textproto

**6,126 metadata blocks across 11 output-type sections. 2.6 MB.**

Origin: `src/alphagenome_research/model/metadata/` in
`github.com/google-deepmind/alphagenome_research`. This is the authoritative
manifest shipped with the research code.

**563 of the 6,126 rows are placeholders literally named "Padding".** They pad
each modality's output head to a round size and carry no biosample. Excluding
them reconciles this file to the CSV above exactly, modality by modality:

| Modality | Real | Padding | Head size |
|---|---|---|---|
| CHIP_TF | 1,617 | 47 | 1,664 |
| CHIP_HISTONE | 1,116 | 36 | 1,152 |
| RNA_SEQ | 667 | 101 | 768 |
| SPLICE_SITE_USAGE | 734 | 0 | 734 |
| CAGE | 546 | 94 | 640 |
| SPLICE_JUNCTIONS | 367 | 0 | 367 |
| DNASE | 305 | 79 | 384 |
| ATAC | 167 | 89 | 256 |
| CONTACT_MAPS | 28 | 0 | 28 |
| PROCAP | 12 | 116 | 128 |
| SPLICE_SITES | 4 | 1 | 5 |
| **Total** | **5,563** | **563** | **6,126** |

Licence: Apache 2.0, as above.

## eqtl_catalogue_dataset_metadata_r7.tsv

Origin: `github.com/eQTL-Catalogue/eQTL-Catalogue-resources`, path
`data_tables/dataset_metadata_r7.tsv`. Release 7, June 2024. 758 datasets across
42 studies.

Kidney content is five rows, all the same study: GTEx kidney cortex at n=73
(the v8 vintage), study QTS000015, datasets QTD000261 through QTD000265,
covering gene expression, exon, transcript, txrev and leafcutter quantifications.
No glomerulus, no tubule, no NephQTL, no Susztak.

Licence: data CC BY 4.0, code Apache 2.0. Commercial use permitted.

Operational note: the REST API is retired and returns HTTP 410. Use the FTP at
`ftp.ebi.ac.uk/pub/databases/spot/eQTL`, and throttle tabix requests or the
server treats them as a denial-of-service attempt and blacklists the address.

## eqtl_catalogue_dataset_metadata_r8_beta.tsv

Origin: same repository, r8 beta metadata, January 2026, full release planned
Q3 2026. **Contains zero kidney datasets**; it is blood, brain and lung cell
types. Kept here so that the claim "r8 adds no kidney" is checkable rather than
asserted.

Licence: as above.

---

## Not yet retrieved, needed for legs 2 and 3

| Source | What | Access | Licence |
|---|---|---|---|
| ENCODE kidney | 570 GRCh38 bigWigs, 592 human kidney-organ experiments | Open, no account, 10 GET/sec limit | "Freely download, analyze and publish... without restrictions" |
| GTEx v10 kidney cortex | eQTL tarball, n=103, 2,718 eGenes | `storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/`, anonymous | No commercial restriction |
| GTEx v10 nominal all-pairs | Non-significant pairs, needed as negatives for calibration | **Requester-pays** bucket `gs://gtex-resources`; anonymous listing returns 401 | Cost, not permission |
| NephQTL2 | 240 glomerular + 311 tubulointerstitial, 332 individuals with WGS | nephqtl2.org, download offered, no account mentioned | **Not stated anywhere** |
| Susztak kidney atlas | eQTL n=686 meta, tubule 356 / glomerulus 303; meQTL n=443 | figshare, open | **CONFLICTED, see below** |

**Susztak licence conflict, unresolved.** The download links point at figshare
article 15183495, which the figshare API reports as CC BY 4.0. The lab's own
user agreement at `susztaklab.com/agree.php` grants a licence "specifically
excluding commercial use or sales" and states "You may not publish the results
of any work resulting from use of the Data without the express written consent
of the Principal Investigator". A publication veto is not a normal data
condition. Two live instruments, opposite effect, same bytes. Resolve in writing
before this atlas appears in any figure.

**Build trap.** Susztak snATAC deposits GSE200547 and GSE172008 are **hg19**;
GSE211785 is GRCh38. Anything compared against AlphaGenome or Borzoi needs
liftOver first.
