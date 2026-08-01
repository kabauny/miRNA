# Data provenance, methodology, and governance

A record of **where the data came from, how it was retrieved, how it was analysed,
and how each result was produced** — for research-governance and reproducibility
review. Nothing in this project generated primary human data; it is a **secondary,
computational analysis of public, de-identified TCGA data**.

## 1. Data provenance — origin of the data

### 1.1 Primary source: TCGA
All molecular and clinical data originate from **The Cancer Genome Atlas (TCGA)**, a
public research program of the US National Cancer Institute (NCI) and National Human
Genome Research Institute (NHGRI). TCGA collected surgically **resected primary
tumour** (and matched normal) biospecimens from consented patients, under each
contributing institution's IRB, and generated standardized molecular assays:

- **mRNA expression** — RNA-seq (used here as RSEM-normalized gene-level values);
- **copy number** — SNP-array segmentation summarized by **GISTIC** (discrete −2…+2);
- **mutations** — whole-exome sequencing (MAF calls);
- **clinical** — abstracted by the Biospecimen Core Resource (BCR).

Two cancer types were used: **lung adenocarcinoma (LUAD)** and **lung squamous-cell
carcinoma (LUSC)**, PanCancer Atlas freeze. **This project did not collect, consent,
sequence, or otherwise generate any of this data** — it reuses TCGA's public output.

### 1.2 Data tier and access terms (governance-relevant)
Only **open-access, de-identified, processed (summary-level) data** was used:
gene-level expression, GISTIC copy-number calls, mutation calls, and clinical
tables. **No controlled-access data** (raw sequence reads, germline variants) was
accessed, so no dbGaP authorization was required. All patient identifiers are TCGA
**barcodes** (e.g. `TCGA-05-4384`); no direct identifiers are present, and time
fields are relative "days-to" offsets as released by TCGA.

## 2. Methodology — how the data was retrieved

### 2.1 cBioPortal REST API (expression, copy number, mutations, clinical)
Retrieved programmatically via the public cBioPortal API (`https://www.cbioportal.org/api`;
client in `src/mirna_tcga/cbioportal.py`). Exact identifiers (in `config.yaml`):

- studies: `luad_tcga_pan_can_atlas_2018`, `lusc_tcga_pan_can_atlas_2018`;
- profiles: `*_rna_seq_v2_mrna` (expression), `*_gistic` (discrete CNA), `*_mutations`;
- sample lists: `*_all`, `*_cna`, `*_sequenced`.

Expression was streamed in gene chunks; sparse copy-number/mutation *events* were
fetched where possible to minimize transfer. No credentials are involved (public API).

### 2.2 BCR Biotab clinical supplements (local GDC download)
The richer distant-metastasis annotation (clinical M stage, follow-up
`new_tumor_event_type`, metastatic sites) is **not** in cBioPortal, so it was taken
from the TCGA **BCR Biotab** clinical supplements — the open-access files fetched by
`TCGAbiolinks::GDCdownload(data.type = "Clinical Supplement", data.format = "BCR
Biotab")` from the NCI **Genomic Data Commons (GDC)**. These were provided locally
under `TCGA/GDCdata/` and parsed by `src/mirna_tcga/biotab.py` (which handles the
three-row Biotab header and per-study column reordering). This local folder is
**large and kept out of version control** (git-ignored); it contains only public,
de-identified clinical text.

### 2.3 Third-party reference data
Pathway gene sets (KEGG, MSigDB Hallmark) were loaded from public GMT mirrors on
GitHub for over-representation tests only; no patient data leaves the local
environment. miRNA reference (mature.fa) and gene-set files are similarly public.

## 3. Methodology — how the data was analysed

### 3.1 Cohort, endpoints, and comparison groups
- **NSCLC cohort** = LUAD + LUSC, keyed by patient barcode across layers.
- **Distant metastasis** = pathologic M1 **∪** clinical M1 **∪** follow-up "Distant
  Metastasis" event (Biotab-enriched: ~33 → 126 patients).
- **Nodal metastasis** = N1/N2/N3 vs N0.
- **"True stage I"** (indolent reference) = stage I at diagnosis, N0/M0, ≥ 2-year
  follow-up free of distant *and* locoregional recurrence.

### 3.2 Statistical methods
- **Copy-number depletion/enrichment**: Cochran–Mantel–Haenszel test **stratified by
  subtype** (`associate.cmh_depletion_screen`).
- **Differential expression**: subtype-stratified **van Elteren rank-sum**
  (`associate.ranksum_screen`); z > 0 = higher in the index group.
- **Survival**: vectorized **Cox score test** (`screen.cox_score_screen`, validated
  against `lifelines`), plus **Kaplan–Meier** estimation and a two-group **log-rank**
  test implemented directly (`scripts/26_stc1_km_curve.py`).
- **Multiple testing**: Benjamini–Hochberg FDR (q-values) within each screen's gene
  universe.

### 3.3 Confounder handling (central to the conclusions)
Metastatic tumours are more proliferative and more advanced, so three explicit
adjustments were applied and reported:
1. **Subtype stratification** (LUAD vs LUSC) in every test — the groups have
   different subtype mixes.
2. **Proliferation adjustment** — strata of subtype × proliferation-tertile
   (`panels.PROLIFERATION_PANEL`); genes/associations that vanish here were
   proliferation-driven (MKI67 was used as a positive control that correctly vanishes).
3. **Dosage-independent test** — expression associations re-tested only in tumours
   where the gene is **not** deleted, to separate cis-dosage from regulation.
4. **Clean reference group** — "true stage I" instead of the heterogeneous "all M0".

### 3.4 Software
Python 3.11; `numpy`, `pandas`, `scipy`, `scikit-learn` for statistics;
`matplotlib` for the figure; `reportlab` for the PDF. Analysis code is the packaged
`mirna_tcga` library with an offline test suite (`pytest`, mocked API / synthetic
data — 69 tests). No proprietary software.

## 4. How the results were produced (reproducibility)

Every number and figure maps to a **versioned script** run against the public data
above; results are regenerated, not hand-curated:

| Result | Script | Commit theme |
|---|---|---|
| No positive bulk signature (DE / CNV) | `17`–`22` | proliferation dominates |
| 8p deleted in true-stage-I, retained in metastatic (q≈0.027) | `23_true_stage_i_lost_genes.py` | negative selection |
| Arm resolved to STC1 (dosage + proliferation adjusted) | `24_chr8p23_which_gene.py` | which gene |
| STC1 → worse OS, proliferation-independent | `25_stc1_survival_check.py` | survival check |
| Kaplan–Meier figure | `26_stc1_km_curve.py` | figure |

Reproducibility controls: all parameters and identifiers live in `config.yaml`;
random seeds are fixed where models are used; the full history is in git (each result
is a labelled commit, pushed to the `master` branch). Derived tables/figures are
regenerated by re-running the scripts — the repository stores the **code and prose
write-ups**, not large derived data.

## 5. Governance, limitations, and appropriate use

- **Nature of the work**: secondary, **observational/associational** computational
  analysis. It generates hypotheses; it does not establish causation.
- **Not clinical-grade**: no analytical/clinical validation, not for diagnosis,
  prognosis, or patient management. The STC1 finding is an internal-consistency
  result awaiting external validation (e.g. MET500 / MSK-MET).
- **Statistical caveats**: exploratory multiple testing; the headline 8p result rests
  partly on a sparse contingency cell; "higher in aggressive disease" associations may
  reflect the tumour micro-environment (e.g. hypoxia) rather than a driver.
- **Data-use compliance**: TCGA data are open-access and, under current NIH policy,
  free of publication embargo; use should acknowledge TCGA/NCI-NHGRI and cite the
  relevant marker papers. No attempt was made — or is possible from these data — to
  re-identify individuals.
- **Privacy**: no PHI processed; no data transmitted to third parties beyond the
  read-only public APIs named above.
