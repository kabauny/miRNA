# miRNA / mRNA / mutation analysis (TCGA)

A small, modern Python toolkit for exploring TCGA lung-cancer data (LUAD &
LUSC): expression signatures, mutation-associated signatures, and clinical
(survival) outcomes.

This is a rewrite of an older R project. The original R scripts are kept under
[`legacy/`](legacy/) for reference; the analysis method (nearest **shrunken
centroids**, a.k.a. PAM) is preserved, now via scikit-learn.

## What changed from the R version

| Concern | Old (R) | New (Python) |
|---|---|---|
| Data access (mRNA/mutations) | manual `MANIFEST.txt` parsing of local files | **cBioPortal REST API** (`mirna_tcga.cbioportal`) |
| Data access (miRNA) | local files | **UCSC Xena** TCGA hub (`mirna_tcga.xena`) |
| Classification | `pamr` | `sklearn.neighbors.NearestCentroid` + CV threshold |
| ID conversion | `biomaRt` | `mygene` / `pybiomart` (`mirna_tcga.idmap`) |
| Normalization | `MASS::boxcox`, custom NZV | `scipy.stats.boxcox`, vectorized filters |
| Survival | (none) | `lifelines` (`mirna_tcga.survival`) |
| Reproducibility | hard-coded `/Users/...` paths | `config.yaml`, packaged, tested |

## Layout

```
src/mirna_tcga/      # the package
  config.py          # load config.yaml
  cbioportal.py      # cBioPortal REST client (mRNA expression, mutations, clinical)
  xena.py            # UCSC Xena loader (TCGA miRNA expression matrices)
  preprocess.py      # variance filter, Box-Cox, standardize
  classify.py        # PAM (nearest shrunken centroids) + CV + signature
  idmap.py           # gene/miRNA id conversion (optional deps)
  survival.py        # Kaplan-Meier / Cox PH (optional deps)
  screen.py          # genome-wide survival screen (vectorized Cox score test)
  associate.py       # binary-endpoint screens (stratified rank-sum, Fisher)
  endpoints.py       # metastasis endpoints (distant M1, nodal N+) from clinical
  biotab.py          # TCGA BCR Biotab clinical-supplement loader (richer distant-met label)
  enrich.py          # pathway over-representation (GMT + hypergeometric)
  layers.py          # multi-omic loaders (expression, deletions, mutations)
  mirbase.py         # MIMAT accession -> hsa-miR name lookup (mature.fa)
  panels.py          # small demo gene panel
scripts/             # runnable pipeline examples
  01_mirna_subtype_signature.py  # LUAD vs LUSC miRNA signature (Xena)
  02_subtype_signature.py        # LUAD vs LUSC mRNA signature (cBioPortal)
  03_mutation_analysis.py        # TP53-mutant vs wild-type signature
  04_survival.py                 # overall survival by stage
  05_signature_survival.py       # subtype signature -> survival (within LUAD)
  06_signature_overlap.py        # overlap between mRNA & miRNA signatures
  07_nsclc_clinical.py           # NSCLC (LUAD+LUSC) clinical summary
  08_nsclc_expression_survival.py # OS in NSCLC stratified by one gene/miRNA
  09_survival_screen.py          # genome-wide mRNA OS screen + pathway enrichment
  10_mirna_survival_screen.py    # genome-wide miRNA OS screen (Xena + cBioPortal OS)
  11_cnv_mutation_screen.py      # deep-deletion & mutation screens vs OS + metastasis
  12_constellation_model.py      # cross-validated multi-omic OS / metastasis models
  13_constellation_sparse.py     # driver-focused L1 model + miRNA layer (constellation)
  14_metastasis_spared_deletions.py # genes never deleted in metastasis (spared machinery)
  15_pancancer_spared_deletions.py  # does the spared-deletion signal generalize? (LUAD-specific)
  16_never_deleted_stage_iv.py      # strict "never deleted in stage IV" re-run on the 126-patient cohort
tests/               # offline tests (synthetic data + mocked API)
config.yaml          # studies, profiles, parameters
legacy/              # original R scripts (reference only)
```

## Setup

```bash
pip install -e .              # core
pip install -e '.[idmap]'     # + mygene / pybiomart for ID conversion
pip install -e '.[survival]'  # + lifelines for survival analysis
pip install -e '.[dev]'       # + pytest / ruff
```

## Usage

```bash
python scripts/01_mirna_subtype_signature.py      # LUAD vs LUSC (miRNA, Xena)
python scripts/02_subtype_signature.py            # LUAD vs LUSC (mRNA, cBioPortal)
python scripts/03_mutation_analysis.py --gene TP53 --study luad
python scripts/04_survival.py --study luad        # needs [survival] extra

# Overall survival in the NSCLC cohort (LUAD + LUSC) split by expression of one
# feature -- Kaplan-Meier medians, log-rank, and Cox PH (needs [survival] extra):
python scripts/08_nsclc_expression_survival.py --gene EGFR        # mRNA (cBioPortal)
python scripts/08_nsclc_expression_survival.py --gene MKI67 --by-study-median
python scripts/08_nsclc_expression_survival.py --mirna hsa-mir-21 # miRNA (Xena)

# Genome-wide screen: every protein-coding gene tested for an OS association in
# NSCLC (subtype-stratified Cox score test, FDR-controlled), then pathway
# over-representation of the significant genes (needs [survival] extra):
python scripts/09_survival_screen.py                 # mRNA: full transcriptome (~20k genes)
python scripts/09_survival_screen.py --max-genes 2000 --save-dir results  # quick
python scripts/10_mirna_survival_screen.py           # miRNA screen (needs Xena egress)

# Copy-number deletions + mutations vs survival AND metastasis, then a
# cross-validated multi-omic model of the two outcomes:
python scripts/11_cnv_mutation_screen.py --save-dir results
python scripts/12_constellation_model.py --save-dir results
python scripts/13_constellation_sparse.py --save-dir results  # + miRNA layer, L1

# Negative-selection screen: genes deletable in non-metastatic tumours but
# spared from deletion in metastatic disease (candidate metastasis machinery).
# Auto-enriches the distant-metastasis label from local BCR Biotab clinical
# supplements (config biotab.root / --biotab-root) -- ~33 -> 126 metastatic in
# NSCLC -- falling back to pathologic M1 if no Biotab folder is present:
python scripts/14_metastasis_spared_deletions.py --save-dir results
python scripts/15_pancancer_spared_deletions.py --save-dir results  # generalize across cancers
python scripts/15_pancancer_spared_deletions.py --protection        # del OR truncating mutation
```

Worked results are written up in
[`docs/nsclc_survival_screen.md`](docs/nsclc_survival_screen.md) (mRNA/miRNA
survival + pathways) and
[`docs/nsclc_multiomics_metastasis.md`](docs/nsclc_multiomics_metastasis.md)
(copy-number, mutations, metastasis, and the constellation model), and
[`docs/nsclc_metastasis_spared_deletions.md`](docs/nsclc_metastasis_spared_deletions.md)
(the negative-selection screen for metastasis-required machinery).

Or from Python:

```python
from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.classify import fit_pam_cv, signature
from mirna_tcga.preprocess import standardize, drop_near_zero_variance

cfg = load_config()
client = CBioPortalClient(**cfg.cbioportal)
mat = client.expression_matrix(
    cfg.mrna_profile("luad"), ["TP53", "EGFR", "KRAS"], cfg.all_samples_list("luad")
)
```

The scripts use a small marker gene panel (`panels.py`) so they run quickly.
For a real analysis, pass the full gene universe instead.

## Tests

```bash
pytest            # all offline; no network required
```

The cBioPortal client is tested against a mocked HTTP session, so the suite
runs without internet access.

## Network access note

The pipeline calls external hosts: `https://www.cbioportal.org/api` (mRNA,
mutations, clinical), `https://*.xenahubs.net` (miRNA), and
`https://raw.githubusercontent.com` (pathway gene-set GMTs for the survival
screen). In a sandboxed/remote environment you may need to **add these hosts to
the network egress allowlist** — otherwise requests return
`403 Host not in allowlist`. See
https://code.claude.com/docs/en/claude-code-on-the-web.

> **miRNA survival note.** The *miRNA* survival screen pairs UCSC Xena miRNA
> matrices with cBioPortal OS (matched by patient barcode); it needs egress to
> `tcga.xenahubs.net`. Xena serves each matrix **gzipped at `<dataset>.gz`** (the
> bare key returns S3 `AccessDenied`), which `mirna_tcga.xena` handles. The miRNA
> profiles cBioPortal hosts for lung (`lusc_tcga_pub`, CPTAC LUAD/LUSC) lack OS
> fields, so Xena is required for the miRNA arm.

## Data sources

- **mRNA / mutations / clinical** → cBioPortal (`mirna_tcga.cbioportal`).
- **miRNA expression** → UCSC Xena TCGA hub (`mirna_tcga.xena`), because
  cBioPortal does not carry miRNA for the PanCancer Atlas studies. Both Xena
  hubs used here (`tcga`, `gdc`) serve genuine TCGA data.
- **Richer clinical (metastasis)** → TCGA **BCR Biotab** clinical supplements
  (`mirna_tcga.biotab`), a local GDC download (`TCGAbiolinks::GDCdownload`,
  `data.format = "BCR Biotab"`) kept under the git-ignored `TCGA/` folder and
  pointed at by `config biotab.root`. Adds clinical M stage + follow-up
  `Distant Metastasis` events that cBioPortal's pathologic-M field lacks.

## Caveats

- Xena miRNA matrices (`miRNA_HiSeq_gene`) are already log2-normalized, so the
  miRNA script skips Box-Cox and goes straight to NZV filter + standardize.
- cBioPortal molecular-profile and sample-list IDs vary by data release;
  discover the exact IDs with `client.get_molecular_profiles(study_id)` and
  `client.get_sample_lists(study_id)`. Xena dataset ids likewise vary by hub.
