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
```

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
mutations, clinical) and `https://*.xenahubs.net` (miRNA). In a sandboxed/remote
environment you may need to **add these hosts to the network egress allowlist**
— otherwise requests return `403 Host not in allowlist`. See
https://code.claude.com/docs/en/claude-code-on-the-web.

## Data sources

- **mRNA / mutations / clinical** → cBioPortal (`mirna_tcga.cbioportal`).
- **miRNA expression** → UCSC Xena TCGA hub (`mirna_tcga.xena`), because
  cBioPortal does not carry miRNA for the PanCancer Atlas studies. Both Xena
  hubs used here (`tcga`, `gdc`) serve genuine TCGA data.

## Caveats

- Xena miRNA matrices (`miRNA_HiSeq_gene`) are already log2-normalized, so the
  miRNA script skips Box-Cox and goes straight to NZV filter + standardize.
- cBioPortal molecular-profile and sample-list IDs vary by data release;
  discover the exact IDs with `client.get_molecular_profiles(study_id)` and
  `client.get_sample_lists(study_id)`. Xena dataset ids likewise vary by hub.
