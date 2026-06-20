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
| Data access | manual `MANIFEST.txt` parsing of local files | **cBioPortal REST API** (`mirna_tcga.cbioportal`) |
| Classification | `pamr` | `sklearn.neighbors.NearestCentroid` + CV threshold |
| ID conversion | `biomaRt` | `mygene` / `pybiomart` (`mirna_tcga.idmap`) |
| Normalization | `MASS::boxcox`, custom NZV | `scipy.stats.boxcox`, vectorized filters |
| Survival | (none) | `lifelines` (`mirna_tcga.survival`) |
| Reproducibility | hard-coded `/Users/...` paths | `config.yaml`, packaged, tested |

## Layout

```
src/mirna_tcga/      # the package
  config.py          # load config.yaml
  cbioportal.py      # REST API client (expression, mutations, clinical)
  preprocess.py      # variance filter, Box-Cox, standardize
  classify.py        # PAM (nearest shrunken centroids) + CV + signature
  idmap.py           # gene/miRNA id conversion (optional deps)
  survival.py        # Kaplan-Meier / Cox PH (optional deps)
  panels.py          # small demo gene panel
scripts/             # runnable pipeline examples
  02_subtype_signature.py   # LUAD vs LUSC signature
  03_mutation_analysis.py   # TP53-mutant vs wild-type signature
  04_survival.py            # overall survival by stage
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
python scripts/02_subtype_signature.py            # LUAD vs LUSC
python scripts/03_mutation_analysis.py --gene TP53 --study luad
python scripts/04_survival.py --study luad        # needs [survival] extra
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

The pipeline calls `https://www.cbioportal.org/api`. In a sandboxed/remote
environment you may need to **add `www.cbioportal.org` to the network egress
allowlist** — otherwise requests return `403 Host not in allowlist`. See
https://code.claude.com/docs/en/claude-code-on-the-web.

## Caveats

- **miRNA expression** is not available for most cBioPortal PanCancer studies;
  those studies cover mRNA, mutations, CNA, and clinical data. For miRNA-level
  analysis, point the data layer at a Firehose-legacy study that has a miRNA
  profile, or swap in UCSC Xena / the GDC API.
- Molecular-profile and sample-list IDs vary by data release; discover the
  exact IDs with `client.get_molecular_profiles(study_id)` and
  `client.get_sample_lists(study_id)`.
