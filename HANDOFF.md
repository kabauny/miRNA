# Session handoff — TCGA NSCLC multi-omic analysis

Continuation notes for a fresh session. Everything described here is already
committed and pushed to **`master`**. Read this first, then the two result docs
in `docs/`.

## Immediate pending task (start here)

The user is providing a **`TCGA` folder** (likely on the connected **Google
Drive**, MCP server `Google_Drive`) containing richer TCGA clinical data — the
**BCR Biotab clinical supplements** (what `TCGAbiolinks::GDCdownload(... data.type
= "Clinical Supplement", data.format = "BCR Biotab")` fetches).

**Why it matters:** our whole metastasis analysis defines "metastatic" as
*pathologic* M1 / stage IV, which **undercounts** — only ~33 M1 patients in lung.
The Biotab has the fuller annotation the user (correctly) wants:
`new_neoplasm_event_type = 'Distant Metastasis'`, clinical M stage (cM), metastatic
sites, days-to-new-tumor-event. Using those would substantially enlarge the
metastatic set and let the metastasis screens run with real power.

**What to do with it:**
1. Search Drive: `mcp__Google_Drive__search_files` (query `title contains 'TCGA'`
   or `fullText contains 'new_tumor_event'`); download with
   `mcp__Google_Drive__download_file_content` / `read_file_content`. (A first
   `search_files` attempt in the prior session was declined by the user — just
   ask them to confirm access or point to the exact file.)
2. Add a loader (e.g. `src/mirna_tcga/biotab.py`) that parses the biotab TSVs and
   extracts a distant-metastasis label per patient barcode.
3. Extend `src/mirna_tcga/endpoints.py` with a `distant_metastasis_biotab()` (or
   augment `distant_metastasis`) that unions pathologic M1 + clinical M1 +
   `new_neoplasm_event_type == 'Distant Metastasis'`.
4. Re-run `scripts/14` and `scripts/15` with the richer endpoint. `endpoints.py`
   is the single integration point — the screens take any 0/1 patient Series.

## Environment gotchas (these cost real time — don't relearn them)

- **Network egress is very restricted.** Reachable: `www.cbioportal.org/api`,
  `raw.githubusercontent.com`, PyPI. **Blocked**: GDC (`api.gdc.cancer.gov`),
  Broad Firehose (`gdac.broadinstitute.org`), MSigDB (`data.broadinstitute.org`),
  Enrichr (`maayanlab.cloud`), g:Profiler, and the Xena *query* hosts
  (`xenabrowser.net`, bare `xenahubs.net`). So **TCGAbiolinks/GDCdownload cannot
  run here** — that's why the user's file must be supplied manually.
- **UCSC Xena download DOES work** at `https://tcga.xenahubs.net/download/<dataset>.gz`
  — the `.gz` suffix is mandatory (bare key → S3 AccessDenied). `mirna_tcga.xena`
  handles this; it's how miRNA expression is fetched.
- **lifelines will NOT install into the system Python** (Debian setuptools breaks
  the `autograd-gamma` build). Use a clean venv:
  ```bash
  python -m venv .venv && . .venv/bin/activate
  pip install -e '.[dev]' && pip install lifelines
  ```
  The venv is in the (ephemeral) scratchpad in the prior session — **recreate it**.
  Offline tests (`pytest`) do NOT need lifelines; the analysis scripts do.
- **cBioPortal has no gene-level methylation** for PanCancer studies (fetch
  returns 0 records), so the 3rd inactivation mode (methylation) is unavailable
  via the API. Deletion + truncating mutation are the two usable modes.
- **cBioPortal expression uses SUMMARY projection** which returns `entrezGeneId`
  but not `hugoGeneSymbol` (already handled in `cbioportal.molecular_matrix`).
- Sparse fetches are cheap: `discrete_cna_events` (HOMDEL) and `mutation_events`.
  Full expression is ~7 GB / ~4 min genome-wide (streamed in `layers.stream_expression`).

## Git / workflow conventions

- **Work on `master`** (the user explicitly asked to consolidate here). Push with
  `git push -u origin master`.
- Gitignored: `results/`, `data/`, `*.csv`, `*.tsv`, `genesets_cache/`, `*.gmt` —
  derived data is NOT committed; scripts regenerate it. Results are written up as
  prose+tables in `docs/`.
- Commit trailers used: `Co-Authored-By:` and `Claude-Session:` (see git log).

## Package layout (`src/mirna_tcga/`)

| Module | Purpose |
|---|---|
| `cbioportal.py` | REST client: expression/CNA matrices, `discrete_cna_events` (HOMDEL), `mutation_events`, `sample_list_ids`, clinical |
| `xena.py` | UCSC Xena miRNA loader (`.gz` fetch) |
| `cohorts.py` | NSCLC = LUAD+LUSC; combined clinical/expression |
| `endpoints.py` | binary metastasis endpoints: `distant_metastasis` (M1/stage IV), `nodal_metastasis` (N+) ← **extend for biotab** |
| `layers.py` | shared loaders: `stream_expression`, `deletion_matrix`, `mutation_matrix` (`truncating_only`), `protein_coding_map` |
| `screen.py` | fast vectorized subtype-stratified Cox score test + BH FDR (validated vs lifelines) |
| `associate.py` | `ranksum_screen` (van Elteren), `fisher_screen`, `cmh_depletion_screen` (subtype-adjusted depletion) |
| `enrich.py` | GMT parsing, `load_gene_sets`, hypergeometric over-representation |
| `survival.py` | lifelines KM / Cox / log-rank helpers |
| `mirbase.py` | MIMAT accession → hsa-miR name (miRBase mature.fa) |
| `panels.py` | `LUNG_MARKER_PANEL`, `NSCLC_DRIVERS`, `DELETION_TSGS`, `CHR8P23_BLOCK` |
| `classify.py`, `preprocess.py`, `integrate.py`, `idmap.py`, `config.py` | earlier pipeline pieces |

## Scripts (`scripts/`) — all runnable, `--save-dir results`

- `09_survival_screen.py` — genome-wide mRNA OS screen + pathway ORA
- `10_mirna_survival_screen.py` — genome-wide miRNA OS screen (Xena)
- `11_cnv_mutation_screen.py` — deletion & mutation vs OS / distant / nodal
- `12_constellation_model.py` — CV multi-omic OS (C-index) + metastasis (AUC), dense L2
- `13_constellation_sparse.py` — driver-focused L1 model + miRNA layer; reports the selected constellation
- `14_metastasis_spared_deletions.py` — negative-selection: genes spared from deletion in metastasis (CMH + pathway burden)
- `15_pancancer_spared_deletions.py` — does §14 generalize across cancers? `--protection` = deletion OR truncating mutation

## Findings so far (see `docs/` for full tables)

`docs/nsclc_survival_screen.md`, `docs/nsclc_multiomics_metastasis.md`,
`docs/nsclc_metastasis_spared_deletions.md`.

1. **mRNA & survival**: 88 genes at FDR<0.05; risk genes enrich for **EMT**
   (q=2e-4), focal adhesion, hypoxia, mitotic spindle, TGF-β. Expression is the
   dominant prognostic layer.
2. **miRNA & survival**: no miRNA clears FDR, but the **miR-200 family** (EMT
   suppressors) + let-7c are the top protective candidates — mirrors the mRNA EMT
   result.
3. **CNV/mutation**: nothing clears FDR individually; nominal STK11→worse OS,
   EGFR/KEAP1→distant met, 9p21(CDKN2A/B)→nodal.
4. **Constellation (CV)**: expression signature carries the signal (OS C-index
   ~0.55–0.57; nodal-met AUC ~0.59). Dense mutation features hurt; **driver-focus
   + L1** fixes that; miRNA adds a small OS gain. Constellation = EMT signature +
   PIK3CA/EGFR/KRAS + 9p21 del + miR-200/miR-155.
5. **Spared-deletion hypothesis** (user's idea: metastasis needs intact
   machinery). Lung: **8p23** block nominally spared in metastatic (M1 0/33; N+ OR
   ~0.5, p~0.03) but **q~0.22, underpowered**. **Pan-cancer test REFUTES it** —
   8p23 does not replicate and *reverses* (pooled z=+2.94, p=0.003, deleted MORE
   in metastatic; only LUAD shows sparing). Likely an underpowered false signal.
   The **method is sound**; the biotab clinical (pending task) is the way to give
   it real power in lung.

## Suggested next steps

1. **Ingest the user's TCGA biotab folder** → richer distant-metastasis label →
   re-run `14`/`15` (the headline follow-up).
2. If a metastatic-enriched cohort becomes reachable (MSK-MET / MET500), rerun the
   spared-deletion + constellation analyses there.
3. Optional: nested-CV hyperparameter search on the L1 penalties in `13`.

## Test status

57 offline tests pass (`pytest`), `ruff check .` clean. Tests mock the HTTP layer
so they need no network and no lifelines.
