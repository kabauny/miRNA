# NSCLC overall-survival screen: genes & pathways

Genome-wide screen for genes whose expression is associated with **overall
survival (OS)** in the TCGA NSCLC cohort (LUAD + LUSC PanCancer Atlas), followed
by pathway over-representation of the significant genes.

Reproduce (writes the full gene table and pathway CSVs to `results/`, which is
git-ignored):

```bash
python scripts/09_survival_screen.py --save-dir results
```

## Method

- **Cohort:** 979 NSCLC patients with OS data (LUAD + LUSC), 384 deaths.
- **Expression:** all 18,484 protein-coding genes with RSEM data (cBioPortal),
  log2-transformed.
- **Per-gene test:** univariate Cox **score test**, *stratified by subtype*
  (LUAD/LUSC) so associations are not driven by the subtype mix. `z > 0` means
  higher expression → higher hazard (worse survival). Validated to match
  lifelines' Cox Wald z. Benjamini–Hochberg FDR across all genes.
- **Hazard ratios:** the top hits were refit with a full Cox model
  (`expr_z + subtype`, lifelines); HR is per 1 SD of log2 expression.
- **Pathways:** hypergeometric over-representation of the FDR-significant genes
  against MSigDB **Hallmark** (50 sets) and **KEGG 2016** (293 sets), with the
  18,484 screened genes as the background universe.

## 1. Genes with a real impact on survival

**88 genes at FDR q < 0.05** — 67 where high expression predicts *worse* survival
("risk"), 21 where it predicts *better* survival ("protective").

### Strongest risk genes (high expression → worse OS)

| Gene | HR / SD | 95% CI | q | Note |
|---|---|---|---|---|
| CIDEC | 1.31 | 1.20–1.43 | <0.001 | |
| IGFBP1 | 1.29 | 1.17–1.42 | 0.002 | hypoxia-responsive |
| GPR78 | 1.24 | 1.14–1.35 | 0.002 | |
| ITGB1 | 1.29 | 1.17–1.42 | 0.003 | integrin, focal adhesion / EMT |
| FLNC | 1.29 | 1.16–1.42 | 0.003 | cytoskeleton / EMT |
| AKAP12 | 1.28 | 1.16–1.42 | 0.003 | |
| NTSR1 | 1.23 | 1.12–1.35 | 0.013 | |
| FSTL3 | 1.26 | 1.14–1.39 | 0.014 | TGF-β / activin |
| VEGFC | 1.26 | 1.14–1.39 | 0.014 | angiogenesis |
| LDHA | 1.26 | 1.13–1.41 | 0.017 | glycolysis (Warburg) |
| SLC16A3 | 1.27 | 1.14–1.43 | 0.019 | MCT4 lactate export |
| ANLN | 1.29 | 1.14–1.46 | 0.017 | mitotic / proliferation |

### Strongest protective genes (high expression → better OS)

| Gene | HR / SD | 95% CI | q |
|---|---|---|---|
| ANKRD65 | 0.77 | 0.70–0.85 | 0.002 |
| MYLIP | 0.79 | 0.72–0.87 | 0.003 |
| GLS2 | 0.79 | 0.72–0.87 | 0.007 |
| ZNF589 | 0.78 | 0.70–0.86 | 0.007 |
| GNMT | 0.78 | 0.70–0.88 | 0.013 |
| RPS6KA5 | 0.78 | 0.70–0.87 | 0.016 |
| HLF | 0.80 | 0.73–0.88 | 0.018 |

## 2. Pathways that are particularly impactful

Over-representation among the 67 **risk** genes (background = all screened genes):

| Pathway | Genes | Fold enrichment | q |
|---|---|---|---|
| **Hallmark: Epithelial–Mesenchymal Transition** | 8 | 11.1× | **0.0002** |
| KEGG: Focal adhesion | 6 | 8.4× | 0.014 |
| Hallmark: Mitotic spindle | 5 | 7.1× | 0.044 |
| KEGG: Proteoglycans in cancer | 5 | 7.0× | 0.044 |
| Hallmark: Apical junction | 5 | 7.0× | 0.044 |
| Hallmark: Hypoxia | 5 | 6.9× | 0.044 |
| Hallmark: TGF-β signaling | 3 | 15.3× | 0.048 |

**Interpretation.** Worse overall survival in NSCLC concentrates in a coherent
program: **epithelial–mesenchymal transition and cell–matrix adhesion / invasion**
(EMT, focal adhesion, apical junction, TGF-β, proteoglycans), **hypoxia and
glycolytic metabolism** (Hypoxia hallmark, plus individual hits LDHA, SLC16A3/MCT4,
PYGB, GPR78), and **proliferation** (mitotic spindle, ANLN). The protective genes
skew toward metabolic/differentiation programs (GLS2, GNMT, HLF) and did not reach
pathway-level significance on their own.

## 3. miRNA screen — not run here (data access)

The analogous *miRNA* survival screen could **not** be run in this environment.
A miRNA screen needs a TCGA miRNA matrix paired with OS; that pairing exists on
UCSC Xena (`tcga.xenahubs.net`), which is **egress-blocked** here. The miRNA
profiles cBioPortal hosts for lung (`lusc_tcga_pub`, CPTAC LUAD/LUSC) carry **no
OS fields**, so they cannot substitute. With egress to `tcga.xenahubs.net`
enabled, the same screen machinery (`mirna_tcga.screen` + `mirna_tcga.xena`)
applies directly to miRNA features.
