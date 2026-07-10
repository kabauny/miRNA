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

## 3. miRNA screen — run, and essentially null

miRNA expression (UCSC Xena, already log2) was matched to the same cBioPortal OS
outcomes by patient barcode and screened identically. Reproduce with:

```bash
python scripts/10_mirna_survival_screen.py --save-dir results
```

- **Cohort:** 768 NSCLC patients with both a tumour miRNA profile and OS
  (291 deaths); 567 mature miRNAs tested (detected in ≥80% of samples,
  remainder per-subtype median-imputed).
- **Result:** **no miRNA reaches FDR q < 0.05** — the best is hsa-miR-1468-5p
  (protective, p = 2e-4, **q = 0.11**). The strongest nominal candidates:

  | miRNA | accession | direction | z | p | q |
  |---|---|---|---|---|---|
  | hsa-miR-1468-5p | MIMAT0006789 | protective | −3.73 | 0.0002 | 0.11 |
  | hsa-miR-132-3p | MIMAT0000426 | risk | +3.47 | 0.0005 | 0.13 |
  | hsa-miR-29c-3p | MIMAT0000681 | protective | −3.40 | 0.0007 | 0.13 |
  | hsa-miR-148a-3p | MIMAT0000243 | protective | −3.25 | 0.0012 | 0.13 |
  | hsa-miR-31-3p | MIMAT0004504 | risk | +3.12 | 0.0018 | 0.13 |
  | hsa-let-7c-5p | MIMAT0000064 | protective | −3.11 | 0.0019 | 0.13 |

**Takeaway.** In this cohort, *individual* mRNA expression carries far stronger
univariate prognostic signal than *individual* miRNA expression: 88 genes clear
genome-wide FDR, but no single miRNA does. This is consistent with miRNAs acting
diffusely (each tunes many targets) rather than as strong stand-alone survival
markers.

But the *named* miRNA ranking is biologically coherent and reinforces the mRNA
story rather than contradicting it:

- The **miR-200 family** (miR-200a-3p/-5p, miR-200b-3p, miR-429) and **let-7c**
  (both arms) cluster among the **protective** miRNAs. The miR-200 family are the
  canonical *suppressors* of epithelial–mesenchymal transition — so higher miR-200
  → less EMT → better survival, the mirror image of the mRNA result where the
  **EMT program is the top risk pathway**.
- Both arms of **miR-31** (5p and 3p) rank on the **risk** side, an internal
  consistency check (a single hairpin agreeing with itself) that a pure-noise
  ranking would not produce.

So the honest summary is: no single miRNA is a stand-alone survival marker at
genome-wide FDR here, but the miRNA signal points at the **same EMT axis** the
mRNA screen pins down with significance.

> Note: this screen needs egress to `tcga.xenahubs.net` (the TCGA miRNA + OS
> pairing). That host is reachable; the matrices are served gzipped at
> `<dataset>.gz` (`mirna_tcga.xena` requests that object).
