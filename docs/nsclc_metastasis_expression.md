# Differential expression: metastatic vs non-metastatic NSCLC

What genes' mRNA expression differs most between metastatic and non-metastatic
NSCLC? The deletion screens ask what copy-number *loss* distinguishes metastasis;
this asks the broader transcriptional question.

Reproduce:

```bash
python scripts/17_metastasis_expression_diff.py --save-dir results
```

## Method

- Genome-wide mRNA expression (18,484 protein-coding genes, log2), TCGA LUAD +
  LUSC, streamed from cBioPortal.
- **Subtype-stratified van Elteren rank-sum** (`associate.ranksum_screen`) of every
  gene vs a binary metastasis endpoint; `z > 0` = higher in the metastatic group.
  Stratifying by LUAD/LUSC is essential — the metastatic groups are LUAD-heavy, so
  an unadjusted test would rediscover subtype differences.
- **Four contrasts, differing only in the comparison groups:**
  1. distant metastasis (Biotab-enriched, 125) vs **all M0** (670);
  2. distant metastasis (122) vs a clean **"true stage I"** control (206): stage I
     / N0 / M0, followed ≥2 yr with no distant *or* locoregional recurrence
     (`biotab.true_stage_i_vs_distant`);
  3. **nodal** metastasis N+ (341) vs N0 (634);
  4. nodal N+ vs N0 **restricted to stage II** (holds stage constant).

Caveat throughout: TCGA samples are resected **primary** tumours, so this
contrasts the primary-tumour transcriptome of patients with vs without metastatic
disease — not metastases themselves.

## Result 1 — distant metastasis has no bulk-primary signature

Neither reference yields a single gene at FDR q < 0.05:

| Distant-met contrast | cases vs controls | genes at q<0.05 | best q |
|---|---|---|---|
| vs all M0 | 125 vs 670 | **0** | 0.44 |
| vs true stage I (clean) | 122 vs 206 | **0** | 0.41 |

The first M0 reference is heterogeneous (spans stage I–III, ~36% node-positive),
so a null there could be dilution. But cleaning the reference to *truly indolent*
stage I tumours (stage I / N0 / M0, ≥2 yr recurrence-free) **does not** rescue a
signal — still zero genes, best q ≈ 0.41, top AUCs ~0.55–0.65 (a faint
proliferation echo: `TYMS`, `KNL1`, `ARHGAP11B`, but nowhere near significant).

**Distant-metastatic capacity is not encoded in the bulk primary-tumour mRNA** of
TCGA NSCLC — robust to the reference. Consistent with distant spread depending on
rare subclones, the microenvironment, or later stochastic events that bulk primary
RNA does not capture.

## Result 2 — nodal metastasis is a strong proliferation signal…

N+ vs N0 is the opposite: **2,186 genes at q < 0.05** (779 up, 1,407 down). The
up-in-N+ genes are textbook cell cycle — `CCNA2`, `TK1`, `KPNA2`, `MAD2L1`,
`CEP55`, `BIRC5`, `EXO1`, `PLK1` — and the pathway enrichment is overwhelming:

| Pathway (up in N+) | overlap | fold | q |
|---|---|---|---|
| Hallmark E2F targets | 89 | 10.7 | 1e-66 |
| Hallmark G2M checkpoint | 75 | 9.3 | 1e-50 |
| Hallmark MYC targets V1 | 60 | 7.3 | 2e-33 |
| Hallmark mTORC1 signaling | 46 | 5.6 | 5e-20 |
| KEGG cell cycle | 34 | 6.5 | 6e-17 |
| KEGG DNA replication | 17 | 11.2 | 7e-13 |

The 1,407 down-in-N+ genes carry no pathway enrichment (diffuse; likely
stromal/immune dilution).

## Result 3 — …but that signal is really tumour stage, not nodal spread

The catch: N0 tumours are ~78% stage I, while N+ tumours are stage II–III. So
"N+ vs N0" across the whole cohort is confounded with stage, and the proliferation
signal could simply be "higher-stage tumours proliferate more." The clean test
holds stage constant — N+ vs N0 **within stage II** (175 vs 105, a balanced split;
stage I is ~all N0, stage III ~all N+):

| Nodal contrast | N+ vs N0 | genes at q<0.05 | proliferation pathways |
|---|---|---|---|
| whole cohort | 341 vs 634 | 2,186 | E2F/G2M/MYC, q down to 1e-66 |
| **stage II only** | 175 vs 105 | **2** | **none (0 sets at q<0.05)** |

Within stage II the signature **collapses**: just 2 genes clear FDR (`ANKRD42`,
`LYRM1` — not cell-cycle genes) and **no pathway** is enriched. The proliferation
genes still lean up (`ASF1B`, `CENPP`, q ≈ 0.10) but as a whisper, not a signal.
So the whole-cohort nodal proliferation signature was **tumour stage, not nodal
spread**: hold stage constant and it disappears.

## Conclusion

Across four contrasts the story is consistent:

- **Distant metastasis** — no bulk-primary-tumour expression signature, against a
  heterogeneous M0 reference *or* a clean true-stage-I indolent control (0 genes,
  best q ≈ 0.41).
- **Nodal metastasis** — a large proliferation signature whole-cohort, but it is
  **explained by tumour stage**: it collapses to nothing within stage II.

Net: in TCGA NSCLC, **metastatic status leaves no transcriptional imprint on the
bulk primary tumour beyond what tumour stage / proliferation already captures.**
The single expression axis that separates "aggressive" from "indolent" here is
proliferation, and it tracks stage. A genuine metastasis-specific program, if one
exists, would need single-cell / subclonal resolution or actual metastasis samples
(MSK-MET / MET500) — bulk primary RNA does not carry it. (Note this is the
metastasis contrast; the *survival* screen in `nsclc_survival_screen.md` is
EMT-dominated — what predicts death differs from what marks tumour stage.)

