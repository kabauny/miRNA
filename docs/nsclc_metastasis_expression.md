# Metastatic vs non-metastatic NSCLC: expression and copy number

What genes distinguish metastatic from non-metastatic NSCLC — in mRNA expression
(script 17) and in copy number (script 18)? Both use the same set of clean,
stage-aware contrasts, and both reach the same conclusion.

Reproduce:

```bash
python scripts/17_metastasis_expression_diff.py --save-dir results  # expression
python scripts/18_metastasis_cnv_diff.py --save-dir results         # copy number
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

## Result 4 — copy number tells the same story (script 18)

The same question for **copy number** — are genes differentially **deleted or
amplified**? — using a two-sided, subtype-adjusted CMH
(`associate.cmh_two_sided_screen`) over both deep deletions (HOMDEL) and high
amplifications (AMP), across the same contrasts:

| Contrast | deep deletions (q<0.05) | amplifications |
|---|---|---|
| distant met vs true stage I | **0** (best q≈0.054) | 0 |
| distant met vs all M0 | 106 *enriched* | 0 |
| nodal N+ vs N0 | **0** | 0 |
| nodal N+ vs N0, stage II | **0** | 0 |

The only contrast with hits is the **dirty** one: distant-met vs all-M0 gives 106
genes with *more* deep deletions in metastatic tumours (0 depleted, 0 amplified).
But they are low-frequency (~0.3–1.3% in M0 → 3–5% in metastatic), **scattered
genome-wide** (`NOTCH2`, `RHOC`, `PHGDH`, `OPCML`, 1q `NBPF`…), and all
one-directional — the signature of **general genomic instability / higher deletion
burden in more-aggressive tumours**, not specific metastasis loci. Against the
clean true-stage-I control it **evaporates to 0** (best q≈0.054). Nodal status
shows *no* CNV signal at all — even whole-cohort — and amplifications show nothing
anywhere. The 8p23 "spared deletion" (scripts 14–16) reappears here at q≈0.054 but
is a borderline minority effect, swamped by the opposite general "more deletions"
trend.

## Conclusion

Across expression and copy number, the same result:

- **Distant metastasis** — no bulk-primary signature in *either* omic, against a
  heterogeneous M0 reference *or* a clean true-stage-I indolent control (expression
  0 genes best q≈0.41; CNV 0 genes best q≈0.054).
- **Nodal metastasis** — a large proliferation *expression* signature whole-cohort
  that is **explained by tumour stage** (collapses within stage II); *no* CNV
  signal at all.

Net: in TCGA NSCLC, **metastatic status leaves no clean molecular imprint on the
bulk primary tumour beyond markers of tumour stage / aggressiveness** —
proliferation (expression) and diffuse deletion burden (CNV), both stage-linked,
neither a metastasis-specific program or locus. A genuine metastasis program, if
one exists, would need single-cell / subclonal resolution or actual metastasis
samples (MSK-MET / MET500). (This is the metastasis contrast; the *survival* screen
in `nsclc_survival_screen.md` is EMT-dominated — what predicts death differs from
what marks tumour stage.)

