# Metastatic vs non-metastatic NSCLC: expression and copy number

What genes distinguish metastatic from non-metastatic NSCLC — in mRNA expression
(script 17) and in copy number (script 18)? Both use the same set of clean,
stage-aware contrasts, and both reach the same conclusion.

Reproduce:

```bash
python scripts/17_metastasis_expression_diff.py --save-dir results  # expression
python scripts/18_metastasis_cnv_diff.py --save-dir results         # copy number (gene-level GISTIC)
python scripts/19_metastasis_segment_cnv.py --save-dir results      # copy number (segments: FGA + arm-level)
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

> **Subtype-specific check (and a cautionary tale).** M1 is 72 % LUAD while the
> true-stage-I control is LUSC-rich (49 % LUAD), and LUSC is intrinsically more
> proliferative — so a *pooled* median fold-change reverses sign (Simpson's
> paradox). A curated proliferation **panel** run **LUAD-only** looked promising
> (FOXM1/MKI67/TYMS up, 8/14 at panel-FDR). But that was selection bias (a
> co-regulated module tested against its own FDR): the **unbiased genome-wide
> LUAD-only** DE (`--studies luad`) is still **null — 0 genes at q<0.05, best
> q ≈ 0.36, no proliferation pathway**. The subtype confound is real, but removing
> it does not uncover a distant-met signature; the LUAD-only lean is nominal at
> best. Panel-FDR over hand-picked genes is not genome-wide evidence.

**Distant-metastatic capacity is not encoded in the bulk primary-tumour mRNA** of
TCGA NSCLC — robust to the reference *and to per-subtype analysis*. Consistent with
distant spread depending on rare subclones, the microenvironment, or later
stochastic events that bulk primary RNA does not capture.

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

**Present in both subtypes, stronger in LUAD (not LUAD-exclusive).** Running the
screen per subtype at nearly matched sample sizes (`--studies luad` / `lusc`;
LUAD 171 N+ vs 328 N0, LUSC 170 vs 306) shows the E2F/G2M program is significant
in *both* — E2F targets q≈5.6e-28 (LUAD) and q≈5.1e-22 (LUSC), G2M q≈6.5e-19 and
1.5e-10. But it is broader in LUAD (**1,461** differential genes vs **306** in
LUSC; MYC/cell-cycle clear FDR only in LUAD). The likely reason is a ceiling
effect: LUSC is intrinsically more proliferative at baseline, so N+ vs N0 has less
dynamic range there. So the nodal proliferation signature is a shared NSCLC
feature, LUAD-enriched — unlike the *distant-met* null, which held in both subtypes.

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

## Result 5 — segment-level copy number: the one arm-level exception (script 19)

The gene-level GISTIC screen (Result 4) tests *focal* deep events and misses broad
copy-number change. Script 19 uses the raw DNAcopy **segments** (continuous log2
ratios, hg19) for two summaries a per-gene binary test cannot give: **fraction
genome altered (FGA)** — one genomic-instability score per sample — and
**arm-level** length-weighted mean log2 per chromosome arm (39 autosomal arms),
each tested subtype-adjusted.

| Contrast | FGA (instability) | arm-level (q<0.05) |
|---|---|---|
| distant met vs true stage I | 0.314 vs 0.323, p=0.33 (ns) | **0 arms** |
| distant met vs all M0 | 0.314 vs 0.346, p=0.39 (ns) | **0 arms** |
| nodal N+ vs N0 | 0.355 vs 0.328, **p=0.009** | 3 gained: **17q** (q=9e-4), 7q, 2p |
| nodal N+ vs N0, stage II | 0.374 vs 0.351, p=0.17 (ns) | **17q only** (q=0.010) |

Two things emerge. First, **distant-metastatic primaries are not more unstable**
and carry no differential arm — FGA is flat (if anything lower) and 0 arms move,
against either reference. Second, the nodal instability *burden* is again **tumour
stage**: N+ tumours have higher FGA whole-cohort (p=0.009) but that collapses
within stage II (p=0.17), and two of the three arm gains (7q, 2p) drop out.

The exception is **17q gain**, the one copy-number feature associated with nodal
status **independent of stage** — significant whole-cohort (q=9e-4) *and* within
stage II (q=0.010, AUC≈0.62).

> **17q is LUSC-specific, not LUAD (`--studies` split).** Splitting the segment
> analysis by subtype: in **LUSC** 17q is the **top arm** (rank 1/39) both
> whole-cohort (z=4.61, q=2e-4) and — crucially — **within stage II** (z=3.64,
> q=0.011, AUC=0.67), so it survives stage matching. In **LUAD** 17q is *not*
> differential (q=0.46 whole, q=0.61 stage II, rank ~13/39) despite comparable
> sample size — a real subtype difference, not power. 17p (TP53's arm) stays
> copy-neutral in both. This is biologically coherent — 17q/3q gains are hallmark
> **squamous** copy-number events — and it means the two durable nodal signals
> *dissociate by subtype*: the proliferation *expression* program is LUAD-enriched,
> while the stage-independent 17q *copy-number* gain is LUSC-specific.

## 17q resolved to genes — arm dosage, not ERBB2 (script 20)

Drilling the 17q arm gain down to genes (`scripts/20_arm_gene_drilldown.py`, which
maps the script-17 DE tables to chromosome arms via MyGene cytobands) shows the DNA
gain is **functional as a coordinated cis-dosage effect** but pins down what it is
— and isn't:

- **Coordinated over-expression.** Of the 782 measured 17q genes, **71 % are shifted
  up** in N+ stage II (vs 56 % genome-wide; whole-cohort 58 % vs 40 %). The arm gain
  raises 17q transcription arm-wide. Pulling the top-hit genes' actual expression
  confirms it: all 14 are up in N+ by **log2FC ≈ 0.1–0.5** (~1.1–1.4×, the modest
  ratio expected from a single-copy dosage gain), and they are **co-expressed as a
  block** (mean pairwise Spearman r ≈ 0.50) — one coordinated program, not
  scattered noise.
- **ERBB2 / HER2 is *not* the driver.** Despite being the famous 17q oncogene, ERBB2
  is unchanged (z = 0.18 stage II, z = −0.02 whole-cohort; rank ~520/782 on the arm).
  Focal HER2 amplification is not what this is.
- **17p / TP53 is not co-lost.** 17p is copy-neutral (arm z = −0.18, p = 0.86), so
  this is a **selective 17q gain, not isochromosome 17q** — not a TP53-deletion event.
- **What's up is proliferation / biogenesis**: `TK1`, `BIRC5`, `KPNA2`, `SKA2`,
  mitochondrial ribosomal proteins (`MRPL12/45/58`, `MRPS7`), proteasome (`PSMB3`,
  `PSMC5`), `NME1`, `PYCR1` — the same program as the whole-cohort proliferation
  signature, here carried by arm dosage.
- **No single 17q gene clears FDR at matched stage** (best q ≈ 0.10) — confirming
  this is a genuine *arm-level* dosage effect, not one driver gene. It is why the
  arm-level CNV test (q ≈ 0.01) sees it while the per-gene screens do not.
- **The cis-dosage is LUSC-specific — and it separates dosage from program.**
  Re-running the drill-down per subtype: in **LUSC** the 17q arm is coordinately
  up — **66 %** of 17q genes up whole-cohort, **70 %** within stage II (vs 40 % /
  53 % genome-wide), matching the LUSC-specific DNA gain. In **LUAD** there is **no
  arm-level skew** (48 % ≈ 44 % baseline). Yet a few 17q genes (`TK1`, `BIRC5`,
  `KPNA2`) *are* up in LUAD N+ — because they are E2F/proliferation targets riding
  the LUAD *program*, not the arm. So the same gene (e.g. `TK1`) is elevated for
  two different reasons: **copy-number dosage in LUSC, proliferation program in
  LUAD.** ERBB2 remains a non-driver in both; the more plausible LUSC 17q-amplicon
  targets are chromatin/p53-pathway genes on the arm (`SMARCD2`, `SUZ12`, and
  `PPM1D` at 17q23, a recurrently amplified p53-suppressing oncogene).

## Conclusion

Across expression and copy number, the same result — with one small exception:

- **Distant metastasis** — no bulk-primary signature in *any* layer (expression,
  gene-level CNV, or segment FGA / arm-level), against a heterogeneous M0 reference
  *or* a clean true-stage-I indolent control. Distant-metastatic primaries are not
  even more genomically unstable.
- **Nodal metastasis** — a proliferation *expression* signature (LUAD-enriched,
  present in both subtypes) and higher CNV instability (FGA) whole-cohort, both
  **explained by tumour stage** (they collapse within stage II). The lone survivor
  is **17q gain**, associated with nodal status even at matched stage (q≈0.01) —
  modest, arm-level, and **LUSC-specific** (top arm in LUSC, absent in LUAD). The
  two durable nodal signals thus dissociate by subtype: LUAD proliferation
  (expression) vs LUSC 17q gain (copy number).

Net: in TCGA NSCLC, **metastatic status leaves essentially no clean molecular
imprint on the bulk primary tumour beyond markers of tumour stage / aggressiveness**
— proliferation (expression) and instability/deletion burden (CNV), both
stage-linked — the sole stage-independent hint being 17q gain in nodal disease.
Nothing marks *distant* metastasis. A genuine metastasis program, if one exists,
would need single-cell / subclonal resolution or actual metastasis samples
(MSK-MET / MET500). (This is the metastasis contrast; the *survival* screen in
`nsclc_survival_screen.md` is EMT-dominated — what predicts death differs from what
marks tumour stage.)

