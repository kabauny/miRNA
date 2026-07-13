# NSCLC multi-omic screen: metastasis, copy-number, mutations, and a constellation model

Extends the [survival screen](nsclc_survival_screen.md) with three DNA/clinical
layers and asks: **which alterations predict survival and metastasis, and does
combining layers beat any single one?**

Reproduce:

```bash
python scripts/11_cnv_mutation_screen.py --save-dir results   # per-alteration screens
python scripts/12_constellation_model.py --save-dir results   # cross-validated models
```

## Endpoints

| Endpoint | Positive / total | Notes |
|---|---|---|
| Overall survival | 384 deaths / 979 | well-powered |
| **Nodal metastasis** (N+ vs N0) | 344 / 982 | best-powered metastasis endpoint |
| **Distant metastasis** (M1 / stage IV) | 34 / 774 | rare — under-powered, exploratory only |

Distant M1 is what was asked for; because only ~34 patients are M1 in this
resected cohort, it is reported but treated as exploratory, with nodal N+ as the
powered anchor.

## 1. Deep-deletion screen (GISTIC −2)

398 genes are homozygously deleted in ≥2% of tumours. **No deletion survives FDR**
for any endpoint, but the **9p21 locus** (co-deleted as one block) is the
consistent nominal signal for **nodal metastasis**:

| Gene (9p21) | deleted n | OR (N+) | p |
|---|---|---|---|
| CDKN2B | 206 | 1.54 | 0.008 |
| CDKN2A | 209 | 1.53 | 0.009 |
| DMRTA1 | 143 | 1.58 | 0.013 |
| MTAP | 152 | 1.50 | 0.026 |

So CDKN2A/B (9p21) deletion trends toward more nodal spread (OR ~1.5), consistent
with its known role in aggressive NSCLC — but it is not significant after
multiple-testing correction, and shows no OS association.

## 2. Mutation screen

1,731 genes are mutated in ≥3% of tumours. Again **nothing clears FDR**, but the
top nominal hits are biologically sensible:

| Gene | endpoint | effect | p |
|---|---|---|---|
| **STK11** | OS | z = +3.06 (worse) | 0.002 |
| **EGFR** | distant met | OR = 3.95 | 0.004 |
| **KEAP1** | distant met | OR = 2.94 | 0.017 |
| **KDM6A** | nodal met | OR = 2.89 | 0.003 |

STK11-mutant → worse survival and EGFR/KEAP1-mutant → distant metastasis are all
established NSCLC biology; they simply don't reach genome-wide significance at
this cohort size / mutation frequency.

## 3. The constellation model

Patient-level multi-omic features (expression signature selected in-fold +
recurrent deep deletions + recurrent mutations + subtype) were combined and
scored by **cross-validation** (no feature leakage): C-index for survival,
ROC-AUC for metastasis. Nested feature sets show each layer's incremental value.

| Model | OS (C-index) | Nodal met (AUC) | Distant met (AUC) |
|---|---|---|---|
| subtype only | 0.525 | 0.488 | 0.66 ± 0.06 * |
| **expression** | **0.556** | **0.585 ± 0.04** | 0.57 ± 0.11 |
| + deletions | 0.551 | 0.583 ± 0.04 | 0.58 ± 0.11 |
| + mutations | 0.504 | 0.553 ± 0.04 | 0.58 ± 0.11 |
| all + subtype | 0.501 | 0.554 ± 0.04 | 0.57 ± 0.11 |

\* distant-met AUCs are unstable (only 33 events, ±0.11); the subtype-only value
is an artifact of the tiny positive set, not a real signal.

**What predicts survival and metastasis — the honest answer.**

- The **transcriptomic program is the carrier of signal.** The expression
  signature is the best single predictor of both overall survival (C-index 0.556)
  and nodal metastasis (AUC 0.585). Absolute performance is modest — these are
  hard endpoints — but expression clearly beats subtype and DNA-only models.
- **Adding copy-number deletions gives no cross-validated gain, and adding the
  sparse mutation matrix actively hurts** (C-index 0.556 → 0.50; AUC 0.585 →
  0.55). At the whole-cohort level the DNA-alteration layers do not add
  predictive value on top of expression — the individually-suggestive hits
  (9p21, STK11, EGFR/KEAP1) are too sparse and weak to help a joint linear model.
- **Distant metastasis cannot be modelled** reliably here (33 events).

So the "constellation" that predicts outcome in NSCLC is essentially the
**EMT / hypoxia / proliferation expression program** identified in the
[survival screen](nsclc_survival_screen.md); deep deletions and mutations
sharpen the mechanistic picture (9p21 loss, STK11/EGFR/KEAP1 mutations) but do
not improve cross-validated prediction at this sample size.

## 4. Sparse, driver-focused model + a fifth (miRNA) layer

`scripts/13_constellation_sparse.py` addresses the two weaknesses above:
mutations are restricted to a **curated NSCLC driver panel** and deletions to
recurrently-deleted **tumour suppressors** (instead of every recurrent gene);
models are **L1-penalized** (embedded feature selection); and a **miRNA layer**
(curated panel incl. the miR-200 / let-7 EMT-suppressor families) is added.

| Endpoint | driver + L1 (DNA + expr) | + miRNA |
|---|---|---|
| OS (C-index) | 0.552 | **0.572** |
| Nodal met (AUC) | **0.592** | 0.583 |
| Distant met (AUC) | 0.52 ± 0.09 | 0.55 ± 0.12 |

Two things improve over the dense L2 model in §3:

- **Focusing the DNA layers + L1 stops them hurting.** Nodal-met AUC is 0.592
  (vs 0.55 when all recurrent mutations were dumped into an L2 model), i.e. the
  DNA features now at least do no harm and slightly help.
- **miRNA adds a small but real gain to survival** (C-index 0.552 → 0.572);
  on nodal metastasis it is roughly neutral.

### The nodal-metastasis constellation (L1-selected, full cohort)

The strongest non-zero features (coefficients in-sample, for interpretation):

| Feature | dir. | reading |
|---|---|---|
| expression signature | ↑ → more met | EMT / hypoxia program dominates |
| PIK3CA mut, EGFR mut, KRAS mut | ↑ → more met | oncogenic driver activation |
| CDKN2B deletion (9p21) | ↑ → more met | echoes the §1 9p21 → nodal trend |
| miR-155-5p | ↑ → more met | oncomiR, expected direction |
| miR-200b-3p, miR-29c-3p | ↓ → more met | EMT-suppressor loss → metastasis |
| subtype (LUSC) | ↓ | LUAD more node-positive than LUSC here |

The signal is still **carried by the expression signature**, but with driver
focus + L1 the model assembles a coherent multi-omic constellation around it:
oncogenic driver mutations (PIK3CA/EGFR/KRAS), 9p21 loss, and the miR-200 /
miR-155 EMT axis all contribute in biologically sensible directions. Cross-
validated gains over expression alone remain modest — the honest conclusion is
that at TCGA-cohort size the transcriptome is the dominant predictor and the
other layers add mechanism plus a small amount of orthogonal signal.
