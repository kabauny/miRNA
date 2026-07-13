# What is never deleted in metastatic NSCLC? A negative-selection screen

**Hypothesis (user).** Metastasis requires *intact* machinery. If deep deletion
of a gene abolishes a tumour's ability to metastasize, metastatic tumours should
almost never carry that deletion — even where the gene is readily deleted in
non-metastatic tumours. So genes that are **deletable but spared in stage IV**
are candidate metastasis-required machinery.

Reproduce:

```bash
python scripts/14_metastasis_spared_deletions.py --save-dir results
```

## Method

- **Deep deletions** (GISTIC −2) genome-wide, 998 NSCLC patients.
- **Endpoints:** distant metastasis (M1 / stage IV, 33 vs 739) as asked, and
  nodal metastasis (N+, 343 vs 637) as a better-powered check.
- **Only deletable genes are testable** — a gene must be deleted in ≥4% of the
  non-metastatic group to be informatively "spared" (a gene never deleted
  anywhere is trivially absent in M1).
- **Subtype adjustment is essential**: M1 is 79% LUAD vs M0 55% LUSC, so an
  unadjusted test would just rediscover LUAD-vs-LUSC deletion differences. Every
  test is a **Cochran–Mantel–Haenszel** depletion test stratified by subtype
  (`mirna_tcga.associate.cmh_depletion_screen`), robust to the zero cells that
  "never deleted in M1" produces.
- Two levels: **gene** (per-gene depletion) and **pathway** (is a gene set's
  per-patient *deletion burden* lower in metastatic tumours? — rank-sum).

## Result: metastatic tumours retain chromosome 8p23

Nothing survives genome-wide FDR (33 M1 patients is too few), **but a single,
specific, biologically coherent signal recurs at both endpoints**: the genes that
are deletable yet spared in metastatic disease **cluster on chromosome 8p23**.

| Gene | cytoband | del freq M0 | del in metastatic | MH OR |
|---|---|---|---|---|
| ARHGEF10 | 8p23.3 | 7.7% | 0 / 33 (M1) | 0.0 |
| MSRA | 8p23.1 | 6.1% | ↓ (N+) | 0.56 |
| PINX1 | 8p23.1 | 6.3% | ↓ (N+) | 0.54 |
| SOX7 | 8p23.1 | 6.3% | ↓ (N+) | 0.54 |
| XKR6 | 8p23.1 | 6.6% | 0.51 (N+) | 0.51 |
| DLGAP2, MFHAS1, FBXO25, CLN8, C8orf74, … | 8p23 | 6–8% | ↓ | ~0.5 |

These genes form one **co-deleted 8p23 block**. In node-positive (N+) tumours the
block is deep-deleted about **half as often** as in N0 (MH OR ≈ 0.5,
nominal p ≈ 0.03–0.05; **q ≈ 0.22**, not significant). In distant-M1 tumours the
block is deleted in **0 / 33** patients despite ~7% deletion in M0 (nominal
p ≈ 0.12). Both endpoints point the same way: **metastatic NSCLC preferentially
retains 8p23** — consistent with the hypothesis that intact 8p23 carries
machinery permissive for (or required by) metastasis.

### Pathway level

The most-spared gene sets (lower deletion burden in metastatic tumours), while
also not FDR-significant, are consistent across endpoints and mechanistically
suggestive:

- **Base-excision (DNA) repair** — lower deletion burden in both M1 and N+
  (nominal p ≈ 0.13 / 0.17): metastatic cells appear to keep DNA-repair capacity.
- **Antigen processing & presentation** and **fatty-acid biosynthesis** — lower
  deletion burden in M1.

## Honest interpretation

- The hypothesis is directly testable and the screen implements it correctly
  (deletable-gene filter + subtype-adjusted depletion at gene and pathway level).
- **The signal is nominal, not statistically significant** — TCGA has only ~33
  distant-metastatic patients with copy-number data. The screen is
  hypothesis-generating.
- The clearest candidate is **8p23**: metastatic tumours retain it ~2× more often
  than expected. 8p loss is a common early NSCLC event, so its *depletion* in
  metastatic tumours is a genuine negative-selection signal worth following.
- As you noted, the same machinery could be silenced by **methylation or
  mutation** instead — absence of deletion is one axis of evidence, not proof. A
  natural extension is to ask whether 8p23 genes are also protected from
  damaging mutation / promoter methylation in metastatic tumours.
- **To power this properly**, repeat in a metastasis-enriched cohort (e.g.
  MSK-MET / MET500), which contains far more true metastatic samples than TCGA's
  resected primaries. The machinery here (`endpoints`, `cmh_depletion_screen`,
  `pathway_burden`) transfers directly.

## Does it generalize? Pan-cancer test (it does not)

`scripts/15_pancancer_spared_deletions.py` runs the same negative-selection
screen across the TCGA cohorts that actually contain metastatic patients —
bladder (BLCA, **135** M1), colorectal (COADREAD, 84), kidney (KIRC, 84), stomach
(STAD, 45), plus lung and melanoma — pooled and stratified by cancer type
(403 metastatic / 2,815 patients). These cohorts have far more metastatic disease
than lung's ~33, so this is the properly-powered test.

**The 8p23-sparing signal does not replicate — it reverses.** Per-cancer
8p23-block deep-deletion rate (M0 → M1):

| Cancer | M1 n | M0 → M1 8p23 deletion |
|---|---|---|
| LUAD | 26 | 8.4% → **3.8%** (spared — the original signal) |
| COADREAD | 84 | 8.4% → **19.0%** (more deleted) |
| BLCA | 135 | 7.5% → **11.9%** (more) |
| STAD | 45 | 3.8% → **11.1%** (more) |
| SKCM | 22 | 0.9% → **4.5%** (more) |
| KIRC | 84 | 2.0% → 1.2% (~equal) |

Pooled (stratified by cancer type), 8p23 is deep-deleted **more** in metastatic
tumours, not less: block-level z = **+2.94, p = 0.003**. Adding truncating
mutations (protection mode) gives no sparing either (z = 0.91, p = 0.36). At the
gene level, no gene is FDR-significantly spared pan-cancer; the top nominal hits
are common fragile sites (MACROD2, WWOX, FHIT, PTPRD, LRP1B, IMMP2L) deleted at
similar rates regardless of metastasis.

**Conclusion.** Only lung adenocarcinoma shows 8p23 sparing, and only nominally.
Across cancers 8p23 is, if anything, *more* deleted in metastatic disease. The
honest reading is that the lung 8p23 result was an **underpowered, non-replicating
signal** — exactly the outcome the pan-cancer check is designed to catch. The
negative-selection *method* is sound and worth applying, but 8p23 is not a
general metastasis-required locus.

### A caveat on the metastasis label (data access)

This all defines "metastatic" as **pathologic M1 / stage IV**, which *undercounts*
— it misses patients who develop distant metastasis over follow-up. The fuller
annotation (clinical M stage, `new_neoplasm_event_type = Distant Metastasis`,
metastatic sites) lives in the TCGA **BCR Biotab** clinical supplements on GDC
(what `TCGAbiolinks::GDCdownload` fetches). In this environment GDC, the Broad
Firehose mirror, and cBioPortal's legacy clinical-M field are all unavailable /
empty, so cBioPortal's pathologic M is the ceiling. Supplying the Biotab clinical
(or enabling GDC egress) would enlarge the metastatic set — but note the
pan-cancer cohorts above (BLCA/COADREAD/KIRC) already carry enough true M1
patients that the non-replication conclusion is not a power artifact.
