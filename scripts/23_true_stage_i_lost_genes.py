"""Genes *lost* in true stage I -- candidate metastasis-permissive machinery.

A different way to hunt metastasis genes, since the bulk primary carries no
positive signature (script 22): look for **negatives**. A "true stage I" tumour
(stage I at diagnosis, N0/M0, and confirmed indolent -- no distant or locoregional
recurrence for >= 2 years) demonstrably *never* acquired the ability to spread. If
metastasis requires some machinery, indolent tumours may have **lost** it -- e.g.
deep-deleted a gene they never needed. So: which genes are **deletable in true
stage I yet retained (not deleted) in tumours that did progress**?

This is the script-14 negative-selection test with a **cleaner reference**: the
comparator is confirmed-indolent *true stage I* rather than the heterogeneous
"all M0" (which secretly contains tumours that later metastasized). Mechanically,
with the metastatic group as the index, a gene "spared from deletion in metastatic"
(CMH or_mh < 1) is exactly a gene "deleted more in true stage I" -- lost there.

Two progression-capable indices vs the true-stage-I reference:
  * distant metastasis (M1 / distant new-tumor event)  -- the strongest contrast;
  * nodal N+                                            -- better powered.

Run:  python scripts/23_true_stage_i_lost_genes.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import cmh_depletion_screen, ranksum_screen
from mirna_tcga.biotab import true_stage_i_vs_distant
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, combined_clinical
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import nodal_metastasis
from mirna_tcga.enrich import load_gene_sets, over_representation
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import deletion_matrix, protein_coding_map

DEFAULT_GENE_SETS = {
    "KEGG_2016": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/enrichr.KEGG_2016.gmt",
    "Hallmark": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/h.all.v7.0.symbols.gmt",
}


def to_patient(df, sub_map):
    df = df.copy()
    pats = sample_to_patient(df.index)
    sub = pd.Series([sub_map.get(s) for s in df.index], index=pats)
    df.index = pats
    keep = ~df.index.duplicated()
    return df[keep], sub[keep]


def pathway_burden(delP, gene_sets, min_set=5):
    cols, present = {}, set(delP.columns)
    for name, genes in gene_sets.items():
        g = list(genes & present)
        if len(g) >= min_set:
            cols[name] = delP[g].sum(axis=1)
    return pd.DataFrame(cols, index=delP.index)


def run(name, delP, subP, outcome, gene_sets, min_ref_freq, save_dir):
    """outcome: 1 = progression-capable (index), 0 = true stage I (reference)."""
    common = delP.index.intersection(outcome.index)
    D, y, strata = delP.loc[common], outcome.loc[common].astype(int), subP.loc[common]
    n1, n0 = int(y.sum()), int((y == 0).sum())
    print(f"\n{'='*70}\n{name}: {n1} progressed vs {n0} true-stage-I "
          f"(subtype-adjusted)\n{'='*70}")
    if n1 < 8 or n0 < 8:
        print("  too few in a group -- skipped")
        return None

    genes = cmh_depletion_screen(D, y, strata, min_ref_freq=min_ref_freq)
    # or_mh < 1 = depleted in the progressed index = deleted MORE in true stage I.
    lost = genes[genes["or_mh"] < 1].copy()
    lost = lost.rename(columns={"freq_ref": "del_in_TSI", "freq_idx": "del_in_progressed"})
    print(f"\nGenes deletable in true stage I (>= {min_ref_freq:.0%}) yet RETAINED in "
          f"progressed tumours -> 'lost' in true stage I (top 15 of {len(genes)}):")
    print(lost.head(15)[["del_in_TSI", "del_in_progressed", "n_idx", "or_mh", "p", "q"]]
          .round(4).to_string())

    burden = pathway_burden(D, gene_sets)
    if not burden.empty:
        pw = ranksum_screen(burden, y, strata)     # z<0 => lower burden in progressed = lost in TSI
        lost_pw = pw[pw["z"] < 0]
        print(f"\nPathways with HIGHER deletion burden in true stage I "
              f"(lost machinery; top 10 of {burden.shape[1]}):")
        print(lost_pw.head(10)[["z", "auc", "p", "q"]].round(4).to_string())

    if save_dir:
        out = Path(save_dir); out.mkdir(parents=True, exist_ok=True)
        tag = name.split()[0]
        genes.to_csv(out / f"lostTSI_{tag}_genes.csv")
    return genes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--studies", nargs="+", default=None)
    ap.add_argument("--min-ref-freq", type=float, default=0.04,
                    help="min deep-deletion freq in true stage I to call a gene 'deletable/lost' there")
    ap.add_argument("--biotab-root", default=None)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = args.studies or cohort_study_keys(cfg, "nsclc")
    root = resolve_path(args.biotab_root or (cfg.raw.get("biotab") or {}).get("root"))
    print(f"Cohort: {', '.join(k.upper() for k in keys)} | biotab {root}")

    clin = combined_clinical(client, cfg, keys, patient_level=True)
    ts = true_stage_i_vs_distant(root, keys)          # 0 = true stage I, 1 = distant met
    tsi = set(ts[ts == 0].index)
    print(f"true stage I (indolent, >=2yr): {len(tsi)} | distant-met: {int((ts == 1).sum())}")

    # N+ vs true-stage-I contrast (true stage I is N0 by definition -> no overlap)
    nmet = nodal_metastasis(clin)
    nplus = set(nmet[nmet == 1].index)
    idx = sorted(tsi | nplus)
    nplus_vs_tsi = pd.Series([1 if p in nplus else 0 for p in idx], index=idx)

    id2sym = protein_coding_map(client)
    print("Fetching deep deletions (HOMDEL) genome-wide ...")
    delB, sub = deletion_matrix(client, cfg, id2sym, keys)
    delP, subP = to_patient(delB, sub.to_dict())
    print(f"  deletion matrix: {delP.shape[0]} patients x {delP.shape[1]} genes")

    cache = Path(args.save_dir or ".") / "genesets_cache"
    gene_sets = load_gene_sets(cfg.raw.get("gene_sets") or DEFAULT_GENE_SETS, cache)

    run("distant-met vs true-stage-I", delP, subP, ts, gene_sets, args.min_ref_freq, args.save_dir)
    run("N+ vs true-stage-I", delP, subP, nplus_vs_tsi, gene_sets, args.min_ref_freq, args.save_dir)
    if args.save_dir:
        print(f"\nWrote results -> {args.save_dir}/")


if __name__ == "__main__":
    main()
