"""Which genes are 'spared' from deletion in metastatic disease?

Hypothesis: metastasis requires *intact* machinery. If deep deletion of a gene
abolishes a tumour's ability to metastasize, then metastatic (stage IV / M1)
tumours should almost never carry that deletion -- even when the gene is readily
deleted in non-metastatic (M0) tumours. Such "deletable but spared in M1" genes
are candidate metastasis-required machinery.

Two levels, both adjusted for tumour subtype (M1 is heavily LUAD-skewed, so an
unadjusted test would just rediscover LUAD-vs-LUSC deletion differences):

1. **Gene level** -- Cochran-Mantel-Haenszel depletion test: deep-deleted in >=
   min-freq of M0 but depleted in M1 (:func:`mirna_tcga.associate.cmh_depletion_screen`).
2. **Pathway level** (better powered) -- is a gene set's *deletion burden* lower
   in M1 than M0? Per-patient set-deletion counts vs M1, subtype-stratified
   rank-sum. This is where the 'machinery' shows up.

Caveat: only ~33 M1 patients have copy-number data, so this is hypothesis-
generating. Nodal metastasis (N+) is run as a better-powered sensitivity check.
The machinery could equally be silenced by methylation/mutation -- absence of
deletion is one signal, not proof.

Run:  python scripts/14_metastasis_spared_deletions.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import cmh_depletion_screen, ranksum_screen
from mirna_tcga.biotab import distant_metastasis_labels
from mirna_tcga.config import resolve_path
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.endpoints import (
    distant_metastasis,
    distant_metastasis_biotab,
    nodal_metastasis,
)
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
    """patients x gene-sets matrix of deep-deletion counts within each set."""
    cols = {}
    present = set(delP.columns)
    for name, genes in gene_sets.items():
        g = list(genes & present)
        if len(g) >= min_set:
            cols[name] = delP[g].sum(axis=1)
    return pd.DataFrame(cols, index=delP.index)


def run_endpoint(delP, subP, endpoint, ep_name, gene_sets, min_ref_freq, save_dir):
    common = delP.index.intersection(endpoint.index)
    D = delP.loc[common]
    y = endpoint.loc[common].astype(int)
    strata = subP.loc[common]
    n1 = int(y.sum())
    print(f"\n{'='*70}\n{ep_name}: {n1} metastatic vs {len(y)-n1} non-metastatic "
          f"(subtype-adjusted)\n{'='*70}")
    if n1 < 8:
        print("  too few metastatic cases -- skipped")
        return

    # 1. gene-level depletion
    genes = cmh_depletion_screen(D, y, strata, min_ref_freq=min_ref_freq)
    spared = genes[(genes["or_mh"] < 1)]
    print(f"\nGenes deletable (>= {min_ref_freq:.0%} of M0) yet depleted in {ep_name} "
          f"(top 15 of {len(genes)} tested):")
    print(spared.head(15).round(4).to_string())

    # 2. pathway-level deletion-burden depletion (the 'machinery')
    burden = pathway_burden(D, gene_sets)
    if not burden.empty:
        pw = ranksum_screen(burden, y, strata)          # z<0 => lower burden in metastatic
        spared_pw = pw[pw["z"] < 0].copy()
        print(f"\nGene sets with LOWER deletion burden in {ep_name} "
              f"(spared machinery; top 12 of {burden.shape[1]}):")
        print(spared_pw.head(12)[["z", "auc", "p", "q"]].round(4).to_string())

    # 3. ORA of the never-deleted-in-metastatic candidates
    cand = genes[(genes["n_idx"] == 0)].index
    universe = list(genes.index)
    ora = over_representation(cand, universe, gene_sets)
    if not ora.empty:
        print(f"\nPathways over-represented among genes NEVER deleted in {ep_name} "
              f"({len(cand)} genes, universe={len(universe)} deletable):")
        print(ora.head(10)[["gene_set", "set_size", "overlap", "fold_enrichment", "q"]]
              .to_string(index=False))

    if save_dir:
        out = Path(save_dir)
        out.mkdir(parents=True, exist_ok=True)
        tag = ep_name.split()[0]
        genes.to_csv(out / f"spared_{tag}_genes.csv")
        if not burden.empty:
            pw.to_csv(out / f"spared_{tag}_pathways.csv")
        if not ora.empty:
            ora.to_csv(out / f"spared_{tag}_ora.csv", index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--min-ref-freq", type=float, default=0.04,
                    help="min deep-deletion frequency in M0 to call a gene 'deletable'")
    ap.add_argument("--biotab-root", default=None,
                    help="path to local BCR Biotab GDCdata dir; enriches the distant-met "
                         "endpoint (default: config biotab.root if it exists on disk)")
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    dmet, nmet = distant_metastasis(clin), nodal_metastasis(clin)

    # Enrich distant metastasis with local BCR Biotab clinical supplements if
    # available -- lifts the confirmed-metastatic set from ~33 (pathologic M1) to
    # ~127 (clinical M1 + follow-up 'Distant Metastasis'). Paths may be given
    # relative to the repo root.
    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    dmet_name = "distant met (M1)"
    if biotab_root:
        root = resolve_path(biotab_root)
        if root.exists():
            bt = distant_metastasis_labels(root, keys)
            n_before = int(dmet.sum())
            dmet = distant_metastasis_biotab(clin, bt)
            print(f"Biotab-enriched distant met: {n_before} -> {int(dmet.sum())} "
                  f"positive of {dmet.size} patients (from {root})")
            dmet_name = "distant met (biotab)"
        else:
            print(f"[warn] biotab root {biotab_root!r} not found; using pathologic-M1 endpoint")

    id2sym = protein_coding_map(client)
    print("Fetching deep deletions (HOMDEL) genome-wide ...")
    delB, sub = deletion_matrix(client, cfg, id2sym, keys)
    delP, subP = to_patient(delB, sub.to_dict())
    print(f"  deletion matrix: {delP.shape[0]} patients x {delP.shape[1]} genes")

    cache = Path(args.save_dir or ".") / "genesets_cache"
    gene_sets = load_gene_sets(cfg.raw.get("gene_sets") or DEFAULT_GENE_SETS, cache)

    run_endpoint(delP, subP, dmet, dmet_name, gene_sets, args.min_ref_freq, args.save_dir)
    run_endpoint(delP, subP, nmet, "nodal met (N+)", gene_sets, args.min_ref_freq, args.save_dir)
    if args.save_dir:
        print(f"\nWrote results -> {args.save_dir}/")


if __name__ == "__main__":
    main()
