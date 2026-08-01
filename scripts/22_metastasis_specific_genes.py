"""Metastasis-*specific* genes: differential expression adjusted for proliferation.

Every metastasis contrast so far is dominated by **proliferation** -- metastatic
TCGA primaries are simply more proliferative (scripts 17: G2M/E2F/MYC up). That is
tumour aggression, not metastasis machinery. This script separates the two by
asking: within tumours of *similar proliferation*, which genes still differ
between metastatic and non-metastatic?

Method: score proliferation per sample (`panels.PROLIFERATION_PANEL`), bin into
tertiles **within each subtype**, and run the subtype-stratified van Elteren
rank-sum (`associate.ranksum_screen`) with strata = subtype x proliferation-tertile.
Genes significant unadjusted but not after adjustment were proliferation-driven;
genes that survive adjustment are candidate **metastasis-specific** markers. Each
endpoint is run both ways and the two hit-lists are diffed.

Endpoints: nodal N+ (best powered), nodal N+ vs N0 within stage II (stage-matched),
and the Biotab distant-met endpoint.

Run (full):   python scripts/22_metastasis_specific_genes.py --save-dir results
Run (quick):  python scripts/22_metastasis_specific_genes.py --max-genes 3000
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import ranksum_screen
from mirna_tcga.biotab import distant_metastasis_labels
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, combined_clinical
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import (
    distant_metastasis,
    distant_metastasis_biotab,
    nodal_metastasis,
    nodal_metastasis_within_stage,
)
from mirna_tcga.enrich import load_gene_sets, over_representation
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import protein_coding_map, stream_expression
from mirna_tcga.panels import PROLIFERATION_PANEL

DEFAULT_GENE_SETS = {
    "KEGG_2016": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/enrichr.KEGG_2016.gmt",
    "Hallmark": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/h.all.v7.0.symbols.gmt",
}


def proliferation_tertile(Xg: pd.DataFrame, subtype: pd.Series) -> pd.Series:
    """Per-sample proliferation tertile (low/mid/high), computed within subtype."""
    panel = [g for g in PROLIFERATION_PANEL if g in Xg.columns]
    z = (Xg[panel] - Xg[panel].mean()) / Xg[panel].std()
    score = z.mean(axis=1)
    tert = pd.Series(index=Xg.index, dtype=object)
    for sub in subtype.dropna().unique():
        m = subtype == sub
        s = score[m]
        if len(s) < 6:
            tert[m] = f"{sub}_all"
            continue
        q1, q2 = s.quantile(1 / 3), s.quantile(2 / 3)
        tert[m] = np.where(s <= q1, f"{sub}_lo", np.where(s >= q2, f"{sub}_hi", f"{sub}_mid"))
    return tert, score


def run(name, Xg, y, subtype, prolif_strata, gene_sets, fdr, save_dir):
    keep = y.notna() & subtype.notna()
    y, sub, strat = y[keep].astype(int), subtype[keep], prolif_strata[keep]
    X = Xg.loc[keep]
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    print(f"\n{'#'*72}\n{name}: {n1} metastatic vs {n0} non-metastatic\n{'#'*72}")

    unadj = ranksum_screen(X, y, sub)                      # subtype only
    adj = ranksum_screen(X, y, strat)                      # subtype x prolif tertile
    hu = set(unadj[unadj["q"] < fdr].index)
    ha = set(adj[adj["q"] < fdr].index)
    prolif = set(PROLIFERATION_PANEL)
    # how much of the unadjusted signal is proliferation biology?
    u_top = unadj.head(200)
    print(f"hits @ q<{fdr}: unadjusted {len(hu)}  |  proliferation-adjusted {len(ha)}")
    print(f"  of the top 200 unadjusted hits, {u_top.index.isin(prolif).sum()}/25 "
          f"proliferation-panel genes present, median |z|={u_top['z'].abs().median():.2f}")
    print(f"  proliferation-driven (lost on adjustment): {len(hu - ha)}")
    print(f"  metastasis-specific (survive adjustment):  {len(ha)}  "
          f"({len(ha & hu)} also unadjusted, {len(ha - hu)} new)")

    adj_hits = adj[adj["q"] < fdr].copy()
    adj_up = adj_hits[adj_hits["z"] > 0]
    adj_dn = adj_hits[adj_hits["z"] < 0]
    # Tertile stratification is conservative, so always show the strongest RESIDUAL
    # candidates by adjusted p (even if q>fdr) -- the honest "is anything left?" view.
    show = adj.head(20).assign(
        dir=np.where(adj.head(20)["z"] > 0, "UP", "dn"),
        prolif=np.where(adj.head(20).index.isin(prolif), "*prolif*", ""))
    label = "metastasis-specific hits" if len(adj_hits) else "strongest residual candidates (none pass FDR)"
    print(f"\nTop 20 {label} (proliferation-adjusted; z>0 = up in metastatic):")
    print(show[["z", "auc", "p", "q", "dir", "prolif"]].round(4).to_string())

    universe = list(X.columns)
    ora_tables = {}
    for label, genes in (("up in metastatic (adj)", adj_up.index),
                         ("down in metastatic (adj)", adj_dn.index)):
        ora = over_representation(genes, universe, gene_sets)
        ora_tables[label] = ora
        sig = ora[ora["q"] < 0.05] if not ora.empty else ora
        print(f"\n--- pathways {label}: {len(genes)} genes, {len(sig)} sets q<0.05 ---")
        if not ora.empty:
            print(ora.head(6)[["gene_set", "set_size", "overlap", "fold_enrichment", "q"]]
                  .to_string(index=False))

    if save_dir:
        out = Path(save_dir); out.mkdir(parents=True, exist_ok=True)
        tag = name.split("(")[0].strip().replace(" ", "_").replace(",", "")
        adj.to_csv(out / f"metspecific_{tag}_adjusted.csv")
        unadj.to_csv(out / f"metspecific_{tag}_unadjusted.csv")
    return adj


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--max-genes", type=int, default=None)
    ap.add_argument("--studies", nargs="+", default=None)
    ap.add_argument("--biotab-root", default=None)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = args.studies or cohort_study_keys(cfg, "nsclc")
    clin = combined_clinical(client, cfg, keys, patient_level=True)
    subtype_pat = clin.set_index("patientId")["cohort"]
    print(f"Cohort: {', '.join(k.upper() for k in keys)}")

    dmet = distant_metastasis(clin)
    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    if biotab_root and resolve_path(biotab_root).exists():
        bt = distant_metastasis_labels(resolve_path(biotab_root), keys)
        dmet = distant_metastasis_biotab(clin, bt)
    endpoints = [
        ("nodal N+", nodal_metastasis(clin)),
        ("nodal N+ vs N0, stage II", nodal_metastasis_within_stage(clin, "II")),
        ("distant met (biotab)", dmet),
    ]

    id2sym = protein_coding_map(client)
    if args.max_genes:
        id2sym = dict(list(id2sym.items())[: args.max_genes])
    print(f"Streaming expression for {len(id2sym)} genes ...")
    expr = stream_expression(client, cfg, list(id2sym), id2sym, keys)
    X = expr.T
    pats = pd.Series(sample_to_patient(X.index), index=X.index)
    Xg = X.dropna(axis=1)
    subtype = pd.Series(subtype_pat.reindex(pats).to_numpy(), index=Xg.index)
    strata, score = proliferation_tertile(Xg, subtype)
    print(f"  matrix {Xg.shape[0]} samples x {Xg.shape[1]} genes; "
          f"proliferation strata: {sorted(strata.dropna().unique())}")

    cache = Path(args.save_dir or ".") / "genesets_cache"
    gene_sets = load_gene_sets(cfg.raw.get("gene_sets") or DEFAULT_GENE_SETS, cache)

    for name, ep in endpoints:
        y = pd.Series(ep.reindex(pats).to_numpy(), index=Xg.index)
        run(name, Xg, y, subtype, strata, gene_sets, args.fdr, args.save_dir)

    if args.save_dir:
        print(f"\nWrote results -> {args.save_dir}/")


if __name__ == "__main__":
    main()
