"""Genome-wide differential expression: metastatic vs non-metastatic NSCLC.

The deletion screens (scripts 14-16) ask what copy-number *loss* distinguishes
metastatic disease. This asks the broader, more direct question: **which genes'
mRNA expression differs most between metastatic and non-metastatic tumours?**

Pipeline (mirrors the survival screen, scripts/09, but for a metastasis endpoint):
1. NSCLC clinical + the **Biotab-enriched** distant-metastasis endpoint
   (126 vs 674); the better-powered nodal endpoint (N+, ~343 vs ~637) is run as a
   sensitivity check.
2. Stream mRNA expression for every protein-coding gene (log2), cBioPortal.
3. Subtype-stratified **van Elteren rank-sum** (`associate.ranksum_screen`) of
   every gene vs metastatic status -- z > 0 means *higher* in metastatic tumours.
   Stratifying by LUAD/LUSC is essential: the metastatic set is LUAD-heavy, so an
   unadjusted test would rediscover LUAD-vs-LUSC expression differences.
4. Pathway over-representation of the FDR hits, split by direction (up vs down in
   metastatic) -- e.g. does the up-in-metastatic set enrich for EMT?

Caveat: TCGA samples are resected **primary** tumours, so this contrasts the
primary-tumour transcriptome of patients *with* vs *without* metastatic disease,
not metastases themselves.

Run (full):   python scripts/17_metastasis_expression_diff.py --save-dir results
Run (quick):  python scripts/17_metastasis_expression_diff.py --max-genes 3000
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import ranksum_screen
from mirna_tcga.biotab import distant_metastasis_labels, true_stage_i_vs_distant
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import (
    distant_metastasis,
    distant_metastasis_biotab,
    nodal_metastasis,
)
from mirna_tcga.enrich import load_gene_sets, over_representation
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import protein_coding_map, stream_expression

DEFAULT_GENE_SETS = {
    "KEGG_2016": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/enrichr.KEGG_2016.gmt",
    "Hallmark": "https://raw.githubusercontent.com/zqfang/GSEApy/master/tests/extdata/h.all.v7.0.symbols.gmt",
}
_N_POS = {"N1", "N2", "N3"}


def _simplify_stage(v: object) -> str | None:
    s = str(v).upper().replace("STAGE ", "").strip()
    for r in ("IV", "III", "II", "I"):
        if s.startswith(r):
            return r
    return None


def stage_restricted_nodal(clin: pd.DataFrame, stage_code: str) -> pd.Series:
    """N+ vs N0 restricted to one pathologic stage -- holds stage constant.

    Isolates nodal spread from tumour stage/progression: within a single stage,
    does node-positive vs node-negative still carry the proliferation signature,
    or was that signal really just higher-stage tumours proliferating more?
    """
    st = clin.set_index("patientId")["AJCC_PATHOLOGIC_TUMOR_STAGE"].map(_simplify_stage)
    n = clin.set_index("patientId")["PATH_N_STAGE"].astype(str).str.upper().str[:2]
    lab = pd.Series(pd.NA, index=st.index, dtype="Int64")
    in_stage = st.eq(stage_code)
    lab[in_stage & n.eq("N0")] = 0
    lab[in_stage & n.isin(_N_POS)] = 1
    return lab.dropna().astype(int)


def run_endpoint(name, Xg, pats, endpoint, subtype, gene_sets, fdr):
    """van Elteren DE screen of Xg (samples x genes) vs one 0/1 endpoint."""
    y = pd.Series(endpoint.reindex(pats).to_numpy(), index=Xg.index)
    strata = pd.Series(subtype.reindex(pats).to_numpy(), index=Xg.index)
    keep = y.notna() & strata.notna()
    y, strata = y[keep].astype(int), strata[keep]
    X = Xg.loc[keep]
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    print(f"\n{'='*70}\n{name}: {n1} metastatic vs {n0} non-metastatic "
          f"(subtype-adjusted, {X.shape[1]} genes)\n{'='*70}")

    res = ranksum_screen(X, y, strata)
    hits = res[res["q"] < fdr]
    up = hits[hits["z"] > 0]      # higher in metastatic
    down = hits[hits["z"] < 0]    # lower in metastatic
    print(f"Genes differential at BH q < {fdr}: {len(hits)}  "
          f"(up in metastatic: {len(up)}, down: {len(down)})")
    print("\nTop 20 by significance (z>0 = higher in metastatic):")
    show = res.head(20).assign(dir=np.where(res.head(20)["z"] > 0, "UP", "down"))
    print(show[["z", "auc", "p", "q", "dir"]].round(4).to_string())

    universe = list(X.columns)
    ora_tables = {}
    for label, genes in (("up in metastatic", up.index), ("down in metastatic", down.index)):
        ora = over_representation(genes, universe, gene_sets)
        ora_tables[label] = ora
        sig = ora[ora["q"] < 0.05] if not ora.empty else ora
        print(f"\n--- pathways {label}: {len(genes)} genes, {len(sig)} sets at q<0.05 ---")
        if not ora.empty:
            print(ora.head(8)[["gene_set", "set_size", "overlap", "fold_enrichment", "q"]]
                  .to_string(index=False))
    return res, ora_tables


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--fdr", type=float, default=0.05, help="BH q threshold")
    ap.add_argument("--max-genes", type=int, default=None, help="cap universe (quick run)")
    ap.add_argument("--biotab-root", default=None,
                    help="local BCR Biotab GDCdata dir (default: config biotab.root if present)")
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    subtype = clin.set_index("patientId")["cohort"]
    dmet = distant_metastasis(clin)
    nmet = nodal_metastasis(clin)

    # Contrasts to screen. The default distant/nodal references are "all M0" /
    # "all N0"; the third uses a clean **true stage I** reference (stage I + N0 +
    # M0, indolent >= 2 years) so the distant-met contrast is not diluted by the
    # heterogeneous, part-node-positive M0 group.
    # Stage-II-only nodal contrast: N+ vs N0 with stage held constant (stage II has
    # a balanced 177 vs 105 split; stage I is ~all N0 and stage III ~all N+).
    nmet_s2 = stage_restricted_nodal(clin, "II")
    contrasts = [
        ("distant metastasis (biotab)", dmet),
        ("nodal metastasis (N+)", nmet),
        ("nodal N+ vs N0, stage II only", nmet_s2),
    ]
    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    if biotab_root and resolve_path(biotab_root).exists():
        root = resolve_path(biotab_root)
        bt = distant_metastasis_labels(root, keys)
        n0 = int(dmet.sum())
        dmet = distant_metastasis_biotab(clin, bt)
        contrasts[0] = ("distant metastasis (biotab)", dmet)
        print(f"Distant-met label: Biotab-enriched {n0} -> {int(dmet.sum())} positive")
        ts1 = true_stage_i_vs_distant(root, keys)
        print(f"True-stage-I contrast: {int((ts1 == 1).sum())} distant-met cases vs "
              f"{int((ts1 == 0).sum())} true-stage-I controls")
        contrasts.append(("distant met vs true stage I", ts1))

    id2sym = protein_coding_map(client)
    if args.max_genes:
        id2sym = dict(list(id2sym.items())[: args.max_genes])
    print(f"Streaming expression for {len(id2sym)} protein-coding genes (LUAD + LUSC) ...")
    expr = stream_expression(client, cfg, list(id2sym), id2sym, keys)
    print(f"  expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    X = expr.T                                   # samples x genes
    pats = pd.Series(sample_to_patient(X.index), index=X.index)
    Xg = X.dropna(axis=1)                         # genes with complete data

    cache = Path(args.save_dir or ".") / "genesets_cache"
    gene_sets = load_gene_sets(cfg.raw.get("gene_sets") or DEFAULT_GENE_SETS, cache)

    results = {}
    for name, ep in contrasts:
        res, ora = run_endpoint(name, Xg, pats, ep, subtype, gene_sets, args.fdr)
        results[name] = (res, ora)

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        for name, (res, ora) in results.items():
            tag = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
            res.to_csv(out / f"de_{tag}_genes.csv")
            for label, tbl in ora.items():
                if not tbl.empty:
                    tbl.to_csv(out / f"de_{tag}_pathways_{label.split()[0]}.csv", index=False)
        print(f"\nWrote results -> {out}/")


if __name__ == "__main__":
    main()
