"""Which 8p23 gene carries the "lost in true stage I" signal?

Script 23 found the 8p23 block deleted in indolent true-stage-I tumours but
retained in ones that progressed (q ~ 0.027). But the block is **co-deleted as one
segment**, so by deletion the genes are collinear -- you cannot tell which one
matters from deletion frequency (they all move together). This script tries to
break the tie with **expression**:

1. Define the block **empirically**: genes whose deep deletion co-occurs with the
   8p23 anchors (XKR6/SOX7/PINX1) across patients (deletion-vector correlation).
2. Groups: progressed (distant-met OR N+) vs true-stage-I indolent.
3. Per block gene, three reads:
   * **cis-dosage** -- expression in deep-deleted vs copy-neutral tumours (sanity:
     deletion should lower expression);
   * **z_all** -- subtype-adjusted rank-sum of expression, progressed vs indolent,
     over *all* tumours (dominated by dosage: indolent delete -> lower);
   * **z_neutral** -- the same test restricted to tumours where *that gene is not
     deleted* (GISTIC >= 0). This removes the shared dosage; a gene whose
     expression still tracks progression here is a **dosage-independent** candidate,
     i.e. the one that plausibly matters -- not just a passenger of the segment.

If no gene shows a z_neutral signal, the honest read is that it is the *segmental
deletion* that matters, not any single gene (which this data cannot resolve).

Run:  python scripts/24_chr8p23_which_gene.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import ranksum_screen
from mirna_tcga.biotab import true_stage_i_vs_distant
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, combined_clinical, combined_expression
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import nodal_metastasis
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import deletion_matrix, protein_coding_map
from mirna_tcga.panels import PROLIFERATION_PANEL

ANCHORS = ["XKR6", "SOX7", "PINX1"]

# Plausible metastasis-relevant function of 8p23 members (for annotating hits).
FUNCTION = {
    "STC1": "stanniocalcin-1, secreted; hypoxia-induced EMT/invasion + angiogenesis (8p21)",
    "LOXL2": "lysyl oxidase-like 2, ECM collagen crosslinking; classic EMT/metastasis driver (8p21)",
    "ANGPT2": "angiopoietin-2, angiogenesis/vascular remodeling",
    "ARHGEF10": "Rho guanine-exchange factor, actin/cell migration",
    "BLK": "SRC-family tyrosine kinase, signaling",
    "SOX7": "TF, Wnt antagonist (usually tumor suppressor)",
    "PINX1": "telomerase inhibitor (usually tumor suppressor)",
    "CSMD1": "complement regulator (tumor suppressor)",
    "MSRA": "methionine sulfoxide reductase, oxidative-stress defense",
    "MTMR9": "myotubularin phosphatase, PI3K/autophagy",
    "PPP1R3B": "glycogen metabolism (PP1 regulatory subunit)",
    "DEFA6": "alpha-defensin, innate immunity",
}


def to_patient(df, sub_map):
    df = df.copy()
    pats = sample_to_patient(df.index)
    sub = pd.Series([sub_map.get(s) for s in df.index], index=pats)
    df.index = pats
    keep = ~df.index.duplicated()
    return df[keep], sub[keep]


def one_gene_z(expr_g: pd.Series, y: pd.Series, strata: pd.Series):
    """Subtype-adjusted rank-sum z for a single gene's expression vs group y."""
    d = pd.concat([expr_g.rename("x"), y.rename("y"), strata.rename("s")], axis=1).dropna()
    if d["y"].nunique() < 2 or (d["y"] == 1).sum() < 8 or (d["y"] == 0).sum() < 8:
        return np.nan, np.nan, len(d)
    res = ranksum_screen(d[["x"]], d["y"], d["s"])
    if res.empty:                       # zero-variance / all-tied gene
        return np.nan, np.nan, len(d)
    return float(res["z"].iloc[0]), float(res["p"].iloc[0]), len(d)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--studies", nargs="+", default=None)
    ap.add_argument("--min-del-freq", type=float, default=0.04)
    ap.add_argument("--corr", type=float, default=0.5, help="deletion-vector corr with anchor to join block")
    ap.add_argument("--biotab-root", default=None)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = args.studies or cohort_study_keys(cfg, "nsclc")
    root = resolve_path(args.biotab_root or (cfg.raw.get("biotab") or {}).get("root"))
    clin = combined_clinical(client, cfg, keys, patient_level=True)
    subtype = clin.set_index("patientId")["cohort"]

    # progressed (distant OR N+) vs true-stage-I indolent
    ts = true_stage_i_vs_distant(root, keys)          # 0 = true stage I, 1 = distant
    nmet = nodal_metastasis(clin)
    tsi = set(ts[ts == 0].index)
    progressed = set(ts[ts == 1].index) | set(nmet[nmet == 1].index)
    pats_all = sorted(tsi | progressed)
    y = pd.Series([1 if p in progressed else 0 for p in pats_all], index=pats_all)
    print(f"Cohort {', '.join(k.upper() for k in keys)}: "
          f"progressed(distant|N+)={int(y.sum())} vs true-stage-I={int((y==0).sum())}")

    # 1. deletion matrix -> define the co-deleted block empirically
    id2sym = protein_coding_map(client)
    print("Fetching deep deletions ...")
    delB, sub = deletion_matrix(client, cfg, id2sym, keys)
    delP, _ = to_patient(delB, sub.to_dict())
    anchors = [a for a in ANCHORS if a in delP.columns]
    anchor_score = delP[anchors].max(axis=1)
    deletable = delP.columns[(delP.mean() >= args.min_del_freq)]
    corr = delP[deletable].apply(lambda c: c.corr(anchor_score))
    block = sorted(corr[corr >= args.corr].index)
    print(f"Empirical 8p23 co-deletion block ({len(block)} genes, "
          f"deletion-corr >= {args.corr} with {'/'.join(anchors)}):\n  {', '.join(block)}")
    if not block:
        raise SystemExit("No block genes found.")

    # 2. expression + discrete CNA for the block (+ proliferation panel for adjustment)
    expr = combined_expression(client, cfg, block + PROLIFERATION_PANEL, keys)   # genes x samples
    cna_frames = [client.cna_matrix(cfg.cna_profile(k), block, cfg.cna_samples_list(k)) for k in keys]
    cna = pd.concat([c for c in cna_frames if not c.empty], axis=1)
    print(f"  expression {expr.shape}, CNA {cna.shape}")

    X = expr.T
    pats = pd.Series(sample_to_patient(X.index), index=X.index)
    ysamp = pd.Series(y.reindex(pats).to_numpy(), index=X.index)
    strata = pd.Series(subtype.reindex(pats).to_numpy(), index=X.index)

    # proliferation tertile within subtype -> strata for the proliferation-adjusted test.
    # z_all is confounded by proliferation (progressed tumours proliferate more); this
    # separates genuine metastasis-associated genes from proliferation passengers.
    panel = [g for g in PROLIFERATION_PANEL if g in X.columns]
    pz = (X[panel] - X[panel].mean()) / X[panel].std()
    pscore = pz.mean(axis=1)
    pstrata = pd.Series(index=X.index, dtype=object)
    for s in strata.dropna().unique():
        m = strata == s
        q1, q2 = pscore[m].quantile(1 / 3), pscore[m].quantile(2 / 3)
        pstrata[m] = np.where(pscore[m] <= q1, f"{s}_lo",
                              np.where(pscore[m] >= q2, f"{s}_hi", f"{s}_mid"))

    rows = []
    for g in block:
        if g not in X.columns:
            continue
        eg = X[g]
        # cis-dosage: expr in deep-deleted vs neutral. CNA is indexed by SAMPLE id
        # (same as expression), so align directly to X.index -- not via patient.
        cn_s = cna.loc[g].reindex(X.index) if g in cna.index else pd.Series(np.nan, index=X.index)
        deleted = cn_s <= -2
        neutral = cn_s >= 0
        dose = eg[deleted].median() - eg[neutral].median() if deleted.sum() >= 5 else np.nan
        z_all, p_all, _ = one_gene_z(eg, ysamp, strata)
        z_neu, p_neu, _ = one_gene_z(eg[neutral], ysamp[neutral], strata[neutral])
        z_pro, p_pro, _ = one_gene_z(eg, ysamp, pstrata)      # proliferation-adjusted
        rows.append((g, delP[g].mean(), dose, z_all, p_all, z_neu, p_neu, z_pro, p_pro,
                     int(neutral.sum())))

    tab = pd.DataFrame(rows, columns=["gene", "del_freq", "cis_dose", "z_all", "p_all",
                                      "z_neutral", "p_neutral", "z_prolif", "p_prolif",
                                      "n_neutral"]).set_index("gene")
    # A real metastasis-associated 8p gene is UP in progressed and survives BOTH
    # confounders: dosage (z_neutral) and proliferation (z_prolif).
    tab = tab.sort_values("p_prolif")
    print("\nPer block gene, UP in progressed, sorted by proliferation-adjusted p "
          "(z>0 = higher in progressed):")
    up = tab[tab["z_all"] > 0]
    print(up[["del_freq", "cis_dose", "z_all", "p_all", "z_neutral", "p_neutral",
              "z_prolif", "p_prolif"]].round(3).head(20).to_string())
    print("\nCandidates surviving BOTH dosage and proliferation adjustment "
          "(z_all>0, p_neutral<0.05, p_prolif<0.05):")
    cand = tab[(tab["z_all"] > 0) & (tab["p_neutral"] < 0.05) & (tab["p_prolif"] < 0.05)]
    if cand.empty:
        print("  NONE -- no single gene separates from the co-deleted segment.")
    else:
        for g in cand.sort_values("p_prolif").index:
            print(f"  {g}: z_prolif={cand.loc[g,'z_prolif']:+.2f} (p={cand.loc[g,'p_prolif']:.3g}), "
                  f"z_neutral={cand.loc[g,'z_neutral']:+.2f}  [{FUNCTION.get(g, 'function n/a')}]")

    if args.save_dir:
        out = Path(args.save_dir); out.mkdir(parents=True, exist_ok=True)
        tab.to_csv(out / "chr8p23_which_gene.csv")
        print(f"\nWrote -> {out}/chr8p23_which_gene.csv")


if __name__ == "__main__":
    main()
