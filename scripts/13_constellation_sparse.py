"""Sparse multi-omic 'constellation' of survival and metastasis (NSCLC).

Refines ``12_constellation_model.py`` on three axes the earlier run motivated:

* **Driver-focused DNA layers** -- mutations restricted to a curated NSCLC driver
  panel and deep deletions to recurrently-deleted tumour suppressors, instead of
  every recurrently-altered gene (which added noise, not signal).
* **Sparse (L1) models** -- L1-penalized logistic / Cox do embedded feature
  selection, so the fitted model *is* a short constellation of named features.
* **A fifth layer: miRNA** (UCSC Xena) -- a curated panel incl. the miR-200 /
  let-7 EMT-suppressor families flagged by the miRNA screen.

Expression enters as an outcome-specific signature score selected *inside* every
CV fold (Cox score test for survival, rank-sum for metastasis) to avoid leakage.
Reported: cross-validated performance, and the full-cohort L1 constellation
(non-zero features; the coefficients are in-sample and for interpretation).

Requires the survival extra + Xena egress.
Run:  python scripts/13_constellation_sparse.py --save-dir results
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

from mirna_tcga import load_config
from mirna_tcga.associate import ranksum_screen
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.endpoints import distant_metastasis, nodal_metastasis
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import deletion_matrix, mutation_matrix, protein_coding_map, stream_expression
from mirna_tcga.mirbase import load_mimat_map
from mirna_tcga.panels import DELETION_TSGS, NSCLC_DRIVERS
from mirna_tcga.screen import cox_score_screen
from mirna_tcga.survival import coerce_clinical
from mirna_tcga.xena import XenaClient

warnings.filterwarnings("ignore")

MATURE_FA = "https://raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/mature.fa"
MIRNA_PANEL = [
    "hsa-miR-200a-3p", "hsa-miR-200a-5p", "hsa-miR-200b-3p", "hsa-miR-200c-3p", "hsa-miR-429",
    "hsa-let-7c-5p", "hsa-miR-31-5p", "hsa-miR-31-3p", "hsa-miR-1468-5p", "hsa-miR-132-3p",
    "hsa-miR-148a-3p", "hsa-miR-29c-3p", "hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-miR-210-3p",
]


def _tumor(s):
    p = str(s).split("-")
    return len(p) < 4 or (p[3][:2].isdigit() and int(p[3][:2]) < 10)


def to_patient(df):
    df = df.loc[[s for s in df.index if _tumor(s)]].copy()
    df.index = sample_to_patient(df.index)
    return df[~df.index.duplicated()]


def zscore(train, test):
    mu, sd = train.mean(), train.std().replace(0, 1)
    return (train - mu) / sd, (test - mu) / sd


def os_signature(expr, tr, te, dur, evt, sub, k=100):
    sc = cox_score_screen(expr.iloc[tr], dur.iloc[tr], evt.iloc[tr], sub.iloc[tr])
    g = sc.head(k).index
    w = np.sign(sc.loc[g, "z"])
    mu, sd = expr[g].iloc[tr].mean(), expr[g].iloc[tr].std().replace(0, 1)
    sig = ((expr[g] - mu) / sd).mul(w, axis=1).mean(axis=1)
    return sig.iloc[tr], sig.iloc[te]


def met_signature(expr, tr, te, y, sub, k=100):
    sc = ranksum_screen(expr.iloc[tr], y.iloc[tr], sub.iloc[tr])
    g = sc.head(k).index
    w = np.sign(sc.loc[g, "z"])
    mu, sd = expr[g].iloc[tr].mean(), expr[g].iloc[tr].std().replace(0, 1)
    sig = ((expr[g] - mu) / sd).mul(w, axis=1).mean(axis=1)
    return sig.iloc[tr], sig.iloc[te]


def build_design(cont_tr, cont_te, bin_blocks, tr, te):
    """z-score continuous cols on train, append 0/1 blocks; return train/test frames."""
    ztr, zte = zscore(cont_tr, cont_te)
    parts_tr, parts_te = [ztr], [zte]
    for b in bin_blocks:
        parts_tr.append(b.iloc[tr])
        parts_te.append(b.iloc[te])
    return pd.concat(parts_tr, axis=1), pd.concat(parts_te, axis=1)


def cv_auc_sparse(expr, y, sub, bin_blocks, C=0.3, folds=5, repeats=3, seed=42):
    rkf = RepeatedStratifiedKFold(n_splits=folds, n_repeats=repeats, random_state=seed)
    idx = np.arange(len(y))
    aucs = []
    for tr, te in rkf.split(idx, y.to_numpy()):
        sig_tr, sig_te = met_signature(expr, tr, te, y, sub)
        cont_tr = pd.DataFrame({"expr_sig": sig_tr})
        cont_te = pd.DataFrame({"expr_sig": sig_te})
        for name, blk in bin_blocks.items():
            if name == "mir":  # continuous miRNA: fold into the z-scored block
                cont_tr[blk.columns.tolist()] = blk.iloc[tr].to_numpy()
                cont_te[blk.columns.tolist()] = blk.iloc[te].to_numpy()
        Xtr, Xte = build_design(cont_tr, cont_te,
                                [b for n, b in bin_blocks.items() if n != "mir"], tr, te)
        clf = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=2000)
        clf.fit(Xtr, y.iloc[tr])
        aucs.append(roc_auc_score(y.iloc[te], clf.predict_proba(Xte)[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


def cv_cindex_sparse(expr, dur, evt, sub, bin_blocks, penalizer=0.05, folds=5, seed=42):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    kf = StratifiedKFold(folds, shuffle=True, random_state=seed)
    cidx = []
    for tr, te in kf.split(np.arange(len(dur)), evt.to_numpy()):
        sig_tr, sig_te = os_signature(expr, tr, te, dur, evt, sub)
        cont_tr = pd.DataFrame({"expr_sig": sig_tr})
        cont_te = pd.DataFrame({"expr_sig": sig_te})
        bins = []
        for name, blk in bin_blocks.items():
            if name == "mir":
                cont_tr[blk.columns.tolist()] = blk.iloc[tr].to_numpy()
                cont_te[blk.columns.tolist()] = blk.iloc[te].to_numpy()
            else:
                bins.append(blk)
        Xtr, Xte = build_design(cont_tr, cont_te, bins, tr, te)
        Xtr = Xtr.assign(_t=dur.iloc[tr].to_numpy(), _e=evt.iloc[tr].to_numpy())
        try:
            cph = CoxPHFitter(penalizer=penalizer, l1_ratio=1.0).fit(Xtr, "_t", "_e")
            risk = cph.predict_partial_hazard(Xte)
            cidx.append(concordance_index(dur.iloc[te], -risk, evt.iloc[te]))
        except Exception:
            continue
    return float(np.mean(cidx)) if cidx else float("nan")


def final_logistic(expr, y, sub, bin_blocks, C=0.3):
    """Full-cohort L1 fit -> non-zero named features (the constellation)."""
    sig, _ = met_signature(expr, np.arange(len(y)), np.arange(len(y)), y, sub)
    cont = pd.DataFrame({"expr_sig": sig})
    bins = []
    for name, blk in bin_blocks.items():
        if name == "mir":
            cont[blk.columns.tolist()] = blk.to_numpy()
        else:
            bins.append(blk)
    z = (cont - cont.mean()) / cont.std().replace(0, 1)
    X = pd.concat([z] + bins, axis=1)
    clf = LogisticRegression(penalty="l1", solver="liblinear", C=C, max_iter=3000).fit(X, y)
    coef = pd.Series(clf.coef_[0], index=X.columns)
    return coef[coef.abs() > 1e-6].sort_values(key=lambda s: s.abs(), ascending=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--max-var-genes", type=int, default=3000)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clin).set_index("patientId")
    nmet, dmet = nodal_metastasis(clin), distant_metastasis(clin)

    id2sym = protein_coding_map(client)
    print(f"Streaming expression ({len(id2sym)} genes) ...")
    expr = to_patient(stream_expression(client, cfg, list(id2sym), id2sym, keys).T).dropna(axis=1)
    v = expr.var().sort_values(ascending=False)
    expr = expr[v.head(args.max_var_genes).index]

    print("Fetching driver mutations / TSG deletions ...")
    mutB, _ = mutation_matrix(client, cfg, id2sym, keys)
    delB, _ = deletion_matrix(client, cfg, id2sym, keys)
    mutP = to_patient(mutB).reindex(columns=[g for g in NSCLC_DRIVERS if g in mutB.columns], fill_value=0)
    delP = to_patient(delB).reindex(columns=[g for g in DELETION_TSGS if g in delB.columns], fill_value=0)
    mutP.columns = [f"mut:{g}" for g in mutP.columns]
    delP.columns = [f"del:{g}" for g in delP.columns]

    print("Fetching miRNA panel (Xena) ...")
    names = load_mimat_map(cfg.raw.get("mirbase_mature_fa", MATURE_FA),
                           Path(args.save_dir or ".") / "genesets_cache" / "mature.fa")
    name2mimat = {v: k for k, v in names.items()}
    want = {name2mimat[n]: n for n in MIRNA_PANEL if n in name2mimat}
    xena = XenaClient(**cfg.xena)
    mir_frames = []
    for key in keys:
        m = xena.expression_matrix(cfg.xena_mirna[key], list(want))
        if not m.empty:
            mir_frames.append(m.rename(index=want).T)
    mirP = to_patient(pd.concat(mir_frames, axis=0)) if mir_frames else pd.DataFrame()
    mirP = mirP.reindex(columns=[n for n in MIRNA_PANEL if n in mirP.columns])
    mirP.columns = [f"mir:{c}" for c in mirP.columns]
    mirP = mirP.astype(float)

    is_lusc = pd.Series({p: int(p in surv.index and surv.loc[p, "cohort"] == "LUSC")
                         for p in expr.index}, name="is_lusc").to_frame()

    # patients present in the DNA + expression layers
    base = expr.index.intersection(mutP.index).intersection(delP.index)
    withmir = base.intersection(mirP.dropna(how="all").index)
    print(f"\nPatients: {len(base)} (DNA+expr) | {len(withmir)} also have miRNA\n")

    def blocks(pat, with_mir):
        b = {"mut": mutP.loc[pat], "del": delP.loc[pat], "subtype": is_lusc.loc[pat]}
        if with_mir:
            b["mir"] = mirP.loc[pat].fillna(mirP.loc[pat].median())
        return b

    results = []

    # ---- survival ----
    for tag, pat, wm in [("DNA+expr", base, False), ("+miRNA", withmir, True)]:
        sp = [p for p in pat if p in surv.index and pd.notna(surv.loc[p, "duration"])]
        dur = surv.loc[sp, "duration"].astype(float)
        evt = surv.loc[sp, "event"].astype(int)
        sub = is_lusc.loc[sp, "is_lusc"].rename("s")
        c = cv_cindex_sparse(expr.loc[sp], dur, evt, sub, blocks(sp, wm))
        print(f"OS  [{tag:9s} n={len(sp)}]  sparse C-index = {c:.3f}")
        results.append(("OS", tag, round(c, 3), None))

    # ---- metastasis ----
    for ep_name, ep in [("nodal_met", nmet), ("distant_met", dmet)]:
        for tag, pat, wm in [("DNA+expr", base, False), ("+miRNA", withmir, True)]:
            mp = [p for p in pat if p in ep.index]
            y = ep.loc[mp].astype(int)
            if y.sum() < 12 or (len(y) - y.sum()) < 12:
                print(f"{ep_name} [{tag}]: too few events ({int(y.sum())}+) -- skipped")
                continue
            sub = is_lusc.loc[mp, "is_lusc"].rename("s")
            m, sd = cv_auc_sparse(expr.loc[mp], y, sub, blocks(mp, wm))
            print(f"{ep_name} [{tag:9s} n={len(mp)}, {int(y.sum())}+]  sparse AUC = {m:.3f} +/- {sd:.3f}")
            results.append((ep_name, tag, round(m, 3), round(sd, 3)))

    # ---- the constellation (full-cohort L1 fit, nodal metastasis) ----
    mp = [p for p in withmir if p in nmet.index]
    y = nmet.loc[mp].astype(int)
    sub = is_lusc.loc[mp, "is_lusc"].rename("s")
    con = final_logistic(expr.loc[mp], y, sub, blocks(mp, True))
    print(f"\nNodal-metastasis constellation (L1-selected, n={len(mp)}, {int(y.sum())}+):")
    for feat, coef in con.items():
        arrow = "up->more met" if coef > 0 else "down->more met"
        print(f"  {feat:24s} {coef:+.3f}  ({arrow})")

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(results, columns=["endpoint", "features", "score", "sd"]).to_csv(
            out / "constellation_sparse_performance.csv", index=False)
        con.rename("coef").to_csv(out / "constellation_nodal_features.csv")
        print(f"\nWrote results -> {out}/")


if __name__ == "__main__":
    main()
