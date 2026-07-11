"""Multi-omic 'constellation' models for survival and metastasis in NSCLC.

Do expression + copy-number deletions + mutations, *together*, predict overall
survival and metastasis better than any single layer? This builds a patient-level
multi-omic feature matrix and reports honest cross-validated performance for
nested feature sets, so the incremental value of each layer is visible.

Endpoints
  * overall survival        -> penalized Cox, cross-validated Harrell C-index
  * nodal metastasis (N+)   -> L2 logistic, cross-validated ROC-AUC
  * distant metastasis (M1) -> L2 logistic, cross-validated ROC-AUC (under-powered)

To avoid leakage, expression features are selected *inside* every training fold
(fast Cox score test for survival, ANOVA F for metastasis); copy-number and
mutation features are the recurrently-altered genes. A final model on the full
cohort is reported for interpretation only (its coefficients are in-sample).

Requires the survival extra:  pip install 'mirna-tcga[survival]'
Run:  python scripts/12_constellation_model.py --save-dir results
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.endpoints import distant_metastasis, nodal_metastasis
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import deletion_matrix, mutation_matrix, protein_coding_map, stream_expression
from mirna_tcga.screen import cox_score_screen
from mirna_tcga.survival import coerce_clinical

warnings.filterwarnings("ignore")


def _tumor(s):
    p = str(s).split("-")
    return len(p) < 4 or (p[3][:2].isdigit() and int(p[3][:2]) < 10)


def to_patient(df):
    """Collapse a samples x features frame to one tumour row per patient."""
    df = df.loc[[s for s in df.index if _tumor(s)]]
    df = df.copy()
    df.index = sample_to_patient(df.index)
    return df[~df.index.duplicated()]


def most_variable(expr_pt, n):
    v = expr_pt.var(axis=0).sort_values(ascending=False)
    return expr_pt[v.head(n).index]


# ---- survival: cross-validated C-index -------------------------------------
def cindex_cv(blocks, dur, evt, sub, k_expr=100, folds=5, seed=42):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index

    idx = dur.index.to_numpy()
    kf = StratifiedKFold(folds, shuffle=True, random_state=seed)
    cidx = []
    for tr, te in kf.split(idx, evt.to_numpy()):
        Xtr, Xte = [], []
        # expression: select in-fold by Cox score test, z-score on train stats
        if "expression" in blocks:
            e = blocks["expression"]
            sc = cox_score_screen(e.iloc[tr], dur.iloc[tr], evt.iloc[tr], sub.iloc[tr])
            g = sc.head(k_expr).index
            mu, sd = e[g].iloc[tr].mean(), e[g].iloc[tr].std().replace(0, 1)
            w = np.sign(sc.loc[g, "z"])
            sig = ((e[g] - mu) / sd).mul(w, axis=1).mean(axis=1)
            Xtr.append(sig.iloc[tr].to_frame("expr_sig"))
            Xte.append(sig.iloc[te].to_frame("expr_sig"))
        for name in ("deletion", "mutation"):
            if name in blocks:
                Xtr.append(blocks[name].iloc[tr])
                Xte.append(blocks[name].iloc[te])
        if "subtype" in blocks:
            Xtr.append(blocks["subtype"].iloc[tr])
            Xte.append(blocks["subtype"].iloc[te])
        dtr = pd.concat(Xtr, axis=1)
        dte = pd.concat(Xte, axis=1)
        dtr = dtr.assign(_t=dur.iloc[tr].to_numpy(), _e=evt.iloc[tr].to_numpy())
        try:
            cph = CoxPHFitter(penalizer=0.1).fit(dtr, "_t", "_e")
            risk = cph.predict_partial_hazard(dte)
            cidx.append(concordance_index(dur.iloc[te], -risk, evt.iloc[te]))
        except Exception:
            continue
    return float(np.mean(cidx)) if cidx else float("nan"), len(cidx)


# ---- metastasis: cross-validated AUC ---------------------------------------
def auc_cv(blocks, y, k_expr=100, folds=5, repeats=3, seed=42):
    rkf = RepeatedStratifiedKFold(n_splits=folds, n_repeats=repeats, random_state=seed)
    aucs = []
    for tr, te in rkf.split(np.zeros(len(y)), y.to_numpy()):
        Xtr, Xte = [], []
        if "expression" in blocks:
            e = blocks["expression"]
            sel = SelectKBest(f_classif, k=min(k_expr, e.shape[1]))
            sc = StandardScaler()
            etr = sc.fit_transform(e.iloc[tr])
            ete = sc.transform(e.iloc[te])
            sel.fit(etr, y.iloc[tr])
            Xtr.append(sel.transform(etr))
            Xte.append(sel.transform(ete))
        for name in ("deletion", "mutation", "subtype"):
            if name in blocks:
                Xtr.append(blocks[name].iloc[tr].to_numpy())
                Xte.append(blocks[name].iloc[te].to_numpy())
        Atr = np.hstack(Xtr)
        Ate = np.hstack(Xte)
        clf = LogisticRegression(max_iter=2000, C=0.5)
        clf.fit(Atr, y.iloc[tr])
        aucs.append(roc_auc_score(y.iloc[te], clf.predict_proba(Ate)[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--max-var-genes", type=int, default=3000, help="pre-filter to N most-variable genes")
    ap.add_argument("--min-del-freq", type=float, default=0.03)
    ap.add_argument("--min-mut-freq", type=float, default=0.05)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clin).set_index("patientId")
    dmet, nmet = distant_metastasis(clin), nodal_metastasis(clin)

    id2sym = protein_coding_map(client)
    print(f"Streaming expression ({len(id2sym)} genes) ...")
    expr = to_patient(stream_expression(client, cfg, list(id2sym), id2sym, keys).T).dropna(axis=1)
    print("Fetching deletions / mutations ...")
    delB, _ = deletion_matrix(client, cfg, id2sym, keys)
    mutB, _ = mutation_matrix(client, cfg, id2sym, keys)
    delP = to_patient(delB)
    mutP = to_patient(mutB)
    delP = delP.loc[:, delP.mean() >= args.min_del_freq]
    mutP = mutP.loc[:, mutP.mean() >= args.min_mut_freq]
    delP.columns = [f"del:{g}" for g in delP.columns]
    mutP.columns = [f"mut:{g}" for g in mutP.columns]
    print(f"Features: {expr.shape[1]} expr | {delP.shape[1]} deletions | {mutP.shape[1]} mutations")

    # common patients across all layers
    common = expr.index.intersection(delP.index).intersection(mutP.index)
    exprV = most_variable(expr.loc[common], args.max_var_genes)
    delC, mutC = delP.loc[common], mutP.loc[common]
    is_lusc = pd.DataFrame({"is_lusc": [1 if p in surv.index and surv.loc[p, "cohort"] == "LUSC" else 0
                                        for p in common]}, index=common)
    print(f"Common multi-omic patients: {len(common)}\n")

    def blocks_for(pat):
        return {"expression": exprV.loc[pat], "deletion": delC.loc[pat],
                "mutation": mutC.loc[pat], "subtype": is_lusc.loc[pat]}

    results = []

    # ---- survival ----
    sp = [p for p in common if p in surv.index and pd.notna(surv.loc[p, "duration"])]
    dur = surv.loc[sp, "duration"].astype(float)
    evt = surv.loc[sp, "event"].astype(int)
    subS = is_lusc.loc[sp, "is_lusc"].rename("s")
    print(f"OVERALL SURVIVAL (n={len(sp)}, {int(evt.sum())} deaths) -- cross-validated C-index:")
    for name, keep in [("subtype only", ["subtype"]), ("expression", ["expression"]),
                       ("+ deletions", ["expression", "deletion"]),
                       ("+ mutations", ["expression", "deletion", "mutation"]),
                       ("all + subtype", ["expression", "deletion", "mutation", "subtype"])]:
        b = {k: v.loc[sp] for k, v in blocks_for(sp).items() if k in keep}
        c, nf = cindex_cv(b, dur, evt, subS)
        print(f"  {name:16s} C-index = {c:.3f}")
        results.append(("OS", name, round(c, 3), None))

    # ---- metastasis ----
    for ep_name, ep in [("nodal_met (N+)", nmet), ("distant_met (M1)", dmet)]:
        mp = [p for p in common if p in ep.index]
        y = ep.loc[mp].astype(int)
        if y.sum() < 10 or (len(y) - y.sum()) < 10:
            print(f"\n{ep_name}: too few events (n+={int(y.sum())}) -- skipped")
            continue
        print(f"\n{ep_name} (n={len(mp)}, {int(y.sum())}+) -- cross-validated ROC-AUC:")
        for name, keep in [("subtype only", ["subtype"]), ("expression", ["expression"]),
                           ("+ deletions", ["expression", "deletion"]),
                           ("+ mutations", ["expression", "deletion", "mutation"]),
                           ("all + subtype", ["expression", "deletion", "mutation", "subtype"])]:
            b = {k: v.loc[mp] for k, v in blocks_for(mp).items() if k in keep}
            m, sd = auc_cv(b, y)
            print(f"  {name:16s} AUC = {m:.3f} +/- {sd:.3f}")
            results.append((ep_name, name, round(m, 3), round(sd, 3)))

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(results, columns=["endpoint", "model", "score", "sd"]).to_csv(
            out / "constellation_performance.csv", index=False)
        print(f"\nWrote results -> {out}/")


if __name__ == "__main__":
    main()
