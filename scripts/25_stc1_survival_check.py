"""Sanity check: does high STC1 predict worse overall survival in NSCLC?

Internal-consistency test for the 8p / STC1 metastasis lead (scripts 23-24). If
STC1 is metastasis-permissive rather than a generic proliferation marker, high
STC1 should predict **worse OS** -- and that association should survive adjustment
for proliferation (the confound that dominated every metastasis contrast).

Uses the validated vectorized **Cox score test** (`screen.cox_score_screen`, no
lifelines needed): z > 0 means higher expression => worse survival. Each gene is
tested two ways:
  * subtype-stratified (LUAD/LUSC);
  * subtype x proliferation-tertile stratified (`panels.PROLIFERATION_PANEL`).
MKI67 is included as a **proliferation positive control**: its OS association
should *vanish* under proliferation adjustment, confirming the adjustment works
and that any gene surviving it is not merely a proliferation proxy.

Run:  python scripts/25_stc1_survival_check.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, combined_expression, nsclc_clinical
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.panels import PROLIFERATION_PANEL
from mirna_tcga.screen import cox_score_screen
from mirna_tcga.survival import coerce_clinical

DEFAULT_GENES = ["STC1", "LOXL2", "ANGPT2", "MKI67"]   # MKI67 = proliferation control


def proliferation_tertile(X: pd.DataFrame, subtype: pd.Series) -> pd.Series:
    panel = [g for g in PROLIFERATION_PANEL if g in X.columns]
    z = (X[panel] - X[panel].mean()) / X[panel].std()
    score = z.mean(axis=1)
    tert = pd.Series(index=X.index, dtype=object)
    for s in subtype.dropna().unique():
        m = subtype == s
        q1, q2 = score[m].quantile(1 / 3), score[m].quantile(2 / 3)
        tert[m] = np.where(score[m] <= q1, f"{s}_lo",
                           np.where(score[m] >= q2, f"{s}_hi", f"{s}_mid"))
    return tert


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--genes", nargs="+", default=DEFAULT_GENES)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clin).set_index("patientId")[["duration", "event", "cohort"]]
    print(f"NSCLC OS: {len(surv)} patients, {int(surv['event'].sum())} deaths")

    genes = list(dict.fromkeys(args.genes + PROLIFERATION_PANEL))
    expr = combined_expression(client, cfg, genes, keys)
    X = expr.T
    X["patientId"] = sample_to_patient(X.index)
    X = X[X["patientId"].isin(surv.index)]
    meta = surv.loc[X["patientId"]]
    dur = pd.Series(meta["duration"].to_numpy(), index=X.index)
    evt = pd.Series(meta["event"].to_numpy(), index=X.index)
    subtype = pd.Series(meta["cohort"].to_numpy(), index=X.index)
    pstrata = proliferation_tertile(X.drop(columns="patientId"), subtype)

    Xg = X[[g for g in args.genes if g in X.columns]].dropna(axis=1)
    r_sub = cox_score_screen(Xg, dur, evt, subtype)
    r_pro = cox_score_screen(Xg, dur, evt, pstrata)

    print("\n(z > 0  =>  higher expression = WORSE overall survival)")
    print(f"{'gene':8s} {'subtype-adj z':>14s} {'p':>9s}   {'prolif-adj z':>13s} {'p':>9s}  note")
    rows = []
    for g in Xg.columns:
        a, b = r_sub.loc[g], r_pro.loc[g]
        note = "proliferation control" if g == "MKI67" else ""
        print(f"{g:8s} {a['z']:>14.2f} {a['p']:>9.4f}   {b['z']:>13.2f} {b['p']:>9.4f}  {note}")
        rows.append((g, a["z"], a["p"], b["z"], b["p"]))

    if args.save_dir:
        out = Path(args.save_dir); out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows, columns=["gene", "z_subtype", "p_subtype", "z_prolif", "p_prolif"]) \
          .set_index("gene").to_csv(out / "stc1_survival_check.csv")
        print(f"\nWrote -> {out}/stc1_survival_check.csv")


if __name__ == "__main__":
    main()
