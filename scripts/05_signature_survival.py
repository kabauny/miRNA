"""Wire a PAM signature to clinical survival.

Full chain: expression -> PAM subtype signature -> per-sample signature score
-> high/low stratification -> overall-survival log-rank + Cox (within LUAD).

This is illustrative: it asks whether a "squamous-like" expression score (from
the LUAD-vs-LUSC signature) is associated with survival among adenocarcinomas.
Swap in any feature set / cohort to ask a different question.

Requires the survival extra:  pip install 'mirna-tcga[survival]'
Run:  python scripts/05_signature_survival.py
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.classify import fit_pam_cv, signature
from mirna_tcga.integrate import attach_clinical, signature_score, stratify
from mirna_tcga.panels import LUNG_MARKER_PANEL
from mirna_tcga.preprocess import drop_near_zero_variance, standardize
from mirna_tcga.survival import coerce_clinical, cox_model, logrank_by_group


def subtype_frame(client: CBioPortalClient, cfg) -> pd.DataFrame:
    frames = []
    for key, label in (("luad", "LUAD"), ("lusc", "LUSC")):
        mat = client.expression_matrix(
            cfg.mrna_profile(key), LUNG_MARKER_PANEL, cfg.all_samples_list(key)
        )
        s = mat.T
        s["label"] = label
        frames.append(s)
    return pd.concat(frames, axis=0).dropna(axis=1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--study", default="luad", help="cohort to test survival in")
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)

    # 1. PAM subtype signature (LUAD vs LUSC)
    df = subtype_frame(client, cfg)
    y = df.pop("label")
    X = standardize(drop_near_zero_variance(df))
    result = fit_pam_cv(X, y, cv_folds=cfg.analysis["cv_folds"],
                        random_state=cfg.analysis["random_state"])
    sig = signature(result, list(X.columns))
    print(f"Signature size: {len(sig)} genes (CV error {result.cv_error:.3f})")

    # 2. score the target cohort by the signature (weighted toward LUSC)
    weights = sig.set_index("gene")["score.LUSC"]
    cohort_samples = y.index[y == args.study.upper()]
    X_cohort = X.loc[X.index.intersection(cohort_samples)]
    scores = signature_score(X_cohort, sig["gene"].tolist(), weights=weights)

    # 3. stratify and attach clinical outcome
    groups = stratify(scores, method="median")
    clinical = client.get_clinical_data(cfg.studies[args.study], patient_level=True)
    merged = attach_clinical(scores, clinical)
    merged["group"] = merged["sample"].map(groups)
    merged = coerce_clinical(merged)
    print(f"Samples with outcome data: {len(merged)}")
    print(merged["group"].value_counts(), "\n")

    # 4. survival tests
    lr = logrank_by_group(merged, "group")
    print("Log-rank (high vs low signature score):")
    lr.print_summary()

    merged["score_z"] = (merged["score"] - merged["score"].mean()) / merged["score"].std()
    cph = cox_model(merged, ["score_z"])
    print("\nCox PH on continuous signature score:")
    cph.print_summary()


if __name__ == "__main__":
    main()
