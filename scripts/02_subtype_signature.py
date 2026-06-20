"""LUAD vs LUSC miRNA/mRNA signature via nearest shrunken centroids (PAM).

Python port of the original ``acgts``/``centroid.R`` idea: pull expression for
both lung subtypes, preprocess, fit a CV-tuned PAM model, and list the genes
that distinguish the subtypes.

Run:  python scripts/02_subtype_signature.py
(uses a small marker-gene panel by default so it runs quickly; pass --all-genes
to use the full coding-gene universe -- slower.)
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.classify import fit_pam_cv, signature
from mirna_tcga.panels import LUNG_MARKER_PANEL
from mirna_tcga.preprocess import drop_near_zero_variance, standardize


def build_subtype_frame(client: CBioPortalClient, cfg, genes: list[str]) -> pd.DataFrame:
    """Samples x features frame with a `label` column (LUAD / LUSC)."""
    frames = []
    for key, label in (("luad", "LUAD"), ("lusc", "LUSC")):
        mat = client.expression_matrix(
            cfg.mrna_profile(key), genes, cfg.all_samples_list(key)
        )
        samples = mat.T  # samples x genes
        samples["label"] = label
        frames.append(samples)
    return pd.concat(frames, axis=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None, help="path to config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)

    df = build_subtype_frame(client, cfg, LUNG_MARKER_PANEL)
    df = df.dropna(axis=1)  # drop genes missing in either cohort
    y = df.pop("label")

    X = drop_near_zero_variance(df)
    X = standardize(X)

    result = fit_pam_cv(
        X, y,
        cv_folds=cfg.analysis["cv_folds"],
        random_state=cfg.analysis["random_state"],
    )
    print(f"Best shrink threshold: {result.best_threshold:.3f}")
    print(f"CV error: {result.cv_error:.3f}\n")
    print("Confusion (cross-validated):")
    print(result.confusion, "\n")
    print("Subtype signature (top genes):")
    print(signature(result, list(X.columns)).head(20).to_string(index=False))


if __name__ == "__main__":
    main()
