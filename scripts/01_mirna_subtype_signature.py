"""LUAD vs LUSC *miRNA* signature via nearest shrunken centroids (PAM).

This is the miRNA counterpart to ``02_subtype_signature.py`` (which uses mRNA
from cBioPortal). Since cBioPortal lacks miRNA expression for the PanCancer
studies, miRNA comes from UCSC Xena's TCGA hub instead -- everything downstream
(preprocess + PAM) is identical.

Run:  python scripts/01_mirna_subtype_signature.py
Note: Xena miRNA matrices are already log2-normalized, so Box-Cox is skipped.
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.classify import fit_pam_cv, signature
from mirna_tcga.preprocess import drop_near_zero_variance, standardize
from mirna_tcga.xena import XenaClient


def build_mirna_frame(client: XenaClient, cfg) -> pd.DataFrame:
    """Samples x miRNA frame with a `label` column (LUAD / LUSC)."""
    frames = []
    for key, label in (("luad", "LUAD"), ("lusc", "LUSC")):
        mat = client.expression_matrix(cfg.xena_mirna[key])  # miRNAs x samples
        samples = mat.T                                       # samples x miRNAs
        samples["label"] = label
        frames.append(samples)
    # outer-align on the union of miRNAs across cohorts
    return pd.concat(frames, axis=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = XenaClient(**cfg.xena)

    df = build_mirna_frame(client, cfg)
    df = df.dropna(axis=1)            # keep miRNAs measured in both cohorts
    y = df.pop("label")

    X = standardize(drop_near_zero_variance(df))
    result = fit_pam_cv(
        X, y,
        cv_folds=cfg.analysis["cv_folds"],
        random_state=cfg.analysis["random_state"],
    )
    print(f"miRNAs analyzed: {X.shape[1]}  |  samples: {X.shape[0]}")
    print(f"Best shrink threshold: {result.best_threshold:.3f}")
    print(f"CV error: {result.cv_error:.3f}\n")
    print("Confusion (cross-validated):")
    print(result.confusion, "\n")
    print("miRNA subtype signature (top):")
    print(signature(result, list(X.columns)).head(20).to_string(index=False))


if __name__ == "__main__":
    main()
