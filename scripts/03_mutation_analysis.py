"""Mutant-vs-wildtype expression signature (default: TP53 in LUAD).

Python port of ``mutationAnalysis.R``: split tumours by mutation status in a
driver gene, then fit a PAM model to find the expression signature separating
mutant from wild-type samples.

Run:  python scripts/03_mutation_analysis.py --gene TP53 --study luad
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.classify import fit_pam_cv, signature
from mirna_tcga.panels import LUNG_MARKER_PANEL as PANEL
from mirna_tcga.preprocess import drop_near_zero_variance, standardize


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--study", default="luad", help="study key from config (luad/lusc)")
    parser.add_argument("--gene", default=None, help="mutation gene (default from config)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    gene = args.gene or cfg.analysis["mutation_gene"]
    client = CBioPortalClient(**cfg.cbioportal)

    sample_list = cfg.all_samples_list(args.study)
    mutated = client.mutated_samples(
        cfg.mutation_profile(args.study), gene, sample_list
    )
    print(f"{gene}-mutant samples in {args.study.upper()}: {len(mutated)}")

    mat = client.expression_matrix(cfg.mrna_profile(args.study), PANEL, sample_list)
    X = mat.T.dropna(axis=1)
    y = X.index.to_series().apply(lambda s: "mutant" if s in mutated else "wildtype")

    X = standardize(drop_near_zero_variance(X))
    result = fit_pam_cv(
        X, y,
        cv_folds=cfg.analysis["cv_folds"],
        random_state=cfg.analysis["random_state"],
    )
    print(f"\nBest shrink threshold: {result.best_threshold:.3f}")
    print(f"CV error: {result.cv_error:.3f}\n")
    print(result.confusion, "\n")
    print(f"{gene} mutation-associated signature:")
    print(signature(result, list(X.columns)).head(20).to_string(index=False))


if __name__ == "__main__":
    main()
