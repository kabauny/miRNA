"""Overlap between two PAM signatures (the modern 'overall.R').

Computes two signatures and reports the features they share:

  * subtype signature  : LUAD vs LUSC
  * mutation signature : TP53-mutant vs wild-type (within --study)

Both are in HUGO gene-symbol space (mRNA), so the intersection is direct.

Cross-modality (miRNA vs mRNA) overlap is NOT a plain id intersection -- miRNAs
must first be mapped to target genes (see mirna_tcga.idmap.mirbase_to_ensembl).
This script demonstrates the same-id-space case end-to-end.

Run:  python scripts/06_signature_overlap.py --study luad
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.classify import fit_pam_cv, signature
from mirna_tcga.integrate import overlap_signatures
from mirna_tcga.panels import LUNG_MARKER_PANEL
from mirna_tcga.preprocess import drop_near_zero_variance, standardize


def _pam_signature(X: pd.DataFrame, y: pd.Series, cfg) -> pd.DataFrame:
    X = standardize(drop_near_zero_variance(X))
    result = fit_pam_cv(X, y, cv_folds=cfg.analysis["cv_folds"],
                        random_state=cfg.analysis["random_state"])
    return signature(result, list(X.columns))


def subtype_signature(client: CBioPortalClient, cfg) -> pd.DataFrame:
    frames = []
    for key, label in (("luad", "LUAD"), ("lusc", "LUSC")):
        mat = client.expression_matrix(
            cfg.mrna_profile(key), LUNG_MARKER_PANEL, cfg.all_samples_list(key)
        )
        s = mat.T
        s["label"] = label
        frames.append(s)
    df = pd.concat(frames, axis=0).dropna(axis=1)
    return _pam_signature(df.drop(columns="label"), df["label"], cfg)


def mutation_signature(client: CBioPortalClient, cfg, study: str, gene: str) -> pd.DataFrame:
    sample_list = cfg.all_samples_list(study)
    mutated = client.mutated_samples(cfg.mutation_profile(study), gene, sample_list)
    mat = client.expression_matrix(cfg.mrna_profile(study), LUNG_MARKER_PANEL, sample_list)
    X = mat.T.dropna(axis=1)
    y = X.index.to_series().apply(lambda s: "mutant" if s in mutated else "wildtype")
    return _pam_signature(X, y, cfg)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--study", default="luad")
    parser.add_argument("--gene", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    gene = args.gene or cfg.analysis["mutation_gene"]
    client = CBioPortalClient(**cfg.cbioportal)

    sub = subtype_signature(client, cfg)
    mut = mutation_signature(client, cfg, args.study, gene)
    shared = overlap_signatures(sub, mut)

    print(f"Subtype signature genes:  {len(sub)}")
    print(f"{gene} mutation signature genes: {len(mut)}")
    print(f"\nShared ({len(shared)}): {', '.join(shared) if shared else '(none)'}")


if __name__ == "__main__":
    main()
