"""Overall survival in NSCLC (LUAD + LUSC) stratified by a single expression feature.

Asks the direct clinical question the project is built around: *do NSCLC patients
whose tumour highly expresses a given gene (or miRNA) survive differently from
those with low expression?*

Chain:  expression (one feature, whole NSCLC cohort)
        -> high/low split at the median
        -> attach overall-survival clinical data
        -> Kaplan-Meier medians + log-rank + Cox PH.

Two data layers, selected by which flag you pass:
  --gene  SYMBOL   mRNA expression via cBioPortal (RSEM), e.g. --gene EGFR
  --mirna FEATURE  miRNA expression via UCSC Xena,   e.g. --mirna hsa-mir-21

Because LUAD and LUSC have different baseline expression, ``--by-study-median``
splits high/low *within* each study before pooling, controlling for subtype.

Requires the survival extra:  pip install 'mirna-tcga[survival]'
Run:  python scripts/08_nsclc_expression_survival.py --gene EGFR
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.integrate import attach_clinical, stratify
from mirna_tcga.survival import coerce_clinical, cox_model, logrank_by_group
from mirna_tcga.xena import XenaClient


def mrna_expression(cfg, feature: str) -> pd.Series:
    """Per-sample mRNA expression of ``feature`` across the NSCLC cohort (cBioPortal)."""
    from mirna_tcga.cohorts import nsclc_expression

    client = CBioPortalClient(**cfg.cbioportal)
    mat = nsclc_expression(client, cfg, [feature])
    if mat.empty or feature not in mat.index:
        raise SystemExit(f"No mRNA expression returned for {feature!r} (check the symbol / egress).")
    return mat.loc[feature].astype(float).dropna()


def mirna_expression(cfg, feature: str) -> pd.Series:
    """Per-sample miRNA expression of ``feature`` across the NSCLC cohort (Xena)."""
    client = XenaClient(**cfg.xena)
    series = []
    for key in cohort_study_keys(cfg, "nsclc"):
        dataset = cfg.xena_mirna.get(key)
        if not dataset:
            continue
        mat = client.expression_matrix(dataset, [feature])
        if not mat.empty and feature in mat.index:
            series.append(mat.loc[feature])
    if not series:
        raise SystemExit(f"No miRNA expression returned for {feature!r} (check the id / egress).")
    return pd.concat(series).astype(float).dropna()


def split_groups(scores: pd.Series, patient_cohort: pd.Series, by_study: bool) -> pd.Series:
    """High/low label per sample, either pooled or within each study cohort."""
    if not by_study:
        return stratify(scores, method="median")
    parts = []
    for _, idx in patient_cohort.groupby(patient_cohort).groups.items():
        members = scores.index.intersection(idx)
        if len(members):
            parts.append(stratify(scores.loc[members], method="median"))
    return pd.concat(parts) if parts else stratify(scores, method="median")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=None)
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--gene", help="HUGO symbol for mRNA expression (cBioPortal)")
    src.add_argument("--mirna", help="miRNA feature id for Xena expression, e.g. hsa-mir-21")
    parser.add_argument("--by-study-median", action="store_true",
                        help="split high/low within LUAD and LUSC separately, then pool")
    parser.add_argument("--save", default=None, help="optional path to write the merged table as CSV")
    args = parser.parse_args()

    cfg = load_config(args.config)
    feature = args.gene or args.mirna
    layer = "mRNA" if args.gene else "miRNA"

    # 1. expression for the feature across the whole NSCLC cohort
    scores = mrna_expression(cfg, feature) if args.gene else mirna_expression(cfg, feature)
    print(f"{layer} feature {feature!r}: expression for {len(scores)} NSCLC samples")

    # 2. clinical outcomes (patient level), keyed to patient barcodes
    client = CBioPortalClient(**cfg.cbioportal)
    clinical = nsclc_clinical(client, cfg, patient_level=True)
    if clinical.empty:
        raise SystemExit("No clinical data returned (check network egress / study ids).")

    merged = attach_clinical(scores, clinical)
    # 3. high/low split (optionally within each subtype), then survival prep
    cohort_by_sample = merged.set_index("sample")["cohort"]
    groups = split_groups(scores.loc[merged["sample"]], cohort_by_sample, args.by_study_median)
    merged["group"] = merged["sample"].map(groups)
    merged = coerce_clinical(merged)
    print(f"Samples with overall-survival data: {len(merged)}")
    print(merged.groupby(["cohort", "group"]).size().rename("n"), "\n")

    if merged["group"].nunique() < 2 or merged.empty:
        raise SystemExit("Not enough data to compare high vs low groups.")

    # 4a. Kaplan-Meier median survival per group (descriptive)
    from lifelines import KaplanMeierFitter

    print(f"Median overall survival (months) by {feature} expression group:")
    for grp, part in merged.groupby("group"):
        km = KaplanMeierFitter().fit(part["duration"], part["event"])
        print(f"  {grp:>4}: n={len(part):3d}  median={km.median_survival_time_:.1f}")
    print()

    # 4b. log-rank test high vs low
    lr = logrank_by_group(merged, "group")
    print(f"Log-rank (high vs low {feature} {layer} expression):")
    lr.print_summary()

    # 4c. Cox PH on continuous, z-scored expression (+ subtype as covariate)
    merged["expr_z"] = (merged["score"] - merged["score"].mean()) / merged["score"].std()
    merged["is_lusc"] = (merged["cohort"] == "LUSC").astype(int)
    cph = cox_model(merged, ["expr_z", "is_lusc"])
    print(f"\nCox PH: hazard vs continuous {feature} expression (adjusted for subtype):")
    cph.print_summary()

    if args.save:
        merged.to_csv(args.save, index=False)
        print(f"\nWrote {len(merged)} rows -> {args.save}")


if __name__ == "__main__":
    main()
