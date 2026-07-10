"""Fetch and summarize NSCLC (LUAD + LUSC) clinical data from cBioPortal.

Run:  python scripts/07_nsclc_clinical.py
(Needs network egress to www.cbioportal.org.)
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import nsclc_clinical

STAGE_COL_CANDIDATES = ["AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_STAGE", "TUMOR_STAGE"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--save", default=None, help="optional path to write the table as CSV")
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)

    df = nsclc_clinical(client, cfg)
    if df.empty:
        raise SystemExit("No clinical data returned (check network egress / study ids).")

    print(f"NSCLC patients: {len(df)}")
    print(df["cohort"].value_counts(), "\n")

    for col in STAGE_COL_CANDIDATES:
        if col in df.columns:
            print(f"Stage ({col}):")
            print(df[col].value_counts(dropna=False).sort_index(), "\n")
            break

    if "OS_STATUS" in df.columns:
        print("Overall survival status:")
        print(df["OS_STATUS"].value_counts(dropna=False), "\n")

    if args.save:
        df.to_csv(args.save, index=False)
        print(f"Wrote {len(df)} rows -> {args.save}")


if __name__ == "__main__":
    main()
