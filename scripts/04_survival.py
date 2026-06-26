"""Overall survival by tumour stage (the clinical angle the R code lacked).

Pulls patient clinical data from cBioPortal, groups by AJCC stage, and runs a
log-rank test. Requires the survival extra:  pip install 'mirna-tcga[survival]'

Run:  python scripts/04_survival.py --study luad
"""

from __future__ import annotations

import argparse

import _bootstrap  # noqa: F401

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.survival import coerce_clinical, logrank_by_group

STAGE_COL_CANDIDATES = ["AJCC_PATHOLOGIC_TUMOR_STAGE", "PATH_STAGE", "TUMOR_STAGE"]


def _stage_group(stage: str) -> str:
    """Collapse 'Stage IIA' etc. to coarse I/II/III/IV groups."""
    s = str(stage).upper().replace("STAGE", "").strip()
    for roman in ("IV", "III", "II", "I"):
        if s.startswith(roman):
            return roman
    return "NA"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=None)
    parser.add_argument("--study", default="luad")
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)

    clinical = client.get_clinical_data(cfg.studies[args.study], patient_level=True)
    stage_col = next((c for c in STAGE_COL_CANDIDATES if c in clinical.columns), None)
    if stage_col is None:
        raise SystemExit(
            f"No stage column found. Available: {sorted(clinical.columns)}"
        )

    clinical = coerce_clinical(clinical)
    clinical["stage_group"] = clinical[stage_col].map(_stage_group)
    clinical = clinical[clinical["stage_group"] != "NA"]

    print(clinical["stage_group"].value_counts().sort_index(), "\n")
    res = logrank_by_group(clinical, "stage_group")
    res.print_summary()


if __name__ == "__main__":
    main()
