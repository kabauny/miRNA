"""Multi-study cohort helpers (e.g. NSCLC = LUAD + LUSC).

Several study-level tables need to be combined into one analysis cohort. These
helpers concatenate per-study cBioPortal results and tag each row with its
source study via a ``cohort`` column.

They accept any object exposing the relevant ``CBioPortalClient`` methods, which
keeps them trivially testable with a fake client.
"""

from __future__ import annotations

from typing import Protocol, Sequence

import pandas as pd


class _ClinicalSource(Protocol):
    def get_clinical_data(self, study_id: str, patient_level: bool = True) -> pd.DataFrame: ...


def cohort_study_keys(cfg, cohort: str = "nsclc") -> list[str]:
    """Study keys for a named cohort, defaulting NSCLC to LUAD + LUSC."""
    default = ["luad", "lusc"] if cohort == "nsclc" else []
    return list(cfg.cohorts.get(cohort, default))


def combined_clinical(
    client: _ClinicalSource,
    cfg,
    study_keys: Sequence[str],
    patient_level: bool = True,
) -> pd.DataFrame:
    """Concatenate clinical tables across studies, tagging each with ``cohort``.

    Studies that return no data are skipped. Returns an empty frame if none of
    the requested studies yield clinical data.
    """
    frames = []
    for key in study_keys:
        df = client.get_clinical_data(cfg.studies[key], patient_level=patient_level)
        if df is None or df.empty:
            continue
        df = df.copy()
        df["cohort"] = key.upper()
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def nsclc_clinical(client: _ClinicalSource, cfg, patient_level: bool = True) -> pd.DataFrame:
    """Combined patient-level clinical table for NSCLC (LUAD + LUSC)."""
    return combined_clinical(client, cfg, cohort_study_keys(cfg, "nsclc"), patient_level)
