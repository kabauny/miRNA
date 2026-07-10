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


class _ExpressionSource(Protocol):
    def expression_matrix(
        self, molecular_profile_id: str, hugo_symbols: Sequence[str], sample_list_id: str
    ) -> pd.DataFrame: ...


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


def combined_expression(
    client: _ExpressionSource,
    cfg,
    hugo_symbols: Sequence[str],
    study_keys: Sequence[str],
) -> pd.DataFrame:
    """Concatenate per-study mRNA expression matrices into one genes x samples frame.

    Each study's expression is fetched via ``client.expression_matrix`` using that
    study's mRNA profile and ``_all`` sample list, then joined column-wise on the
    shared gene index. Studies returning no data are skipped; an empty frame is
    returned if none yield data. Sample ids are globally unique across TCGA
    studies, so no de-duplication is needed.
    """
    frames = []
    for key in study_keys:
        mat = client.expression_matrix(
            cfg.mrna_profile(key), list(hugo_symbols), cfg.all_samples_list(key)
        )
        if mat is None or mat.empty:
            continue
        frames.append(mat)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=1)


def nsclc_expression(client: _ExpressionSource, cfg, hugo_symbols: Sequence[str]) -> pd.DataFrame:
    """Combined mRNA expression (genes x samples) across the NSCLC cohort."""
    return combined_expression(client, cfg, hugo_symbols, cohort_study_keys(cfg, "nsclc"))
