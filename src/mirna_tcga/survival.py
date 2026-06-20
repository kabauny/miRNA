"""Survival / clinical-outcome helpers (Kaplan-Meier + Cox PH).

This is the natural extension the original R pipeline never had: linking a
signature, stage, or mutation status to patient outcome. Built on **lifelines**
(optional extra: ``pip install 'mirna-tcga[survival]'``).

Expects a tidy clinical table with a duration column (e.g. months) and a binary
event column (1 = event/death, 0 = censored).
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd


def _require_lifelines():
    try:
        import lifelines  # noqa: F401
    except ImportError as exc:  # pragma: no cover - import-guard
        raise ImportError(
            "'lifelines' is required for survival analysis. Install with: "
            "pip install 'mirna-tcga[survival]'"
        ) from exc
    return __import__("lifelines")


def coerce_clinical(
    df: pd.DataFrame,
    duration_col: str = "OS_MONTHS",
    event_col: str = "OS_STATUS",
) -> pd.DataFrame:
    """Return a clean frame with numeric ``duration`` and 0/1 ``event`` columns.

    cBioPortal encodes status as e.g. ``"1:DECEASED"`` / ``"0:LIVING"``; the
    leading integer is used as the event indicator.
    """
    out = df.copy()
    out["duration"] = pd.to_numeric(out[duration_col], errors="coerce")
    status = out[event_col].astype(str).str.split(":").str[0]
    out["event"] = pd.to_numeric(status, errors="coerce")
    return out.dropna(subset=["duration", "event"])


def logrank_by_group(
    df: pd.DataFrame,
    group_col: str,
    duration_col: str = "duration",
    event_col: str = "event",
):
    """Multivariate log-rank test comparing survival across ``group_col``."""
    lifelines = _require_lifelines()
    from lifelines.statistics import multivariate_logrank_test

    return multivariate_logrank_test(
        df[duration_col], df[group_col], df[event_col]
    )


def cox_model(
    df: pd.DataFrame,
    covariates: Sequence[str],
    duration_col: str = "duration",
    event_col: str = "event",
):
    """Fit and return a Cox proportional-hazards model on ``covariates``."""
    lifelines = _require_lifelines()
    cph = lifelines.CoxPHFitter()
    cols = list(covariates) + [duration_col, event_col]
    cph.fit(df[cols].dropna(), duration_col=duration_col, event_col=event_col)
    return cph
