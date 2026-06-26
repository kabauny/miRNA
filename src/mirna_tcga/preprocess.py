"""Feature preprocessing: variance filtering, Box-Cox normalization, scaling.

These mirror the steps in the original R ``preprocessing`` script
(near-zero-variance removal, per-feature Box-Cox, standardization) but use
vectorized scipy/sklearn implementations.

Convention: input ``X`` is a samples x features :class:`pandas.DataFrame`
(the transpose of the genes x samples matrix returned by the API client).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def near_zero_variance(X: pd.DataFrame, freq_cut: float = 20.0, unique_cut: float = 10.0) -> list[str]:
    """Return columns flagged as near-zero-variance (caret-style heuristic).

    A column is flagged when the ratio of the most- to second-most-common value
    exceeds ``freq_cut`` *and* the percentage of distinct values is below
    ``unique_cut``. Constant columns are always flagged. This reproduces the
    intent of the original hand-rolled ``near_zero_var`` function.
    """
    flagged: list[str] = []
    n = len(X)
    for col in X.columns:
        values = X[col].to_numpy()
        counts = pd.Series(values).value_counts()
        if len(counts) <= 1:  # constant
            flagged.append(col)
            continue
        freq_ratio = counts.iloc[0] / counts.iloc[1]
        pct_unique = 100.0 * len(counts) / n
        if freq_ratio > freq_cut and pct_unique < unique_cut:
            flagged.append(col)
    return flagged


def drop_near_zero_variance(X: pd.DataFrame, **kwargs) -> pd.DataFrame:
    """Return ``X`` with near-zero-variance columns removed."""
    flagged = near_zero_variance(X, **kwargs)
    return X.drop(columns=flagged)


def boxcox_normalize(X: pd.DataFrame, offset: float = 1.0) -> pd.DataFrame:
    """Per-feature Box-Cox transform (lambda fit by MLE per column).

    A positive ``offset`` is added first because Box-Cox requires strictly
    positive input (counts contain zeros). Columns that are constant after the
    shift are returned unchanged.
    """
    out = {}
    for col in X.columns:
        shifted = X[col].to_numpy(dtype=float) + offset
        if np.allclose(shifted, shifted[0]):
            out[col] = shifted - offset
            continue
        transformed, _ = stats.boxcox(shifted)
        out[col] = transformed
    return pd.DataFrame(out, index=X.index)


def standardize(X: pd.DataFrame) -> pd.DataFrame:
    """Z-score each feature (mean 0, sd 1). Zero-variance columns become 0."""
    mean = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1).replace(0.0, 1.0)
    return (X - mean) / sd


def skewed_features(X: pd.DataFrame, kurtosis_cut: float = 2.0) -> list[str]:
    """Columns whose absolute excess kurtosis exceeds ``kurtosis_cut``.

    Mirrors the original ``skew()`` helper (which used ``e1071::kurtosis``).
    """
    return [col for col in X.columns if abs(stats.kurtosis(X[col], fisher=True)) >= kurtosis_cut]
