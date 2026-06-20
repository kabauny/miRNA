"""Nearest shrunken centroid classification (PAM) via scikit-learn.

This is the Python equivalent of the original ``pamr`` workflow:

* ``pamr.train`` / ``pamr.cv``  -> :class:`sklearn.neighbors.NearestCentroid`
  with ``shrink_threshold``, tuned by cross-validation.
* ``pamr.listgenes``           -> :func:`signature`, the features whose
  class centroids still differ from the overall centroid after shrinkage.

Convention: ``X`` is samples x features, ``y`` is the label vector.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.neighbors import NearestCentroid


@dataclass
class PamResult:
    """Outcome of a cross-validated PAM fit."""

    model: NearestCentroid
    best_threshold: float
    cv_error: float
    cv_errors: pd.Series          # error rate per candidate threshold
    confusion: pd.DataFrame       # cross-validated confusion matrix
    labels: list[str]


def _default_thresholds(X: pd.DataFrame, n: int = 30) -> np.ndarray:
    """A geometric grid of shrinkage thresholds spanning 0..~max feature std."""
    max_std = float(np.nanmax(X.std(axis=0, ddof=1).to_numpy()))
    if not np.isfinite(max_std) or max_std <= 0:
        max_std = 1.0
    return np.linspace(0.0, max_std, n)


def fit_pam_cv(
    X: pd.DataFrame,
    y: pd.Series | np.ndarray,
    thresholds: np.ndarray | None = None,
    cv_folds: int = 5,
    random_state: int = 42,
) -> PamResult:
    """Fit nearest shrunken centroids, choosing the shrink threshold by CV.

    The threshold minimizing cross-validated misclassification error is kept
    (ties broken toward the *largest* threshold, i.e. the sparsest model, as
    PAM conventionally does), then the model is refit on all data.
    """
    y = np.asarray(y)
    if thresholds is None:
        thresholds = _default_thresholds(X)

    cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
    errors = {}
    for t in thresholds:
        model = NearestCentroid(shrink_threshold=t if t > 0 else None)
        try:
            pred = cross_val_predict(model, X.to_numpy(), y, cv=cv)
        except ValueError:
            # A threshold so large it removes every feature -> skip.
            continue
        errors[float(t)] = float(np.mean(pred != y))

    if not errors:
        raise RuntimeError("No threshold produced a valid model; check the input data.")

    cv_errors = pd.Series(errors).sort_index()
    min_err = cv_errors.min()
    # largest threshold achieving the minimum error (sparsest)
    best_threshold = float(cv_errors[cv_errors == min_err].index.max())

    final = NearestCentroid(shrink_threshold=best_threshold if best_threshold > 0 else None)
    final.fit(X.to_numpy(), y)

    pred = cross_val_predict(
        NearestCentroid(shrink_threshold=best_threshold if best_threshold > 0 else None),
        X.to_numpy(), y, cv=cv,
    )
    labels = sorted(np.unique(y).tolist())
    conf = pd.DataFrame(
        confusion_matrix(y, pred, labels=labels),
        index=pd.Index(labels, name="true"),
        columns=pd.Index(labels, name="predicted"),
    )

    return PamResult(
        model=final,
        best_threshold=best_threshold,
        cv_error=float(min_err),
        cv_errors=cv_errors,
        confusion=conf,
        labels=labels,
    )


def signature(result: PamResult, feature_names: list[str]) -> pd.DataFrame:
    """Features retained by the shrunken-centroid model and their class scores.

    Equivalent to ``pamr.listgenes``: for each class we report the (shrunken)
    centroid's deviation from the overall centroid. Features where every class
    deviation is ~0 have been shrunk away and are excluded.
    """
    model = result.model
    centroids = np.asarray(model.centroids_)          # (n_classes, n_features)
    overall = centroids.mean(axis=0)
    deviation = centroids - overall                   # (n_classes, n_features)

    kept = np.where(np.any(np.abs(deviation) > 1e-12, axis=0))[0]
    cols = {"gene": [feature_names[i] for i in kept]}
    for ci, cls in enumerate(model.classes_):
        cols[f"score.{cls}"] = deviation[ci, kept]
    out = pd.DataFrame(cols)
    # rank by overall importance (max absolute class score)
    out["abs_max"] = out[[c for c in out.columns if c.startswith("score.")]].abs().max(axis=1)
    return out.sort_values("abs_max", ascending=False).drop(columns="abs_max").reset_index(drop=True)
