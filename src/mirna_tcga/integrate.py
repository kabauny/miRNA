"""Wire signatures, expression, and clinical data together.

Bridges the modelling output (a PAM signature) to per-sample scores, patient
stratification, clinical survival, and cross-signature overlap -- the
integrative step the original R project gestured at (``overall.R``) but never
completed.

Conventions:
* ``X`` is a samples x features :class:`pandas.DataFrame` (ideally standardized).
* a "signature" is the DataFrame returned by :func:`mirna_tcga.classify.signature`
  (a ``gene`` column plus per-class ``score.<class>`` columns).
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd


def signature_score(
    X: pd.DataFrame,
    genes: Sequence[str],
    weights: pd.Series | None = None,
) -> pd.Series:
    """Per-sample score summarizing a signature's expression.

    With no weights, the score is the mean (standardized) expression across the
    signature features present in ``X``. With ``weights`` (e.g. a class centroid
    deviation), it is a sign-aware weighted mean, so a high score means "more
    like" the weighted class.
    """
    cols = [g for g in genes if g in X.columns]
    if not cols:
        raise ValueError("None of the signature features are present in X.")
    sub = X[cols]
    if weights is None:
        return sub.mean(axis=1)
    w = weights.reindex(cols).fillna(0.0)
    denom = w.abs().sum()
    if denom == 0:
        return sub.mean(axis=1)
    return sub.mul(w, axis=1).sum(axis=1) / denom


def stratify(scores: pd.Series, method: str = "median", threshold: float | None = None) -> pd.Series:
    """Split a continuous score into "high"/"low" groups.

    ``method`` is "median" (default) or "mean"; pass an explicit ``threshold``
    to override. Samples equal to the cut-point fall in "low".
    """
    if threshold is None:
        threshold = scores.median() if method == "median" else scores.mean()
    return scores.apply(lambda v: "high" if v > threshold else "low").rename("group")


def overlap_signatures(*signatures: pd.DataFrame, key: str = "gene") -> list[str]:
    """Features common to all provided signatures (intersection on ``key``)."""
    if not signatures:
        return []
    common = set(signatures[0][key])
    for sig in signatures[1:]:
        common &= set(sig[key])
    return sorted(common)


def sample_to_patient(sample_ids: Sequence[str]) -> list[str]:
    """Derive TCGA patient barcodes from sample barcodes.

    ``TCGA-05-4244-01`` (sample) -> ``TCGA-05-4244`` (patient). Ids that do not
    look like TCGA barcodes are returned unchanged.
    """
    out = []
    for s in sample_ids:
        parts = str(s).split("-")
        out.append("-".join(parts[:3]) if len(parts) >= 3 else str(s))
    return out


def attach_clinical(
    scores: pd.Series,
    clinical: pd.DataFrame,
    patient_col: str = "patientId",
) -> pd.DataFrame:
    """Join per-sample scores to a patient-level clinical table.

    Sample ids in ``scores`` are collapsed to patient barcodes and merged with
    ``clinical``. Returns one row per sample that matched a patient.
    """
    df = pd.DataFrame({"sample": scores.index, "score": scores.to_numpy()})
    df["patientId"] = sample_to_patient(df["sample"])
    merged = df.merge(clinical, left_on="patientId", right_on=patient_col, how="inner")
    return merged
