"""Fast, genome-wide survival screening via a vectorized Cox score test.

Fitting a full Cox model per gene (e.g. with lifelines) across ~20k genes is far
too slow for a genome-wide screen. For a *univariate* screen we only need the
Cox partial-likelihood **score test** statistic evaluated at ``beta = 0`` -- the
same statistic as the log-rank test -- which has a closed form and can be
computed for every gene at once with cumulative sums over the risk sets.

The test is optionally **stratified** by a categorical confounder (e.g. tumour
subtype LUAD vs LUSC): score and variance contributions are accumulated within
each stratum and summed, exactly as a stratified Cox model does. Breslow's
approximation is used for tied event times.

The signed z-score (``z > 0`` => higher expression associates with higher hazard
/ worse survival) matches lifelines' Wald z closely. Use this to rank and
FDR-filter genes, then refit a full Cox model on the survivors
(:func:`mirna_tcga.survival.cox_model`) for reportable hazard ratios.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR q-values for a 1-D array of p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def _score_stat_one_stratum(X: np.ndarray, time: np.ndarray, event: np.ndarray):
    """Return score ``U`` and variance ``V`` (per gene) for one stratum."""
    order = np.argsort(time, kind="mergesort")
    Xo = X[order]
    to = time[order]
    eo = event[order]
    n, g = Xo.shape

    # Suffix sums == risk-set sums: S1[i] = sum_{j>=i} Xo[j] (samples with time >= to[i]).
    S1 = np.cumsum(Xo[::-1], axis=0)[::-1]
    S2 = np.cumsum((Xo[::-1]) ** 2, axis=0)[::-1]
    cnt = np.arange(n, 0, -1)

    ev_idx = np.where(eo == 1)[0]
    if ev_idx.size == 0:
        return np.zeros(g), np.zeros(g)

    ev_t = to[ev_idx]
    uniq, first = np.unique(ev_t, return_index=True)
    lo = np.searchsorted(to, uniq, side="left")  # first index of the risk set per event time

    meanR = S1[lo] / cnt[lo][:, None]
    meanR2 = S2[lo] / cnt[lo][:, None]
    varR = meanR2 - meanR**2

    nd = np.diff(np.append(first, ev_idx.size))  # tied deaths per unique event time
    sumX = np.add.reduceat(Xo[ev_idx], first, axis=0)  # sum of X over deaths per event time

    U = (sumX - nd[:, None] * meanR).sum(0)
    V = (nd[:, None] * varR).sum(0)
    return U, V


def cox_score_screen(
    X: pd.DataFrame,
    durations: pd.Series,
    events: pd.Series,
    strata: pd.Series | None = None,
) -> pd.DataFrame:
    """Univariate Cox score test for every column (gene) of ``X``.

    Parameters
    ----------
    X:
        samples x genes expression matrix (rows aligned to the survival data).
    durations:
        follow-up time per sample.
    events:
        1 = event/death, 0 = censored.
    strata:
        optional per-sample categorical label; the test is stratified by it.

    Returns
    -------
    DataFrame indexed by gene with columns ``z`` (signed: >0 = higher expression
    -> worse survival), ``p`` (two-sided), and ``q`` (Benjamini-Hochberg FDR),
    sorted by ascending ``p``. Genes with no informative variance are dropped.
    """
    genes = list(X.columns)
    Xv = X.to_numpy(dtype=float)
    time = durations.to_numpy(dtype=float)
    event = events.to_numpy(dtype=int)
    if strata is None:
        labels = np.zeros(len(time))
    else:
        labels = strata.to_numpy()

    U = np.zeros(Xv.shape[1])
    V = np.zeros(Xv.shape[1])
    for st in np.unique(labels):
        m = labels == st
        if m.sum() < 2:
            continue
        u, v = _score_stat_one_stratum(Xv[m], time[m], event[m])
        U += u
        V += v

    with np.errstate(divide="ignore", invalid="ignore"):
        z = U / np.sqrt(V)
    out = pd.DataFrame({"gene": genes, "z": z})
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["z"])
    out["p"] = 2 * norm.sf(np.abs(out["z"].to_numpy()))
    out["q"] = benjamini_hochberg(out["p"].to_numpy())
    return out.sort_values("p").set_index("gene")
