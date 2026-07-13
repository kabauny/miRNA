"""Association screens against a binary endpoint (e.g. metastasis).

Complements :mod:`mirna_tcga.screen` (which tests features against *survival*):

* :func:`ranksum_screen` -- continuous features (gene / miRNA expression) vs a
  binary outcome, via a subtype-**stratified** Wilcoxon rank-sum (van Elteren)
  test. Fully vectorized across features.
* :func:`fisher_screen` -- binary features (copy-number deletions, mutations) vs
  a binary outcome, via Fisher's exact test with odds ratios.

Both return a tidy DataFrame with a signed effect, p-value, and BH FDR ``q``.
A binary feature vs *survival* needs no new code -- pass the 0/1 matrix straight
to :func:`mirna_tcga.screen.cox_score_screen` (the score test is the log-rank).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2_dist
from scipy.stats import fisher_exact, norm, rankdata

from .screen import benjamini_hochberg


def ranksum_screen(
    X: pd.DataFrame,
    outcome: pd.Series,
    strata: pd.Series | None = None,
) -> pd.DataFrame:
    """Stratified Wilcoxon rank-sum test of every column of ``X`` vs ``outcome``.

    ``outcome`` is 0/1 (aligned to ``X``'s rows). Within each stratum the rank
    sum of the positive group is compared to its null; contributions are pooled
    across strata (van Elteren). ``z > 0`` means the feature is *higher* in the
    positive (outcome == 1) group. Returns columns ``z``, ``auc`` (P(feature
    higher in a positive than a negative), pooled), ``p``, ``q`` sorted by ``p``.
    """
    y = outcome.reindex(X.index).to_numpy()
    labels = np.zeros(len(y)) if strata is None else strata.reindex(X.index).to_numpy()
    Xv = X.to_numpy(dtype=float)

    num = np.zeros(Xv.shape[1])  # sum_k (R_k - E_k)
    var = np.zeros(Xv.shape[1])
    for st in np.unique(labels):
        m = labels == st
        yk = y[m]
        n1, n0 = int((yk == 1).sum()), int((yk == 0).sum())
        if n1 == 0 or n0 == 0:
            continue
        ranks = rankdata(Xv[m], axis=0, method="average")  # (n_k, g)
        N = n1 + n0
        rbar = (N + 1) / 2
        R1 = ranks[yk == 1].sum(0)
        E1 = n1 * rbar
        # tie-corrected variance via the finite-population formula on actual ranks
        var_k = n1 * n0 / (N * (N - 1)) * ((ranks**2).sum(0) - N * rbar**2)
        num += R1 - E1
        var += var_k

    with np.errstate(divide="ignore", invalid="ignore"):
        z = num / np.sqrt(var)

    # pooled AUC (Mann-Whitney U1 / (n1 n0)): P(feature higher in a positive)
    yb = y == 1
    n1, n0 = int(yb.sum()), int((y == 0).sum())
    ranks_all = rankdata(Xv, axis=0, method="average")
    auc = (ranks_all[yb].sum(0) - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else np.full(z.shape, np.nan)

    out = pd.DataFrame({"gene": list(X.columns), "z": z, "auc": auc})
    out = out.replace([np.inf, -np.inf], np.nan).dropna(subset=["z"])
    out["p"] = 2 * norm.sf(np.abs(out["z"].to_numpy()))
    out["q"] = benjamini_hochberg(out["p"].to_numpy())
    return out[["gene", "z", "auc", "p", "q"]].sort_values("p").set_index("gene")


def cmh_depletion_screen(
    B: pd.DataFrame,
    group: pd.Series,
    strata: pd.Series,
    min_ref_freq: float = 0.03,
) -> pd.DataFrame:
    """Subtype-adjusted test for a binary feature being *depleted* in a group.

    Motivating question: which genes are recurrently deep-deleted in the
    reference group (``group == 0``, e.g. non-metastatic) yet **spared** from
    deletion in the index group (``group == 1``, e.g. stage IV)? A gene whose
    intact copy is required for metastasis should show this depletion.

    Uses a Cochran-Mantel-Haenszel test stratified by ``strata`` (e.g. tumour
    subtype), which is essential when the index group is unbalanced across strata
    and is robust to the zero cells produced by "never deleted in the index
    group". Only genes altered in >= ``min_ref_freq`` of the reference group are
    tested (a gene must be *deletable* to be informatively spared).

    Returns per gene: ``freq_ref`` / ``freq_idx`` (deletion frequency in each
    group), ``n_idx`` (deletions in the index group), ``or_mh`` (MH odds ratio;
    < 1 = depleted in the index group), one-sided depletion ``p``, and BH ``q``,
    sorted by ``p``.
    """
    y = group.reindex(B.index).to_numpy()
    st = strata.reindex(B.index).to_numpy()
    Bv = B.to_numpy(dtype=float)
    g1, g0 = (y == 1), (y == 0)

    freq_ref = Bv[g0].mean(0)
    freq_idx = Bv[g1].mean(0) if g1.sum() else np.zeros(Bv.shape[1])
    n_idx = Bv[g1].sum(0)

    num = np.zeros(Bv.shape[1])   # sum_k (a - E[a])
    var = np.zeros(Bv.shape[1])
    or_num = np.zeros(Bv.shape[1])
    or_den = np.zeros(Bv.shape[1])
    for s in np.unique(st):
        m = st == s
        n1, n0 = int((m & g1).sum()), int((m & g0).sum())
        n = n1 + n0
        if n1 == 0 or n0 == 0 or n < 2:
            continue
        a = Bv[m & g1].sum(0)          # deleted & index
        b = Bv[m & g0].sum(0)          # deleted & reference
        c, d = n1 - a, n0 - b
        row1 = a + b                    # deleted (either group)
        num += a - row1 * n1 / n
        var += row1 * (c + d) * n1 * n0 / (n * n * (n - 1))
        or_num += a * d / n
        or_den += b * c / n

    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = (np.abs(num) - 0.5).clip(min=0) ** 2 / var
        p_two = chi2_dist.sf(chi2, 1)
        or_mh = or_num / or_den
    # one-sided depletion (fewer deletions in the index group than expected)
    p_dep = np.where(num < 0, p_two / 2, 1 - p_two / 2)

    out = pd.DataFrame({
        "gene": list(B.columns), "freq_ref": freq_ref, "freq_idx": freq_idx,
        "n_idx": n_idx.astype(int), "or_mh": or_mh, "p": p_dep,
    })
    out = out[out["freq_ref"] >= min_ref_freq].replace([np.inf, -np.inf], np.nan).dropna(subset=["p"])
    if out.empty:
        return out.set_index("gene")
    out["q"] = benjamini_hochberg(out["p"].to_numpy())
    return out.sort_values("p").set_index("gene")


def fisher_screen(B: pd.DataFrame, outcome: pd.Series) -> pd.DataFrame:
    """Fisher's exact test of every binary feature in ``B`` vs ``outcome``.

    ``B`` is a samples x features 0/1 matrix (e.g. deep-deletion or mutation
    flags); ``outcome`` is 0/1. Returns per-feature ``n_alt`` (feature-positive
    count), ``odds_ratio`` (feature vs outcome), ``p``, and BH ``q``, sorted by
    ``p``. A Haldane-Anscombe 0.5 correction keeps the odds ratio finite when a
    cell is zero.
    """
    y = outcome.reindex(B.index).to_numpy()
    rows = []
    for feat in B.columns:
        b = B[feat].to_numpy()
        a11 = int(((b == 1) & (y == 1)).sum())
        a10 = int(((b == 1) & (y == 0)).sum())
        a01 = int(((b == 0) & (y == 1)).sum())
        a00 = int(((b == 0) & (y == 0)).sum())
        _, p = fisher_exact([[a11, a10], [a01, a00]])
        orr = ((a11 + 0.5) * (a00 + 0.5)) / ((a10 + 0.5) * (a01 + 0.5))
        rows.append((feat, a11 + a10, orr, p))
    out = pd.DataFrame(rows, columns=["gene", "n_alt", "odds_ratio", "p"])
    if out.empty:
        return out.set_index("gene")
    out["q"] = benjamini_hochberg(out["p"].to_numpy())
    return out.sort_values("p").set_index("gene")
