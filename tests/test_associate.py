import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from mirna_tcga.associate import fisher_screen, ranksum_screen


def test_ranksum_matches_scipy_unstratified():
    rng = np.random.default_rng(0)
    n = 200
    y = rng.integers(0, 2, n)
    signal = rng.standard_normal(n) + 1.5 * y   # higher in positives
    noise = rng.standard_normal(n)
    X = pd.DataFrame({"SIGNAL": signal, "NOISE": noise})
    res = ranksum_screen(X, pd.Series(y))
    # scipy two-sided p for the signal gene
    u, p = mannwhitneyu(signal[y == 1], signal[y == 0], alternative="two-sided")
    assert abs(res.loc["SIGNAL", "p"] - p) < 1e-6
    assert res.loc["SIGNAL", "z"] > 0          # higher in positives
    assert res.loc["SIGNAL", "auc"] > 0.6
    assert res.loc["NOISE", "p"] > 0.05


def test_ranksum_stratified_controls_confounder():
    # feature differs by stratum, but NOT by outcome within stratum -> null
    rng = np.random.default_rng(1)
    strata = np.array(["A"] * 150 + ["B"] * 150)
    y = rng.integers(0, 2, 300)
    x = np.where(strata == "A", 0.0, 5.0) + rng.standard_normal(300)  # driven by stratum only
    res = ranksum_screen(pd.DataFrame({"F": x}), pd.Series(y), pd.Series(strata))
    assert res.loc["F", "p"] > 0.05  # stratification removes the confound


def test_fisher_screen_detects_and_ranks():
    # DEL co-occurs with the outcome; CLEAN does not
    y = np.array([1] * 20 + [0] * 20)
    dele = np.array([1] * 16 + [0] * 4 + [0] * 4 + [1] * 16)   # enriched in y=1... see below
    dele = np.array([1] * 15 + [0] * 5 + [0] * 17 + [1] * 3)   # 15/20 pos vs 3/20 neg
    clean = np.array(([0, 1] * 20))
    B = pd.DataFrame({"DEL": dele, "CLEAN": clean})
    res = fisher_screen(B, pd.Series(y))
    assert res.loc["DEL", "odds_ratio"] > 1
    assert res.loc["DEL", "p"] < 0.01
    assert res.loc["DEL", "p"] < res.loc["CLEAN", "p"]
    assert res.loc["DEL", "n_alt"] == 18
    assert "q" in res.columns
