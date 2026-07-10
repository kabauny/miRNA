import numpy as np
import pandas as pd

from mirna_tcga.screen import benjamini_hochberg, cox_score_screen


def test_benjamini_hochberg_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.03, 0.5, 0.9])
    q = benjamini_hochberg(p)
    assert np.all((q >= 0) & (q <= 1))
    assert np.all(np.diff(q) >= -1e-12)  # non-decreasing with p (already sorted)
    assert q[0] >= p[0]  # q >= p for the smallest


def _survival_data(seed=0, n=400, beta=0.8):
    rng = np.random.default_rng(seed)
    strata = rng.choice(["A", "B"], n)
    risky = rng.standard_normal(n)     # hazardous
    protective = rng.standard_normal(n)
    noise = rng.standard_normal(n)
    lp = beta * risky - beta * protective + (strata == "B") * 0.3
    t = rng.exponential(np.exp(-lp))
    c = rng.exponential(1.5, n)
    event = (t <= c).astype(int)
    time = np.minimum(t, c)
    X = pd.DataFrame({"RISKY": risky, "PROTECTIVE": protective, "NOISE": noise})
    return X, pd.Series(time), pd.Series(event), pd.Series(strata)


def test_score_screen_recovers_direction_and_signal():
    X, t, e, s = _survival_data()
    res = cox_score_screen(X, t, e, strata=s)
    # both true genes are far more significant than noise
    assert res.loc["RISKY", "p"] < 1e-3
    assert res.loc["PROTECTIVE", "p"] < 1e-3
    assert res.loc["NOISE", "p"] > 0.05
    # sign convention: higher expression -> higher hazard is z > 0
    assert res.loc["RISKY", "z"] > 0
    assert res.loc["PROTECTIVE", "z"] < 0


def test_score_screen_unstratified_runs_and_sorts_by_p():
    X, t, e, _ = _survival_data(seed=3)
    res = cox_score_screen(X, t, e)
    assert list(res.columns) == ["z", "p", "q"]
    assert res["p"].is_monotonic_increasing
    assert set(res.index) == {"RISKY", "PROTECTIVE", "NOISE"}
