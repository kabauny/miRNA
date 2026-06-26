import numpy as np
import pandas as pd

from mirna_tcga.classify import fit_pam_cv, signature


def _separable_data(n=60, p=8, seed=0):
    """Two classes separable on the first two features; rest are noise."""
    rng = np.random.default_rng(seed)
    half = n // 2
    X = rng.normal(0, 1, (n, p))
    y = np.array(["A"] * half + ["B"] * half)
    X[:half, 0] += 4.0   # informative
    X[:half, 1] -= 4.0   # informative
    cols = [f"g{i}" for i in range(p)]
    return pd.DataFrame(X, columns=cols), pd.Series(y)


def test_pam_learns_separable_classes():
    X, y = _separable_data()
    result = fit_pam_cv(X, y, cv_folds=5, random_state=0)
    assert result.cv_error < 0.1
    assert result.confusion.shape == (2, 2)
    assert result.labels == ["A", "B"]


def test_signature_ranks_informative_features_first():
    X, y = _separable_data()
    result = fit_pam_cv(X, y, cv_folds=5, random_state=0)
    sig = signature(result, list(X.columns))
    assert not sig.empty
    # the two injected features should rank at the top
    top2 = set(sig["gene"].head(2))
    assert {"g0", "g1"}.issubset(top2)


def test_cv_errors_indexed_by_threshold():
    X, y = _separable_data()
    result = fit_pam_cv(X, y, cv_folds=5, random_state=0)
    assert result.best_threshold in set(result.cv_errors.index)
    assert (result.cv_errors >= 0).all()
