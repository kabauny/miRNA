import numpy as np
import pandas as pd

from mirna_tcga.preprocess import (
    boxcox_normalize,
    drop_near_zero_variance,
    near_zero_variance,
    skewed_features,
    standardize,
)


def _frame():
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "good": rng.normal(10, 3, 100),
        "constant": np.ones(100),
        "near_zero": [1.0] * 98 + [2.0, 3.0],  # dominant value, few uniques
    })


def test_near_zero_variance_flags_constant_and_nzv():
    flagged = near_zero_variance(_frame())
    assert "constant" in flagged
    assert "near_zero" in flagged
    assert "good" not in flagged


def test_drop_near_zero_variance_keeps_informative():
    out = drop_near_zero_variance(_frame())
    assert list(out.columns) == ["good"]


def test_standardize_zero_mean_unit_sd():
    out = standardize(_frame()[["good"]])
    assert abs(out["good"].mean()) < 1e-9
    assert abs(out["good"].std(ddof=1) - 1.0) < 1e-9


def test_standardize_handles_constant_column():
    out = standardize(_frame()[["constant"]])
    assert (out["constant"] == 0).all()  # no div-by-zero blow-up


def test_boxcox_returns_finite_same_shape():
    X = _frame()[["good"]].abs()
    out = boxcox_normalize(X)
    assert out.shape == X.shape
    assert np.isfinite(out.to_numpy()).all()


def test_skewed_features_detects_heavy_tail():
    rng = np.random.default_rng(1)
    df = pd.DataFrame({
        "normal": rng.normal(0, 1, 2000),
        "heavy": rng.standard_t(2, 2000),  # heavy-tailed -> high kurtosis
    })
    flagged = skewed_features(df)
    assert "heavy" in flagged
