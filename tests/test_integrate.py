import numpy as np
import pandas as pd
import pytest

from mirna_tcga.integrate import (
    attach_clinical,
    overlap_signatures,
    sample_to_patient,
    signature_score,
    stratify,
)


def _X():
    return pd.DataFrame(
        {"g1": [1.0, 2.0, 3.0, 4.0], "g2": [4.0, 3.0, 2.0, 1.0], "noise": [0.0, 0.0, 0.0, 0.0]},
        index=["TCGA-AA-0001-01", "TCGA-AA-0002-01", "TCGA-AA-0003-01", "TCGA-AA-0004-01"],
    )


def test_signature_score_unweighted_is_mean():
    s = signature_score(_X(), ["g1", "g2"])
    assert s.tolist() == [2.5, 2.5, 2.5, 2.5]


def test_signature_score_weighted_is_directional():
    w = pd.Series({"g1": 1.0, "g2": -1.0})  # reward g1, penalize g2
    s = signature_score(_X(), ["g1", "g2"], weights=w)
    assert s.is_monotonic_increasing  # g1 up + g2 down => increasing score


def test_signature_score_ignores_missing_and_raises_when_none():
    s = signature_score(_X(), ["g1", "absent"])
    assert s.tolist() == [1.0, 2.0, 3.0, 4.0]
    with pytest.raises(ValueError):
        signature_score(_X(), ["absent"])


def test_stratify_median_split():
    scores = pd.Series([1.0, 2.0, 3.0, 4.0])
    g = stratify(scores, method="median")
    assert g.tolist() == ["low", "low", "high", "high"]


def test_overlap_signatures_intersects_on_gene():
    a = pd.DataFrame({"gene": ["x", "y", "z"]})
    b = pd.DataFrame({"gene": ["y", "z", "w"]})
    c = pd.DataFrame({"gene": ["z", "y"]})
    assert overlap_signatures(a, b, c) == ["y", "z"]


def test_sample_to_patient_collapses_barcode():
    assert sample_to_patient(["TCGA-05-4244-01", "TCGA-05-4244-11"]) == [
        "TCGA-05-4244", "TCGA-05-4244",
    ]
    assert sample_to_patient(["weird_id"]) == ["weird_id"]


def test_attach_clinical_joins_on_patient():
    scores = pd.Series([0.5, 1.5], index=["TCGA-AA-0001-01", "TCGA-AA-0002-01"])
    clinical = pd.DataFrame(
        {"patientId": ["TCGA-AA-0001", "TCGA-AA-0099"], "OS_MONTHS": ["12", "30"]}
    )
    merged = attach_clinical(scores, clinical)
    assert len(merged) == 1                       # only patient 0001 matched
    assert merged.iloc[0]["OS_MONTHS"] == "12"
    assert np.isclose(merged.iloc[0]["score"], 0.5)
