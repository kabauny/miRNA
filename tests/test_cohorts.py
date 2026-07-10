import pandas as pd

from mirna_tcga.cohorts import (
    cohort_study_keys,
    combined_clinical,
    combined_expression,
    nsclc_clinical,
    nsclc_expression,
)


class FakeConfig:
    studies = {"luad": "luad_study", "lusc": "lusc_study"}
    cohorts = {"nsclc": ["luad", "lusc"]}

    def mrna_profile(self, key):
        return self.studies[key] + "_mrna"

    def all_samples_list(self, key):
        return self.studies[key] + "_all"


class FakeClient:
    """Returns a canned clinical table per study id."""

    def __init__(self, tables):
        self.tables = tables
        self.calls = []

    def get_clinical_data(self, study_id, patient_level=True):
        self.calls.append((study_id, patient_level))
        return self.tables.get(study_id, pd.DataFrame())


def test_cohort_study_keys_default_nsclc():
    assert cohort_study_keys(FakeConfig(), "nsclc") == ["luad", "lusc"]


def test_combined_clinical_tags_cohort_and_concatenates():
    client = FakeClient({
        "luad_study": pd.DataFrame({"patientId": ["P1"], "OS_MONTHS": ["10"]}),
        "lusc_study": pd.DataFrame({"patientId": ["P2"], "OS_MONTHS": ["20"]}),
    })
    df = nsclc_clinical(client, FakeConfig())
    assert len(df) == 2
    assert set(df["cohort"]) == {"LUAD", "LUSC"}
    assert df.loc[df["patientId"] == "P1", "cohort"].item() == "LUAD"


def test_combined_clinical_skips_empty_studies():
    client = FakeClient({
        "luad_study": pd.DataFrame({"patientId": ["P1"]}),
        "lusc_study": pd.DataFrame(),  # no data for this study
    })
    df = combined_clinical(client, FakeConfig(), ["luad", "lusc"])
    assert list(df["cohort"]) == ["LUAD"]


def test_combined_clinical_all_empty_returns_empty():
    client = FakeClient({})
    df = combined_clinical(client, FakeConfig(), ["luad", "lusc"])
    assert df.empty


class FakeExprClient:
    """Returns a canned genes x samples matrix per mRNA profile id."""

    def __init__(self, mats):
        self.mats = mats
        self.calls = []

    def expression_matrix(self, molecular_profile_id, hugo_symbols, sample_list_id):
        self.calls.append((molecular_profile_id, tuple(hugo_symbols), sample_list_id))
        return self.mats.get(molecular_profile_id, pd.DataFrame())


def _mat(samples, value):
    return pd.DataFrame({s: [value] for s in samples}, index=["EGFR"])


def test_combined_expression_concatenates_samples_across_studies():
    client = FakeExprClient({
        "luad_study_mrna": _mat(["TCGA-05-0001-01", "TCGA-05-0002-01"], 3.0),
        "lusc_study_mrna": _mat(["TCGA-18-0003-01"], 7.0),
    })
    mat = nsclc_expression(client, FakeConfig(), ["EGFR"])
    assert list(mat.index) == ["EGFR"]
    assert mat.shape == (1, 3)
    assert mat.loc["EGFR", "TCGA-18-0003-01"] == 7.0
    # correct profile + sample-list ids were requested per study
    assert ("luad_study_mrna", ("EGFR",), "luad_study_all") in client.calls


def test_combined_expression_skips_empty_and_all_empty_returns_empty():
    client = FakeExprClient({"luad_study_mrna": _mat(["TCGA-05-0001-01"], 1.0)})
    mat = combined_expression(client, FakeConfig(), ["EGFR"], ["luad", "lusc"])
    assert mat.shape == (1, 1)
    empty = combined_expression(FakeExprClient({}), FakeConfig(), ["EGFR"], ["luad", "lusc"])
    assert empty.empty
