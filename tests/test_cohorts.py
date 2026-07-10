import pandas as pd

from mirna_tcga.cohorts import cohort_study_keys, combined_clinical, nsclc_clinical


class FakeConfig:
    studies = {"luad": "luad_study", "lusc": "lusc_study"}
    cohorts = {"nsclc": ["luad", "lusc"]}


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
