import numpy as np
import pandas as pd

from mirna_tcga.layers import deletion_matrix, mutation_matrix, stream_expression


class FakeCfg:
    studies = {"luad": "luad", "lusc": "lusc"}

    def mrna_profile(self, k):
        return f"{k}_mrna"

    def all_samples_list(self, k):
        return f"{k}_all"

    def cna_profile(self, k):
        return f"{k}_gistic"

    def cna_samples_list(self, k):
        return f"{k}_cna"

    def sequenced_samples_list(self, k):
        return f"{k}_seq"

    def mutation_profile(self, k):
        return f"{k}_mut"


ID2SYM = {1: "AAA", 2: "BBB"}


class FakeClient:
    def __init__(self):
        self.samples = {"luad_cna": ["TCGA-1-01", "TCGA-2-01"], "luad_seq": ["TCGA-1-01", "TCGA-2-01"],
                        "lusc_cna": ["TCGA-3-01"], "lusc_seq": ["TCGA-3-01"]}

    def sample_list_ids(self, sid):
        return self.samples[sid]

    def discrete_cna_events(self, profile, sample_list, event_type="HOMDEL"):
        if profile == "luad_gistic":
            # sample 1 has AAA deleted; sample 2 nothing
            return pd.DataFrame({"sampleId": ["TCGA-1-01"], "entrezGeneId": [1], "alteration": [-2]})
        return pd.DataFrame()  # lusc: no deletions

    def mutation_events(self, profile, entrez, sample_list):
        if profile == "luad_mut":
            return pd.DataFrame({
                "sampleId": ["TCGA-1-01", "TCGA-2-01"],
                "entrezGeneId": [2, 2],
                "mutationType": ["Missense_Mutation", "Silent"],  # silent dropped
            })
        return pd.DataFrame()

    def fetch_molecular_data(self, profile, entrez, sample_list_id):
        s = "S1" if "luad" in profile else "S2"
        return pd.DataFrame({"entrezGeneId": [1, 2], "sampleId": [s, s], "value": [3.0, 7.0]})


def test_deletion_matrix_flags_and_tags_subtype():
    B, sub = deletion_matrix(FakeClient(), FakeCfg(), ID2SYM, ["luad", "lusc"])
    assert B.loc["TCGA-1-01", "AAA"] == 1
    assert B.loc["TCGA-2-01", "AAA"] == 0     # profiled but not deleted
    assert "TCGA-3-01" in B.index and B.loc["TCGA-3-01"].sum() == 0
    assert sub["TCGA-1-01"] == "LUAD" and sub["TCGA-3-01"] == "LUSC"


def test_mutation_matrix_excludes_silent():
    B, sub = mutation_matrix(FakeClient(), FakeCfg(), ID2SYM, ["luad", "lusc"])
    assert B.loc["TCGA-1-01", "BBB"] == 1     # missense kept
    assert B.loc["TCGA-2-01", "BBB"] == 0     # silent dropped -> not flagged


def test_stream_expression_log2_and_concatenates_studies():
    expr = stream_expression(FakeClient(), FakeCfg(), [1, 2], ID2SYM, ["luad", "lusc"])
    assert expr.shape == (2, 2)               # 2 genes x 2 study-samples
    assert np.isclose(expr.loc["BBB", "S1"], np.log2(8.0))  # log2(7+1)
