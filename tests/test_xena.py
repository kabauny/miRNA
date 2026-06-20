"""Offline tests for the Xena loader using a fake HTTP session."""

import gzip

import pandas as pd

from mirna_tcga.xena import HOSTS, XenaClient

EXPR_TSV = (
    "sample\tS1\tS2\tS3\n"
    "hsa-miR-21\t5.0\t6.0\t7.0\n"
    "hsa-miR-205\t1.0\t0.5\t0.2\n"
)

PHENO_TSV = (
    "sample\tsample_type\tstage\n"
    "S1\tTumor\tI\n"
    "S2\tTumor\tII\n"
)


class FakeResponse:
    def __init__(self, content: bytes):
        self.content = content

    def raise_for_status(self):
        pass


class FakeSession:
    def __init__(self, content: bytes):
        self._content = content
        self.last_url = None

    def get(self, url, timeout=None):
        self.last_url = url
        return FakeResponse(self._content)


def test_host_alias_resolves():
    c = XenaClient(host="tcga")
    assert c.host == HOSTS["tcga"]
    assert XenaClient(host="custom.host.net").host == "custom.host.net"


def test_expression_matrix_plain_tsv():
    client = XenaClient(session=FakeSession(EXPR_TSV.encode("utf-8")))
    mat = client.expression_matrix("TCGA.LUAD.sampleMap/miRNA_HiSeq_gene")
    assert mat.shape == (2, 3)              # 2 miRNAs x 3 samples
    assert mat.loc["hsa-miR-21", "S2"] == 6.0
    assert list(mat.columns) == ["S1", "S2", "S3"]


def test_expression_matrix_handles_gzip():
    gz = gzip.compress(EXPR_TSV.encode("utf-8"))
    client = XenaClient(session=FakeSession(gz))
    mat = client.expression_matrix("ds")
    assert mat.loc["hsa-miR-205", "S3"] == 0.2


def test_expression_matrix_gene_subset_ignores_missing():
    client = XenaClient(session=FakeSession(EXPR_TSV.encode("utf-8")))
    mat = client.expression_matrix("ds", genes=["hsa-miR-21", "not-present"])
    assert list(mat.index) == ["hsa-miR-21"]


def test_phenotype_orientation_samples_as_rows():
    client = XenaClient(session=FakeSession(PHENO_TSV.encode("utf-8")))
    pheno = client.phenotype("ds")
    assert pheno.loc["S1", "sample_type"] == "Tumor"
    assert pheno.index.name == "sample"


def test_download_url_built_from_host():
    sess = FakeSession(EXPR_TSV.encode("utf-8"))
    XenaClient(host="tcga", session=sess).expression_matrix("TCGA.LUAD.sampleMap/miRNA_HiSeq_gene")
    assert sess.last_url == "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/miRNA_HiSeq_gene"
