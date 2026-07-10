"""Offline tests for the Xena loader using a fake HTTP session."""

import gzip

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
    def __init__(self, content: bytes, status_code: int = 200):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise AssertionError(f"HTTP {self.status_code}")


class FakeSession:
    """Serves ``content`` for the ``.gz`` URL; 404s everything else by default."""

    def __init__(self, content: bytes, gz_only: bool = True):
        self._content = content
        self.gz_only = gz_only
        self.last_url = None
        self.urls = []

    def get(self, url, timeout=None):
        self.urls.append(url)
        self.last_url = url
        if self.gz_only and not url.endswith(".gz"):
            return FakeResponse(b"", status_code=403)
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


def test_download_requests_gz_object():
    # Xena stores matrices gzipped at <dataset>.gz; the loader must request that.
    sess = FakeSession(EXPR_TSV.encode("utf-8"))
    XenaClient(host="tcga", session=sess).expression_matrix("TCGA.LUAD.sampleMap/miRNA_HiSeq_gene")
    assert sess.last_url == (
        "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/miRNA_HiSeq_gene.gz"
    )


def test_download_falls_back_to_plain_when_gz_absent():
    # A dataset stored uncompressed: .gz 404s, loader retries the bare key.
    class PlainOnly(FakeSession):
        def get(self, url, timeout=None):
            self.urls.append(url)
            self.last_url = url
            if url.endswith(".gz"):
                return FakeResponse(b"", status_code=404)
            return FakeResponse(self._content)

    sess = PlainOnly(EXPR_TSV.encode("utf-8"))
    mat = XenaClient(host="tcga", session=sess).expression_matrix("ds")
    assert mat.loc["hsa-miR-21", "S1"] == 5.0
    assert sess.urls[0].endswith(".gz") and sess.last_url.endswith("/ds")
