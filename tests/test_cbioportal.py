"""Offline tests for the cBioPortal client using a fake HTTP session."""

from mirna_tcga.cbioportal import CBioPortalClient


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        pass

    def json(self):
        return self._payload


class FakeSession:
    """Records calls and replays queued responses."""

    def __init__(self):
        self.headers = {}
        self.get_queue = []
        self.post_queue = []
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append(("GET", url, params))
        return FakeResponse(self.get_queue.pop(0))

    def post(self, url, json=None, params=None, timeout=None):
        self.calls.append(("POST", url, json))
        return FakeResponse(self.post_queue.pop(0))


def _client(session):
    return CBioPortalClient(base_url="https://example/api", page_size=2, session=session)


def test_get_paginates_until_short_page():
    sess = FakeSession()
    # page 0 full (==page_size), page 1 short -> stop
    sess.get_queue = [[{"a": 1}, {"a": 2}], [{"a": 3}]]
    client = _client(sess)
    out = client._get("studies")
    assert [r["a"] for r in out] == [1, 2, 3]
    # two GETs issued with incrementing pageNumber
    assert sess.calls[0][2]["pageNumber"] == 0
    assert sess.calls[1][2]["pageNumber"] == 1


def test_entrez_ids_from_symbols():
    sess = FakeSession()
    sess.post_queue = [[{"hugoGeneSymbol": "TP53", "entrezGeneId": 7157}]]
    client = _client(sess)
    assert client.entrez_ids(["TP53"]) == [7157]


def test_expression_matrix_pivots_long_to_wide():
    sess = FakeSession()
    # 1) genes/fetch  2) molecular-data/fetch
    sess.post_queue = [
        [{"hugoGeneSymbol": "TP53", "entrezGeneId": 7157},
         {"hugoGeneSymbol": "EGFR", "entrezGeneId": 1956}],
        [
            {"hugoGeneSymbol": "TP53", "sampleId": "S1", "value": 1.0},
            {"hugoGeneSymbol": "TP53", "sampleId": "S2", "value": 2.0},
            {"hugoGeneSymbol": "EGFR", "sampleId": "S1", "value": 3.0},
            {"hugoGeneSymbol": "EGFR", "sampleId": "S2", "value": 4.0},
        ],
    ]
    client = _client(sess)
    mat = client.expression_matrix("prof", ["TP53", "EGFR"], "list_all")
    assert mat.shape == (2, 2)
    assert mat.loc["TP53", "S2"] == 2.0
    assert mat.loc["EGFR", "S1"] == 3.0


def test_expression_matrix_maps_entrez_when_symbol_absent():
    # The real SUMMARY projection returns entrezGeneId but NOT hugoGeneSymbol;
    # the matrix must still be keyed by HUGO symbol via the gene lookup.
    sess = FakeSession()
    sess.post_queue = [
        [{"hugoGeneSymbol": "TP53", "entrezGeneId": 7157},
         {"hugoGeneSymbol": "EGFR", "entrezGeneId": 1956}],
        [
            {"entrezGeneId": 7157, "sampleId": "S1", "value": 1.0},
            {"entrezGeneId": 1956, "sampleId": "S1", "value": 3.0},
            {"entrezGeneId": 1956, "sampleId": "S2", "value": 4.0},
        ],
    ]
    client = _client(sess)
    mat = client.expression_matrix("prof", ["TP53", "EGFR"], "list_all")
    assert mat.loc["TP53", "S1"] == 1.0
    assert mat.loc["EGFR", "S2"] == 4.0


def test_mutated_samples_excludes_silent():
    sess = FakeSession()
    sess.post_queue = [
        [{"hugoGeneSymbol": "TP53", "entrezGeneId": 7157}],  # genes/fetch
        [
            {"sampleId": "S1", "mutationType": "Missense_Mutation"},
            {"sampleId": "S2", "mutationType": "Silent"},
            {"sampleId": "S3", "mutationType": "Nonsense_Mutation"},
        ],
    ]
    client = _client(sess)
    muts = client.mutated_samples("mutprof", "TP53", "list_all")
    assert muts == {"S1", "S3"}


def test_clinical_data_pivots_to_wide():
    sess = FakeSession()
    # page_size is 2, so a 3-record result spills onto a second (empty) page
    sess.get_queue = [
        [
            {"patientId": "P1", "clinicalAttributeId": "OS_MONTHS", "value": "12"},
            {"patientId": "P1", "clinicalAttributeId": "OS_STATUS", "value": "1:DECEASED"},
        ],
        [
            {"patientId": "P2", "clinicalAttributeId": "OS_MONTHS", "value": "30"},
        ],
    ]
    client = _client(sess)
    df = client.get_clinical_data("study", patient_level=True)
    row = df.set_index("patientId").loc["P1"]
    assert row["OS_MONTHS"] == "12"
    assert row["OS_STATUS"] == "1:DECEASED"
