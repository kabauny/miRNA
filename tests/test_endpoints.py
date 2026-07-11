import pandas as pd

from mirna_tcga.endpoints import distant_metastasis, nodal_metastasis

CLIN = pd.DataFrame({
    "patientId": ["P1", "P2", "P3", "P4", "P5", "P6"],
    "PATH_M_STAGE": ["M0", "M1", "MX", "M0", "M1B", None],
    "AJCC_PATHOLOGIC_TUMOR_STAGE": ["STAGE IB", "STAGE IV", "STAGE IV", "STAGE IIA", "STAGE IV", "STAGE IA"],
    "PATH_N_STAGE": ["N0", "N2", "NX", "N1", "N0", None],
})


def test_distant_metastasis_pools_m1_and_stage_iv_and_drops_unknown():
    d = distant_metastasis(CLIN)
    assert d.loc["P1"] == 0            # M0
    assert d.loc["P2"] == 1            # M1
    assert d.loc["P3"] == 1            # MX but STAGE IV -> positive
    assert d.loc["P4"] == 0            # M0
    assert d.loc["P5"] == 1            # M1B + stage IV
    assert "P6" not in d.index         # no M call, not stage IV -> dropped
    assert set(d.unique()) <= {0, 1}


def test_nodal_metastasis_binarizes_and_drops_nx():
    n = nodal_metastasis(CLIN)
    assert n.loc["P1"] == 0
    assert n.loc["P2"] == 1            # N2
    assert n.loc["P4"] == 1            # N1
    assert n.loc["P5"] == 0            # N0
    assert "P3" not in n.index         # NX dropped
    assert "P6" not in n.index         # missing dropped
