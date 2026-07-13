import pandas as pd

from mirna_tcga.endpoints import (
    distant_metastasis,
    distant_metastasis_biotab,
    nodal_metastasis,
    nodal_metastasis_within_stage,
)

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


def test_distant_metastasis_biotab_unions_both_sources():
    # Biotab adds P6 (a follow-up distant met the cBioPortal frame missed) and
    # upgrades P4 (M0 in clin) to positive on a metachronous distant met.
    biotab = pd.Series({"P4": 1, "P6": 1, "P1": 0}, dtype="Int64")
    d = distant_metastasis_biotab(CLIN, biotab)
    assert d.loc["P1"] == 0            # M0 in both sources
    assert d.loc["P2"] == 1            # M1 from clin, absent in biotab
    assert d.loc["P4"] == 1            # M0 in clin but biotab distant met -> positive wins
    assert d.loc["P6"] == 1            # only biotab has it (clin dropped it)
    assert set(d.unique()) <= {0, 1}
    # Never fewer patients than the cBioPortal-only endpoint.
    assert set(distant_metastasis(CLIN).index) <= set(d.index)


def test_nodal_metastasis_within_stage_restricts_to_one_stage():
    clin = pd.DataFrame({
        "patientId": ["Q1", "Q2", "Q3", "Q4", "Q5"],
        "PATH_N_STAGE": ["N0", "N2", "N0", "N1", "N0"],
        "AJCC_PATHOLOGIC_TUMOR_STAGE":
            ["STAGE IIA", "STAGE IIB", "STAGE IB", "STAGE IIIA", "STAGE IIA"],
    })
    d = nodal_metastasis_within_stage(clin, "II")
    assert d.loc["Q1"] == 0            # stage II, N0
    assert d.loc["Q2"] == 1            # stage II, N2
    assert d.loc["Q5"] == 0            # stage II, N0
    assert "Q3" not in d.index         # stage I -> excluded
    assert "Q4" not in d.index         # stage III -> excluded
    assert set(d.unique()) <= {0, 1}


def test_nodal_metastasis_binarizes_and_drops_nx():
    n = nodal_metastasis(CLIN)
    assert n.loc["P1"] == 0
    assert n.loc["P2"] == 1            # N2
    assert n.loc["P4"] == 1            # N1
    assert n.loc["P5"] == 0            # N0
    assert "P3" not in n.index         # NX dropped
    assert "P6" not in n.index         # missing dropped
