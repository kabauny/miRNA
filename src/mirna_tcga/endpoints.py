"""Derive binary clinical endpoints (metastasis) from a cBioPortal clinical table.

Each helper returns a per-patient ``pandas.Series`` of {0, 1} indexed by
``patientId``, with patients whose status is unknown/ambiguous dropped (so the
positive/negative counts are honest). These feed the association screens and the
predictive models the same way survival duration/event feed the Cox tests.
"""

from __future__ import annotations

import pandas as pd

_M_POS = {"M1", "M1A", "M1B"}
_N_POS = {"N1", "N2", "N3"}


def _series(clin: pd.DataFrame, values: pd.Series, name: str) -> pd.Series:
    out = pd.Series(values.to_numpy(), index=clin["patientId"], name=name)
    return out.dropna().astype(int)


def distant_metastasis(clin: pd.DataFrame) -> pd.Series:
    """1 = distant metastasis (pathologic M1 or clinical stage IV), 0 = M0.

    Pools ``PATH_M_STAGE`` in {M1, M1A, M1B} with ``AJCC_PATHOLOGIC_TUMOR_STAGE``
    == "STAGE IV". Patients coded MX (unassessed) with no stage-IV call and those
    with neither field are dropped -- only confident M0 vs M1 patients are kept.
    Note: distant metastasis is rare in resected TCGA lung cohorts, so this
    endpoint is under-powered for genome-wide discovery (report effect sizes).
    """
    m = clin.get("PATH_M_STAGE", pd.Series(index=clin.index, dtype=object)).astype(str)
    stage = clin.get("AJCC_PATHOLOGIC_TUMOR_STAGE", pd.Series(index=clin.index, dtype=object)).astype(str)
    pos = m.isin(_M_POS) | stage.eq("STAGE IV")
    neg = m.eq("M0")
    val = pd.Series(pd.NA, index=clin.index, dtype="Float64")
    val[neg] = 0
    val[pos] = 1  # positive overrides an M0 call (stage IV without M1 code)
    return _series(clin, val, "distant_met")


def nodal_metastasis(clin: pd.DataFrame) -> pd.Series:
    """1 = regional lymph-node metastasis (N1/N2/N3), 0 = N0.

    ``PATH_N_STAGE`` NX (unassessed) / missing are dropped. This is the
    best-powered metastasis endpoint in TCGA lung.
    """
    n = clin.get("PATH_N_STAGE", pd.Series(index=clin.index, dtype=object)).astype(str)
    val = pd.Series(pd.NA, index=clin.index, dtype="Float64")
    val[n.eq("N0")] = 0
    val[n.isin(_N_POS)] = 1
    return _series(clin, val, "nodal_met")
