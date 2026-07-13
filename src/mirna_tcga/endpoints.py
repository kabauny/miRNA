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


def distant_metastasis_biotab(
    clin: pd.DataFrame, biotab_label: pd.Series
) -> pd.Series:
    """Distant metastasis, enriched with BCR Biotab annotation.

    Unions the cBioPortal-derived :func:`distant_metastasis` call (pathologic M1 /
    stage IV) with a ``biotab_label`` from
    :func:`mirna_tcga.biotab.distant_metastasis_labels` (clinical M1 + metachronous
    'Distant Metastasis' new-tumor events). The Biotab label is keyed by
    ``bcr_patient_barcode``, which equals cBioPortal ``patientId``.

    Resolution rule per patient (positive dominates, then negative, then whichever
    source has a call):

    * ``1`` if **either** source calls distant metastasis;
    * ``0`` if **either** source confidently calls M0 and neither calls positive;
    * dropped if neither source has a confident call.

    This lifts the confirmed-metastatic lung set from ~33 (pathologic-M1 only) to
    ~120, giving the metastasis screens real power. ``clin`` alone (no Biotab)
    still works via :func:`distant_metastasis`; this is the enriched variant.
    """
    base = distant_metastasis(clin)  # indexed by patientId, {0,1}
    bt = biotab_label.astype("Int64")
    idx = base.index.union(bt.index)
    b = base.reindex(idx).astype("Int64")
    t = bt.reindex(idx)
    val = pd.Series(pd.NA, index=idx, dtype="Int64")
    val[(b == 0) | (t == 0)] = 0
    val[(b == 1) | (t == 1)] = 1  # positive from either source wins
    out = val.dropna().astype(int)
    out.name = "distant_met"
    return out


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


def _simplify_stage(v: object) -> str | None:
    s = str(v).upper().replace("STAGE ", "").strip()
    for r in ("IV", "III", "II", "I"):
        if s.startswith(r):
            return r
    return None


def nodal_metastasis_within_stage(clin: pd.DataFrame, stage_code: str = "II") -> pd.Series:
    """1 = N+, 0 = N0, restricted to a single pathologic stage (e.g. "II").

    Holds tumour stage constant so a nodal contrast is not confounded by stage:
    across the whole cohort N0 tumours are ~78% stage I while N+ tumours are stage
    II-III, so plain ``nodal_metastasis`` partly measures stage/progression. Within
    one stage it measures nodal status per se. Stage II is the balanced choice in
    TCGA NSCLC (stage I is nearly all N0, stage III nearly all N+).
    """
    n = clin.get("PATH_N_STAGE", pd.Series(index=clin.index, dtype=object)).astype(str).str.upper().str[:2]
    stage = clin.get("AJCC_PATHOLOGIC_TUMOR_STAGE", pd.Series(index=clin.index, dtype=object)).map(_simplify_stage)
    in_stage = stage.eq(stage_code)
    val = pd.Series(pd.NA, index=clin.index, dtype="Float64")
    val[in_stage & n.eq("N0")] = 0
    val[in_stage & n.isin(_N_POS)] = 1
    return _series(clin, val, "nodal_met_stage")
