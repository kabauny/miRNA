"""Loader for TCGA **BCR Biotab** clinical supplements (local GDC download).

cBioPortal's clinical table only carries *pathologic* M stage, so "metastatic"
is defined as pathologic M1 / stage IV -- which in resected TCGA lung cohorts is
only ~33 patients and badly under-powers the metastasis screens. The BCR Biotab
clinical supplements (``TCGAbiolinks::GDCdownload(data.type = "Clinical
Supplement", data.format = "BCR Biotab")``) carry the fuller annotation:

* patient file  ``nationwidechildrens.org_clinical_patient_<study>.txt``
    - ``ajcc_metastasis_pathologic_pm``  (pathologic M: M0 / M1 / M1a / M1b / MX)
    - ``ajcc_metastasis_clinical_cm``    (clinical M)
    - ``metastatic_site_at_diagnosis``   (distant site => M1 at diagnosis)
* follow-up file ``nationwidechildrens.org_clinical_follow_up_v*_<study>.txt``
    - ``new_tumor_event_type`` == "Distant Metastasis" (possibly compound, e.g.
      "Locoregional Recurrence|Distant Metastasis") -- metachronous distant met.

Unioning these lifts the confirmed-metastatic lung set from ~33 to ~120.

Format notes
------------
Biotab TSVs carry **three header rows**: row 1 = machine column names, row 2 = a
human-readable alias, row 3 = ``CDE_ID:...``. Only row 1 is a real header; rows 2
and 3 are skipped. **Column order is NOT stable across studies** (LUSC's clinical
M sits four columns over from LUAD's), so every field is looked up by name, never
by position. Missing values use the sentinels ``[Not Available]``,
``[Not Applicable]``, ``[Unknown]``, ``[Not Evaluated]``.

The download lives outside the repo (it is large and not committed); point
``root`` at the ``GDCdata`` directory of the local TCGA folder.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Missing-value sentinels used throughout BCR Biotab files.
_NA = {
    "",
    "[not available]",
    "[not applicable]",
    "[unknown]",
    "[not evaluated]",
    "[discrepancy]",
    "[pending]",
    "nan",
}

# Pathologic / clinical M codes that denote distant metastasis.
_M_POS = {"m1", "m1a", "m1b", "cm1", "cm1a", "cm1b", "pm1", "pm1a", "pm1b"}
_M_POS_UPPER = {m.upper() for m in _M_POS}


def _clean(s: pd.Series) -> pd.Series:
    """Lower-case, strip, and map Biotab NA sentinels to ``pd.NA``."""
    out = s.astype(str).str.strip()
    return out.mask(out.str.lower().isin(_NA))


def _read_biotab(path: Path) -> pd.DataFrame:
    """Read one Biotab TSV, dropping the two sub-header rows (alias + CDE_ID).

    ``read_csv`` consumes line 1 as the header, so the alias and ``CDE_ID`` rows
    land at positions 0 and 1 -- ``iloc[2:]`` drops them and is safe on files with
    fewer than two data rows.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return df.iloc[2:].reset_index(drop=True)


def _find(root: Path, study: str, kind: str) -> list[Path]:
    """Glob the Biotab file(s) for a study under the GDC layout.

    ``kind`` is ``"patient"`` or ``"follow_up"``. UUID sub-directories vary, and a
    study may have several follow-up versions (v1.0, v4.0, ``_nte_`` ...), so this
    returns every match.
    """
    pat = (
        f"TCGA-{study.upper()}/**/Clinical_Supplement/**/"
        f"nationwidechildrens.org_clinical_{kind}_*_{study.lower()}.txt"
        if kind == "follow_up"
        else
        f"TCGA-{study.upper()}/**/Clinical_Supplement/**/"
        f"nationwidechildrens.org_clinical_{kind}_{study.lower()}.txt"
    )
    return sorted(root.glob(pat))


def _patient_positive(root: Path, study: str) -> tuple[set[str], set[str]]:
    """Barcodes with distant-met evidence and confirmed-M0 barcodes (patient file)."""
    pos: set[str] = set()
    m0: set[str] = set()
    for path in _find(root, study, "patient"):
        df = _read_biotab(path)
        bc = _clean(df["bcr_patient_barcode"])
        pm = _clean(df.get("ajcc_metastasis_pathologic_pm", pd.Series(index=df.index, dtype=object)))
        cm = _clean(df.get("ajcc_metastasis_clinical_cm", pd.Series(index=df.index, dtype=object)))
        site = _clean(df.get("metastatic_site_at_diagnosis", pd.Series(index=df.index, dtype=object)))
        pm_l, cm_l = pm.str.lower(), cm.str.lower()
        is_pos = pm_l.isin(_M_POS) | cm_l.isin(_M_POS) | site.notna()
        is_m0 = pm_l.eq("m0") | cm_l.eq("m0")
        pos |= set(bc[is_pos].dropna())
        m0 |= set(bc[is_m0].dropna())
    return pos, m0


def _followup_positive(root: Path, study: str) -> set[str]:
    """Barcodes whose follow-up records a 'Distant Metastasis' new-tumor event."""
    pos: set[str] = set()
    for path in _find(root, study, "follow_up"):
        df = _read_biotab(path)
        if "new_tumor_event_type" not in df.columns:
            continue
        bc = _clean(df["bcr_patient_barcode"])
        nte = _clean(df["new_tumor_event_type"]).str.lower()
        # Substring catches compound codes ("Locoregional Recurrence|Distant Metastasis").
        pos |= set(bc[nte.str.contains("distant metastasis", na=False)].dropna())
    return pos


def distant_metastasis_label(root: str | Path, study: str) -> pd.Series:
    """Per-patient distant-metastasis label for one study from local Biotab files.

    Returns a ``pandas.Series`` of {0, 1} indexed by ``bcr_patient_barcode``
    (e.g. ``TCGA-05-4384`` -- the same key as cBioPortal ``patientId``). A patient
    is:

    * ``1`` if the patient file codes pathologic/clinical **M1** (any subclass),
      records a metastatic site at diagnosis, or any follow-up records a
      **Distant Metastasis** new-tumor event;
    * ``0`` if coded **M0** with no positive evidence;
    * **dropped** otherwise (MX / missing with no distant-met signal), so the
      0/1 counts stay honest.
    """
    root = Path(root)
    pos_pat, m0 = _patient_positive(root, study)
    pos = pos_pat | _followup_positive(root, study)
    neg = m0 - pos  # a later distant-met event overrides an initial M0 call
    label = {bc: 1 for bc in pos}
    label.update({bc: 0 for bc in neg})
    return pd.Series(label, name="distant_met_biotab", dtype=int).sort_index()


def distant_metastasis_labels(root: str | Path, studies: list[str]) -> pd.Series:
    """Concatenate :func:`distant_metastasis_label` across studies (e.g. NSCLC)."""
    parts = [distant_metastasis_label(root, s) for s in studies]
    return pd.concat(parts).sort_index()


def _num(s: pd.Series) -> pd.Series:
    """Parse a Biotab day-count column to numeric (sentinels -> NaN)."""
    return pd.to_numeric(_clean(s), errors="coerce")


def _simplify_stage(v: object) -> str | None:
    """'Stage IIIA' -> 'III'; unknown/sentinel -> None."""
    s = str(v).upper().replace("STAGE ", "").strip()
    for r in ("IV", "III", "II", "I"):
        if s.startswith(r):
            return r
    return None


def patient_followup(root: str | Path, studies: list[str]) -> pd.DataFrame:
    """Per-patient staging + longitudinal follow-up from the BCR Biotab files.

    Aggregates the patient file and every follow-up file (a patient has one
    patient row and zero or more follow-up rows) into one record per
    ``bcr_patient_barcode``, indexed by barcode, with columns:

    * ``stage`` -- pathologic stage collapsed to I / II / III / IV;
    * ``pN`` / ``pM`` -- pathologic nodal / metastasis codes at diagnosis;
    * ``followup_days`` -- longest observed follow-up (max of every
      ``last_contact_days_to`` / ``death_days_to`` across patient + follow-up rows);
    * ``distant_day`` / ``loco_day`` -- earliest days-to-new-tumor-event for a
      **Distant Metastasis** / **Locoregional Recurrence** event (``inf`` if none).
    """
    root = Path(root)
    rec: dict[str, dict] = {}

    def _slot(bc: str) -> dict:
        return rec.setdefault(bc, {"stage": None, "pN": None, "pM": None,
                                   "followup_days": float("-inf"),
                                   "distant_day": float("inf"), "loco_day": float("inf")})

    for study in studies:
        for path in _find(root, study, "patient"):
            df = _read_biotab(path)
            bc = _clean(df["bcr_patient_barcode"])
            stage = df.get("ajcc_pathologic_tumor_stage", pd.Series(index=df.index)).map(_simplify_stage)
            pN = _clean(df.get("ajcc_nodes_pathologic_pn", pd.Series(index=df.index))).str.upper().str[:2]
            pM = _clean(df.get("ajcc_metastasis_pathologic_pm", pd.Series(index=df.index))).str.upper()
            fu = pd.concat([_num(df.get("last_contact_days_to", pd.Series(index=df.index))),
                            _num(df.get("death_days_to", pd.Series(index=df.index)))], axis=1).max(axis=1)
            for i, b in bc.items():
                if pd.isna(b):
                    continue
                r = _slot(b)
                r["stage"], r["pN"], r["pM"] = stage[i], pN[i], pM[i]
                if pd.notna(fu[i]):
                    r["followup_days"] = max(r["followup_days"], float(fu[i]))

        for path in _find(root, study, "follow_up"):
            df = _read_biotab(path)
            bc = _clean(df["bcr_patient_barcode"])
            nte = _clean(df.get("new_tumor_event_type", pd.Series(index=df.index))).str.lower()
            day = _num(df.get("new_tumor_event_dx_days_to", pd.Series(index=df.index)))
            fu = pd.concat([_num(df.get("last_contact_days_to", pd.Series(index=df.index))),
                            _num(df.get("death_days_to", pd.Series(index=df.index)))], axis=1).max(axis=1)
            for i, b in bc.items():
                if pd.isna(b):
                    continue
                r = _slot(b)
                if pd.notna(fu[i]):
                    r["followup_days"] = max(r["followup_days"], float(fu[i]))
                t = nte[i] if isinstance(nte[i], str) else ""
                d = float(day[i]) if pd.notna(day[i]) else float("inf")
                if "distant metastasis" in t:
                    r["distant_day"] = min(r["distant_day"], d)
                if "locoregional" in t:
                    r["loco_day"] = min(r["loco_day"], d)

    out = pd.DataFrame.from_dict(rec, orient="index")
    out["followup_days"] = out["followup_days"].replace(float("-inf"), pd.NA)
    return out.sort_index()


def true_stage_i_vs_distant(
    root: str | Path, studies: list[str], window_days: int = 730
) -> pd.Series:
    """Distant-metastatic cases vs a **"true stage I"** indolent control group.

    A cleaner reference than "all M0": the standard endpoint contrasts distant
    metastasis against every non-metastatic patient, but that reference spans
    stage I-III and includes node-positive, locally-advanced disease. Here the
    reference is restricted to patients who were **stage I at diagnosis** (hence
    N0 / M0) and stayed **free of distant *and* locoregional recurrence for at
    least ``window_days``** (default 2 years), with follow-up long enough to
    confirm it.

    Returns a ``pandas.Series`` of {0, 1} indexed by ``bcr_patient_barcode``:

    * ``1`` -- distant metastasis ever (pathologic M1 at diagnosis or a Distant
      Metastasis new-tumor event);
    * ``0`` -- "true stage I": stage I + N0 + M0, ``followup_days >= window_days``,
      no distant metastasis ever, and no locoregional recurrence within the
      window;
    * **dropped** -- everyone else (M0 but not confirmed indolent, stage II-III
      non-metastatic, or too little follow-up to confirm the 2-year window).
    """
    fu = patient_followup(root, studies)
    days = pd.to_numeric(fu["followup_days"], errors="coerce")
    pM = fu["pM"].fillna("")
    distant_ever = (fu["distant_day"] < float("inf")) | pM.isin(_M_POS_UPPER)

    control = (
        fu["stage"].eq("I")
        & fu["pN"].eq("N0")
        & (pM.isin({"M0", "MX"}) | pM.eq(""))
        & (days >= window_days)
        & ~distant_ever
        & ~(fu["loco_day"] <= window_days)
    )
    label = pd.Series(pd.NA, index=fu.index, dtype="Int64")
    label[control] = 0
    label[distant_ever] = 1          # a case is never also a control
    return label.dropna().astype(int).sort_index()
