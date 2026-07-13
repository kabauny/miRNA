"""Offline tests for the BCR Biotab clinical-supplement loader.

Builds a tiny synthetic GDC directory tree (patient + follow-up TSVs, including
the three-row Biotab header) so the loader is exercised with no download.
"""

from mirna_tcga.biotab import (
    distant_metastasis_label,
    distant_metastasis_labels,
    true_stage_i_vs_distant,
)

# Row 1 = real header; rows 2-3 = the alias + CDE_ID sub-headers the loader skips.
_ALIAS = "alias"
_CDE = "CDE_ID:0"


def _write_biotab(path, columns, rows):
    lines = ["\t".join(columns), "\t".join([_ALIAS] * len(columns)), "\t".join([_CDE] * len(columns))]
    lines += ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _make_study(root, study, patient_rows, followup_rows):
    base = root / f"TCGA-{study.upper()}" / "harmonized" / "Clinical" / "Clinical_Supplement"
    pdir = base / "uuid-patient"
    fdir = base / "uuid-followup"
    pdir.mkdir(parents=True)
    fdir.mkdir(parents=True)
    _write_biotab(
        pdir / f"nationwidechildrens.org_clinical_patient_{study}.txt",
        ["bcr_patient_barcode", "ajcc_metastasis_pathologic_pm",
         "ajcc_metastasis_clinical_cm", "metastatic_site_at_diagnosis"],
        patient_rows,
    )
    _write_biotab(
        fdir / f"nationwidechildrens.org_clinical_follow_up_v1.0_{study}.txt",
        ["bcr_patient_barcode", "new_tumor_event_type"],
        followup_rows,
    )


def test_label_unions_pathologic_clinical_site_and_followup(tmp_path):
    _make_study(
        tmp_path, "luad",
        patient_rows=[
            ["TCGA-AA-0001", "M0", "[Not Applicable]", "[Not Available]"],   # M0
            ["TCGA-AA-0002", "M1", "[Not Applicable]", "[Not Available]"],   # pathologic M1
            ["TCGA-AA-0003", "M1b", "[Not Applicable]", "[Not Available]"],  # M1 subclass
            ["TCGA-AA-0004", "MX", "[Not Applicable]", "Bone"],              # site => positive
            ["TCGA-AA-0005", "MX", "[Not Applicable]", "[Not Available]"],   # MX, no evidence -> dropped
            ["TCGA-AA-0006", "M0", "[Not Applicable]", "[Not Available]"],   # M0 now, met on follow-up
        ],
        followup_rows=[
            ["TCGA-AA-0006", "Locoregional Recurrence|Distant Metastasis"],  # compound code
            ["TCGA-AA-0001", "[Not Available]"],
        ],
    )
    lab = distant_metastasis_label(tmp_path, "luad")
    assert lab.loc["TCGA-AA-0001"] == 0
    assert lab.loc["TCGA-AA-0002"] == 1
    assert lab.loc["TCGA-AA-0003"] == 1
    assert lab.loc["TCGA-AA-0004"] == 1          # metastatic site at diagnosis
    assert "TCGA-AA-0005" not in lab.index       # MX with no evidence dropped
    assert lab.loc["TCGA-AA-0006"] == 1          # follow-up distant met overrides M0
    assert set(lab.unique()) <= {0, 1}


def _make_study_full(root, study, patient_rows, followup_rows):
    """Patient + follow-up files carrying the staging / timing columns."""
    base = root / f"TCGA-{study.upper()}" / "harmonized" / "Clinical" / "Clinical_Supplement"
    pdir, fdir = base / "uuid-patient", base / "uuid-followup"
    pdir.mkdir(parents=True)
    fdir.mkdir(parents=True)
    _write_biotab(
        pdir / f"nationwidechildrens.org_clinical_patient_{study}.txt",
        ["bcr_patient_barcode", "ajcc_pathologic_tumor_stage", "ajcc_nodes_pathologic_pn",
         "ajcc_metastasis_pathologic_pm", "last_contact_days_to", "death_days_to"],
        patient_rows,
    )
    _write_biotab(
        fdir / f"nationwidechildrens.org_clinical_follow_up_v1.0_{study}.txt",
        ["bcr_patient_barcode", "new_tumor_event_type", "new_tumor_event_dx_days_to",
         "last_contact_days_to", "death_days_to"],
        followup_rows,
    )


def test_true_stage_i_vs_distant(tmp_path):
    _make_study_full(
        tmp_path, "luad",
        patient_rows=[
            # true stage I: stage I, N0, M0, followed >2yr, no events -> control (0)
            ["TCGA-AA-0001", "Stage IA", "N0", "M0", "1500", "[Not Available]"],
            # stage IV / M1 -> case (1)
            ["TCGA-AA-0002", "Stage IV", "N2", "M1", "400", "[Not Available]"],
            # stage I but followed <2yr, no events -> dropped (can't confirm indolent)
            ["TCGA-AA-0003", "Stage IB", "N0", "M0", "300", "[Not Available]"],
            # stage IIIA, node-positive, non-metastatic -> dropped (not a clean control)
            ["TCGA-AA-0004", "Stage IIIA", "N2", "M0", "1500", "[Not Available]"],
            # stage I, but distant recurrence on follow-up -> case (1)
            ["TCGA-AA-0005", "Stage IA", "N0", "M0", "1600", "[Not Available]"],
            # stage I, followed >2yr, but locoregional recurrence within window -> dropped
            ["TCGA-AA-0006", "Stage IB", "N0", "M0", "1500", "[Not Available]"],
        ],
        followup_rows=[
            ["TCGA-AA-0005", "Distant Metastasis", "500", "1600", "[Not Available]"],
            ["TCGA-AA-0006", "Locoregional Recurrence", "300", "1500", "[Not Available]"],
        ],
    )
    lab = true_stage_i_vs_distant(tmp_path, ["luad"], window_days=730)
    assert lab.loc["TCGA-AA-0001"] == 0        # true stage I control
    assert lab.loc["TCGA-AA-0002"] == 1        # M1 at diagnosis
    assert "TCGA-AA-0003" not in lab.index     # too little follow-up
    assert "TCGA-AA-0004" not in lab.index     # stage III, not a clean control
    assert lab.loc["TCGA-AA-0005"] == 1        # distant recurrence -> case
    assert "TCGA-AA-0006" not in lab.index     # locoregional recurrence in window
    assert set(lab.unique()) <= {0, 1}


def test_labels_concatenates_studies(tmp_path):
    _make_study(tmp_path, "luad",
                patient_rows=[["TCGA-AA-0001", "M1", "x", "x"]], followup_rows=[])
    _make_study(tmp_path, "lusc",
                patient_rows=[["TCGA-BB-0001", "M0", "x", "[Not Available]"]], followup_rows=[])
    both = distant_metastasis_labels(tmp_path, ["luad", "lusc"])
    assert both.loc["TCGA-AA-0001"] == 1
    assert both.loc["TCGA-BB-0001"] == 0
    assert len(both) == 2
