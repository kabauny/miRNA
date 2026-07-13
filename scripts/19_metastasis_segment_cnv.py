"""Segment-level copy number: metastatic vs non-metastatic NSCLC.

The gene-level GISTIC screens (14-18) test *focal* deep events and miss broad,
whole-arm copy-number change. This uses the raw DNAcopy **segments** (continuous
log2 ratios, hg19) to ask the burden / arm-level questions those cannot:

1. **Fraction genome altered (FGA)** -- do metastatic tumours have more of their
   genome copy-number-altered overall (genomic instability)? One score per sample,
   tested with a subtype-stratified rank-sum.
2. **Arm-level** -- is any of the 39 autosomal arms differentially gained or lost?
   Length-weighted mean log2 per arm, one rank-sum per arm, BH-FDR across arms.

Same clean contrasts as scripts 17/18: distant met vs a true-stage-I indolent
control, and nodal N+ vs N0 both whole-cohort and within stage II (stage held
constant). Primary tumours only (-01).

Run:  python scripts/19_metastasis_segment_cnv.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import ranksum_screen
from mirna_tcga.biotab import distant_metastasis_labels, true_stage_i_vs_distant
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, combined_clinical, nsclc_clinical  # noqa: F401
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import (
    distant_metastasis,
    distant_metastasis_biotab,
    nodal_metastasis,
    nodal_metastasis_within_stage,
)
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.segments import arm_level_means, fraction_genome_altered


def to_patient(obj):
    """Reindex a sampleId-keyed Series/DataFrame to patientId (drop dup patients)."""
    obj = obj.copy()
    obj.index = sample_to_patient(obj.index)
    return obj[~obj.index.duplicated()]


def run_fga(name, fga, subP, label):
    common = fga.index.intersection(label.index)
    y = label.loc[common].astype(int)
    strata = subP.loc[common]
    res = ranksum_screen(fga.loc[common].to_frame("FGA"), y, strata)
    r = res.loc["FGA"]
    m1 = fga.loc[common][y == 1].median()
    m0 = fga.loc[common][y == 0].median()
    print(f"  FGA: median metastatic {m1:.3f} vs non-metastatic {m0:.3f}  "
          f"(z={r['z']:.2f}, AUC={r['auc']:.3f}, p={r['p']:.2g})")
    return r


def run_arms(name, armP, subP, label, fdr):
    common = armP.index.intersection(label.index)
    y = label.loc[common].astype(int)
    strata = subP.loc[common]
    X = armP.loc[common].fillna(0.0)          # unassayed arm -> neutral
    res = ranksum_screen(X, y, strata)
    hits = res[res["q"] < fdr]
    gained = hits[hits["z"] > 0]
    lost = hits[hits["z"] < 0]
    print(f"  arms differential at q<{fdr}: {len(hits)}  "
          f"(gained in metastatic: {len(gained)}, lost: {len(lost)})")
    if not res.empty:
        print(res.head(6)[["z", "auc", "p", "q"]].round(4).to_string())
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--fga-threshold", type=float, default=0.2,
                    help="|log2 ratio| above which a segment counts as altered")
    ap.add_argument("--studies", nargs="+", default=None,
                    help="study keys (default: nsclc = luad+lusc); pass one subtype "
                         "to test whether an arm signal is subtype-specific")
    ap.add_argument("--biotab-root", default=None)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = args.studies or cohort_study_keys(cfg, "nsclc")
    print(f"Cohort: {', '.join(k.upper() for k in keys)}")

    clin = combined_clinical(client, cfg, keys, patient_level=True)
    dmet = distant_metastasis(clin)
    contrasts = [("nodal N+ vs N0", nodal_metastasis(clin)),
                 ("nodal N+ vs N0, stage II", nodal_metastasis_within_stage(clin, "II"))]
    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    if biotab_root and resolve_path(biotab_root).exists():
        root = resolve_path(biotab_root)
        dmet = distant_metastasis_biotab(clin, distant_metastasis_labels(root, keys))
        ts1 = true_stage_i_vs_distant(root, keys)
        contrasts = [("distant met vs true stage I", ts1),
                     ("distant met vs all M0", dmet)] + contrasts
    else:
        contrasts = [("distant met vs all M0", dmet)] + contrasts

    # Fetch segments for primary-tumour samples of both studies.
    seg_frames, subtype = [], {}
    for key in keys:
        sids = [s for s in client.sample_list_ids(cfg.all_samples_list(key))
                if s.split("-")[-1].startswith("01")]           # primary tumours
        print(f"Fetching segments for {len(sids)} {key.upper()} samples ...")
        seg_frames.append(client.copy_number_segments(cfg.studies[key], sids))
        subtype.update({s: key.upper() for s in sids})
    seg = pd.concat(seg_frames, ignore_index=True)
    print(f"  {len(seg)} segments across {seg['sampleId'].nunique()} samples")

    fga = to_patient(fraction_genome_altered(seg, args.fga_threshold))
    arm = to_patient(arm_level_means(seg))
    subP = to_patient(pd.Series(subtype))
    print(f"  FGA + {arm.shape[1]}-arm matrix for {len(fga)} patients "
          f"(median FGA {fga.median():.3f})")

    fga_rows, arm_res = {}, {}
    for name, label in contrasts:
        print(f"\n{'='*70}\n{name}\n{'='*70}")
        fga_rows[name] = run_fga(name, fga, subP, label)
        arm_res[name] = run_arms(name, arm, subP, label, args.fdr)

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(fga_rows).T.to_csv(out / "segment_fga_by_contrast.csv")
        for name, res in arm_res.items():
            tag = name.lower().replace(" ", "_").replace(",", "").replace("+", "p")
            res.to_csv(out / f"segment_arms_{tag}.csv")
        print(f"\nWrote results -> {out}/")


if __name__ == "__main__":
    main()
