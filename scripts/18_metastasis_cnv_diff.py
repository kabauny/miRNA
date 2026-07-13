"""Copy-number differences: metastatic vs non-metastatic NSCLC (two-sided).

The expression screen (script 17) found no distant-met transcriptional signature
and a nodal signal that was really tumour stage. This asks the same question of
**copy number**: are any genes differentially **deleted or amplified** between
metastatic and non-metastatic tumours?

Distinct from the spared-deletion screens (14-16), which were one-sided (only
deep deletions *depleted* in metastatic disease). Here:

* both focal events -- deep deletion (HOMDEL, GISTIC -2) **and** high amplification
  (AMP, +2);
* **two-sided** -- an alteration enriched *or* depleted in the metastatic group
  (`associate.cmh_two_sided_screen`), subtype-adjusted;
* the same clean contrasts as script 17 -- distant met vs a **true stage I**
  indolent control, and nodal N+ vs N0 both whole-cohort and **within stage II**
  (to separate nodal spread from tumour stage).

Run:  python scripts/18_metastasis_cnv_diff.py --save-dir results
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import cmh_two_sided_screen
from mirna_tcga.biotab import distant_metastasis_labels, true_stage_i_vs_distant
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import (
    distant_metastasis,
    distant_metastasis_biotab,
    nodal_metastasis,
    nodal_metastasis_within_stage,
)
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import cna_matrix, protein_coding_map


def to_patient(df, sub_map):
    """samples x genes -> patients x genes (drop duplicate patients) + subtype."""
    df = df.copy()
    pats = sample_to_patient(df.index)
    sub = pd.Series([sub_map.get(s) for s in df.index], index=pats)
    df.index = pats
    keep = ~df.index.duplicated()
    return df[keep], sub[keep]


def run(name, alt, B, sub, label, min_freq, fdr):
    common = B.index.intersection(label.index)
    y = label.loc[common].astype(int)
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    res = cmh_two_sided_screen(B.loc[common], y, sub.loc[common], min_freq=min_freq)
    hits = res[res["q"] < fdr]
    enr = hits[hits["direction"] == "enriched"]
    dep = hits[hits["direction"] == "depleted"]
    print(f"\n{'-'*70}\n{name}  |  {alt}   ({n1} metastatic vs {n0} non-metastatic)\n{'-'*70}")
    print(f"genes differential at q<{fdr}: {len(hits)}  "
          f"(enriched in metastatic: {len(enr)}, depleted: {len(dep)})")
    if not res.empty:
        cols = ["freq_ref", "freq_idx", "n_idx", "or_mh", "direction", "p", "q"]
        print("top 12 by significance:")
        print(res.head(12)[cols].round(4).to_string())
    return res


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-freq", type=float, default=0.03,
                    help="min alteration freq in either group for a gene to be tested")
    ap.add_argument("--biotab-root", default=None)
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    dmet = distant_metastasis(clin)
    contrasts = [("nodal N+ vs N0", nodal_metastasis(clin)),
                 ("nodal N+ vs N0, stage II", nodal_metastasis_within_stage(clin, "II"))]
    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    if biotab_root and resolve_path(biotab_root).exists():
        root = resolve_path(biotab_root)
        dmet = distant_metastasis_biotab(clin, distant_metastasis_labels(root, keys))
        ts1 = true_stage_i_vs_distant(root, keys)
        print(f"distant met (biotab): {int(dmet.sum())} positive | "
              f"true-stage-I contrast: {int((ts1==1).sum())} vs {int((ts1==0).sum())}")
        contrasts = [("distant met vs true stage I", ts1),
                     ("distant met vs all M0", dmet)] + contrasts
    else:
        contrasts = [("distant met vs all M0", dmet)] + contrasts

    id2sym = protein_coding_map(client)
    print("Fetching deep deletions (HOMDEL) and amplifications (AMP) genome-wide ...")
    delB, sub = cna_matrix(client, cfg, id2sym, keys, "HOMDEL")
    ampB, _ = cna_matrix(client, cfg, id2sym, keys, "AMP")
    delP, subP = to_patient(delB, sub.to_dict())
    ampP, _ = to_patient(ampB, sub.to_dict())
    print(f"  deletions: {delP.shape[1]} genes, amplifications: {ampP.shape[1]} genes, "
          f"{delP.shape[0]} patients")

    out_rows = {}
    for name, label in contrasts:
        for alt, mat in (("deletion", delP), ("amplification", ampP)):
            res = run(name, alt, mat, subP, label, args.min_freq, args.fdr)
            out_rows[(name, alt)] = res

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        for (name, alt), res in out_rows.items():
            tag = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
            res.to_csv(out / f"cnv_{tag}_{alt}.csv")
        print(f"\nWrote results -> {out}/")


if __name__ == "__main__":
    main()
