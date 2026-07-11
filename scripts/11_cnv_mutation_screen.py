"""Copy-number-deletion and mutation screens vs survival and metastasis (NSCLC).

For every recurrently **deep-deleted** gene (GISTIC -2) and every recurrently
**mutated** gene, test association with three endpoints across the NSCLC cohort
(LUAD + LUSC):

* overall survival        -- subtype-stratified Cox score / log-rank
  (:func:`mirna_tcga.screen.cox_score_screen` on the 0/1 alteration matrix)
* distant metastasis (M1) -- Fisher exact (:func:`mirna_tcga.associate.fisher_screen`)
* nodal metastasis (N+)   -- Fisher exact

Deep deletions and mutations are pulled from cBioPortal as sparse event lists
(fast), then turned into 0/1 sample x gene matrices over the profiled samples.

Requires the survival extra:  pip install 'mirna-tcga[survival]'
Run:  python scripts/11_cnv_mutation_screen.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import fisher_screen
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.endpoints import distant_metastasis, nodal_metastasis
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import deletion_matrix, mutation_matrix, protein_coding_map
from mirna_tcga.screen import cox_score_screen
from mirna_tcga.survival import coerce_clinical


def screen_layer(B, subtype, surv, endpoints, min_freq, label, top, save_dir):
    """Run OS + metastasis screens for one 0/1 alteration matrix."""
    freq = B.mean(axis=0)
    B = B.loc[:, freq >= min_freq]
    print(f"\n=== {label}: {B.shape[1]} genes altered in >= {min_freq:.0%} of samples ===")
    patient = pd.Series(sample_to_patient(B.index), index=B.index)

    # survival (subtype-stratified log-rank on the 0/1 matrix)
    s = surv.reindex(patient).dropna(subset=["duration", "event"])
    idx = patient.index[patient.isin(s.index)]
    idx = idx[~patient[idx].duplicated().to_numpy()]
    dur = pd.Series(s.loc[patient[idx], "duration"].to_numpy(), index=idx)
    evt = pd.Series(s.loc[patient[idx], "event"].to_numpy(), index=idx)
    os_res = cox_score_screen(B.loc[idx], dur, evt, pd.Series(subtype.reindex(idx).to_numpy(), index=idx))
    os_res = os_res.rename(columns={"z": "os_z", "p": "os_p", "q": "os_q"})
    print(f"\n[{label}] top alterations vs OVERALL SURVIVAL (z>0 = worse):")
    print(os_res.head(10).round(4).to_string())

    tables = {"survival": os_res}
    for ep_name, ep in endpoints.items():
        yy = ep.reindex(patient[idx]).dropna()
        sub = idx[patient[idx].isin(yy.index)]
        y = pd.Series(ep.reindex(patient[sub]).to_numpy(), index=sub)
        fr = fisher_screen(B.loc[sub], y)
        tables[ep_name] = fr
        pos = int(y.sum())
        print(f"\n[{label}] top alterations vs {ep_name.upper()} ({pos}+/{len(y)}; OR>1 = more metastasis):")
        print(fr.head(10).round(4).to_string())

    if save_dir:
        out = Path(save_dir)
        out.mkdir(parents=True, exist_ok=True)
        for name, t in tables.items():
            t.to_csv(out / f"{label}_{name}.csv")
    return tables


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=None)
    parser.add_argument("--min-del-freq", type=float, default=0.02, help="min deep-deletion frequency")
    parser.add_argument("--min-mut-freq", type=float, default=0.03, help="min mutation frequency")
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--save-dir", default=None)
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)

    clin = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clin).set_index("patientId")[["duration", "event"]]
    endpoints = {"distant_met": distant_metastasis(clin), "nodal_met": nodal_metastasis(clin)}
    print(f"Endpoints: OS={len(surv)} patients | "
          f"distant_met={int(endpoints['distant_met'].sum())}+/{len(endpoints['distant_met'])} | "
          f"nodal_met={int(endpoints['nodal_met'].sum())}+/{len(endpoints['nodal_met'])}")

    id2sym = protein_coding_map(client)
    keys = cohort_study_keys(cfg, "nsclc")

    print("\nFetching deep deletions (HOMDEL) ...")
    delB, del_sub = deletion_matrix(client, cfg, id2sym, keys)
    print(f"  deep-deletion matrix: {delB.shape[0]} samples x {delB.shape[1]} genes")
    print("Fetching mutations ...")
    mutB, mut_sub = mutation_matrix(client, cfg, id2sym, keys)
    print(f"  mutation matrix: {mutB.shape[0]} samples x {mutB.shape[1]} genes")

    screen_layer(delB, del_sub, surv, endpoints, args.min_del_freq, "deletion", args.top, args.save_dir)
    screen_layer(mutB, mut_sub, surv, endpoints, args.min_mut_freq, "mutation", args.top, args.save_dir)
    if args.save_dir:
        print(f"\nWrote results -> {args.save_dir}/")


if __name__ == "__main__":
    main()
