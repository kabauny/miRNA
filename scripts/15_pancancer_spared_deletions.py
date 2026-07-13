"""Is 'spared from deletion in metastasis' lung-specific? A pan-cancer test.

Runs the negative-selection screen of ``scripts/14`` across the TCGA PanCancer
cohorts that actually contain metastatic (stage IV / M1) patients -- bladder
(BLCA, 136 M1), colorectal (COADREAD, 84), kidney (KIRC, 84), stomach (STAD, 45)
and lung -- which are far better powered than lung alone (~34 M1).

Two outputs:
  * **Pooled** Cochran-Mantel-Haenszel depletion test **stratified by cancer
    type**: which genes are deep-deleted in non-metastatic tumours but spared in
    metastatic disease *consistently across cancers*?
  * The **chr8p23 block** (the lung hit) tracked per cancer, to see whether it
    replicates or is lung-specific.

Optionally counts a gene as inactivated by deep deletion **or** truncating
mutation (``--protection``), to test whether spared deletion is simply replaced
by loss-of-function mutation.

Run:  python scripts/15_pancancer_spared_deletions.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import cmh_depletion_screen, ranksum_screen
from mirna_tcga.biotab import distant_metastasis_labels
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import distant_metastasis, distant_metastasis_biotab
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import TRUNCATING_CLASSES, _binary_matrix, protein_coding_map
from mirna_tcga.panels import CHR8P23_BLOCK

DEFAULT_STUDIES = ["blca", "coadread", "kirc", "stad", "luad", "lusc", "skcm"]
SUFFIX = "_tcga_pan_can_atlas_2018"

# cBioPortal study key -> local BCR Biotab folder code(s) (TCGA-<CODE>). Most map
# 1:1; the merged colorectal study maps to COAD (the local download has no READ),
# and cohorts with no local Biotab (e.g. skcm) fall back to the cBioPortal M1 call.
BIOTAB_ALIASES = {"coadread": ["coad"]}


def _tumor(s):
    p = str(s).split("-")
    return len(p) < 4 or (p[3][:2].isdigit() and int(p[3][:2]) < 10)


def to_patient(df):
    df = df.loc[[s for s in df.index if _tumor(s)]].copy()
    df.index = sample_to_patient(df.index)
    return df[~df.index.duplicated()]


def study_inactivation(client, study, id2sym, protection):
    """patient x gene 0/1 deep-deletion (optionally OR truncating-mutation) matrix."""
    cna_sl = f"{study}_cna"
    samples = client.sample_list_ids(cna_sl)
    if not samples:
        return pd.DataFrame()
    ev = client.discrete_cna_events(f"{study}_gistic", cna_sl, "HOMDEL")
    inact = _binary_matrix(ev, samples, id2sym)
    if protection:
        seq_sl = f"{study}_sequenced"
        muts = client.mutation_events(f"{study}_mutations", list(id2sym), seq_sl)
        if not muts.empty and "mutationType" in muts:
            muts = muts[muts["mutationType"].astype(str).str.lower().isin(TRUNCATING_CLASSES)]
        mut_mat = _binary_matrix(muts, samples, id2sym)  # same CNA samples as denominator
        inact = inact.add(mut_mat.reindex(index=inact.index, columns=inact.columns, fill_value=0),
                          fill_value=0).clip(upper=1).astype("int8")
    return to_patient(inact)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--studies", nargs="+", default=DEFAULT_STUDIES)
    ap.add_argument("--min-ref-freq", type=float, default=0.03)
    ap.add_argument("--protection", action="store_true",
                    help="count deletion OR truncating mutation as inactivation")
    ap.add_argument("--biotab-root", default=None,
                    help="path to local BCR Biotab GDCdata dir; enriches each cohort's "
                         "distant-met endpoint where a matching download exists "
                         "(default: config biotab.root if it exists on disk)")
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config()
    client = CBioPortalClient(**cfg.cbioportal)
    id2sym = protein_coding_map(client)
    mode = "inactivation (del + truncating mut)" if args.protection else "deep deletion"
    print(f"Pan-cancer spared-{mode} screen across: {', '.join(args.studies)}\n")

    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    biotab_root = resolve_path(biotab_root) if biotab_root else None
    if biotab_root and biotab_root.exists():
        print(f"Biotab enrichment: on ({biotab_root})\n")
    else:
        biotab_root = None

    mats, ys, cancer, block_rows = [], [], [], []
    for st in args.studies:
        study = st + SUFFIX
        clin = client.get_clinical_data(study, patient_level=True)
        if clin.empty:
            print(f"  {st}: no clinical -- skipped")
            continue
        dm = distant_metastasis(clin)
        if biotab_root is not None:
            bt = distant_metastasis_labels(biotab_root, BIOTAB_ALIASES.get(st, [st]))
            if not bt.empty:
                n0 = int(dm.sum())
                dm = distant_metastasis_biotab(clin, bt)
                print(f"  {st.upper():9s} biotab: M1 {n0} -> {int(dm.sum())}")
        B = study_inactivation(client, study, id2sym, args.protection)
        common = B.index.intersection(dm.index)
        if len(common) < 20 or dm.loc[common].sum() < 5:
            print(f"  {st}: too few metastatic ({int(dm.loc[common].sum())}) -- skipped")
            continue
        B, y = B.loc[common], dm.loc[common].astype(int)
        mats.append(B)
        ys.append(y)
        cancer += [st.upper()] * len(common)
        # chr8p23 block rate in this cancer
        blk = [g for g in CHR8P23_BLOCK if g in B.columns]
        score = (B[blk].sum(axis=1) > 0).astype(int) if blk else pd.Series(0, index=B.index)
        f1 = score[y == 1].mean() if (y == 1).any() else np.nan
        f0 = score[y == 0].mean() if (y == 0).any() else np.nan
        block_rows.append((st.upper(), int(y.sum()), int((y == 0).sum()), round(f0, 3), round(f1, 3)))
        print(f"  {st.upper():9s} n={len(common):4d}  M1={int(y.sum()):3d}  "
              f"8p23-block deleted: M0={f0:.1%}  M1={f1:.1%}")

    if not mats:
        raise SystemExit("No usable cohorts.")

    # pooled matrix (union of genes), stratified by cancer type
    allgenes = sorted(set().union(*[m.columns for m in mats]))
    Bpool = pd.concat([m.reindex(columns=allgenes, fill_value=0) for m in mats], axis=0).astype("int8")
    ypool = pd.concat(ys)
    strata = pd.Series(cancer, index=Bpool.index)
    print(f"\nPooled: {Bpool.shape[0]} patients ({int(ypool.sum())} metastatic) x "
          f"{Bpool.shape[1]} genes, stratified by cancer type")

    res = cmh_depletion_screen(Bpool, ypool, strata, min_ref_freq=args.min_ref_freq)
    spared = res[res["or_mh"] < 1]
    print(f"\nGenes spared from {mode} in metastatic disease, PAN-CANCER "
          f"(top 20 of {len(res)} tested):")
    print(spared.head(20).round(4).to_string())

    # chr8p23 block, pooled block-level test (stratified rank-sum on block count)
    blk = [g for g in CHR8P23_BLOCK if g in Bpool.columns]
    block_count = Bpool[blk].sum(axis=1).to_frame("chr8p23_block")
    bt = ranksum_screen(block_count, ypool, strata)
    print("\nchr8p23 block, pooled (z<0 => less deleted in metastatic):")
    print(bt.round(4).to_string())
    print("\nPer-cancer 8p23-block deletion rate (M0 vs M1):")
    print(pd.DataFrame(block_rows, columns=["cancer", "M1_n", "M0_n", "M0_rate", "M1_rate"]).to_string(index=False))

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        res.to_csv(out / f"pancancer_spared_{'protection' if args.protection else 'deletion'}.csv")
        print(f"\nWrote results -> {out}/")


if __name__ == "__main__":
    main()
