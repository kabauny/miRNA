"""Genes *never* deep-deleted in stage IV NSCLC -- Biotab-enriched re-run.

The original spared-deletion screen (``scripts/14``) found genes deletable in
non-metastatic tumours but seen in **0 / 33** stage IV patients (the chr8p23
block). That "never deleted" call rested on only ~33 metastatic patients. This
script re-asks the strict question -- *which deletable genes have ZERO deep
deletions among stage IV patients?* -- on the **Biotab-enriched** distant-met
endpoint (126 stage IV vs 674 non-metastatic), and is explicit about two things
the small cohort hid:

1. **Power.** With 126 stage IV patients a gene deleted in fraction ``f`` of
   non-metastatic tumours is expected to be hit ``126 x f`` times by chance, so
   "never deleted" is only surprising for genes with an appreciable ``f``. The
   script sweeps the deletability floor and reports, at each, how many genes are
   strictly never-deleted -- the count collapses as the floor (and thus the
   expected hits) rises.
2. **Subtype.** Stage IV NSCLC is LUAD-heavy while non-metastatic tumours are
   LUSC-richer, so a gene deleted mostly in LUSC looks "spared" in the LUAD-heavy
   stage IV group for reasons that have nothing to do with metastasis. Every
   never-deleted call therefore carries its subtype-adjusted CMH ``q`` and a raw
   LUAD/LUSC breakdown of where its non-metastatic deletions fall.

Run:  python scripts/16_never_deleted_stage_iv.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.associate import cmh_depletion_screen
from mirna_tcga.biotab import distant_metastasis_labels
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.config import resolve_path
from mirna_tcga.endpoints import distant_metastasis, distant_metastasis_biotab
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.layers import deletion_matrix, protein_coding_map
from mirna_tcga.panels import CHR8P23_BLOCK

FLOOR_SWEEP = [0.02, 0.03, 0.04, 0.05, 0.06, 0.08]


def to_patient(df, sub_map):
    """samples x genes -> patients x genes (drop duplicate patients), + subtype."""
    df = df.copy()
    pats = sample_to_patient(df.index)
    sub = pd.Series([sub_map.get(s) for s in df.index], index=pats)
    df.index = pats
    keep = ~df.index.duplicated()
    return df[keep], sub[keep]


def subtype_split(delP, subP, y, genes):
    """Per gene: where do the non-metastatic (M0) deletions fall by subtype?"""
    m0 = y.index[y == 0]
    rows = []
    for g in genes:
        col = delP[g]
        del_m0 = col.loc[m0]
        sub_m0 = subP.loc[m0]
        n_luad = int(((del_m0 == 1) & (sub_m0 == "LUAD")).sum())
        n_lusc = int(((del_m0 == 1) & (sub_m0 == "LUSC")).sum())
        rows.append((g, n_luad, n_lusc))
    return pd.DataFrame(rows, columns=["gene", "m0_del_LUAD", "m0_del_LUSC"]).set_index("gene")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--min-ref-freq", type=float, default=0.03,
                    help="deletability floor: min deep-deletion freq in non-metastatic tumours")
    ap.add_argument("--biotab-root", default=None,
                    help="local BCR Biotab GDCdata dir (default: config biotab.root if present)")
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    dmet = distant_metastasis(clin)
    biotab_root = args.biotab_root or (cfg.raw.get("biotab") or {}).get("root")
    if biotab_root and resolve_path(biotab_root).exists():
        bt = distant_metastasis_labels(resolve_path(biotab_root), keys)
        n0 = int(dmet.sum())
        dmet = distant_metastasis_biotab(clin, bt)
        print(f"Stage IV label: Biotab-enriched {n0} -> {int(dmet.sum())} positive")
    else:
        print(f"Stage IV label: pathologic M1 only ({int(dmet.sum())} positive) -- no Biotab")

    id2sym = protein_coding_map(client)
    print("Fetching deep deletions (HOMDEL) genome-wide ...")
    delB, sub = deletion_matrix(client, cfg, id2sym, keys)
    delP, subP = to_patient(delB, sub.to_dict())

    common = delP.index.intersection(dmet.index)
    delP, subP, y = delP.loc[common], subP.loc[common], dmet.loc[common].astype(int)
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    luad_iv = float((subP[y == 1] == "LUAD").mean())
    luad_m0 = float((subP[y == 0] == "LUAD").mean())
    print(f"\nCohort: {n1} stage IV vs {n0} non-metastatic, {delP.shape[1]} genes.")
    print(f"Subtype skew: stage IV is {luad_iv:.0%} LUAD vs non-metastatic {luad_m0:.0%} LUAD "
          f"(the confound the CMH q adjusts for).")

    # Full genome-wide depletion screen (no floor), subtype-stratified.
    res = cmh_depletion_screen(delP, y, subP, min_ref_freq=0.0)

    # 1) Power sweep: how many genes are strictly never-deleted at each floor?
    print("\nStrictly never-deleted-in-stage-IV genes vs deletability floor:")
    print(f"  {'floor':>6}  {'deletable':>9}  {'never-del':>9}  {'expected hits @floor':>20}")
    for f in FLOOR_SWEEP:
        sub_res = res[res["freq_ref"] >= f]
        never = sub_res[sub_res["n_idx"] == 0]
        print(f"  {f:>6.0%}  {len(sub_res):>9}  {len(never):>9}  {f * n1:>17.1f} hits")

    # 2) The never-deleted list at the chosen floor, ranked by expected hits.
    tested = res[res["freq_ref"] >= args.min_ref_freq].copy()
    never = tested[tested["n_idx"] == 0].copy()
    never["expected_iv"] = (never["freq_ref"] * n1).round(2)
    split = subtype_split(delP, subP, y, list(never.index))
    never = never.join(split)
    never = never.sort_values("expected_iv", ascending=False)
    cols = ["freq_ref", "expected_iv", "m0_del_LUAD", "m0_del_LUSC", "or_mh", "p", "q"]
    print(f"\nGenes NEVER deep-deleted in stage IV (0/{n1}) yet deletable in "
          f">= {args.min_ref_freq:.0%} of non-metastatic tumours "
          f"({len(never)} genes, ranked by expected hits):")
    print(never[cols].head(30).round(4).to_string() if not never.empty
          else "  (none -- every deletable gene is hit at least once in stage IV)")

    # 3) The original chr8p23 block: what happened to the 0/33 call?
    blk = [g for g in CHR8P23_BLOCK if g in res.index]
    if blk:
        b = res.loc[blk, ["freq_ref", "freq_idx", "n_idx", "or_mh", "q"]].copy()
        b = b.join(subtype_split(delP, subP, y, blk)).sort_values("freq_ref", ascending=False)
        print(f"\nchr8p23 block (the original 0/33 hit) in the {n1}-patient cohort:")
        print(b.round(4).to_string())

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        never[cols + ["expected_iv"]].to_csv(out / "never_deleted_stage_iv.csv")
        print(f"\nWrote results -> {out}/never_deleted_stage_iv.csv")


if __name__ == "__main__":
    main()
