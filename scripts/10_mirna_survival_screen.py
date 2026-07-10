"""Genome-wide miRNA survival screen in NSCLC (LUAD + LUSC).

The miRNA analogue of ``09_survival_screen.py``. miRNA expression is not on
cBioPortal for these cohorts, so it is loaded from **UCSC Xena** (already
log2-normalized), while overall-survival outcomes come from cBioPortal and are
matched to Xena samples by TCGA patient barcode.

Every mature miRNA is tested for an OS association with the same fast,
subtype-stratified Cox score test used for genes, then FDR-corrected; the top
hits are refit with a full Cox model for hazard ratios. (No pathway
over-representation: miRNAs are not members of the gene-level pathway sets.)

Requires the survival extra + Xena egress:
    pip install 'mirna-tcga[survival]'
Run:  python scripts/10_mirna_survival_screen.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.cohorts import cohort_study_keys, nsclc_clinical
from mirna_tcga.integrate import sample_to_patient
from mirna_tcga.mirbase import load_mimat_map
from mirna_tcga.screen import cox_score_screen
from mirna_tcga.survival import coerce_clinical, cox_model
from mirna_tcga.xena import XenaClient

DEFAULT_MATURE_FA = (
    "https://raw.githubusercontent.com/nf-core/test-datasets/smrnaseq/reference/mature.fa"
)


def _is_tumor(sample_id: str) -> bool:
    """True for TCGA primary/tumour samples (sample-type code < 10)."""
    parts = str(sample_id).split("-")
    if len(parts) < 4:
        return True
    code = parts[3][:2]
    return code.isdigit() and int(code) < 10


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", default=None)
    parser.add_argument("--fdr", type=float, default=0.05, help="BH q threshold for 'real impact'")
    parser.add_argument("--top", type=int, default=25, help="hits to refit with full Cox")
    parser.add_argument("--min-call-rate", type=float, default=0.8,
                        help="keep miRNAs detected in >= this fraction of samples (rest imputed)")
    parser.add_argument("--save-dir", default=None, help="directory to write result CSVs")
    args = parser.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    xena = XenaClient(**cfg.xena)

    # 1. survival outcomes for the NSCLC cohort (cBioPortal), keyed by patient
    clinical = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clinical).set_index("patientId")[["duration", "event", "cohort"]]
    print(f"NSCLC patients with overall-survival data: {len(surv)}")

    # 2. miRNA expression (Xena, already log2); tag each sample with its subtype
    frames = []
    for key in cohort_study_keys(cfg, "nsclc"):
        mat = xena.expression_matrix(cfg.xena_mirna[key])  # miRNAs x samples
        mat = mat.loc[:, [s for s in mat.columns if _is_tumor(s)]]
        m = mat.T  # samples x miRNAs
        m["_subtype"] = key.upper()
        frames.append(m)
        print(f"  {key.upper()}: {mat.shape[0]} miRNAs x {mat.shape[1]} tumour samples")
    X = pd.concat(frames, axis=0)
    subtype = X.pop("_subtype")

    # 3. align samples -> patients -> survival
    patient = pd.Series(sample_to_patient(X.index), index=X.index)
    keep = patient.isin(surv.index)
    X, patient, subtype = X[keep], patient[keep], subtype[keep]
    # one (tumour) sample per patient
    X = X[~patient.duplicated().to_numpy()]
    patient = patient[X.index]
    subtype = subtype[X.index]
    meta = surv.loc[patient]
    dur = pd.Series(meta["duration"].to_numpy(), index=X.index)
    evt = pd.Series(meta["event"].to_numpy(), index=X.index)

    # miRNA-seq has heavy dropout: keep miRNAs detected in >= min_call_rate of
    # samples, then impute the remaining missing values with the per-subtype
    # median (undetected == low expression) so the score test uses every sample.
    call_rate = X.notna().mean(axis=0)
    Xm = X.loc[:, call_rate >= args.min_call_rate]
    Xm = Xm.groupby(subtype.to_numpy()).transform(lambda c: c.fillna(c.median()))
    Xm = Xm.fillna(Xm.median()).dropna(axis=1)
    print(f"Analysis matrix: {Xm.shape[0]} patients x {Xm.shape[1]} miRNAs "
          f"(call-rate >= {args.min_call_rate}; {int(evt.sum())} deaths)\n")

    # 4. score-test screen across all miRNAs
    res = cox_score_screen(Xm, dur, evt, pd.Series(subtype.to_numpy(), index=X.index))
    # attach readable hsa-miR names (falls back to the accession if unavailable)
    try:
        cache = Path(args.save_dir or ".") / "genesets_cache" / "mature.fa"
        names = load_mimat_map(cfg.raw.get("mirbase_mature_fa", DEFAULT_MATURE_FA), cache)
    except Exception:  # pragma: no cover - network/offline fallback
        names = {}
    res.insert(0, "miRNA", [names.get(a, a) for a in res.index])

    hits = res[res["q"] < args.fdr]
    print(f"miRNAs with real survival impact (BH q < {args.fdr}): {len(hits)}")
    print(f"  higher expression -> WORSE survival: {int((hits['z'] > 0).sum())}")
    print(f"  higher expression -> BETTER survival: {int((hits['z'] < 0).sum())}\n")
    print("Top 20 by significance:")
    show = res.head(20).assign(dir=np.where(res.head(20)["z"] > 0, "risk", "prot"))
    print(show.round(4).to_string(), "\n")

    # 5. reportable hazard ratios via full Cox on the strongest hits
    print(f"Full Cox hazard ratios (per SD, subtype-adjusted) for top {args.top} hits:")
    rows = []
    for feat in hits.head(args.top).index:
        d = pd.DataFrame({
            "expr_z": (Xm[feat] - Xm[feat].mean()) / Xm[feat].std(),
            "is_lusc": (subtype.to_numpy() == "LUSC").astype(int),
            "duration": dur.to_numpy(),
            "event": evt.to_numpy(),
        }).dropna()
        try:
            s = cox_model(d, ["expr_z", "is_lusc"]).summary.loc["expr_z"]
            rows.append((names.get(feat, feat), feat, s["exp(coef)"],
                         s["exp(coef) lower 95%"], s["exp(coef) upper 95%"], s["p"]))
        except Exception:  # pragma: no cover - numerical edge cases
            continue
    hr = pd.DataFrame(
        rows, columns=["miRNA", "accession", "HR", "HR_low", "HR_high", "cox_p"]
    ).set_index("miRNA")
    print(hr.round(3), "\n")

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        res.to_csv(out / "mirna_survival_screen.csv")
        hr.to_csv(out / "mirna_top_hits_hazard_ratios.csv")
        print(f"Wrote results -> {out}/")


if __name__ == "__main__":
    main()
