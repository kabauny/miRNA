"""Kaplan-Meier: STC1-high vs STC1-low overall survival in NSCLC (writeup figure).

Produces the KM figure for the STC1 metastasis lead. High/low is split at the
**within-subtype median** (LUAD/LUSC separately, then pooled) so the contrast is
not just LUAD-vs-LUSC baseline. Kaplan-Meier estimator and the two-group log-rank
test are implemented directly (no lifelines dependency); figure via matplotlib.

Run:  python scripts/26_stc1_km_curve.py --gene STC1 --save-dir docs/figures
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import chi2  # noqa: E402

from mirna_tcga import load_config  # noqa: E402
from mirna_tcga.cbioportal import CBioPortalClient  # noqa: E402
from mirna_tcga.cohorts import cohort_study_keys, combined_expression, nsclc_clinical  # noqa: E402
from mirna_tcga.integrate import sample_to_patient  # noqa: E402
from mirna_tcga.survival import coerce_clinical  # noqa: E402


def km_curve(dur: np.ndarray, evt: np.ndarray):
    """Kaplan-Meier step function: returns (times, survival, at_risk_at_time)."""
    order = np.argsort(dur)
    dur, evt = dur[order], evt[order]
    times = np.unique(dur[evt == 1])
    surv, s, at_risk = [], 1.0, []
    for t in times:
        n = int((dur >= t).sum())
        d = int(((dur == t) & (evt == 1)).sum())
        s *= (1 - d / n) if n else 1.0
        surv.append(s)
        at_risk.append(n)
    return np.asarray(times), np.asarray(surv), np.asarray(at_risk)


def median_survival(times, surv):
    below = np.where(surv <= 0.5)[0]
    return float(times[below[0]]) if len(below) else np.nan


def logrank(dur, evt, grp):
    """Two-group log-rank test; returns (chi2, p, O1, E1)."""
    times = np.unique(dur[evt == 1])
    O1 = E1 = V = 0.0
    g1 = grp == 1
    for t in times:
        at = dur >= t
        n, n1 = int(at.sum()), int((at & g1).sum())
        d = int(((dur == t) & (evt == 1)).sum())
        d1 = int(((dur == t) & (evt == 1) & g1).sum())
        if n <= 1:
            continue
        O1 += d1
        E1 += d * n1 / n
        V += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    x2 = (O1 - E1) ** 2 / V if V > 0 else 0.0
    return x2, float(chi2.sf(x2, 1)), O1, E1


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--gene", default="STC1")
    ap.add_argument("--save-dir", default="docs/figures")
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    keys = cohort_study_keys(cfg, "nsclc")

    clin = nsclc_clinical(client, cfg, patient_level=True)
    surv = coerce_clinical(clin).set_index("patientId")[["duration", "event", "cohort"]]
    expr = combined_expression(client, cfg, [args.gene], keys)
    X = expr.T
    X["patientId"] = sample_to_patient(X.index)
    X = X[X["patientId"].isin(surv.index)].drop_duplicates("patientId")
    df = X.set_index("patientId").join(surv, how="inner")
    df = df.dropna(subset=[args.gene, "duration", "event"])

    # high/low split at the within-subtype median, then pool
    df["high"] = 0
    for c in df["cohort"].unique():
        m = df["cohort"] == c
        df.loc[m, "high"] = (df.loc[m, args.gene] > df.loc[m, args.gene].median()).astype(int)

    dur = df["duration"].to_numpy(float)
    evt = df["event"].to_numpy(int)
    grp = df["high"].to_numpy(int)
    x2, p, O1, E1 = logrank(dur, evt, grp)
    print(f"{args.gene}: {len(df)} patients, {int(evt.sum())} deaths | "
          f"high={int(grp.sum())} low={int((grp==0).sum())}")
    print(f"log-rank chi2={x2:.2f}, p={p:.4g}")

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    colors = {1: "#c0392b", 0: "#2c7fb8"}
    labels = {1: f"{args.gene}-high (n={int(grp.sum())})", 0: f"{args.gene}-low (n={int((grp==0).sum())})"}
    for g in (0, 1):
        m = grp == g
        t, s, _ = km_curve(dur[m], evt[m])
        t = np.concatenate([[0], t])
        s = np.concatenate([[1.0], s])
        ax.step(t, s, where="post", color=colors[g], lw=2, label=labels[g])
        med = median_survival(t, s)
        print(f"  {labels[g]}: median OS = {med:.1f} months" if not np.isnan(med)
              else f"  {labels[g]}: median OS not reached")

    ax.set_xlabel("Overall survival (months)")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.02)
    ax.set_title(f"NSCLC overall survival by {args.gene} expression\n"
                 f"(within-subtype median split; log-rank p = {p:.3g})")
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out = Path(args.save_dir); out.mkdir(parents=True, exist_ok=True)
    path = out / f"{args.gene.lower()}_km_nsclc.png"
    fig.savefig(path, dpi=150)
    print(f"Wrote figure -> {path}")


if __name__ == "__main__":
    main()
