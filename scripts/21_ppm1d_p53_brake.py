"""Is PPM1D (WIP1) the 17q-amplicon p53 brake? A functional test.

The LUSC nodal-metastasis 17q gain over-expresses the whole arm (script 20). ERBB2
is not the driver; a better candidate is **PPM1D** (WIP1, 17q23.2) -- a phosphatase
that dephosphorylates and **inactivates p53**. If PPM1D is the functional target of
the gain, then in tumours with **wild-type TP53** (where p53 can still act) higher
PPM1D should mean a **suppressed p53 transcriptional program** -- its direct
targets (p21/CDKN1A, MDM2, GADD45A, BBC3/PUMA, ...) should go *down*. Where TP53 is
already mutated, p53 is dead and PPM1D cannot matter -- so the effect must be
**TP53-wild-type-specific**. That stratification is the whole test: LUSC is ~80%
TP53-mutant, so an unstratified correlation would be badly confounded.

Design (per subtype, since PPM1D amplification is LUSC-enriched but TP53-WT is
LUAD-enriched -- both angles are informative):
  1. Expression of PPM1D + 24 chr17-free direct p53 targets (`panels.P53_TARGETS`)
     + FDXR (a p53 target that sits on 17q -- a dosage positive control) + TP53.
  2. TP53 mutation status; PPM1D discrete copy number (amp/gain frequency).
  3. p53-target signature = mean of per-gene z-scored expression.
  4. Within TP53-WT vs TP53-mut separately: Spearman(PPM1D expr, p53 signature) and
     a PPM1D-high vs -low rank-sum. Prediction: negative in TP53-WT, ~null in mut.
  5. FDXR control: Spearman(PPM1D, FDXR) should be *positive* (both ride 17q
     dosage) -- proving the arm is co-expressed while the p53 targets move opposite.
  6. Same test using PPM1D **copy number** (amp/gain vs neutral) as the exposure.

Run:  python scripts/21_ppm1d_p53_brake.py --save-dir results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from mirna_tcga import load_config
from mirna_tcga.cbioportal import CBioPortalClient
from mirna_tcga.panels import P53_TARGETS, P53_TARGETS_ON_17Q

REGULATOR = "PPM1D"


def zscore(df: pd.DataFrame) -> pd.DataFrame:
    """Per-gene (row) z-score across samples."""
    mu = df.mean(axis=1)
    sd = df.std(axis=1).replace(0, np.nan)
    return df.sub(mu, axis=0).div(sd, axis=0)


def _corr(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    d = pd.concat([x, y], axis=1).dropna()
    if len(d) < 10:
        return float("nan"), float("nan"), len(d)
    r, p = spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    return r, p, len(d)


def _highlow(reg: pd.Series, sig: pd.Series) -> tuple[float, float, float]:
    """Rank-sum of the signature in top vs bottom PPM1D tertile; returns (hi, lo, p)."""
    d = pd.concat([reg.rename("r"), sig.rename("s")], axis=1).dropna()
    if len(d) < 15:
        return float("nan"), float("nan"), float("nan")
    lo_c, hi_c = d["r"].quantile(1 / 3), d["r"].quantile(2 / 3)
    lo, hi = d[d["r"] <= lo_c]["s"], d[d["r"] >= hi_c]["s"]
    if len(lo) < 5 or len(hi) < 5:
        return float("nan"), float("nan"), float("nan")
    p = mannwhitneyu(hi, lo, alternative="two-sided").pvalue
    return hi.median(), lo.median(), p


def analyse_subtype(client, cfg, key: str) -> dict:
    study = cfg.studies[key]
    genes = [REGULATOR, "TP53", *P53_TARGETS, *P53_TARGETS_ON_17Q]
    expr = client.expression_matrix(cfg.mrna_profile(key), genes, cfg.all_samples_list(key))
    if expr.empty:
        print(f"  {key.upper()}: no expression -- skipped")
        return {}
    tp53_mut = client.mutated_samples(cfg.mutation_profile(key), "TP53", cfg.sequenced_samples_list(key))
    cna = client.cna_matrix(cfg.cna_profile(key), [REGULATOR], cfg.cna_samples_list(key))

    z = zscore(expr)
    targets = [g for g in P53_TARGETS if g in z.index]
    sig = z.loc[targets].mean(axis=0)                 # per-sample p53-target score
    reg = expr.loc[REGULATOR]                          # PPM1D expression (log2)
    samples = expr.columns

    # copy-number summary + amp/gain vs neutral exposure
    cn = cna.loc[REGULATOR] if not cna.empty and REGULATOR in cna.index else pd.Series(dtype=float)
    amp = float((cn >= 2).mean()) if len(cn) else float("nan")
    gain = float((cn >= 1).mean()) if len(cn) else float("nan")

    seq_samples = set(client.sample_list_ids(cfg.sequenced_samples_list(key)))
    wt = [s for s in samples if s in seq_samples and s not in tp53_mut]
    mt = [s for s in samples if s in tp53_mut]

    print(f"\n{'='*70}\n{key.upper()}  (study {study})\n{'='*70}")
    print(f"samples: {len(samples)} expr | TP53 sequenced {len(seq_samples)} "
          f"-> mut {len(mt)} ({len(mt)/max(len(seq_samples),1):.0%}), WT {len(wt)}")
    print(f"PPM1D copy number: amp(+2) {amp:.1%}, gain(>=+1) {gain:.1%}")

    out = {"key": key, "amp": amp, "gain": gain, "n_wt": len(wt), "n_mut": len(mt)}
    print(f"\np53-target signature ({len(targets)} chr17-free targets) vs PPM1D expression:")
    for grp, cols in (("TP53-WT", wt), ("TP53-mut", mt)):
        r, p, n = _corr(reg[cols], sig[cols])
        hi, lo, php = _highlow(reg[cols], sig[cols])
        arrow = "  <-- predicted r<0" if grp == "TP53-WT" else ""
        print(f"  {grp:9s} n={n:3d}: Spearman r={r:+.3f} p={p:.2g}   "
              f"tertile hi={hi:+.2f} lo={lo:+.2f} rank-sum p={php:.2g}{arrow}")
        out[f"r_{grp}"] = r
        out[f"p_{grp}"] = p

    # positive control: FDXR (p53 target ON 17q) should track PPM1D UP (co-dosage)
    if P53_TARGETS_ON_17Q[0] in z.index:
        rc, pc, _ = _corr(reg[wt], expr.loc[P53_TARGETS_ON_17Q[0], wt])
        print(f"  FDXR (17q dosage control), TP53-WT: Spearman r={rc:+.3f} p={pc:.2g} "
              f"(expected r>0: 17q co-expression)")

    # CDKN1A / p21 alone -- the canonical readout
    if "CDKN1A" in z.index:
        r21, p21v, _ = _corr(reg[wt], expr.loc["CDKN1A", wt])
        print(f"  CDKN1A/p21 alone, TP53-WT: Spearman r={r21:+.3f} p={p21v:.2g}")

    # copy-number exposure: PPM1D amp/gain vs neutral, p53 signature (TP53-WT)
    if len(cn):
        cn_wt = cn.reindex(wt).dropna()
        hi = sig.reindex(cn_wt[cn_wt >= 1].index).dropna()
        neu = sig.reindex(cn_wt[cn_wt <= 0].index).dropna()
        if len(hi) >= 5 and len(neu) >= 5:
            pcn = mannwhitneyu(hi, neu, alternative="two-sided").pvalue
            print(f"  PPM1D gain/amp vs neutral (TP53-WT): p53 sig "
                  f"{hi.median():+.2f} (n={len(hi)}) vs {neu.median():+.2f} (n={len(neu)}), "
                  f"rank-sum p={pcn:.2g}")

    # per-target detail (TP53-WT): which p53 targets move with PPM1D
    rows = []
    for g in targets:
        r, p, _ = _corr(reg[wt], expr.loc[g, wt])
        rows.append((g, r, p))
    det = pd.DataFrame(rows, columns=["target", "r", "p"]).set_index("target").sort_values("r")
    out["detail"] = det
    print(f"\n  most-suppressed p53 targets as PPM1D rises (TP53-WT, top 8 by r<0):")
    print(det.head(8).round(3).to_string())
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--studies", nargs="+", default=["lusc", "luad"],
                    help="subtypes to test separately (default: lusc luad)")
    ap.add_argument("--save-dir", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    client = CBioPortalClient(**cfg.cbioportal)
    print(f"Testing {REGULATOR} (WIP1, 17q23) as a p53 brake in: "
          f"{', '.join(k.upper() for k in args.studies)}")

    results = [analyse_subtype(client, cfg, k) for k in args.studies]

    if args.save_dir:
        out = Path(args.save_dir)
        out.mkdir(parents=True, exist_ok=True)
        for r in results:
            if r and "detail" in r:
                r["detail"].to_csv(out / f"ppm1d_p53_targets_{r['key']}.csv")
        print(f"\nWrote per-target tables -> {out}/")


if __name__ == "__main__":
    main()
