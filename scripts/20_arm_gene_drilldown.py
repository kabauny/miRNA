"""Drill an arm-level CNV hit down to genes: is it coordinated cis-dosage?

Script 19 found chromosome-arm **17q gain** associated with nodal metastasis even
at matched stage. This asks whether that DNA gain shows up as coordinated
over-expression of 17q genes (an arm-level *cis*-dosage effect) and which genes
carry it -- reading the differential-expression tables written by script 17.

For a chosen arm it reports, per DE table: the share of the arm's genes shifted up
vs the genome-wide share (the dosage skew), the top up-regulated genes on the arm,
and a named gene of interest (default ERBB2, the famous 17q oncogene -- worth
checking whether it is actually the driver).

Needs the idmap extra (mygene) for gene->cytoband; the mapping is cached to
``<save-dir>/gene_cytoband.csv`` so reruns are offline.

Run:  python scripts/20_arm_gene_drilldown.py --arm 17q --gene ERBB2 --save-dir results
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import _bootstrap  # noqa: F401
import pandas as pd

from mirna_tcga.idmap import cytoband_to_arm, symbols_to_arm


def load_cytoband_arm(genes, cache: Path) -> pd.Series:
    """symbol -> arm, cached. Falls back to mygene on a cache miss.

    Reads either an ``arm`` column (this script's cache) or a ``cytoband`` column
    (e.g. a raw MyGene ``map_location`` dump), reducing the latter with
    :func:`cytoband_to_arm`.
    """
    if cache.exists():
        c = pd.read_csv(cache, index_col=0)
        arm = c["arm"] if "arm" in c else c["cytoband"].map(cytoband_to_arm)
        if not [g for g in genes if g not in arm.index]:
            return arm.dropna()
    arm = symbols_to_arm(genes)
    arm.rename("arm").to_frame().to_csv(cache)
    return arm


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", default="17q", help="chromosome arm to drill into")
    ap.add_argument("--gene", default="ERBB2", help="named gene to highlight")
    ap.add_argument("--save-dir", default="results", help="dir holding script-17 de_*_genes.csv")
    ap.add_argument("--pattern", default="de_nodal*_genes.csv",
                    help="glob of DE tables to drill (relative to save-dir)")
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()

    save = Path(args.save_dir)
    csvs = sorted(save.glob(args.pattern))
    if not csvs:
        raise SystemExit(f"No DE tables matching {args.pattern} in {save}/ -- run script 17 first.")

    tables = {p.stem.replace("de_", "").replace("_genes", ""): pd.read_csv(p, index_col=0) for p in csvs}
    all_genes = sorted(set().union(*[t.index for t in tables.values()]))
    arm = load_cytoband_arm(all_genes, save / "gene_cytoband.csv")

    for name, de in tables.items():
        d = de.join(arm.rename("arm"))
        on_arm = d[d["arm"] == args.arm]
        if on_arm.empty:
            print(f"\n[{name}] no {args.arm} genes mapped -- skipped")
            continue
        up_arm = (on_arm["z"] > 0).mean()
        up_all = (d["z"] > 0).mean()
        print(f"\n{'='*70}\n{name}\n{'='*70}")
        print(f"{args.arm}: {len(on_arm)} genes | median z {on_arm['z'].median():.2f} | "
              f"up (z>0) {up_arm:.0%} vs genome-wide {up_all:.0%}  "
              f"-> {'coordinated cis-dosage up' if up_arm - up_all > 0.05 else 'no arm-level skew'}")
        top = on_arm.sort_values("z", ascending=False).head(args.top)
        print(top[["z", "auc", "p", "q"]].round(4).to_string())
        if args.gene in on_arm.index:
            g = on_arm.loc[args.gene]
            rank = list(on_arm.sort_values("z", ascending=False).index).index(args.gene) + 1
            verdict = "IS elevated" if g["z"] > 2 else "is NOT the driver (unchanged)"
            print(f"{args.gene}: z={g['z']:.2f} auc={g['auc']:.3f} q={g['q']:.3g}  "
                  f"(rank {rank}/{len(on_arm)} on {args.arm}) -> {verdict}")


if __name__ == "__main__":
    main()
