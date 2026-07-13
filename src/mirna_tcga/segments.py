"""Segment-level copy-number summaries from DNAcopy segments (hg19).

The discrete GISTIC screens (scripts 14-18) test *focal* deep events per gene and
lose segment length + continuous magnitude. Raw copy-number segments
(:meth:`mirna_tcga.cbioportal.CBioPortalClient.copy_number_segments`) keep both,
so they support the two broad-CNV summaries a per-gene binary test cannot:

* :func:`fraction_genome_altered` -- one genomic-instability score per sample
  (fraction of the assayed genome whose |log2 ratio| exceeds a threshold);
* :func:`arm_level_means` -- a sample x chromosome-arm matrix of length-weighted
  mean log2 ratio (whole-arm gains / losses / aneuploidy).

Coordinates are hg19 (cBioPortal PanCancer). Each segment is assigned to the arm
containing its midpoint via approximate hg19 centromere midpoints; the five
acrocentric short arms (13p/14p/15p/21p/22p) carry no assayable material and are
excluded, leaving 39 autosomal arms.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Approximate hg19 centromere midpoints (bp), autosomes 1-22 + X ('23').
_CENTROMERE = {
    "1": 125_000_000, "2": 93_300_000, "3": 91_000_000, "4": 50_400_000,
    "5": 48_400_000, "6": 61_000_000, "7": 59_900_000, "8": 45_600_000,
    "9": 49_000_000, "10": 40_200_000, "11": 53_700_000, "12": 35_800_000,
    "13": 17_900_000, "14": 17_600_000, "15": 19_000_000, "16": 36_600_000,
    "17": 24_000_000, "18": 17_200_000, "19": 26_500_000, "20": 27_500_000,
    "21": 13_200_000, "22": 14_700_000, "23": 60_600_000,
}
# Acrocentric short arms with no assayable material -- excluded from arm analysis.
_ACROCENTRIC_P = {"13p", "14p", "15p", "21p", "22p"}
# The 39 standard autosomal arms, ordered.
ARMS = [c + a for c in [str(i) for i in range(1, 23)] for a in ("p", "q")
        if c + a not in _ACROCENTRIC_P]


def _prep(seg: pd.DataFrame) -> pd.DataFrame:
    s = seg.copy()
    s["chr"] = s["chromosome"].astype(str)
    s["len"] = (s["end"] - s["start"]).clip(lower=0)
    s["mean"] = pd.to_numeric(s["segmentMean"], errors="coerce")
    return s.dropna(subset=["mean"])


def fraction_genome_altered(seg: pd.DataFrame, threshold: float = 0.2) -> pd.Series:
    """Per-sample fraction of the assayed genome with |log2 ratio| > ``threshold``.

    A standard genomic-instability / aneuploidy score in [0, 1]: total length of
    altered segments divided by total assayed length. ``threshold`` = 0.2 (log2)
    is the common cBioPortal/TCGA cutoff. Autosomes only.
    """
    s = _prep(seg)
    s = s[s["chr"].isin(_CENTROMERE) & (s["chr"] != "23")]
    s["alt_len"] = np.where(s["mean"].abs() > threshold, s["len"], 0.0)
    g = s.groupby("sampleId")[["len", "alt_len"]].sum()
    return (g["alt_len"] / g["len"].replace(0, np.nan)).rename("fga")


def arm_level_means(seg: pd.DataFrame) -> pd.DataFrame:
    """Sample x arm matrix of length-weighted mean log2 ratio (39 autosomal arms).

    Each segment is assigned to the chromosome arm containing its midpoint; within
    each (sample, arm) the log2 ratios are averaged weighted by segment length.
    A positive value = net gain of that arm, negative = net loss. Arms with no
    covered segment for a sample are NaN.
    """
    s = _prep(seg)
    s = s[s["chr"].isin(_CENTROMERE)]
    mid = (s["start"] + s["end"]) / 2
    cen = s["chr"].map(_CENTROMERE)
    s["arm"] = s["chr"] + np.where(mid < cen, "p", "q")
    s = s[s["arm"].isin(ARMS)]
    s["wm"] = s["len"] * s["mean"]
    g = s.groupby(["sampleId", "arm"]).agg(wm=("wm", "sum"), length=("len", "sum"))
    g["arm_mean"] = g["wm"] / g["length"].replace(0, np.nan)
    return g["arm_mean"].unstack("arm").reindex(columns=ARMS)
