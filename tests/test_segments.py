"""Offline tests for segment-level copy-number summaries."""

import numpy as np
import pandas as pd

from mirna_tcga.segments import ARMS, arm_level_means, fraction_genome_altered

# One sample: chr1p neutral (100 Mb), chr1q gain +0.5 (100 Mb), chr8p loss -1.0
# (40 Mb), chr8q neutral (100 Mb). Midpoints: 1p=50M<125M, 1q=150M>125M,
# 8p=20M<45.6M, 8q=96M>45.6M.
SEG = pd.DataFrame([
    {"sampleId": "S1", "chromosome": "1", "start": 0, "end": 100_000_000,
     "numberOfProbes": 100, "segmentMean": 0.0},
    {"sampleId": "S1", "chromosome": "1", "start": 100_000_000, "end": 200_000_000,
     "numberOfProbes": 100, "segmentMean": 0.5},
    {"sampleId": "S1", "chromosome": "8", "start": 0, "end": 40_000_000,
     "numberOfProbes": 50, "segmentMean": -1.0},
    {"sampleId": "S1", "chromosome": "8", "start": 46_000_000, "end": 146_000_000,
     "numberOfProbes": 100, "segmentMean": 0.0},
])


def test_fraction_genome_altered():
    fga = fraction_genome_altered(SEG, threshold=0.2)
    # altered = 1q gain (100 Mb) + 8p loss (40 Mb) = 140 Mb of 340 Mb total.
    assert abs(fga.loc["S1"] - 140 / 340) < 1e-6


def test_fraction_genome_altered_threshold():
    # With a threshold above 0.5, only the 8p loss (|-1.0|) counts: 40/340.
    fga = fraction_genome_altered(SEG, threshold=0.6)
    assert abs(fga.loc["S1"] - 40 / 340) < 1e-6


def test_arm_level_means_assigns_and_weights():
    arm = arm_level_means(SEG)
    assert list(arm.columns) == ARMS and len(ARMS) == 39
    assert abs(arm.loc["S1", "1p"] - 0.0) < 1e-9      # neutral p-arm
    assert abs(arm.loc["S1", "1q"] - 0.5) < 1e-9      # gained q-arm
    assert abs(arm.loc["S1", "8p"] - (-1.0)) < 1e-9   # lost 8p
    assert abs(arm.loc["S1", "8q"] - 0.0) < 1e-9      # neutral 8q
    assert np.isnan(arm.loc["S1", "2p"])              # unassayed arm -> NaN
