"""Tests for the cytoband->arm parser (mygene-free)."""

import numpy as np

from mirna_tcga.idmap import cytoband_to_arm


def test_cytoband_to_arm_basic():
    assert cytoband_to_arm("17q12") == "17q"
    assert cytoband_to_arm("17p13.1") == "17p"
    assert cytoband_to_arm("8q24.21") == "8q"
    assert cytoband_to_arm("Xp22.1") == "Xp"


def test_cytoband_to_arm_spans_and_missing():
    assert cytoband_to_arm("17q11-q21") == "17q"   # span -> first arm
    assert cytoband_to_arm("-") is None
    assert cytoband_to_arm(np.nan) is None
    assert cytoband_to_arm("") is None
    assert cytoband_to_arm("not a band") is None
