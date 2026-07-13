import numpy as np
import pandas as pd

from mirna_tcga.associate import cmh_depletion_screen, cmh_two_sided_screen


def test_cmh_flags_gene_spared_in_index_group():
    # SPARED: deleted in ~30% of reference, never in the index group -> depleted
    # RANDOM: deleted equally in both -> not depleted
    rng = np.random.default_rng(0)
    n_ref, n_idx = 200, 40
    group = pd.Series([0] * n_ref + [1] * n_idx)
    strata = pd.Series(["A"] * n_ref + ["A"] * n_idx)
    spared = np.array([1 if rng.random() < 0.3 else 0 for _ in range(n_ref)] + [0] * n_idx)
    rnd = np.array([1 if rng.random() < 0.3 else 0 for _ in range(n_ref + n_idx)])
    B = pd.DataFrame({"SPARED": spared, "RANDOM": rnd})
    res = cmh_depletion_screen(B, group, strata, min_ref_freq=0.05)
    assert res.loc["SPARED", "n_idx"] == 0
    assert res.loc["SPARED", "or_mh"] < 1
    assert res.loc["SPARED", "p"] < 0.05
    assert res.loc["SPARED", "p"] < res.loc["RANDOM", "p"]


def test_cmh_adjusts_for_stratum_confounding():
    # Gene deleted only in stratum B; index group is mostly stratum A.
    # Pooled it looks "spared in index", but that is pure subtype confounding;
    # stratified CMH should NOT call it depleted.
    idx = (["A"] * 30 + ["B"] * 5)      # index: mostly A
    ref = (["A"] * 30 + ["B"] * 60)     # reference: lots of B
    strata = pd.Series(idx + ref)
    group = pd.Series([1] * len(idx) + [0] * len(ref))
    # deletion happens only in stratum B, at the same rate regardless of group
    dele = [0] * 30 + [1, 1, 0, 1, 1] + [0] * 30 + [1 if i % 2 else 0 for i in range(60)]
    B = pd.DataFrame({"BONLY": dele})
    res = cmh_depletion_screen(B, group, strata, min_ref_freq=0.05)
    # within stratum B the deletion rate is similar in both groups -> not depleted
    assert res.loc["BONLY", "p"] > 0.05


def test_cmh_two_sided_catches_enrichment_and_depletion():
    # ENRICHED: altered in ~40% of the index group but ~5% of reference.
    # DEPLETED: altered in ~30% of reference but never in the index.
    # RANDOM: equal in both.
    rng = np.random.default_rng(1)
    n_ref, n_idx = 200, 120
    group = pd.Series([0] * n_ref + [1] * n_idx)
    strata = pd.Series(["A"] * (n_ref + n_idx))
    enriched = np.array([1 if rng.random() < 0.05 else 0 for _ in range(n_ref)]
                        + [1 if rng.random() < 0.40 else 0 for _ in range(n_idx)])
    depleted = np.array([1 if rng.random() < 0.30 else 0 for _ in range(n_ref)] + [0] * n_idx)
    rnd = np.array([1 if rng.random() < 0.20 else 0 for _ in range(n_ref + n_idx)])
    B = pd.DataFrame({"ENRICHED": enriched, "DEPLETED": depleted, "RANDOM": rnd})
    res = cmh_two_sided_screen(B, group, strata, min_freq=0.05)
    assert res.loc["ENRICHED", "direction"] == "enriched"
    assert res.loc["ENRICHED", "or_mh"] > 1
    assert res.loc["ENRICHED", "p"] < 0.05
    assert res.loc["DEPLETED", "direction"] == "depleted"
    assert res.loc["DEPLETED", "or_mh"] < 1
    assert res.loc["DEPLETED", "p"] < 0.05
    assert res.loc["RANDOM", "p"] > 0.05


def test_cmh_two_sided_keeps_index_only_alteration():
    # Amplification common in the index group but rare in reference must survive
    # the min_freq filter (kept on EITHER group's frequency).
    group = pd.Series([0] * 100 + [1] * 60)
    strata = pd.Series(["A"] * 160)
    amp = np.array([0] * 100 + [1] * 25 + [0] * 35)   # 0% ref, ~42% index
    B = pd.DataFrame({"AMP": amp})
    res = cmh_two_sided_screen(B, group, strata, min_freq=0.05)
    assert "AMP" in res.index                # not dropped despite freq_ref = 0
    assert res.loc["AMP", "direction"] == "enriched"


def test_cmh_requires_deletable_reference():
    group = pd.Series([0] * 50 + [1] * 20)
    strata = pd.Series(["A"] * 70)
    rare = np.zeros(70, dtype=int)
    rare[0] = 1  # 1/50 in reference -> below min_ref_freq
    B = pd.DataFrame({"RARE": rare})
    res = cmh_depletion_screen(B, group, strata, min_ref_freq=0.05)
    assert "RARE" not in res.index  # filtered: not deletable enough to be informative
