import pandas as pd

from mirna_tcga.enrich import over_representation, parse_gmt

GMT = (
    "SET_A\tdesc\tG1\tG2\tG3\tG4\tG5\tG6\n"
    "SET_B\tdesc\tG7\tG8\tG9\tG10\tG11\tG12\n"
    "TOO_SMALL\tdesc\tG1\tG2\n"
)


def test_parse_gmt_reads_sets_and_upcases():
    sets = parse_gmt(GMT)
    assert sets["SET_A"] == {"G1", "G2", "G3", "G4", "G5", "G6"}
    assert "TOO_SMALL" in sets
    assert parse_gmt("g1\tg2\tg3")["g1"] == {"G3"}  # lower-case symbols upper-cased


def test_over_representation_flags_enriched_set():
    universe = [f"G{i}" for i in range(1, 21)]
    hits = ["G1", "G2", "G3", "G4", "G5"]  # all of SET_A
    gene_sets = parse_gmt(GMT)
    ora = over_representation(hits, universe, gene_sets, min_set_size=5)
    top = ora.iloc[0]
    assert top["gene_set"] == "SET_A"
    assert top["overlap"] == 5
    assert top["fold_enrichment"] > 1
    assert top["p"] < ora.set_index("gene_set").loc["SET_B", "p"]
    # small set filtered out
    assert "TOO_SMALL" not in set(ora["gene_set"])
    assert "q" in ora.columns


def test_over_representation_empty_when_no_sets_pass_filters():
    ora = over_representation(["G1"], ["G1", "G2"], parse_gmt(GMT), min_set_size=5)
    assert isinstance(ora, pd.DataFrame)
    assert ora.empty
