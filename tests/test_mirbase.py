from mirna_tcga.mirbase import load_mimat_map, name_or_accession, parse_mature_fa

MATURE_FA = (
    ">hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p\n"
    "UGAGGUAGUAGGUUGUAUAGUU\n"
    ">hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p\n"
    "UAGCUUAUCAGACUGAUGUUGA\n"
    ">cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p\n"
    "UGAGGUAGUAGGUUGUAUAGUU\n"
)


def test_parse_mature_fa_maps_accession_to_name():
    m = parse_mature_fa(MATURE_FA)
    assert m["MIMAT0000062"] == "hsa-let-7a-5p"
    assert m["MIMAT0000076"] == "hsa-miR-21-5p"
    assert m["MIMAT0000001"] == "cel-let-7-5p"  # species-agnostic by accession
    assert len(m) == 3


def test_parse_mature_fa_ignores_malformed_headers():
    assert parse_mature_fa(">no_accession here\nACGU\nplain line\n") == {}


def test_name_or_accession_falls_back():
    m = {"MIMAT0000076": "hsa-miR-21-5p"}
    assert name_or_accession("MIMAT0000076", m) == "hsa-miR-21-5p"
    assert name_or_accession("MIMAT9999999", m) == "MIMAT9999999"


def test_load_mimat_map_from_local_path(tmp_path):
    p = tmp_path / "mature.fa"
    p.write_text(MATURE_FA)
    m = load_mimat_map(str(p))
    assert m["MIMAT0000062"] == "hsa-let-7a-5p"
