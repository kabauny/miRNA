"""Gene / miRNA identifier conversion.

Replaces the original ``biomaRt`` scripts. Two backends:

* :func:`symbols_to_ensembl` / :func:`ensembl_to_symbols` via **mygene**
  (MyGene.info) -- fast and no local install of marts.
* :func:`mirbase_to_ensembl` via **pybiomart** -- the same Ensembl marts the
  original R code queried, for miRBase ids which mygene covers less well.

Both backends are optional extras (``pip install mirna-tcga[idmap]``); the
functions raise a clear error if the dependency is missing.
"""

from __future__ import annotations

from typing import Sequence

import pandas as pd


def _require(module: str):
    try:
        return __import__(module)
    except ImportError as exc:  # pragma: no cover - import-guard
        raise ImportError(
            f"'{module}' is required for ID conversion. Install with: "
            "pip install 'mirna-tcga[idmap]'"
        ) from exc


def symbols_to_ensembl(symbols: Sequence[str], species: str = "human") -> pd.DataFrame:
    """Map HUGO gene symbols to Ensembl gene ids via MyGene.info."""
    mygene = _require("mygene")
    mg = mygene.MyGeneInfo()
    res = mg.querymany(
        list(symbols), scopes="symbol", fields="ensembl.gene,symbol",
        species=species, as_dataframe=True, df_index=False,
    )
    return res


def ensembl_to_symbols(ensembl_ids: Sequence[str], species: str = "human") -> pd.DataFrame:
    """Map Ensembl gene ids to HUGO symbols via MyGene.info.

    Version suffixes (``ENSG00000141510.16``) are stripped before querying.
    """
    mygene = _require("mygene")
    cleaned = [e.split(".")[0] for e in ensembl_ids]
    mg = mygene.MyGeneInfo()
    return mg.querymany(
        cleaned, scopes="ensembl.gene", fields="symbol,name",
        species=species, as_dataframe=True, df_index=False,
    )


def mirbase_to_ensembl(mirbase_ids: Sequence[str]) -> pd.DataFrame:
    """Map miRBase ids to Ensembl gene info via the Ensembl BioMart.

    Mirrors the attributes requested by the original ``miRNA_conversion.R``.
    """
    pybiomart = _require("pybiomart")
    dataset = pybiomart.Dataset(name="hsapiens_gene_ensembl", host="http://www.ensembl.org")
    return dataset.query(
        attributes=[
            "mirbase_id", "ensembl_gene_id", "external_gene_name",
            "chromosome_name", "start_position", "end_position",
        ],
        filters={"mirbase_id": list(mirbase_ids)},
    )
