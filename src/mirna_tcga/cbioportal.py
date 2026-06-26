"""Thin, dependency-light client for the cBioPortal REST API (v3).

Only the handful of endpoints this project needs are implemented. The client
returns tidy :class:`pandas.DataFrame` objects and provides high-level helpers
to build a genes x samples expression matrix and a clinical table.

API reference: https://www.cbioportal.org/api/swagger-ui/index.html
"""

from __future__ import annotations

from typing import Any, Iterable, Sequence

import pandas as pd
import requests


def _chunked(seq: Sequence[Any], size: int) -> Iterable[Sequence[Any]]:
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


class CBioPortalClient:
    """Minimal cBioPortal API client.

    Parameters
    ----------
    base_url:
        API root, e.g. ``https://www.cbioportal.org/api``.
    timeout:
        Per-request timeout in seconds.
    page_size:
        Page size used for paginated GET endpoints.
    session:
        Optional pre-configured :class:`requests.Session` (useful for tests).
    """

    def __init__(
        self,
        base_url: str = "https://www.cbioportal.org/api",
        timeout: int = 60,
        page_size: int = 10000,
        session: requests.Session | None = None,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout
        self.page_size = page_size
        self.session = session or requests.Session()
        self.session.headers.update({"Accept": "application/json"})

    # -- low-level helpers -------------------------------------------------
    def _url(self, path: str) -> str:
        return f"{self.base_url}/{path.lstrip('/')}"

    def _get(self, path: str, params: dict[str, Any] | None = None) -> list[dict]:
        """GET with transparent pagination; returns the concatenated records."""
        base = dict(params or {})
        base.setdefault("pageSize", self.page_size)
        page_size = base["pageSize"]
        results: list[dict] = []
        page = 0
        while True:
            page_params = {**base, "pageNumber": page}
            resp = self.session.get(self._url(path), params=page_params, timeout=self.timeout)
            resp.raise_for_status()
            batch = resp.json()
            results.extend(batch)
            if len(batch) < page_size:
                break
            page += 1
        return results

    def _post(self, path: str, json: Any, params: dict[str, Any] | None = None) -> list[dict]:
        resp = self.session.post(
            self._url(path), json=json, params=params or {}, timeout=self.timeout
        )
        resp.raise_for_status()
        return resp.json()

    # -- discovery ---------------------------------------------------------
    def get_studies(self) -> pd.DataFrame:
        return pd.DataFrame(self._get("studies"))

    def get_molecular_profiles(self, study_id: str) -> pd.DataFrame:
        return pd.DataFrame(self._get(f"studies/{study_id}/molecular-profiles"))

    def get_sample_lists(self, study_id: str) -> pd.DataFrame:
        return pd.DataFrame(self._get(f"studies/{study_id}/sample-lists"))

    # -- gene id mapping ---------------------------------------------------
    def get_genes(self, hugo_symbols: Sequence[str]) -> pd.DataFrame:
        """Map HUGO symbols to Entrez gene ids (cBioPortal molecular data is
        keyed by Entrez id)."""
        records = self._post("genes/fetch", json=list(hugo_symbols), params={"geneIdType": "HUGO_GENE_SYMBOL"})
        return pd.DataFrame(records)

    def entrez_ids(self, hugo_symbols: Sequence[str]) -> list[int]:
        genes = self.get_genes(hugo_symbols)
        if genes.empty:
            return []
        return genes["entrezGeneId"].astype(int).tolist()

    # -- molecular data ----------------------------------------------------
    def fetch_molecular_data(
        self,
        molecular_profile_id: str,
        entrez_gene_ids: Sequence[int],
        sample_list_id: str | None = None,
        sample_ids: Sequence[str] | None = None,
        gene_chunk: int = 2000,
    ) -> pd.DataFrame:
        """Fetch molecular (e.g. expression) values as a tidy long DataFrame.

        Exactly one of ``sample_list_id`` or ``sample_ids`` must be given.
        Gene ids are chunked to keep request bodies reasonable.
        """
        if (sample_list_id is None) == (sample_ids is None):
            raise ValueError("Provide exactly one of sample_list_id or sample_ids.")

        records: list[dict] = []
        for chunk in _chunked(list(entrez_gene_ids), gene_chunk):
            body: dict[str, Any] = {"entrezGeneIds": list(chunk)}
            if sample_list_id is not None:
                body["sampleListId"] = sample_list_id
            else:
                body["sampleIds"] = list(sample_ids)
            records.extend(
                self._post(
                    f"molecular-profiles/{molecular_profile_id}/molecular-data/fetch",
                    json=body,
                    params={"projection": "SUMMARY"},
                )
            )
        return pd.DataFrame(records)

    def expression_matrix(
        self,
        molecular_profile_id: str,
        hugo_symbols: Sequence[str],
        sample_list_id: str,
        gene_chunk: int = 2000,
    ) -> pd.DataFrame:
        """Return a genes (HUGO) x samples expression matrix."""
        entrez = self.entrez_ids(hugo_symbols)
        long = self.fetch_molecular_data(
            molecular_profile_id, entrez, sample_list_id=sample_list_id, gene_chunk=gene_chunk
        )
        if long.empty:
            return pd.DataFrame()
        wide = long.pivot_table(
            index="hugoGeneSymbol", columns="sampleId", values="value", aggfunc="first"
        )
        wide.index.name = "gene"
        wide.columns.name = "sample"
        return wide

    # -- mutations ---------------------------------------------------------
    def fetch_mutations(
        self,
        molecular_profile_id: str,
        hugo_symbols: Sequence[str],
        sample_list_id: str,
    ) -> pd.DataFrame:
        entrez = self.entrez_ids(hugo_symbols)
        records = self._post(
            f"molecular-profiles/{molecular_profile_id}/mutations/fetch",
            json={"entrezGeneIds": entrez, "sampleListId": sample_list_id},
            params={"projection": "DETAILED"},
        )
        return pd.DataFrame(records)

    def mutated_samples(
        self,
        molecular_profile_id: str,
        hugo_symbol: str,
        sample_list_id: str,
    ) -> set[str]:
        """Set of sample ids carrying a (non-silent) mutation in ``hugo_symbol``."""
        muts = self.fetch_mutations(molecular_profile_id, [hugo_symbol], sample_list_id)
        if muts.empty:
            return set()
        if "mutationType" in muts.columns:
            muts = muts[muts["mutationType"].str.lower() != "silent"]
        return set(muts["sampleId"].unique())

    # -- clinical ----------------------------------------------------------
    def get_clinical_data(
        self,
        study_id: str,
        patient_level: bool = True,
    ) -> pd.DataFrame:
        """Return clinical attributes as a wide (sample/patient x attribute) table."""
        clinical_type = "PATIENT" if patient_level else "SAMPLE"
        long = pd.DataFrame(
            self._get(
                f"studies/{study_id}/clinical-data",
                params={"clinicalDataType": clinical_type},
            )
        )
        if long.empty:
            return long
        id_col = "patientId" if patient_level else "sampleId"
        wide = long.pivot_table(
            index=id_col, columns="clinicalAttributeId", values="value", aggfunc="first"
        )
        wide.columns.name = None
        return wide.reset_index()
