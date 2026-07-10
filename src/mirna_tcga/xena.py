"""UCSC Xena loader for TCGA miRNA (and other) expression matrices.

cBioPortal does not carry miRNA expression for the PanCancer Atlas studies, but
UCSC Xena does host **TCGA** miRNA mature-strand expression. This module
downloads Xena dataset matrices and returns them in the same genes x samples
shape as :meth:`mirna_tcga.cbioportal.CBioPortalClient.expression_matrix`, so
the downstream preprocessing / PAM code is unchanged.

Xena data model
---------------
* Genomic matrices (expression) are stored **features x samples**: the first
  column holds the gene/miRNA id, the header row holds sample ids.
* Phenotype/clinical matrices are stored **samples x variables**.

Datasets are served at ``https://<host>/download/<dataset>`` and may be gzipped;
this loader detects and decompresses transparently.
"""

from __future__ import annotations

import gzip
import io
from typing import Sequence

import pandas as pd
import requests

# Common Xena hubs (all serve TCGA data).
HOSTS = {
    "tcga": "tcga.xenahubs.net",   # original TCGA
    "gdc": "gdc.xenahubs.net",     # TCGA reprocessed by the GDC pipeline
    "pancanatlas": "pancanatlas.xenahubs.net",
}


class XenaClient:
    """Download and parse UCSC Xena dataset matrices.

    Parameters
    ----------
    host:
        A key in :data:`HOSTS` or a full hostname.
    timeout:
        Per-request timeout in seconds.
    session:
        Optional pre-configured :class:`requests.Session` (useful for tests).
    """

    def __init__(
        self,
        host: str = "tcga",
        timeout: int = 120,
        session: requests.Session | None = None,
    ) -> None:
        self.host = HOSTS.get(host, host)
        self.timeout = timeout
        self.session = session or requests.Session()

    # -- raw fetch ---------------------------------------------------------
    def _fetch_raw(self, dataset: str) -> str:
        """Download a dataset and return its decoded text (handles gzip).

        Xena stores matrices gzipped, served at ``<dataset>.gz``; the bare key
        does not exist (S3 answers ``AccessDenied``). We therefore request the
        ``.gz`` object first and fall back to the plain key for any dataset that
        happens to be stored uncompressed.
        """
        base = f"https://{self.host}/download/{dataset}"
        resp = self.session.get(base + ".gz", timeout=self.timeout)
        if resp.status_code != 200:
            resp = self.session.get(base, timeout=self.timeout)
        resp.raise_for_status()
        content = resp.content
        if content[:2] == b"\x1f\x8b":  # gzip magic number
            content = gzip.decompress(content)
        return content.decode("utf-8")

    # -- matrices ----------------------------------------------------------
    def expression_matrix(
        self, dataset: str, genes: Sequence[str] | None = None
    ) -> pd.DataFrame:
        """Return a genes (rows) x samples (cols) expression matrix.

        ``genes`` optionally restricts to a subset (missing ids are ignored).
        """
        text = self._fetch_raw(dataset)
        mat = pd.read_csv(io.StringIO(text), sep="\t", index_col=0)
        mat.index.name = "gene"
        mat.columns.name = "sample"
        if genes is not None:
            keep = [g for g in genes if g in mat.index]
            mat = mat.loc[keep]
        return mat

    def phenotype(self, dataset: str) -> pd.DataFrame:
        """Return a samples (rows) x variables (cols) phenotype/clinical table."""
        text = self._fetch_raw(dataset)
        pheno = pd.read_csv(io.StringIO(text), sep="\t", index_col=0)
        pheno.index.name = "sample"
        return pheno
