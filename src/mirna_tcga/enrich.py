"""Pathway over-representation analysis (ORA) for a gene hit list.

Given a set of "interesting" genes (e.g. those significantly associated with
survival) and the background universe they were selected from, test each pathway
/ gene set for enrichment with a one-sided hypergeometric (Fisher exact) test,
then FDR-correct across sets.

Gene sets are read from the MSigDB / Enrichr **GMT** format (one set per line:
``name <tab> description <tab> gene1 <tab> gene2 ...``). This keeps the analysis
self-contained -- no live enrichment web service is required.
"""

from __future__ import annotations

from typing import Iterable, Mapping

import pandas as pd
from scipy.stats import hypergeom

from .screen import benjamini_hochberg


def parse_gmt(text: str) -> dict[str, set[str]]:
    """Parse GMT text into ``{set_name: {genes}}`` (upper-cased symbols)."""
    sets: dict[str, set[str]] = {}
    for line in text.splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            continue
        name = parts[0].strip()
        genes = {g.strip().upper() for g in parts[2:] if g.strip()}
        if genes:
            sets[name] = genes
    return sets


def over_representation(
    hits: Iterable[str],
    universe: Iterable[str],
    gene_sets: Mapping[str, set[str]],
    min_set_size: int = 5,
    max_set_size: int = 500,
) -> pd.DataFrame:
    """One-sided hypergeometric enrichment of ``hits`` within each gene set.

    Every gene set is restricted to the tested ``universe`` before counting, so
    the background is the genes actually screened (not the whole genome). Sets
    with fewer than ``min_set_size`` or more than ``max_set_size`` universe genes
    are skipped.

    Returns a DataFrame (one row per tested set) with the overlap, expected
    count, fold enrichment, hypergeometric ``p``, BH ``q``, and the overlapping
    genes -- sorted by ascending ``p``.
    """
    universe = {g.upper() for g in universe}
    hits = {g.upper() for g in hits} & universe
    N = len(universe)
    n = len(hits)

    rows = []
    for name, genes in gene_sets.items():
        in_uni = genes & universe
        K = len(in_uni)
        if K < min_set_size or K > max_set_size:
            continue
        overlap = hits & in_uni
        k = len(overlap)
        # P(X >= k) with X ~ Hypergeometric(N, K, n)
        p = hypergeom.sf(k - 1, N, K, n) if k > 0 else 1.0
        expected = K * n / N if N else 0.0
        rows.append({
            "gene_set": name,
            "set_size": K,
            "overlap": k,
            "expected": round(expected, 2),
            "fold_enrichment": round(k / expected, 2) if expected > 0 else 0.0,
            "p": p,
            "genes": ",".join(sorted(overlap)),
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["q"] = benjamini_hochberg(df["p"].to_numpy())
    cols = ["gene_set", "set_size", "overlap", "expected", "fold_enrichment", "p", "q", "genes"]
    return df[cols].sort_values("p").reset_index(drop=True)
