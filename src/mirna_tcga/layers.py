"""Assemble multi-omic feature matrices for the NSCLC cohort.

Shared data-loading used by the screen / model scripts, so expression, deep
deletions, and mutations are built one consistent way:

* :func:`stream_expression`  -- genes x samples log2 mRNA matrix (chunked).
* :func:`deletion_matrix`    -- samples x genes 0/1 deep-deletion (HOMDEL) flags.
* :func:`mutation_matrix`    -- samples x genes 0/1 non-silent mutation flags.

The two alteration matrices come back with a companion ``subtype`` Series
(LUAD/LUSC per sample) for stratified testing / covariate adjustment.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Non-coding / silent mutation classes to exclude when flagging "mutated".
SILENT_CLASSES = {
    "silent", "3'utr", "5'utr", "3'flank", "5'flank", "intron", "igr", "rna",
}

# Truncating / clearly loss-of-function mutation classes. Used when "inactivated"
# should mean gene knockout (comparable to a deep deletion), not any coding change.
TRUNCATING_CLASSES = {
    "nonsense_mutation", "frame_shift_del", "frame_shift_ins", "splice_site",
    "nonstop_mutation", "translation_start_site", "splice_region",
}


def protein_coding_map(client) -> dict[int, str]:
    """``{entrezGeneId: hugoGeneSymbol}`` for all protein-coding genes."""
    g = pd.DataFrame(client._get("genes", {"projection": "SUMMARY"}))
    g = g[g["type"] == "protein-coding"]
    return dict(zip(g["entrezGeneId"].astype(int), g["hugoGeneSymbol"]))


def stream_expression(client, cfg, entrez, id2sym, study_keys, chunk: int = 2000) -> pd.DataFrame:
    """Genes (HUGO) x samples log2 mRNA matrix over ``study_keys`` (streamed)."""
    mats = []
    for key in study_keys:
        parts = []
        for i in range(0, len(entrez), chunk):
            sub = entrez[i : i + chunk]
            long = client.fetch_molecular_data(
                cfg.mrna_profile(key), sub, sample_list_id=cfg.all_samples_list(key)
            )
            if long.empty:
                continue
            long["gene"] = long["entrezGeneId"].astype(int).map(id2sym)
            wide = long.pivot_table(index="gene", columns="sampleId", values="value", aggfunc="first")
            parts.append(wide.astype("float32"))
            del long, wide
        if parts:
            mats.append(np.log2(pd.concat(parts).clip(lower=0) + 1.0))
    return pd.concat(mats, axis=1) if mats else pd.DataFrame()


def _binary_matrix(events, samples, id2sym, gene_col="entrezGeneId"):
    mat = pd.DataFrame(index=pd.Index(samples), dtype="int8")
    if events is None or events.empty:
        return mat
    ev = events.assign(gene=events[gene_col].astype(int).map(id2sym)).dropna(subset=["gene"])
    if ev.empty:
        return mat
    wide = ev.assign(v=1).pivot_table(
        index="sampleId", columns="gene", values="v", aggfunc="max", fill_value=0
    )
    return wide.reindex(index=samples).fillna(0).astype("int8")


def deletion_matrix(client, cfg, id2sym, study_keys):
    """samples x genes 0/1 deep-deletion (HOMDEL) matrix + per-sample subtype."""
    flags, subtype = [], {}
    for key in study_keys:
        samples = client.sample_list_ids(cfg.cna_samples_list(key))
        ev = client.discrete_cna_events(cfg.cna_profile(key), cfg.cna_samples_list(key), "HOMDEL")
        flags.append(_binary_matrix(ev, samples, id2sym))
        subtype.update({s: key.upper() for s in samples})
    B = pd.concat(flags, axis=0).fillna(0).astype("int8") if flags else pd.DataFrame()
    return B, pd.Series(subtype)


def mutation_matrix(client, cfg, id2sym, study_keys, truncating_only: bool = False):
    """samples x genes 0/1 mutation matrix + per-sample subtype.

    By default flags any non-silent mutation; with ``truncating_only`` keeps only
    clearly loss-of-function classes (nonsense / frameshift / splice / nonstop),
    i.e. gene knockouts comparable to a deep deletion.
    """
    entrez = list(id2sym)
    flags, subtype = [], {}
    for key in study_keys:
        samples = client.sample_list_ids(cfg.sequenced_samples_list(key))
        muts = client.mutation_events(cfg.mutation_profile(key), entrez, cfg.sequenced_samples_list(key))
        if not muts.empty and "mutationType" in muts:
            mtype = muts["mutationType"].astype(str).str.lower()
            muts = muts[mtype.isin(TRUNCATING_CLASSES)] if truncating_only \
                else muts[~mtype.isin(SILENT_CLASSES)]
        flags.append(_binary_matrix(muts, samples, id2sym))
        subtype.update({s: key.upper() for s in samples})
    B = pd.concat(flags, axis=0).fillna(0).astype("int8") if flags else pd.DataFrame()
    return B, pd.Series(subtype)
