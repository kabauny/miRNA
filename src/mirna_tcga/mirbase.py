"""miRBase mature-miRNA name lookup (MIMAT accession -> ``hsa-miR-*`` name).

UCSC Xena's ``miRNA_HiSeq_gene`` matrices are keyed by stable miRBase **mature
accessions** (``MIMAT...``), which are precise but unreadable. This module maps
them to human-readable names by parsing a miRBase ``mature.fa`` file, whose
headers look like::

    >hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p

Accessions are stable across miRBase releases, so a current ``mature.fa`` names
the older accessions used by the TCGA miRNA-seq pipeline. The file is loaded from
a URL (disk-cached) or a local path; callers should degrade gracefully to raw
accessions if it is unavailable.
"""

from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen


def parse_mature_fa(text: str) -> dict[str, str]:
    """Parse ``mature.fa`` header lines into ``{MIMAT accession: name}``."""
    mapping: dict[str, str] = {}
    for line in text.splitlines():
        if not line.startswith(">"):
            continue
        parts = line[1:].split()
        if len(parts) >= 2 and parts[1].startswith("MIMAT"):
            mapping[parts[1]] = parts[0]
    return mapping


def load_mimat_map(source: str, cache: str | Path | None = None) -> dict[str, str]:
    """Load a MIMAT->name map from a URL (optionally disk-cached) or local path."""
    if source.startswith("http"):
        if cache is not None:
            cache = Path(cache)
            if not cache.exists():
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_bytes(urlopen(source, timeout=60).read())  # noqa: S310
            text = cache.read_text()
        else:
            text = urlopen(source, timeout=60).read().decode("utf-8")  # noqa: S310
    else:
        text = Path(source).read_text()
    return parse_mature_fa(text)


def name_or_accession(accession: str, mapping: dict[str, str]) -> str:
    """``hsa-miR-x`` name if known, else the accession unchanged."""
    return mapping.get(accession, accession)
