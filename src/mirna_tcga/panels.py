"""Small curated gene panels for quick, runnable demonstrations.

Using a panel keeps the example scripts fast. For real analyses, pass the full
coding-gene universe (e.g. all symbols from ``client.get_genes`` / a reference)
instead of a panel.
"""

# Lung adenocarcinoma (LUAD) vs squamous (LUSC) markers + common drivers.
LUNG_MARKER_PANEL = [
    "TP53", "KRAS", "EGFR", "STK11", "KEAP1", "NF1", "BRAF", "PIK3CA",
    "NAPSA", "NKX2-1", "SFTPC", "KRT5", "KRT6A", "TP63", "SOX2",
    "MET", "ALK", "RET", "ROS1", "CDKN2A", "PTEN", "RB1", "MYC",
]
