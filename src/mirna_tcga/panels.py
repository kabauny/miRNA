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

# Curated NSCLC driver genes (recurrently mutated oncogenes + tumour suppressors
# from the TCGA LUAD/LUSC studies and OncoKB). Used to focus the mutation layer
# of the constellation model on biology rather than every recurrently-mutated
# gene (which added noise, not signal).
NSCLC_DRIVERS = [
    # oncogenes
    "KRAS", "EGFR", "BRAF", "PIK3CA", "MET", "ERBB2", "MAP2K1", "NRAS", "HRAS",
    "RIT1", "CTNNB1", "MYC", "ALK", "RET", "ROS1", "FGFR1", "FGFR3", "AKT1",
    # tumour suppressors / other drivers
    "TP53", "STK11", "KEAP1", "NF1", "RB1", "CDKN2A", "PTEN", "SMARCA4",
    "ARID1A", "ARID2", "RBM10", "SETD2", "KMT2D", "KMT2C", "NOTCH1", "FAT1",
    "NFE2L2", "ATM", "U2AF1", "MGA", "SMAD4", "PTPRD", "KDM6A", "TSC1",
]

# Tumour-suppressor genes whose homozygous (deep) deletion is recurrent and
# functionally meaningful in NSCLC -- the focus of the copy-number layer.
DELETION_TSGS = [
    "CDKN2A", "CDKN2B", "MTAP", "PTEN", "RB1", "STK11", "SMAD4", "NF1",
    "ATM", "PTPRD", "FAT1", "LRP1B", "WWOX", "MACROD2", "CSMD1", "PARK2",
]
