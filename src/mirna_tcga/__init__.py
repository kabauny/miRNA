"""TCGA miRNA / mRNA / mutation analysis toolkit.

A small, modern Python rewrite of the original R pipeline. Data comes from the
cBioPortal REST API; the modelling keeps the original nearest-shrunken-centroid
(PAM) method via scikit-learn's ``NearestCentroid``.
"""

__version__ = "0.1.0"

from .config import Config, load_config

__all__ = ["Config", "load_config", "__version__"]
