"""Make ``src/`` importable when running scripts without installing the package.

Import this first in each script:  ``import _bootstrap  # noqa``
(Prefer ``pip install -e .`` for real use; this is just for quick runs.)
"""

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
