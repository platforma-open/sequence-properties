"""Put this regression package dir on sys.path so the byte-compare test can do
`from golden_cases import CASES` (a sibling module), matching regen_golden.py's
import. pytest auto-loads conftest before collecting sibling test modules.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
