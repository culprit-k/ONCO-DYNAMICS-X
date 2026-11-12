"""Pytest path setup for repository-local workflow scripts."""
from __future__ import annotations

import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
root_string = str(REPOSITORY_ROOT)
if root_string not in sys.path:
    sys.path.insert(0, root_string)
