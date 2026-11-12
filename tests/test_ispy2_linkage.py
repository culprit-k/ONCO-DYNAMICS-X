from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).parents[1] / "scripts" / "validate_ispy2_linkage.py"
SPEC = importlib.util.spec_from_file_location("ispy2_linkage", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_one_to_one_join_reports_counts_without_key_values() -> None:
    feature = pd.DataFrame({"subject": ["001", "002", "003"], "x": [1, 2, 3]})
    clinical = pd.DataFrame({"patient": ["001", "002", "004"], "pcr": [0, 1, 0]})
    report, merged = MODULE.validate_linkage(
        feature, clinical, feature_key="subject", clinical_key="patient"
    )
    assert report["matched_feature_rows"] == 2
    assert report["unmatched_feature_rows"] == 1
    assert report["unmatched_clinical_rows"] == 1
    assert len(merged) == 2
    assert "001" not in str(report)


def test_many_to_many_is_prohibited() -> None:
    feature = pd.DataFrame({"subject": ["x"], "x": [1]})
    clinical = pd.DataFrame({"patient": ["x"], "pcr": [0]})
    with pytest.raises(ValueError, match="many-to-many is prohibited"):
        MODULE.validate_linkage(
            feature,
            clinical,
            feature_key="subject",
            clinical_key="patient",
            cardinality="many_to_many",
        )


def test_duplicate_clinical_keys_are_rejected() -> None:
    feature = pd.DataFrame({"subject": ["x"], "x": [1]})
    clinical = pd.DataFrame({"patient": ["x", "x"], "pcr": [0, 1]})
    with pytest.raises(ValueError, match="Clinical linkage key is not unique"):
        MODULE.validate_linkage(
            feature, clinical, feature_key="subject", clinical_key="patient"
        )
