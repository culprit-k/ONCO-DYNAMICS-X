from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).parents[1] / "scripts" / "analyze_ispy2_multifeature.py"
SPEC = importlib.util.spec_from_file_location("ispy2_multifeature", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def make_frame(n: int = 12) -> pd.DataFrame:
    frame: dict[str, object] = {
        MODULE.ID_COLUMN: [f"private-id-{i}" for i in range(n)],
    }
    visit_scale = [1.0, 0.75, 0.48, 0.32]
    for key, cols in MODULE.FEATURES.items():
        baseline = np.linspace(8.0, 20.0, n)
        if key == "sphericity":
            baseline = np.linspace(0.1, 0.35, n)
            visit_scale = [1.0, 0.99, 1.05, 1.18]
        elif key == "longest_diameter":
            baseline = np.linspace(2.0, 8.0, n)
            visit_scale = [1.0, 0.8, 0.6, 0.4]
        elif key == "background_parenchymal_enhancement":
            baseline = np.linspace(10.0, 40.0, n)
            visit_scale = [1.0, 0.9, 0.8, 0.7]
        else:
            visit_scale = [1.0, 0.75, 0.48, 0.32]
        variation = np.linspace(-0.04, 0.04, n)
        matrix = np.column_stack(
            [baseline * (scale + variation * (idx + 1)) for idx, scale in enumerate(visit_scale)]
        )
        for idx, col in enumerate(cols):
            frame[col] = matrix[:, idx]
        for idx, col in enumerate(MODULE.CHANGE_COLUMNS[key], start=1):
            frame[col] = (matrix[:, idx] / matrix[:, 0] - 1.0) * 100.0
    return pd.DataFrame(frame)


def test_analysis_is_aggregate_only_and_checks_all_percentage_fields() -> None:
    report = MODULE.analyze_frame(make_frame(), bootstrap_iterations=200, seed=7)
    assert report["n_participants"] == 12
    assert all(
        value < 1e-6
        for value in report["formula_max_absolute_error_percentage_points"].values()
    )
    rendered = str(report)
    assert "private-id-0" not in rendered
    assert "subject_id" not in rendered
    assert len(report["baseline_contrasts"]) == 12
    assert len(report["t3_change_associations_with_volume_field"]) == 3


def test_analysis_rejects_corrupted_change_column() -> None:
    frame = make_frame()
    frame.loc[0, "LD_pch_T0_T3"] += 1.0
    with pytest.raises(ValueError, match="fails the derived-percentage formula check"):
        MODULE.analyze_frame(frame, bootstrap_iterations=200, seed=7)


def test_analysis_rejects_duplicate_subject_rows() -> None:
    frame = make_frame()
    frame.loc[1, MODULE.ID_COLUMN] = frame.loc[0, MODULE.ID_COLUMN]
    with pytest.raises(ValueError, match="present and unique"):
        MODULE.analyze_frame(frame, bootstrap_iterations=200, seed=7)
