import pandas as pd
import pytest

from scripts.prepare_ispy2_series_manifest import build_series_manifest


def test_manifest_contains_only_reviewed_series_instance_uids():
    inventory = pd.DataFrame(
        {
            "patient_id": ["123456", "123456", "654321"],
            "visit_stage": ["T0", "T1", "T0"],
            "series_instance_uid": ["uid-a", "uid-b", "uid-c"],
        }
    )
    selection = pd.DataFrame(
        {
            "patient_id": ["123456", "123456", "654321"],
            "visit_stage": ["T0", "T1", "T0"],
            "series_instance_uid": ["uid-a", "uid-b", "uid-c"],
            "include": ["1", "0", "yes"],
        }
    )
    manifest, summary = build_series_manifest(inventory, selection)
    assert manifest["SeriesInstanceUID"].tolist() == ["uid-a", "uid-c"]
    assert summary["selected_patient_count"] == 2


def test_manifest_refuses_unmatched_uid_patient_visit_triple():
    inventory = pd.DataFrame(
        {"patient_id": ["123456"], "visit_stage": ["T0"], "series_instance_uid": ["uid-a"]}
    )
    selection = pd.DataFrame(
        {
            "patient_id": ["123456"],
            "visit_stage": ["T1"],
            "series_instance_uid": ["uid-a"],
            "include": ["1"],
        }
    )
    with pytest.raises(ValueError, match="do not match"):
        build_series_manifest(inventory, selection)


def test_manifest_refuses_duplicates_and_unknown_flags():
    inventory = pd.DataFrame(
        {"patient_id": ["123456"], "visit_stage": ["T0"], "series_instance_uid": ["uid-a"]}
    )
    selection = pd.DataFrame(
        {
            "patient_id": ["123456", "123456"],
            "visit_stage": ["T0", "T0"],
            "series_instance_uid": ["uid-a", "uid-a"],
            "include": ["1", "maybe"],
        }
    )
    with pytest.raises(ValueError, match="duplicate"):
        build_series_manifest(inventory, selection)
