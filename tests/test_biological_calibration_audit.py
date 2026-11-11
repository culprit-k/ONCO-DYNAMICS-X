from onco_dynamics_x.data_engine.biological_calibration_audit import (  # noqa: E501
    summarise_calibration_targets,
)


def test_calibration_target_audit_matches_names_and_never_emits_values():
    rows = [
        {
            "attribute_id": "KI67_PERCENT",
            "display_name": "Ki-67 proliferation index",
            "level": "SAMPLE",
            "data_type": "NUMBER",
            "unique_entities_with_value": 27,
            "records_with_value": 27,
        },
        {
            "attribute_id": "ER_STATUS",
            "display_name": "Estrogen receptor status",
            "level": "SAMPLE",
            "data_type": "STRING",
            "unique_entities_with_value": 42,
            "records_with_value": 42,
        },
        {
            "attribute_id": "NOISE",
            "display_name": "Unrelated field",
            "level": "PATIENT",
            "data_type": "STRING",
            "unique_entities_with_value": 3,
            "records_with_value": 3,
        },
    ]
    result = summarise_calibration_targets("study", "Study", rows)

    assert result["matched_target_categories"] == [
        "ki67_or_proliferation_index",
        "er_status",
    ]
    assert result["target_attribute_matches"]["ki67_or_proliferation_index"][0][
        "entities_with_value"
    ] == 27
    assert "patient_id" not in str(result).lower()
    assert "sample_id" not in str(result).lower()
    ki67_field = result["target_attribute_matches"][
        "ki67_or_proliferation_index"
    ][0]
    assert "value" not in ki67_field

def test_calibration_target_audit_does_not_misclassify_cancer_status_as_er():
    rows = [
        {
            "attribute_id": "PERSON_NEOPLASM_CANCER_STATUS",
            "display_name": "Person Neoplasm Cancer Status",
            "level": "PATIENT",
            "data_type": "STRING",
            "unique_entities_with_value": 961,
            "records_with_value": 961,
        },
        {
            "attribute_id": "ER_IHC",
            "display_name": "ER status measured by IHC",
            "level": "PATIENT",
            "data_type": "STRING",
            "unique_entities_with_value": 2426,
            "records_with_value": 2426,
        },
        {
            "attribute_id": "HER2_STATUS",
            "display_name": "HER2 Status",
            "level": "SAMPLE",
            "data_type": "STRING",
            "unique_entities_with_value": 1980,
            "records_with_value": 1980,
        },
    ]
    result = summarise_calibration_targets("study", "Study", rows)
    er_matches = result["target_attribute_matches"]["er_status"]
    her2_matches = result["target_attribute_matches"]["her2_status"]

    assert [item["attribute_id"] for item in er_matches] == ["ER_IHC"]
    assert [item["attribute_id"] for item in her2_matches] == ["HER2_STATUS"]
