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


def test_target_value_audit_reports_exact_expression_overlap_without_identifiers():
    from onco_dynamics_x.data_engine.biological_calibration_audit import (
        summarise_target_value_overlap,
    )

    attributes = [
        {
            "clinicalAttributeId": "ER_STATUS",
            "displayName": "ER Status",
            "patientAttribute": False,
        },
        {
            "clinicalAttributeId": "HER2_SNP6",
            "displayName": "HER2 status measured by SNP6",
            "patientAttribute": True,
        },
    ]
    patients = [
        {"patientId": "P1", "clinicalAttributeId": "HER2_SNP6", "value": "Positive"},
        {"patientId": "P2", "clinicalAttributeId": "HER2_SNP6", "value": "Negative"},
        {"patientId": "P3", "clinicalAttributeId": "HER2_SNP6", "value": "Positive"},
    ]
    samples = [
        {"sampleId": "S1", "patientId": "P1"},
        {"sampleId": "S2", "patientId": "P2"},
        {"sampleId": "S3", "patientId": "P3"},
    ]
    sample_records = [
        {"sampleId": "S1", "clinicalAttributeId": "ER_STATUS", "value": "Positive"},
        {"sampleId": "S2", "clinicalAttributeId": "ER_STATUS", "value": "Negative"},
        {"sampleId": "S3", "clinicalAttributeId": "ER_STATUS", "value": "Positive"},
    ]

    result = summarise_target_value_overlap(
        attributes,
        patients,
        sample_records,
        samples,
        {"S1", "S2"},
    )
    by_attribute = {item["attribute_id"]: item for item in result}

    assert by_attribute["ER_STATUS"]["expression_entities_with_value"] == 2
    assert by_attribute["ER_STATUS"]["value_summary"]["category_counts"] == {
        "Negative": 1,
        "Positive": 1,
    }
    assert by_attribute["HER2_SNP6"]["expression_entities_with_value"] == 2
    assert by_attribute["HER2_SNP6"]["expression_entities_without_unambiguous_value"] == 0
    assert "P1" not in str(result)
    assert "S1" not in str(result)


def test_target_value_audit_counts_ambiguous_entities_separately():
    from onco_dynamics_x.data_engine.biological_calibration_audit import (
        summarise_target_value_overlap,
    )

    attributes = [
        {
            "clinicalAttributeId": "ER_STATUS",
            "displayName": "ER Status",
            "patientAttribute": False,
        }
    ]
    records = [
        {"sampleId": "S1", "clinicalAttributeId": "ER_STATUS", "value": "Positive"},
        {"sampleId": "S1", "clinicalAttributeId": "ER_STATUS", "value": "Negative"},
    ]
    result = summarise_target_value_overlap(
        attributes, [], records, [{"sampleId": "S1", "patientId": "P1"}], {"S1"}
    )
    assert result[0]["expression_entities_with_value"] == 0
    assert result[0]["expression_entities_ambiguous"] == 1
    assert result[0]["expression_entities_without_unambiguous_value"] == 1
