from onco_dynamics_x.data_engine.survival_analysis import (
    build_survival_records,
    normalise_stage,
    parse_os_event,
)


def test_stage_normaliser_maps_common_tcga_and_metabric_values():
    assert normalise_stage("Stage IIA") == 2
    assert normalise_stage("Stage IIIC") == 3
    assert normalise_stage("4") == 4
    assert normalise_stage("Stage 1") == 1
    assert normalise_stage("Stage X") is None
    assert normalise_stage("Stage 0") is None


def test_os_event_parser_respects_cbioportal_status_prefix():
    assert parse_os_event("1:DECEASED") == 1
    assert parse_os_event("0:LIVING") == 0
    assert parse_os_event("unknown") is None


def test_survival_record_builder_excludes_nonpositive_time_and_missing_stage():
    patients = [{"patientId": "P1"}, {"patientId": "P2"}, {"patientId": "P3"}]
    samples = [
        {"sampleId": "S1", "patientId": "P1"},
        {"sampleId": "S2", "patientId": "P2"},
        {"sampleId": "S3", "patientId": "P3"},
    ]
    patient_clinical = [
        {"clinicalAttributeId": "AGE", "patientId": "P1", "value": "50"},
        {"clinicalAttributeId": "AGE", "patientId": "P2", "value": "60"},
        {"clinicalAttributeId": "AGE", "patientId": "P3", "value": "70"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P1", "value": "12"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1", "value": "1:DECEASED"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P2", "value": "0"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P2", "value": "0:LIVING"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P3", "value": "24"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P3", "value": "0:LIVING"},
        {"clinicalAttributeId": "STAGE", "patientId": "P1", "value": "Stage II"},
        {"clinicalAttributeId": "STAGE", "patientId": "P2", "value": "Stage I"},
    ]
    frame, summary = build_survival_records(
        patients=patients,
        samples=samples,
        expression_sample_ids=["S1", "S2", "S3"],
        patient_clinical=patient_clinical,
        sample_clinical=[],
        age_attribute="AGE",
        stage_attribute="STAGE",
        stage_level="PATIENT",
    )
    assert len(frame) == 1
    assert frame[0]["patient_id"] == "P1"
    assert summary["patients_retained_for_model"] == 1
    assert summary["exclusion_counts"]["nonpositive_os_time"] == 1
    assert summary["exclusion_counts"]["missing_or_unrecognized_stage"] == 1


def test_survival_record_builder_can_retain_patients_without_stage_for_age_only_sensitivity():
    patients = [{"patientId": "P1"}, {"patientId": "P2"}]
    samples = [
        {"sampleId": "S1", "patientId": "P1"},
        {"sampleId": "S2", "patientId": "P2"},
    ]
    clinical = [
        {"clinicalAttributeId": "AGE", "patientId": "P1", "value": "50"},
        {"clinicalAttributeId": "AGE", "patientId": "P2", "value": "60"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P1", "value": "12"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1", "value": "1:DECEASED"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P2", "value": "24"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P2", "value": "0:LIVING"},
        {"clinicalAttributeId": "STAGE", "patientId": "P1", "value": "Stage II"},
    ]
    rows, summary = build_survival_records(
        patients=patients,
        samples=samples,
        expression_sample_ids=["S1", "S2"],
        patient_clinical=clinical,
        sample_clinical=[],
        age_attribute="AGE",
        stage_attribute="STAGE",
        stage_level="PATIENT",
        require_stage=False,
    )
    assert len(rows) == 2
    assert rows[0]["stage_group"] == 2
    assert rows[1]["stage_group"] is None
    assert summary["patients_retained_for_model"] == 2
    assert summary["exclusion_counts"]["missing_or_unrecognized_stage"] == 1
