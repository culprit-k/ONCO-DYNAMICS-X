from onco_dynamics_x.data_engine.cbioportal import (
    CBioPortalClient,
    _candidate_endpoint,
    summarize_clinical_coverage,
    summarize_expression_survival_linkage,
    summarize_survival_endpoint_pairs,
)


def test_client_uses_cbioportal_api_root():
    assert CBioPortalClient().base_url == "https://www.cbioportal.org/api"


def test_candidate_endpoint_identifies_overall_survival_fields():
    assert _candidate_endpoint({
        "clinicalAttributeId": "OS_MONTHS",
        "displayName": "Overall Survival (Months)",
    }) == "overall_survival_candidate"


def test_candidate_endpoint_does_not_infer_unrelated_fields():
    assert _candidate_endpoint({
        "clinicalAttributeId": "ER_STATUS",
        "displayName": "ER status",
    }) is None


def test_clinical_coverage_reports_aggregate_counts_only():
    attributes = [
        {"clinicalAttributeId": "OS_STATUS", "displayName": "Overall survival status",
         "datatype": "STRING", "patientAttribute": True},
        {"clinicalAttributeId": "SAMPLE_QUALITY", "displayName": "Sample quality",
         "datatype": "STRING", "patientAttribute": False},
    ]
    patient_data = [
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1", "value": "1:DECEASED"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P2", "value": "0:LIVING"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P3", "value": ""},
    ]
    sample_data = [
        {"clinicalAttributeId": "SAMPLE_QUALITY", "sampleId": "S1", "value": "GOOD"},
        {"clinicalAttributeId": "SAMPLE_QUALITY", "sampleId": "S2", "value": ""},
    ]
    result = summarize_clinical_coverage(attributes, patient_data, sample_data)
    assert len(result) == 2
    survival = next(row for row in result if row["attribute_id"] == "OS_STATUS")
    assert survival["unique_entities_with_value"] == 2
    assert survival["candidate_survival_endpoint"] == "overall_survival_candidate"
    assert "P1" not in str(result)
    assert "DECEASED" not in str(result)


def test_client_pagination_stops_on_short_page():
    client = CBioPortalClient(retries=0)
    calls = []

    def fake_get_json(path, params):
        calls.append(params["pageNumber"])
        return [{"id": "a"}, {"id": "b"}] if params["pageNumber"] == 0 else []

    client.get_json = fake_get_json
    assert client.get_all("fake", page_size=3) == [{"id": "a"}, {"id": "b"}]
    assert calls == [0]




def test_survival_pair_summary_checks_patient_level_overlap_and_time_quality():
    records = [
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P1", "value": "12"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1", "value": "1:DECEASED"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P2", "value": "0"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P2", "value": "0:LIVING"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P3", "value": "not-a-number"},
        {"clinicalAttributeId": "RFS_MONTHS", "patientId": "P1", "value": "9"},
        {"clinicalAttributeId": "RFS_STATUS", "patientId": "P1", "value": "1:Recurred"},
    ]
    result = summarize_survival_endpoint_pairs(records)
    overall = next(row for row in result if row["endpoint"] == "overall_survival")
    assert overall["patients_with_time"] == 3
    assert overall["patients_with_status"] == 2
    assert overall["patients_with_both_fields"] == 2
    assert overall["positive_time_count"] == 1
    assert overall["zero_or_negative_time_count"] == 1
    assert overall["unparseable_time_count"] == 0
    relapse = next(row for row in result if row["endpoint"] == "relapse_free_survival")
    assert relapse["patients_with_both_fields"] == 1


def test_expression_survival_linkage_uses_patient_sample_mapping():
    samples = [
        {"sampleId": "S1", "patientId": "P1"},
        {"sampleId": "S2", "patientId": "P2"},
        {"sampleId": "S3", "patientId": "P3"},
    ]
    expression_ids = ["S1", "S2", "UNKNOWN"]
    clinical = [
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P1", "value": "12"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1", "value": "1:DECEASED"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P2", "value": "24"},
        {"clinicalAttributeId": "RFS_MONTHS", "patientId": "P1", "value": "10"},
        {"clinicalAttributeId": "RFS_STATUS", "patientId": "P1", "value": "1:Recurred"},
    ]
    result = summarize_expression_survival_linkage(samples, expression_ids, clinical)
    assert result["expression_sample_list_sample_count"] == 3
    assert result["expression_samples_linked_to_patient"] == 2
    assert result["unique_expression_patients"] == 2
    assert result["expression_patients_with_complete_overall_survival"] == 1
    assert result["expression_patients_with_complete_relapse_free_survival"] == 1
