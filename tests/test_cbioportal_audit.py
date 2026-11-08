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


def test_generic_study_audit_supports_non_metabric_study_without_expression_list():
    class FakeClient:
        def get_json(self, path, params=None):
            assert path == "studies/example_study"
            return {"studyId": "example_study", "name": "Example"}

        def get_all(self, path, params=None, page_size=1000):
            if path.endswith("/patients"):
                return [{"patientId": "P1"}]
            if path.endswith("/samples"):
                return [{"sampleId": "S1", "patientId": "P1"}]
            if path.endswith("/clinical-attributes"):
                return [
                    {"clinicalAttributeId": "OS_MONTHS", "displayName": "OS months",
                     "datatype": "NUMBER", "patientAttribute": True},
                    {"clinicalAttributeId": "OS_STATUS", "displayName": "OS status",
                     "datatype": "STRING", "patientAttribute": True},
                ]
            if path.endswith("/clinical-data"):
                if params and params.get("clinicalDataType") == "PATIENT":
                    return [
                        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P1", "value": "12"},
                        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1",
                         "value": "0:LIVING"},
                    ]
                return []
            return []

    from onco_dynamics_x.data_engine.cbioportal import audit_cbioportal_study

    result = audit_cbioportal_study(FakeClient(), study_id="example_study")
    assert result["study"]["study_id"] == "example_study"
    assert result["cohort_counts"]["unique_patients"] == 1
    assert result["expression_survival_linkage"] is None


def test_endpoint_fields_are_not_misclassified_as_treatment_exposure():
    from onco_dynamics_x.data_engine.cbioportal import _candidate_treatment_field

    assert not _candidate_treatment_field({
        "clinicalAttributeId": "DFS_MONTHS",
        "displayName": "Disease Free (Months)",
        "description": "Disease free survival since initial treatment",
    })
    assert _candidate_treatment_field({
        "clinicalAttributeId": "HISTORY_NEOADJUVANT_TRTYN",
        "displayName": "Neoadjuvant Therapy Type Administered",
    })
    assert not _candidate_treatment_field({
        "clinicalAttributeId": "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
        "displayName": "New Neoplasm Event Post Initial Therapy Indicator",
    })


def test_survival_summary_includes_disease_specific_and_progression_endpoints():
    records = [
        {"clinicalAttributeId": "DSS_MONTHS", "patientId": "P1", "value": "20"},
        {"clinicalAttributeId": "DSS_STATUS", "patientId": "P1", "value": "0:CENSORED"},
        {"clinicalAttributeId": "DFS_MONTHS", "patientId": "P1", "value": "10"},
        {"clinicalAttributeId": "DFS_STATUS", "patientId": "P1", "value": "1:RECURRED"},
        {"clinicalAttributeId": "PFS_MONTHS", "patientId": "P1", "value": "8"},
        {"clinicalAttributeId": "PFS_STATUS", "patientId": "P1", "value": "1:PROGRESSED"},
    ]
    result = summarize_survival_endpoint_pairs(records)
    by_endpoint = {row["endpoint"]: row for row in result}
    assert by_endpoint["disease_specific_survival"]["patients_with_both_fields"] == 1
    assert by_endpoint["disease_free_survival"]["patients_with_both_fields"] == 1
    assert by_endpoint["progression_free_survival"]["patients_with_both_fields"] == 1
    assert by_endpoint["disease_free_survival"]["status_value_counts_complete_pairs"] == {
        "1:RECURRED": 1
    }


def test_patient_covariate_coverage_reports_complete_cases_without_ids():
    from onco_dynamics_x.data_engine.cbioportal import summarize_patient_covariate_coverage

    records = [
        {"clinicalAttributeId": "AGE", "patientId": "P1", "value": "55"},
        {"clinicalAttributeId": "AGE", "patientId": "P2", "value": "61"},
        {"clinicalAttributeId": "STAGE", "patientId": "P1", "value": "II"},
        {"clinicalAttributeId": "STAGE", "patientId": "P2", "value": ""},
        {"clinicalAttributeId": "HYPOXIA", "patientId": "P1", "value": "0.7"},
    ]
    result = summarize_patient_covariate_coverage(
        records, ("AGE", "STAGE", "HYPOXIA")
    )
    assert result["complete_case_patients_across_all_fields"] == 1
    assert result["fields"]["AGE"]["patients_with_value"] == 2
    assert "P1" not in str(result)


def test_expression_linkage_counts_positive_followup_and_event_categories():
    samples = [{"sampleId": "S1", "patientId": "P1"},
               {"sampleId": "S2", "patientId": "P2"}]
    clinical = [
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P1", "value": "12"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P1", "value": "1:DECEASED"},
        {"clinicalAttributeId": "OS_MONTHS", "patientId": "P2", "value": "0"},
        {"clinicalAttributeId": "OS_STATUS", "patientId": "P2", "value": "0:LIVING"},
    ]
    result = summarize_expression_survival_linkage(
        samples, ["S1", "S2"], clinical
    )
    assert result["expression_patients_with_complete_overall_survival"] == 2
    assert result["expression_patients_with_positive_overall_survival_time"] == 1
    assert result["expression_positive_time_overall_survival_status_counts"] == {
        "1:DECEASED": 1
    }
