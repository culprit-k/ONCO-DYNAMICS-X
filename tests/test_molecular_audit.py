from onco_dynamics_x.data_engine.molecular_audit import (
    CohortProfileSpec,
    _summarise_gene_records,
    audit_breast_cancer_marker_panel,
    audit_molecular_profile,
)


class FakeMolecularClient:
    def __init__(self, cohorts, marker_panel):
        self.cohorts = cohorts
        self.marker_panel = marker_panel

    def get_json(self, path, params=None):
        for cohort in self.cohorts:
            if path == f"studies/{cohort.study_id}":
                return {
                    "studyId": cohort.study_id,
                    "referenceGenome": "hg19",
                }
            if path == f"molecular-profiles/{cohort.molecular_profile_id}":
                return {
                    "molecularProfileId": cohort.molecular_profile_id,
                    "studyId": cohort.study_id,
                    "name": "mRNA expression",
                    "description": "Synthetic test profile",
                    "datatype": "CONTINUOUS",
                    "molecularAlterationType": "MRNA_EXPRESSION",
                }
            if path == f"sample-lists/{cohort.sample_list_id}":
                return {
                    "sampleListId": cohort.sample_list_id,
                    "studyId": cohort.study_id,
                    "sampleIds": ["SAMPLE-1", "SAMPLE-2", "SAMPLE-3"],
                }
        for symbol, entrez_id, _ in self.marker_panel:
            if path == f"genes/{entrez_id}":
                return {"entrezGeneId": entrez_id, "hugoGeneSymbol": symbol}
        if "/molecular-data" in path:
            gene_id = params["entrezGeneId"]
            return [
                {
                    "entrezGeneId": gene_id,
                    "sampleId": "SAMPLE-1",
                    "patientId": "PATIENT-SECRET",
                    "value": "1.5",
                },
                {
                    "entrezGeneId": gene_id,
                    "sampleId": "SAMPLE-2",
                    "patientId": "PATIENT-SECRET-2",
                    "value": "2.5",
                },
            ]
        raise AssertionError(f"Unexpected request: {path}")


def test_gene_summary_reports_coverage_and_values_without_identifiers():
    summary = _summarise_gene_records(
        [
            {"sampleId": "PRIVATE-SAMPLE-1", "patientId": "P1", "value": "1.0"},
            {"sampleId": "PRIVATE-SAMPLE-2", "patientId": "P2", "value": "3.0"},
            {"sampleId": "PRIVATE-SAMPLE-3", "patientId": "P3", "value": None},
        ],
        expected_sample_count=4,
    )

    assert summary["records_returned"] == 3
    assert summary["unique_samples_with_record"] == 3
    assert summary["sample_coverage_fraction"] == 0.75
    assert summary["missing_value_records"] == 1
    assert summary["finite_value_count"] == 2
    assert summary["distribution"]["median"] == 2.0
    assert "PRIVATE-SAMPLE" not in str(summary)
    assert "P1" not in str(summary)


def test_gene_summary_counts_duplicates_nonfinite_and_non_numeric_values():
    summary = _summarise_gene_records(
        [
            {"sampleId": "S1", "value": "1.0"},
            {"sampleId": "S1", "value": "2.0"},
            {"sampleId": "S2", "value": "not-numeric"},
            {"sampleId": "S3", "value": "nan"},
        ],
        expected_sample_count=3,
    )

    assert summary["duplicate_sample_records"] == 1
    assert summary["non_numeric_or_nonfinite_value_records"] == 2
    assert summary["finite_value_count"] == 2


def test_marker_panel_audit_aggregates_two_cohorts_without_row_level_data():
    cohorts = (
        CohortProfileSpec("TCGA test", "study_tcga", "profile_tcga", "list_tcga"),
        CohortProfileSpec(
            "METABRIC test", "study_metabric", "profile_metabric", "list_metabric"
        ),
    )
    panel = (
        ("ESR1", 2099, "Hormone receptor signalling"),
        ("MKI67", 4288, "Proliferation"),
    )
    client = FakeMolecularClient(cohorts, panel)
    report = audit_breast_cancer_marker_panel(
        client=client, cohorts=cohorts, marker_panel=panel
    )

    assert len(report["cohorts"]) == 2
    assert report["cross_cohort_overlap"]["marker_count_present_in_all_profiles"] == 2
    assert all(
        cohort["markers_with_finite_values"] == 2
        for cohort in report["cohorts"]
    )
    assert report["interpretation"]["cross_platform_harmonization_performed"] is False
    assert report["data_governance"]["aggregate_metrics_only"] is True
    assert "PATIENT-SECRET" not in str(report)
    assert "SAMPLE-1" not in str(report)


def test_marker_audit_rejects_mismatched_gene_identity():
    cohort = CohortProfileSpec(
        "TCGA test", "study_tcga", "profile_tcga", "list_tcga"
    )

    class MismatchedGeneClient(FakeMolecularClient):
        def get_json(self, path, params=None):
            if path == "genes/2099":
                return {"entrezGeneId": 2099, "hugoGeneSymbol": "WRONG"}
            return super().get_json(path, params)

    client = MismatchedGeneClient((cohort,), (("ESR1", 2099, "test"),))
    import pytest

    with pytest.raises(ValueError, match="resolves to WRONG"):
        audit_molecular_profile(client, cohort, (("ESR1", 2099, "test"),))
