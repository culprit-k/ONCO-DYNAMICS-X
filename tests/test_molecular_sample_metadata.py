from scripts.run_pam50_molecular_validation import (
    _parse_sample_metadata,
)


def test_parse_cBioPortal_sample_objects_with_nested_patient_and_tissue_type():
    samples = [
        {
            "id": "TCGA-AA-0001-01",
            "type": "Primary Solid Tumor",
            "patient": {"id": "TCGA-AA-0001"},
        },
        {
            "sampleId": "METABRIC-001",
            "sampleType": "Metastatic",
            "patientId": "METABRIC-P001",
        },
        {
            "id": "TCGA-AA-0002-11",
            "type": "Solid Tissues Normal",
            "patient": {"patientId": "TCGA-AA-0002"},
        },
    ]

    result = _parse_sample_metadata(samples)

    assert result["TCGA-AA-0001-01"] == {
        "patientId": "TCGA-AA-0001",
        "type": "Primary Solid Tumor",
    }
    assert result["METABRIC-001"] == {
        "patientId": "METABRIC-P001",
        "type": "Metastatic",
    }
    assert result["TCGA-AA-0002-11"] == {
        "patientId": "TCGA-AA-0002",
        "type": "Solid Tissues Normal",
    }


def test_parse_cBioPortal_sample_objects_ignores_unlinked_records():
    samples = [
        {"type": "Primary Solid Tumor", "patient": {"id": "P1"}},
        {"id": "S2", "type": "Primary Solid Tumor", "patient": {}},
    ]

    result = _parse_sample_metadata(samples)

    assert "S2" in result
    assert result["S2"]["patientId"] == ""
    assert len(result) == 1
