import numpy as np
import pandas as pd
import pytest

from scripts.develop_ispy2_pcr_response_models import (
    analyze,
    build_response_frame,
    reaction_diffusion_identifiability_gate,
)


def synthetic_inputs(n=100):
    rng = np.random.default_rng(17)
    ids = np.arange(1000, 1000 + n)
    baseline = rng.lognormal(mean=3.0, sigma=0.5, size=n)
    reduction_1 = rng.normal(-0.5, 0.35, size=n)
    reduction_2 = rng.normal(-1.1, 0.5, size=n)
    reduction_3 = rng.normal(-1.8, 0.7, size=n)
    features = pd.DataFrame(
        {
            "CLINICAL-TRIAL-SUBJECT-ID": ids,
            "VOLUME_TUM_BLU_V10": baseline,
            "VOLUME_TUM_BLU_V20": baseline * np.exp(reduction_1),
            "VOLUME_TUM_BLU_V30": baseline * np.exp(reduction_2),
            "VOLUME_TUM_BLU_V40": baseline * np.exp(reduction_3),
        }
    )
    clinical = pd.DataFrame(
        {
            "Patient_ID": ids,
            "Arm": rng.choice(["Control", "Arm A", "Arm B"], size=n),
            "HR": rng.integers(0, 2, size=n),
            "HER2": rng.integers(0, 2, size=n),
            "MP": rng.integers(0, 2, size=n),
            "pCR": np.tile([0, 1], n // 2),
            "Age_at_Screening": rng.normal(50, 10, size=n),
            "Race": ["Unknown"] * n,
            "menopausal_status": ["Unknown"] * n,
            "ethnicity": ["Unknown"] * n,
        }
    )
    return features, clinical


def test_response_frame_uses_log_ratios_not_elapsed_time():
    features, clinical = synthetic_inputs()
    frame = build_response_frame(features, clinical)
    expected = np.log(
        features.loc[0, "VOLUME_TUM_BLU_V20"]
        / features.loc[0, "VOLUME_TUM_BLU_V10"]
    )
    assert np.isclose(frame.loc[0, "r1"], expected)
    assert "time_days" not in frame.columns
    assert frame["Patient_ID"].is_unique


def test_response_frame_rejects_nonpositive_volume():
    features, clinical = synthetic_inputs()
    features.loc[0, "VOLUME_TUM_BLU_V20"] = 0
    with pytest.raises(ValueError, match="strictly positive"):
        build_response_frame(features, clinical)


def test_rd_gate_blocks_missing_spatial_and_timing_inputs():
    gate = reaction_diffusion_identifiability_gate({})
    assert gate["status"] == "blocked_nonidentifiable"
    assert "registered_spatial_tumour_fields" in gate["missing_inputs"]
    assert gate["diffusion_identifiable_from_whole_volume_only"] is False


def test_rd_gate_only_opens_when_all_required_inputs_are_declared():
    fields = {
        "registered_spatial_tumour_fields": True,
        "actual_acquisition_times": True,
        "treatment_exposure_timeline": True,
        "validated_mri_observation_operator": True,
        "spatial_microenvironment_fields": True,
    }
    assert (
        reaction_diffusion_identifiability_gate(fields)["status"]
        == "eligible_for_spatial_calibration"
    )


def test_analysis_is_aggregate_only_and_reports_zero_patient_overlap():
    features, clinical = synthetic_inputs()
    result = analyze(features, clinical, n_bootstrap=100)
    assert result["split"]["patient_overlap_n"] == 0
    assert result["split"]["all_linked_patients_assigned_once"] is True
    assert result["reaction_diffusion_comparator"]["status"] == "blocked_nonidentifiable"
    serialized = str(result)
    assert "1000" not in serialized
    assert "Patient_ID" not in result
