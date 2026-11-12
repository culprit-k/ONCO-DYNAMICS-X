import numpy as np
import pytest

from onco_dynamics_x.imaging.observation import (
    ftv_from_analysis_mask,
    ftv_from_dce_maps,
    soft_ftv_from_dce_maps,
)


def test_official_mask_volume_uses_physical_voxel_geometry():
    mask = np.zeros((2, 3, 4), dtype=bool)
    mask.flat[:10] = True
    result = ftv_from_analysis_mask(mask, (2.0, 2.0, 2.5), reference_ftv_cc=0.1)
    assert result.voxel_count == 10
    assert result.voxel_volume_cc == pytest.approx(0.01)
    assert result.ftv_cc == pytest.approx(0.1)
    assert result.difference_from_reference_cc == pytest.approx(0.0)


def test_hard_dce_operator_applies_pe_and_ser_thresholds():
    pe = np.array([[[50.0, 70.0], [90.0, 85.0]]])
    ser = np.array([[[1.2, 0.4], [0.2, 0.95]]])
    eligible = np.ones_like(pe, dtype=bool)
    result, mask = ftv_from_dce_maps(
        pe, ser, eligible, (1.0, 1.0, 1.0),
        pe_threshold_percent=70.0, ser_minimum=0.0,
    )
    assert mask.sum() == 3
    assert result.ftv_cc == pytest.approx(0.003)
    assert result.pe_threshold_percent == 70.0


def test_dce_operator_rejects_misaligned_geometry():
    pe = np.ones((2, 2, 2))
    with pytest.raises(ValueError, match="aligned 3D"):
        ftv_from_dce_maps(
            pe, pe, np.ones((2, 2, 3), dtype=bool), (1, 1, 1),
            pe_threshold_percent=70, ser_minimum=0,
        )


def test_binary_mask_rejects_unknown_fractional_dicom_seg_encoding():
    mask = np.array([[[0, 0.25], [1, 0]]], dtype=float)
    with pytest.raises(ValueError, match="binary 0/1"):
        ftv_from_analysis_mask(mask, (1, 1, 1))


def test_soft_operator_is_finite_and_bounded_by_eligible_voxel_volume():
    pe = np.full((2, 2, 2), 70.0)
    ser = np.full((2, 2, 2), 1.0)
    eligible = np.ones_like(pe, dtype=bool)
    result = soft_ftv_from_dce_maps(
        pe, ser, eligible, (1, 1, 1),
        pe_threshold_percent=70, ser_minimum=0.9,
        pe_temperature_percent=2.0, ser_temperature=0.1,
    )
    assert np.isfinite(result.ftv_cc)
    assert 0 < result.ftv_cc <= 0.008


def test_soft_operator_requires_explicit_positive_smoothing_parameters():
    values = np.ones((2, 2, 2))
    with pytest.raises(ValueError, match="temperatures"):
        soft_ftv_from_dce_maps(
            values, values, values.astype(bool), (1, 1, 1),
            pe_threshold_percent=70, ser_minimum=0.9,
            pe_temperature_percent=0, ser_temperature=0.1,
        )
