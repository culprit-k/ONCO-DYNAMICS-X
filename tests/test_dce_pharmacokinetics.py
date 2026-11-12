import numpy as np
import pytest

from onco_dynamics_x.imaging.dce_pharmacokinetics import (
    constrain_pk_parameter_logits,
    extended_tofts_concentration,
    sample_pe_at_times,
    ser_from_pe_maps,
    spoiled_gre_percent_enhancement,
)


def test_pk_logits_are_constrained_to_valid_parameter_ranges():
    logits = np.array([[0.0, 0.0, 0.0], [2.0, -1.0, 1.0]])
    params = constrain_pk_parameter_logits(logits)
    assert (params["ktrans_per_min"] > 0).all()
    assert ((params["ve_fraction"] > 0) & (params["ve_fraction"] < 1)).all()
    assert (params["vp_fraction"] >= 0).all()
    assert (params["vp_fraction"] + params["ve_fraction"] < 1).all()


def test_extended_tofts_plasma_term_and_extracellular_contribution():
    times = np.array([0.0, 60.0, 120.0])
    cp = np.array([1.0, 1.0, 0.0])
    ct = extended_tofts_concentration(
        times,
        cp,
        ktrans_per_min=0.2,
        ve_fraction=0.4,
        vp_fraction=0.1,
    )
    assert ct.shape == (3,)
    assert ct[0] == pytest.approx(0.1)
    assert ct[1] > ct[0]
    assert ct[2] > 0


def test_extended_tofts_supports_spatial_parameter_maps():
    times = np.array([0.0, 60.0])
    cp = np.array([1.0, 0.5])
    ktrans = np.full((2, 3, 4), 0.2)
    ve = np.full_like(ktrans, 0.4)
    vp = np.full_like(ktrans, 0.1)
    ct = extended_tofts_concentration(
        times, cp, ktrans_per_min=ktrans, ve_fraction=ve, vp_fraction=vp
    )
    assert ct.shape == (2, 2, 3, 4)
    assert np.isfinite(ct).all()


def test_extended_tofts_rejects_invalid_parameter_fractions():
    with pytest.raises(ValueError, match=r"vp \+ ve"):
        extended_tofts_concentration(
            np.array([0.0, 60.0]),
            np.array([1.0, 0.5]),
            ktrans_per_min=0.2,
            ve_fraction=0.8,
            vp_fraction=0.3,
        )


def test_spgr_enhancement_is_zero_without_contrast():
    concentration = np.zeros((3, 2, 2, 2))
    enhancement = spoiled_gre_percent_enhancement(
        concentration,
        baseline_t1_seconds=np.full((2, 2, 2), 1.2),
        flip_angle_degrees=15,
        repetition_time_seconds=0.005,
        longitudinal_relaxivity_s_inv_per_mM=4.0,
    )
    assert np.allclose(enhancement, 0.0)


def test_spgr_enhancement_increases_with_positive_concentration():
    concentration = np.array([0.0, 0.2, 0.5])[:, None]
    enhancement = spoiled_gre_percent_enhancement(
        concentration,
        baseline_t1_seconds=1.2,
        flip_angle_degrees=15,
        repetition_time_seconds=0.005,
        longitudinal_relaxivity_s_inv_per_mM=4.0,
    )
    assert enhancement[0, 0] == pytest.approx(0.0)
    assert enhancement[2, 0] > enhancement[1, 0] > 0


def test_sampling_pe_maps_uses_actual_requested_timepoints():
    times = np.array([0.0, 60.0, 120.0, 180.0])
    pe = np.array([0.0, 10.0, 20.0, 30.0])[:, None, None, None]
    sampled = sample_pe_at_times(times, pe, np.array([90.0, 150.0]))
    assert sampled[:, 0, 0, 0].tolist() == pytest.approx([15.0, 25.0])


def test_pe_sampling_refuses_extrapolation():
    with pytest.raises(ValueError, match="within the measured"):
        sample_pe_at_times(
            np.array([0.0, 60.0]),
            np.array([0.0, 10.0])[:, None],
            np.array([90.0]),
        )


def test_ser_ratio_returns_explicit_denominator_validity_mask():
    early = np.array([[100.0, 50.0]])
    late = np.array([[50.0, 0.0]])
    ser, valid = ser_from_pe_maps(
        early, late, minimum_absolute_late_enhancement_percent=1.0
    )
    assert valid.tolist() == [[True, False]]
    assert ser[0, 0] == pytest.approx(2.0)
    assert np.isnan(ser[0, 1])


def test_spgr_t2star_attenuation_requires_complete_parameter_set():
    with pytest.raises(ValueError, match="supplied together"):
        spoiled_gre_percent_enhancement(
            np.zeros((2, 1)),
            baseline_t1_seconds=1.2,
            flip_angle_degrees=15,
            repetition_time_seconds=0.005,
            longitudinal_relaxivity_s_inv_per_mM=4.0,
            echo_time_seconds=0.002,
        )


def test_spgr_t2star_attenuation_is_accepted_when_fully_specified():
    pe = spoiled_gre_percent_enhancement(
        np.array([0.0, 0.2])[:, None],
        baseline_t1_seconds=1.2,
        flip_angle_degrees=15,
        repetition_time_seconds=0.005,
        longitudinal_relaxivity_s_inv_per_mM=4.0,
        echo_time_seconds=0.002,
        baseline_r2star_s_inv=20.0,
        transverse_relaxivity_s_inv_per_mM=5.0,
    )
    assert pe[0, 0] == pytest.approx(0.0)
    assert np.isfinite(pe).all()
