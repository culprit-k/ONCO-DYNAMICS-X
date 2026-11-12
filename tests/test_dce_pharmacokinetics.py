import numpy as np
import pytest

from onco_dynamics_x.imaging.dce_pharmacokinetics import (
    constrain_pk_parameter_logits,
    extended_tofts_concentration,
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
    with pytest.raises(ValueError, match="vp + ve"):
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
