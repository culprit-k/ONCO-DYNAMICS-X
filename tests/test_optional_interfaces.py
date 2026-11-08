import numpy as np
import pytest

from onco_dynamics_x.calibration.bayesian import PriorSpec, build_gaussian_parameter_model
from onco_dynamics_x.control.mean_field_ocp import optimize_hypothetical_input


def test_prior_spec_validates_distribution_and_scale():
    PriorSpec("normal", 0.0, 1.0).validate()
    with pytest.raises(ValueError):
        PriorSpec("uniform", 0.0, 1.0).validate()
    with pytest.raises(ValueError):
        PriorSpec("normal", 0.0, 0.0).validate()


def test_bayesian_builder_rejects_empty_prior_before_optional_import():
    with pytest.raises(ValueError, match="At least one"):
        build_gaussian_parameter_model(
            np.array([1.0]), lambda _params: 1.0, {}, observation_sd=0.1
        )


def test_control_interface_rejects_invalid_horizon_before_optional_import():
    with pytest.raises(ValueError, match="Horizon"):
        optimize_hypothetical_input(horizon=0.0)
