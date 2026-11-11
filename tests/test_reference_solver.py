import numpy as np
import pytest

pytest.importorskip("scipy.integrate")

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state
from onco_dynamics_x.validation.reference_solver import (
    compare_explicit_euler_to_reference,
    simulate_adaptive_reference,
)


def test_adaptive_reference_returns_finite_state_and_solver_diagnostics():
    state = initial_state(7)
    parameters = ModelParameters()
    dx = 2.0 / 6
    reference, diagnostics = simulate_adaptive_reference(
        state, parameters, dx, final_time=2e-5
    )

    reference.validate()
    assert reference.shape == state.shape
    assert diagnostics["success"] is True
    assert diagnostics["method"] == "DOP853"
    assert diagnostics["function_evaluations"] > 0
    assert diagnostics["post_integration_clipping_applied"] is False


def test_explicit_euler_can_be_compared_to_adaptive_reference():
    state = initial_state(7)
    parameters = ModelParameters()
    dx = 2.0 / 6
    dt = min(1e-5, diffusion_timestep_limit(parameters, dx) * 0.1)
    report = compare_explicit_euler_to_reference(
        state, parameters, dx=dx, dt=dt, steps=2
    )

    assert report["comparison"] == "explicit_euler_vs_adaptive_method_of_lines"
    assert np.isfinite(report["mean_relative_l2_error"])
    assert np.isfinite(report["maximum_relative_l2_error"])
    assert set(report["relative_l2_error_by_field"]) == set(state.FIELD_NAMES)
    assert report["maximum_relative_l2_error"] < 0.01


def test_adaptive_reference_rejects_invalid_grid_spacing():
    with pytest.raises(ValueError, match="dx"):
        simulate_adaptive_reference(
            initial_state(7), ModelParameters(), dx=0.0, final_time=1e-4
        )
