import numpy as np
import pytest

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state
from onco_dynamics_x.validation.diagnostics import relative_l2_error, timestep_refinement


def test_relative_l2_error_is_zero_for_identical_arrays():
    values = np.arange(9, dtype=float).reshape(3, 3)
    assert relative_l2_error(values, values.copy()) == 0.0


def test_relative_l2_error_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        relative_l2_error(np.zeros((2, 2)), np.zeros((3, 3)))


def test_timestep_refinement_returns_all_state_fields():
    params = ModelParameters()
    state = initial_state(9)
    dx = 2.0 / 8
    dt = min(0.0001, diffusion_timestep_limit(params, dx) * 0.25)
    errors = timestep_refinement(state, params, dx, dt, steps=2)
    assert set(errors) == set(state.FIELD_NAMES)
    assert all(np.isfinite(value) and value >= 0 for value in errors.values())
