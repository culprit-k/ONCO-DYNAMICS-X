import numpy as np
import pytest

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state
from onco_dynamics_x.validation.diagnostics import (
    grid_refinement,
    relative_l2_error,
    timestep_refinement,
)


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


@pytest.mark.parametrize("dimensions", [2, 3])
def test_grid_refinement_compares_nested_grids(dimensions):
    params = ModelParameters()
    coarse = initial_state(5, spatial_dimensions=dimensions)
    fine = initial_state(9, spatial_dimensions=dimensions)
    dx = 0.5
    dt_limit = diffusion_timestep_limit(params, dx, dimensions)
    dt = min(0.0001, dt_limit * 0.25)
    errors = grid_refinement(
        coarse, fine, params, dx, dt, coarse_steps=1, refinement_ratio=2
    )
    assert set(errors) == set(coarse.FIELD_NAMES)
    assert all(np.isfinite(value) and value >= 0 for value in errors.values())


def test_grid_refinement_rejects_non_nested_shapes():
    params = ModelParameters()
    with pytest.raises(ValueError, match="Fine grid shape"):
        grid_refinement(
            initial_state(5),
            initial_state(8),
            params,
            coarse_dx=0.5,
            coarse_dt=0.0001,
            coarse_steps=1,
        )
