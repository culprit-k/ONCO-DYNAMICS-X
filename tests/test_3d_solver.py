import numpy as np
import pytest

from onco_dynamics_x.diffusion import laplacian_no_flux
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit, step
from onco_dynamics_x.state import initial_state


def test_laplacian_of_constant_3d_volume_is_zero():
    field = np.full((7, 7, 7), 2.75, dtype=float)
    assert np.allclose(laplacian_no_flux(field, dx=0.2), 0.0)


def test_laplacian_rejects_unsupported_dimension():
    with pytest.raises(ValueError):
        laplacian_no_flux(np.ones((3, 3, 3, 3)), dx=0.1)


def test_3d_initial_state_and_solver_step():
    state = initial_state(grid_size=9, spatial_dimensions=3)
    assert state.shape == (9, 9, 9)
    parameters = ModelParameters()
    dx = 2.0 / 8
    dt = min(0.0001, diffusion_timestep_limit(parameters, dx, 3) * 0.25)
    updated = step(state, parameters, dx=dx, dt=dt)
    updated.validate()
    assert updated.cancer.shape == (9, 9, 9)


def test_3d_diffusion_limit_is_stricter_than_2d():
    parameters = ModelParameters()
    assert diffusion_timestep_limit(parameters, 0.1, 3) < (
        diffusion_timestep_limit(parameters, 0.1, 2)
    )
