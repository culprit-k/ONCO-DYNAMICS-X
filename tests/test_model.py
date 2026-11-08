import numpy as np
import pytest
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.reactions import reaction_terms
from onco_dynamics_x.solver import diffusion_timestep_limit, simulate, step
from onco_dynamics_x.state import initial_state

def test_initial_state_has_eight_finite_fields():
    state = initial_state(15)
    state.validate()
    assert len(state.as_dict()) == 8
    assert state.shape == (15, 15)

def test_state_validation_rejects_shape_mismatch():
    state = initial_state(9)
    state.immune = np.ones((8, 8))
    with pytest.raises(ValueError):
        state.validate()

def test_reaction_terms_return_all_fields_and_finite_values():
    state = initial_state(15)
    reactions = reaction_terms(state, ModelParameters())
    assert set(reactions) == set(state.FIELD_NAMES)
    assert all(np.isfinite(field).all() for field in reactions.values())

def test_solver_preserves_finiteness_and_nonnegative_fields():
    state = initial_state(21)
    params = ModelParameters()
    dx = 2.0 / 20
    dt = min(0.0005, diffusion_timestep_limit(params, dx) * 0.5)
    updated = step(state, params, dx=dx, dt=dt)
    updated.validate()
    assert all(float(np.min(field)) >= 0.0 for field in updated.as_dict().values())

def test_diffusion_stability_guard_rejects_large_step():
    state = initial_state(21)
    params = ModelParameters()
    limit = diffusion_timestep_limit(params, dx=0.1)
    with pytest.raises(ValueError, match="diffusion limit"):
        step(state, params, dx=0.1, dt=limit * 2)

def test_simulation_is_reproducible_for_same_initial_state():
    params = ModelParameters()
    state = initial_state(11)
    dx = 2.0 / 10
    dt = min(0.0002, diffusion_timestep_limit(params, dx) * 0.5)
    times_a, states_a = simulate(state, params, dx, dt, steps=4, save_every=2)
    times_b, states_b = simulate(state, params, dx, dt, steps=4, save_every=2)
    assert np.array_equal(times_a, times_b)
    assert np.allclose(states_a[-1].cancer, states_b[-1].cancer)
