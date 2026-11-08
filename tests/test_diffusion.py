import numpy as np
import pytest
from onco_dynamics_x.diffusion import laplacian_no_flux

def test_laplacian_of_constant_field_is_zero():
    field = np.ones((9, 9), dtype=float) * 3.5
    assert np.allclose(laplacian_no_flux(field, dx=0.1), 0.0)

def test_laplacian_rejects_invalid_spacing():
    with pytest.raises(ValueError):
        laplacian_no_flux(np.ones((5, 5)), dx=0.0)
