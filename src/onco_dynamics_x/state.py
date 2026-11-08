"""State variables and reproducible synthetic initial conditions."""
from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]


@dataclass
class TumourState:
    """Eight dimensionless spatial fields on a shared 2D or 3D grid."""

    cancer: FloatArray
    immune: FloatArray
    fibroblast: FloatArray
    oxygen: FloatArray
    nutrient: FloatArray
    drug: FloatArray
    evasion: FloatArray
    vasculature: FloatArray

    FIELD_NAMES: ClassVar[tuple[str, ...]] = (
        "cancer",
        "immune",
        "fibroblast",
        "oxygen",
        "nutrient",
        "drug",
        "evasion",
        "vasculature",
    )

    def as_dict(self) -> dict[str, FloatArray]:
        return {name: getattr(self, name) for name in self.FIELD_NAMES}

    def copy(self) -> "TumourState":
        return TumourState(**{name: value.copy() for name, value in self.as_dict().items()})

    @property
    def shape(self) -> tuple[int, ...]:
        return self.cancer.shape

    def validate(self) -> None:
        arrays = self.as_dict()
        if len({array.shape for array in arrays.values()}) != 1:
            raise ValueError("All state fields must have the same shape.")
        if self.cancer.ndim not in (2, 3):
            raise ValueError("The model supports two-dimensional or three-dimensional grids.")
        for name, array in arrays.items():
            if not np.isfinite(array).all():
                raise ValueError(f"State field '{name}' contains NaN or infinite values.")


def initial_state(grid_size: int = 51, spatial_dimensions: int = 2) -> TumourState:
    """Create a deterministic synthetic 2D slice or 3D volume.

    This is a numerical demonstration fixture, not a reconstruction of a real tumour.
    """
    if grid_size < 5:
        raise ValueError("grid_size must be at least 5.")
    if spatial_dimensions not in (2, 3):
        raise ValueError("spatial_dimensions must be either 2 or 3.")

    axis = np.linspace(-1.0, 1.0, grid_size)
    if spatial_dimensions == 2:
        xx, yy = np.meshgrid(axis, axis, indexing="xy")
        radius_squared = xx**2 + yy**2
    else:
        xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
        radius_squared = xx**2 + yy**2 + zz**2

    core = np.exp(-radius_squared / 0.055)
    radius = np.sqrt(radius_squared)
    rim = np.exp(-((radius - 0.32) ** 2) / 0.025)
    state = TumourState(
        cancer=0.04 + 0.48 * core,
        immune=0.06 + 0.12 * rim,
        fibroblast=0.08 + 0.22 * rim,
        oxygen=0.98 - 0.30 * core,
        nutrient=0.98 - 0.18 * core,
        drug=np.zeros_like(core),
        evasion=0.04 + 0.20 * core,
        vasculature=1.0 + 0.08 * rim,
    )
    state.validate()
    return state
