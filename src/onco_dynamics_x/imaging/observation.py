"""Measurement-side DCE-MRI observation operator for I-SPY 2 FTV.

The implemented operator maps PE/SER maps and validated eligibility masks to
functional tumour volume. The separate latent-PDE-to-PE/SER map is intentionally
not supplied: it must be calibrated against registered voxel-level observations.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping, Sequence

import numpy as np


def _validate_spacing(spacing_mm: Sequence[float]) -> tuple[float, float, float]:
    values = tuple(float(x) for x in spacing_mm)
    if len(values) != 3 or not np.isfinite(values).all() or any(x <= 0 for x in values):
        raise ValueError("spacing_mm must contain three finite positive voxel dimensions")
    return values  # type: ignore[return-value]


def _validate_binary_mask(mask: np.ndarray, name: str = "mask") -> np.ndarray:
    array = np.asarray(mask)
    if array.ndim != 3:
        raise ValueError(f"{name} must be a 3D array")
    if array.dtype == bool:
        return array
    if not np.isfinite(array).all() or not np.isin(np.unique(array), [0, 1]).all():
        raise ValueError(f"{name} must be boolean or binary 0/1; convert DICOM SEG explicitly")
    return array.astype(bool)


@dataclass(frozen=True)
class FTVMeasurement:
    """Observed/modelled FTV in cc, with physical and threshold provenance."""

    ftv_cc: float
    voxel_count: float
    voxel_volume_cc: float
    pe_threshold_percent: float | None
    ser_minimum: float | None
    method: str
    reference_ftv_cc: float | None = None
    difference_from_reference_cc: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)



def decode_ftv_bitmask(
    encoded_mask: np.ndarray,
    *,
    bit_masks: Mapping[str, int],
) -> dict[str, np.ndarray]:
    """Decode FTV mask layers from explicit, externally verified bit values.

    Required layer names are precontrast_pass, pe_pass, voi_pass and omit.
    Optional connectivity_pass is supported if the dataset-specific mask
    documentation defines such a bit. This function intentionally does not
    hard-code layer bit values.
    """
    required = {"precontrast_pass", "pe_pass", "voi_pass", "omit"}
    missing = required - set(bit_masks)
    if missing:
        raise ValueError(f"Missing verified bit definitions: {sorted(missing)}")
    data = np.asarray(encoded_mask)
    if data.ndim != 3 or not np.issubdtype(data.dtype, np.integer):
        raise ValueError("encoded_mask must be a 3D integer bitmask array")
    if np.any(data < 0):
        raise ValueError("encoded_mask cannot contain negative values")
    names = sorted(required)
    if "connectivity_pass" in bit_masks:
        names.append("connectivity_pass")
    values = [int(bit_masks[name]) for name in names]
    if any(value <= 0 or value & (value - 1) for value in values):
        raise ValueError("Each bit definition must be a positive power of two")
    if len(set(values)) != len(values):
        raise ValueError("Bit definitions must be distinct")
    layers = {name: (data & int(bit_masks[name])) != 0 for name in names}
    eligible = (
        layers["precontrast_pass"]
        & layers["voi_pass"]
        & ~layers["omit"]
    )
    if "connectivity_pass" in layers:
        eligible &= layers["connectivity_pass"]
    layers["eligibility_mask_without_pe_ser_threshold"] = eligible
    layers["pe_pass"] = layers["pe_pass"]
    return layers


def ftv_from_analysis_mask(
    analysis_mask: np.ndarray,
    voxel_spacing_mm: Sequence[float],
    *,
    reference_ftv_cc: float | None = None,
    comparison_tolerance_cc: float = 0.05,
) -> FTVMeasurement:
    """Integrate the official binary FTV analysis mask in physical voxel units.

    The input must be a decoded binary FTV-positive mask. Do not convert the
    raw bit-encoded DICOM analysis mask with >0: that would combine distinct
    filtering layers. Decode it first using bit assignments verified against
    the official Analysis mask files description.
    """
    mask = _validate_binary_mask(analysis_mask, "analysis_mask")
    spacing = _validate_spacing(voxel_spacing_mm)
    if comparison_tolerance_cc < 0:
        raise ValueError("comparison_tolerance_cc cannot be negative")
    voxel_volume_cc = float(np.prod(spacing) / 1000.0)
    voxel_count = int(mask.sum())
    ftv_cc = voxel_count * voxel_volume_cc
    difference = (
        ftv_cc - float(reference_ftv_cc) if reference_ftv_cc is not None else None
    )
    return FTVMeasurement(
        ftv_cc=float(ftv_cc),
        voxel_count=float(voxel_count),
        voxel_volume_cc=voxel_volume_cc,
        pe_threshold_percent=None,
        ser_minimum=None,
        method="binary official FTV analysis mask integrated using DICOM voxel geometry",
        reference_ftv_cc=reference_ftv_cc,
        difference_from_reference_cc=difference,
    )


def ftv_from_dce_maps(
    pe_map_percent: np.ndarray,
    ser_map: np.ndarray,
    eligibility_mask: np.ndarray,
    voxel_spacing_mm: Sequence[float],
    *,
    pe_threshold_percent: float,
    ser_minimum: float,
) -> tuple[FTVMeasurement, np.ndarray]:
    """Apply the published PE/SER thresholds inside a pre-validated mask.

    eligibility_mask must encode the correctly transformed VOI, exclusion
    regions, pre-contrast tissue/background mask and the documented spatial
    connectivity filter. This function does not recreate those manual/derived
    masks from raw images.
    """
    pe = np.asarray(pe_map_percent, dtype=float)
    ser = np.asarray(ser_map, dtype=float)
    eligible = _validate_binary_mask(eligibility_mask, "eligibility_mask")
    if pe.ndim != 3 or ser.ndim != 3 or pe.shape != ser.shape or pe.shape != eligible.shape:
        raise ValueError("PE, SER and eligibility arrays must be aligned 3D arrays")
    if not np.isfinite(pe).all() or not np.isfinite(ser).all():
        raise ValueError("PE and SER maps must be finite")
    if not np.isfinite(pe_threshold_percent) or not np.isfinite(ser_minimum):
        raise ValueError("PE and SER thresholds must be finite")
    mask = eligible & (pe >= float(pe_threshold_percent)) & (ser > float(ser_minimum))
    spacing = _validate_spacing(voxel_spacing_mm)
    voxel_volume_cc = float(np.prod(spacing) / 1000.0)
    count = int(mask.sum())
    measurement = FTVMeasurement(
        ftv_cc=float(count * voxel_volume_cc),
        voxel_count=float(count),
        voxel_volume_cc=voxel_volume_cc,
        pe_threshold_percent=float(pe_threshold_percent),
        ser_minimum=float(ser_minimum),
        method="hard DCE PE/SER threshold operator inside pre-validated eligibility mask",
    )
    return measurement, mask


def soft_ftv_from_dce_maps(
    pe_map_percent: np.ndarray,
    ser_map: np.ndarray,
    eligibility_mask: np.ndarray,
    voxel_spacing_mm: Sequence[float],
    *,
    pe_threshold_percent: float,
    ser_minimum: float,
    pe_temperature_percent: float,
    ser_temperature: float,
) -> FTVMeasurement:
    """Differentiable approximation to the hard PE/SER threshold operator.

    This is suitable as a future optimisation loss only after voxel-wise PE/SER
    predictions have been calibrated. Temperatures are explicit parameters and
    must be selected on training data, not silently assumed.
    """
    pe = np.asarray(pe_map_percent, dtype=float)
    ser = np.asarray(ser_map, dtype=float)
    eligible = _validate_binary_mask(eligibility_mask, "eligibility_mask")
    if pe.ndim != 3 or ser.ndim != 3 or pe.shape != ser.shape or pe.shape != eligible.shape:
        raise ValueError("PE, SER and eligibility arrays must be aligned 3D arrays")
    if not np.isfinite(pe).all() or not np.isfinite(ser).all():
        raise ValueError("PE and SER maps must be finite")
    if pe_temperature_percent <= 0 or ser_temperature <= 0:
        raise ValueError("Smoothing temperatures must be positive and explicitly specified")
    z_pe = np.clip((pe - pe_threshold_percent) / pe_temperature_percent, -60.0, 60.0)
    z_ser = np.clip((ser - ser_minimum) / ser_temperature, -60.0, 60.0)
    probability = (1.0 / (1.0 + np.exp(-z_pe))) * (1.0 / (1.0 + np.exp(-z_ser)))
    probability *= eligible
    voxel_volume_cc = float(np.prod(_validate_spacing(voxel_spacing_mm)) / 1000.0)
    expected_voxels = float(probability.sum())
    return FTVMeasurement(
        ftv_cc=expected_voxels * voxel_volume_cc,
        voxel_count=expected_voxels,
        voxel_volume_cc=voxel_volume_cc,
        pe_threshold_percent=float(pe_threshold_percent),
        ser_minimum=float(ser_minimum),
        method="soft sigmoid approximation to the DCE PE/SER threshold operator",
    )
