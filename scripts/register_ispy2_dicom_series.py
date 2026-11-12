"""Rigidly register one curated I-SPY 2 DICOM series to a T0 reference.

The user must curate a single 3D anatomical reference series per visit from the
metadata inventory. This script does not infer DCE phases or decode DICOM SEG.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np


def _read_series(directory: Path, series_uid: str):
    try:
        import SimpleITK as sitk
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("Install the optional imaging extra: pip install -e '.[imaging]'") from exc
    if not directory.is_dir():
        raise ValueError(f"DICOM series directory does not exist: {directory}")
    available = list(sitk.ImageSeriesReader.GetGDCMSeriesIDs(str(directory)) or [])
    if series_uid not in available:
        raise ValueError(f"Series UID {series_uid} not found in {directory}; found {len(available)} series")
    files = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(str(directory), series_uid)
    if not files:
        raise ValueError(f"No DICOM instances returned for Series UID {series_uid}")
    reader = sitk.ImageSeriesReader()
    reader.SetFileNames(files)
    image = reader.Execute()
    if image.GetDimension() != 3:
        raise ValueError(
            f"Expected one curated 3D volume, got dimension={image.GetDimension()}. "
            "Select a single-volume anatomical series, not a multi-phase DCE series."
        )
    if any(size < 2 for size in image.GetSize()):
        raise ValueError(f"Image dimensions are too small: {image.GetSize()}")
    return sitk, image, len(files)


def _ncc(left: np.ndarray, right: np.ndarray) -> float | None:
    a = np.asarray(left, dtype=float)[::2, ::2, ::2].ravel()
    b = np.asarray(right, dtype=float)[::2, ::2, ::2].ravel()
    valid = np.isfinite(a) & np.isfinite(b) & ((a != 0) | (b != 0))
    if valid.sum() < 10:
        return None
    a, b = a[valid], b[valid]
    a -= a.mean()
    b -= b.mean()
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denominator == 0:
        return None
    return float(np.dot(a, b) / denominator)


def register_series(
    fixed_directory: Path,
    fixed_series_uid: str,
    moving_directory: Path,
    moving_series_uid: str,
    output_image: Path,
    output_transform: Path,
    qc_json: Path,
    *,
    moving_mask_path: Path | None = None,
    fixed_mask_path: Path | None = None,
    seed: int = 20261112,
    iterations: int = 300,
) -> dict[str, Any]:
    import SimpleITK as sitk

    sitk, fixed, fixed_file_count = _read_series(fixed_directory, fixed_series_uid)
    _, moving, moving_file_count = _read_series(moving_directory, moving_series_uid)
    if fixed.GetDimension() != moving.GetDimension():
        raise ValueError("Fixed and moving images must have the same dimension")
    if fixed.GetDimension() != 3:
        raise ValueError("Only 3D registration is supported")
    for name, image in (("fixed", fixed), ("moving", moving)):
        if any(not np.isfinite(x) or x <= 0 for x in image.GetSpacing()):
            raise ValueError(f"{name} image has invalid voxel spacing")

    # Rigid Euler transform is a conservative first registration stage.
    initial = sitk.CenteredTransformInitializer(
        fixed,
        moving,
        sitk.Euler3DTransform(),
        sitk.CenteredTransformInitializerFilter.GEOMETRY,
    )
    registration = sitk.ImageRegistrationMethod()
    registration.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration.SetMetricSamplingStrategy(registration.RANDOM)
    registration.SetMetricSamplingPercentage(0.20, seed)
    registration.SetInterpolator(sitk.sitkLinear)
    registration.SetOptimizerAsGradientDescent(
        learningRate=1.0,
        numberOfIterations=int(iterations),
        convergenceMinimumValue=1e-6,
        convergenceWindowSize=10,
    )
    registration.SetOptimizerScalesFromPhysicalShift()
    registration.SetShrinkFactorsPerLevel([4, 2, 1])
    registration.SetSmoothingSigmasPerLevel([2, 1, 0])
    registration.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
    registration.SetInitialTransform(initial, inPlace=False)
    final_transform = registration.Execute(fixed, moving)

    def resample(image, transform, interpolator, pixel_type=None):
        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(fixed)
        resampler.SetTransform(transform)
        resampler.SetInterpolator(interpolator)
        resampler.SetDefaultPixelValue(0)
        if pixel_type is not None:
            resampler.SetOutputPixelType(pixel_type)
        return resampler.Execute(image)

    identity = sitk.Transform(3, sitk.sitkIdentity)
    moving_before = resample(moving, identity, sitk.sitkLinear, sitk.sitkFloat32)
    moving_after = resample(moving, final_transform, sitk.sitkLinear, sitk.sitkFloat32)
    fixed_array = sitk.GetArrayFromImage(sitk.Cast(fixed, sitk.sitkFloat32))
    before_array = sitk.GetArrayFromImage(moving_before)
    after_array = sitk.GetArrayFromImage(moving_after)

    output_image.parent.mkdir(parents=True, exist_ok=True)
    output_transform.parent.mkdir(parents=True, exist_ok=True)
    qc_json.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(moving_after, str(output_image))
    sitk.WriteTransform(final_transform, str(output_transform))

    moving_mask_out = None
    dice = None
    if moving_mask_path is not None:
        if not moving_mask_path.exists():
            raise ValueError(f"Moving mask not found: {moving_mask_path}")
        moving_mask = sitk.ReadImage(str(moving_mask_path))
        if moving_mask.GetDimension() != 3:
            raise ValueError("Moving mask must be a 3D image")
        moved_mask = resample(
            moving_mask, final_transform, sitk.sitkNearestNeighbor, sitk.sitkUInt8
        )
        moving_mask_out = output_image.with_name(output_image.stem + "_mask.nii.gz")
        sitk.WriteImage(moved_mask, str(moving_mask_out))
        if fixed_mask_path is not None:
            if not fixed_mask_path.exists():
                raise ValueError(f"Fixed mask not found: {fixed_mask_path}")
            fixed_mask = sitk.ReadImage(str(fixed_mask_path))
            if fixed_mask.GetSize() != moved_mask.GetSize() or fixed_mask.GetSpacing() != moved_mask.GetSpacing():
                raise ValueError("Fixed mask geometry must match the fixed image")
            a = sitk.GetArrayFromImage(sitk.Cast(fixed_mask > 0, sitk.sitkUInt8)).astype(bool)
            b = sitk.GetArrayFromImage(sitk.Cast(moved_mask > 0, sitk.sitkUInt8)).astype(bool)
            denominator = int(a.sum() + b.sum())
            dice = float(2 * np.logical_and(a, b).sum() / denominator) if denominator else None

    fixed_array = np.asarray(fixed_array)
    before_array = np.asarray(before_array)
    after_array = np.asarray(after_array)
    qc = {
        "registration": "3D rigid Euler registration, Mattes mutual information",
        "fixed_series_uid": fixed_series_uid,
        "moving_series_uid": moving_series_uid,
        "fixed_instance_count": int(fixed_file_count),
        "moving_instance_count": int(moving_file_count),
        "fixed_size_xyz": [int(x) for x in fixed.GetSize()],
        "moving_size_xyz": [int(x) for x in moving.GetSize()],
        "fixed_spacing_mm": [float(x) for x in fixed.GetSpacing()],
        "moving_spacing_mm": [float(x) for x in moving.GetSpacing()],
        "registered_size_xyz": [int(x) for x in fixed.GetSize()],
        "final_metric_value": float(registration.GetMetricValue()),
        "optimizer_stop_condition": str(registration.GetOptimizerStopConditionDescription()),
        "optimizer_iteration": int(registration.GetOptimizerIteration()),
        "ncc_before": _ncc(fixed_array, before_array),
        "ncc_after": _ncc(fixed_array, after_array),
        "moving_mask_resampled": str(moving_mask_out) if moving_mask_out else None,
        "mask_dice_if_both_masks_supplied": dice,
        "transform_type": "Euler3D rigid",
        "clinical_qc_status": (
            "requires visual review and anatomical landmark/segmentation QC; "
            "numeric registration metrics are not proof of anatomical validity"
        ),
    }
    qc_json.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-dicom-dir", required=True)
    parser.add_argument("--fixed-series-uid", required=True, help="Curated T0 reference series UID")
    parser.add_argument("--moving-dicom-dir", required=True)
    parser.add_argument("--moving-series-uid", required=True, help="Curated T1/T2/T3 series UID")
    parser.add_argument("--output-image", required=True, help="Registered moving volume (.nii.gz)")
    parser.add_argument("--output-transform", required=True, help="Saved rigid transform (.tfm)")
    parser.add_argument("--qc-json", required=True, help="Private local registration QC JSON")
    parser.add_argument("--moving-mask", help="Optional already-converted moving mask in NIfTI")
    parser.add_argument("--fixed-mask", help="Optional fixed/reference mask in NIfTI for Dice QC")
    parser.add_argument("--seed", type=int, default=20261112)
    parser.add_argument("--iterations", type=int, default=300)
    args = parser.parse_args()
    if args.iterations < 10:
        parser.error("--iterations must be at least 10")
    qc = register_series(
        Path(args.fixed_dicom_dir),
        args.fixed_series_uid,
        Path(args.moving_dicom_dir),
        args.moving_series_uid,
        Path(args.output_image),
        Path(args.output_transform),
        Path(args.qc_json),
        moving_mask_path=Path(args.moving_mask) if args.moving_mask else None,
        fixed_mask_path=Path(args.fixed_mask) if args.fixed_mask else None,
        seed=args.seed,
        iterations=args.iterations,
    )
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
