"""Resample a visit DCE map or label mask using a T2-derived rigid transform.

A transform estimated between visit anatomical images is only valid for the
DCE maps if each visit map shares the same DICOM FrameOfReferenceUID as that
visit's anatomical reference. The corresponding baseline map defines output
geometry. This check is mandatory and cannot be bypassed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def resample_map(
    fixed_map_reference: Path,
    moving_map: Path,
    transform_path: Path,
    output_path: Path,
    qc_path: Path,
    *,
    fixed_map_frame_uid: str,
    fixed_anatomical_frame_uid: str,
    moving_map_frame_uid: str,
    moving_anatomical_frame_uid: str,
    kind: str = "map",
) -> dict[str, Any]:
    try:
        import SimpleITK as sitk
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("Install the optional imaging extra: pip install -e '.[imaging]'") from exc

    pairs = [
        ("fixed map", fixed_map_frame_uid, fixed_anatomical_frame_uid),
        ("moving map", moving_map_frame_uid, moving_anatomical_frame_uid),
    ]
    for label, map_uid, anatomy_uid in pairs:
        if not map_uid or not anatomy_uid or map_uid != anatomy_uid:
            raise ValueError(
                f"{label} FrameOfReferenceUID does not match its anatomical reference; "
                "a same-visit map-to-anatomy registration is required first"
            )
    if kind not in {"map", "label"}:
        raise ValueError("kind must be 'map' or 'label'")
    for path in (fixed_map_reference, moving_map, transform_path):
        if not path.exists():
            raise ValueError(f"Input file does not exist: {path}")

    fixed = sitk.ReadImage(str(fixed_map_reference))
    moving = sitk.ReadImage(str(moving_map))
    transform = sitk.ReadTransform(str(transform_path))
    if fixed.GetDimension() != 3 or moving.GetDimension() != 3:
        raise ValueError("Fixed and moving DCE maps must be 3D")
    if transform.GetDimension() != 3:
        raise ValueError("Registration transform must be 3D")
    if any(size < 2 for size in fixed.GetSize() + moving.GetSize()):
        raise ValueError("DCE map dimensions are too small")

    interpolator = sitk.sitkLinear if kind == "map" else sitk.sitkNearestNeighbor
    output_type = sitk.sitkFloat32 if kind == "map" else sitk.sitkUInt8
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fixed)
    resampler.SetTransform(transform)
    resampler.SetInterpolator(interpolator)
    resampler.SetDefaultPixelValue(0)
    resampler.SetOutputPixelType(output_type)
    registered = resampler.Execute(moving)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    qc_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(registered, str(output_path))
    qc = {
        "kind": kind,
        "transform_direction": "moving visit physical coordinates -> baseline visit physical coordinates",
        "interpolator": "linear" if kind == "map" else "nearest_neighbor",
        "fixed_map_size_xyz": [int(x) for x in fixed.GetSize()],
        "moving_map_size_xyz": [int(x) for x in moving.GetSize()],
        "output_size_xyz": [int(x) for x in registered.GetSize()],
        "fixed_map_spacing_mm": [float(x) for x in fixed.GetSpacing()],
        "moving_map_spacing_mm": [float(x) for x in moving.GetSpacing()],
        "fixed_map_frame_uid": fixed_map_frame_uid,
        "moving_map_frame_uid": moving_map_frame_uid,
        "frame_of_reference_checks_passed": True,
        "same_visit_map_anatomy_frame_checks_passed": True,
        "note": (
            "The same-visit FrameOfReferenceUID check does not prove registration "
            "quality or correct DICOM-to-NIfTI conversion. Inspect overlay and "
            "landmark/segmentation QC before biological modelling."
        ),
    }
    qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return qc


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixed-map-reference", required=True, help="Baseline T0 PE/SER/label map defining output geometry")
    parser.add_argument("--moving-map", required=True, help="Follow-up PE/SER/label map")
    parser.add_argument("--transform", required=True, help="Visit anatomical moving-to-fixed transform")
    parser.add_argument("--output", required=True, help="Resampled map or label output")
    parser.add_argument("--qc-json", required=True, help="Private local QC report")
    parser.add_argument("--fixed-map-frame-uid", required=True)
    parser.add_argument("--fixed-anatomical-frame-uid", required=True)
    parser.add_argument("--moving-map-frame-uid", required=True)
    parser.add_argument("--moving-anatomical-frame-uid", required=True)
    parser.add_argument("--kind", choices=["map", "label"], default="map")
    args = parser.parse_args()
    qc = resample_map(
        Path(args.fixed_map_reference),
        Path(args.moving_map),
        Path(args.transform),
        Path(args.output),
        Path(args.qc_json),
        fixed_map_frame_uid=args.fixed_map_frame_uid,
        fixed_anatomical_frame_uid=args.fixed_anatomical_frame_uid,
        moving_map_frame_uid=args.moving_map_frame_uid,
        moving_anatomical_frame_uid=args.moving_anatomical_frame_uid,
        kind=args.kind,
    )
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
