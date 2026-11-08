"""Command-line entry points."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from onco_dynamics_x.data_engine.gdc import (
    GDCClient,
    GDCFile,
    write_sample_manifest,
)
from onco_dynamics_x.data_engine.gdc_expression import build_gdc_expression_matrix
from onco_dynamics_x.data_engine.manifest import audit_manifest
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import initial_state


def run_simulation(args: argparse.Namespace) -> Path:
    parameters = ModelParameters()
    state = initial_state(args.grid_size, spatial_dimensions=args.dimensions)
    dx = 2.0 / (args.grid_size - 1)
    times, snapshots = simulate(
        state,
        parameters,
        dx=dx,
        dt=args.dt,
        steps=args.steps,
        save_every=max(1, args.save_every),
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["time", *state.FIELD_NAMES])
        for time, snapshot in zip(times, snapshots, strict=True):
            writer.writerow(
                [f"{time:.8g}"]
                + [
                    f"{float(np.mean(snapshot.as_dict()[name])):.8g}"
                    for name in state.FIELD_NAMES
                ]
            )
    print(f"Saved {len(times)} snapshots to {output}")
    print(f"Final mean cancer-cell field: {float(np.mean(snapshots[-1].cancer)):.5f}")
    return output


def run_gdc_query(args: argparse.Namespace) -> Path:
    """Query public GDC metadata and create a sample-level audit manifest."""
    client = GDCClient(timeout=args.timeout)
    files = client.query_open_files(
        project_id=args.project_id,
        data_type=args.data_type,
        workflow_type=args.workflow_type,
        page_size=args.page_size,
        max_results=args.max_results,
    )
    metadata_path = Path(args.metadata_output)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "access",
        "data_type",
        "workflow_type",
        "case_ids",
        "patient_submitter_ids",
        "sample_submitter_ids",
    ]
    with metadata_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for item in files:
            writer.writerow(
                {
                    "file_id": item.file_id,
                    "file_name": item.file_name,
                    "md5sum": item.md5sum,
                    "file_size": item.file_size,
                    "access": item.access,
                    "data_type": item.data_type,
                    "workflow_type": item.workflow_type,
                    "case_ids": "|".join(item.case_ids),
                    "patient_submitter_ids": "|".join(item.patient_submitter_ids),
                    "sample_submitter_ids": "|".join(item.sample_submitter_ids),
                }
            )

    manifest_path = write_sample_manifest(
        files, args.manifest_output, args.project_id, args.modality
    )
    audit = audit_manifest(manifest_path)
    print(f"Found {len(files)} open-access file records.")
    print(f"Metadata: {metadata_path}")
    print(f"Sample manifest: {manifest_path}")
    print(
        f"Manifest audit: {'PASSED' if audit.passed else 'FAILED'}; "
        f"{audit.unique_patients} patients, {audit.unique_samples} samples."
    )
    for warning in audit.warnings:
        print(f"WARNING: {warning}")
    for error in audit.errors:
        print(f"ERROR: {error}")
    return metadata_path


def run_gdc_download_open(args: argparse.Namespace) -> Path:
    """Download one explicitly selected open GDC file from a prior metadata CSV."""
    metadata_path = Path(args.metadata_csv)
    if not metadata_path.exists():
        raise FileNotFoundError(metadata_path)
    required = {
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "access",
        "data_type",
        "workflow_type",
    }
    with metadata_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Metadata CSV is missing columns: {', '.join(sorted(missing))}."
            )
        rows = [row for row in reader if row.get("file_id") == args.file_id]
    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one metadata row for file ID {args.file_id}; found {len(rows)}."
        )
    row = rows[0]
    if row["access"].lower() != "open":
        raise PermissionError("Selected GDC file is not marked open-access.")
    metadata = GDCFile(
        file_id=row["file_id"],
        file_name=row["file_name"],
        md5sum=row["md5sum"],
        file_size=int(row["file_size"]),
        access=row["access"].lower(),
        data_type=row["data_type"],
        workflow_type=row["workflow_type"],
        case_ids=tuple(filter(None, row.get("case_ids", "").split("|"))),
        patient_submitter_ids=tuple(
            filter(None, row.get("patient_submitter_ids", "").split("|"))
        ),
        sample_submitter_ids=tuple(
            filter(None, row.get("sample_submitter_ids", "").split("|"))
        ),
        associated_samples=(),
    )
    downloaded = GDCClient(timeout=args.timeout).download_open_file(
        metadata, args.destination
    )
    print(f"Downloaded and checksum-verified open file: {downloaded}")
    return downloaded


def run_gdc_build_expression_matrix(args: argparse.Namespace) -> Path:
    """Assemble downloaded sample-level open GDC expression files."""
    matrix = build_gdc_expression_matrix(
        metadata_csv=args.metadata_csv,
        data_directory=args.data_directory,
        gene_id_column=args.gene_id_column,
        expression_column=args.expression_column,
        output_npz=args.output,
    )
    print(
        f"Saved expression matrix {matrix.values.shape[0]} samples x "
        f"{matrix.values.shape[1]} genes to {args.output}"
    )
    return Path(args.output)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="onco-dynamics",
        description="Run ONCO-DYNAMICS X research data and simulation tools.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sim = subparsers.add_parser("simulate", help="run the reaction-diffusion demo")
    sim.add_argument("--grid-size", type=int, default=51)
    sim.add_argument("--dimensions", type=int, choices=(2, 3), default=2)
    sim.add_argument("--steps", type=int, default=400)
    sim.add_argument("--dt", type=float, default=0.0003)
    sim.add_argument("--save-every", type=int, default=20)
    sim.add_argument("--output", default="results/simulation_summary.csv")
    sim.set_defaults(func=run_simulation)

    gdc = subparsers.add_parser(
        "gdc-query",
        help="query open-access GDC file metadata (does not download files)",
    )
    gdc.add_argument("--project-id", default="TCGA-BRCA")
    gdc.add_argument("--data-type", default="Gene Expression Quantification")
    gdc.add_argument("--workflow-type", default=None)
    gdc.add_argument("--page-size", type=int, default=100)
    gdc.add_argument("--max-results", type=int, default=100)
    gdc.add_argument("--timeout", type=float, default=60.0)
    gdc.add_argument("--modality", default="transcriptomics")
    gdc.add_argument("--metadata-output", default="data/metadata/gdc_open_files.csv")
    gdc.add_argument("--manifest-output", default="data/metadata/gdc_sample_manifest.csv")
    gdc.set_defaults(func=run_gdc_query)

    download = subparsers.add_parser(
        "gdc-download-open",
        help="download one explicitly selected, metadata-verified open GDC file",
    )
    download.add_argument("--metadata-csv", required=True)
    download.add_argument("--file-id", required=True)
    download.add_argument("--destination", default="data/raw/gdc")
    download.add_argument("--timeout", type=float, default=120.0)
    download.set_defaults(func=run_gdc_download_open)

    matrix = subparsers.add_parser(
        "gdc-build-expression-matrix",
        help="assemble downloaded open GDC expression files into a sample matrix",
    )
    matrix.add_argument("--metadata-csv", required=True)
    matrix.add_argument("--data-directory", default="data/raw/gdc")
    matrix.add_argument("--gene-id-column", default="gene_id")
    matrix.add_argument("--expression-column", default="tpm_unstranded")
    matrix.add_argument("--output", default="data/processed/tcga_brca_expression.npz")
    matrix.set_defaults(func=run_gdc_build_expression_matrix)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "simulate":
        if args.grid_size < 5:
            parser.error("--grid-size must be at least 5")
        if args.steps < 1 or args.dt <= 0:
            parser.error("--steps and --dt must be positive")
    args.func(args)


if __name__ == "__main__":
    main()
