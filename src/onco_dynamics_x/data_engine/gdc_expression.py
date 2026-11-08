"""Assemble downloaded open GDC expression files into a linked sample-by-gene matrix."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from onco_dynamics_x.data_engine.expression import ExpressionMatrix


def read_gdc_expression_table(
    path: str | Path,
    gene_id_column: str = "gene_id",
    expression_column: str = "tpm_unstranded",
    skip_summary_rows: bool = True,
) -> dict[str, float]:
    """Read one tab-separated gene-expression table with explicit column names."""
    source = Path(path)
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        required = {gene_id_column, expression_column}
        header_fields: list[str] | None = None
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            candidate = [
                field.strip().lstrip("#").strip()
                for field in line.rstrip("\r\n").split("\t")
            ]
            if required.issubset(candidate):
                header_fields = candidate
                break
        if header_fields is None:
            raise ValueError(
                f"{source.name} has no tab-separated header containing "
                f"{gene_id_column!r} and {expression_column!r}."
            )
        reader = csv.DictReader(handle, fieldnames=header_fields, delimiter="\t")
        values: dict[str, float] = {}
        for line_number, row in enumerate(reader, start=line_number + 1):
            gene_id = (row.get(gene_id_column) or "").strip()
            if not gene_id:
                raise ValueError(f"{source.name}:{line_number}: empty gene identifier.")
            if skip_summary_rows and gene_id.startswith("N_"):
                continue
            if gene_id in values:
                raise ValueError(f"{source.name}:{line_number}: duplicate gene ID '{gene_id}'.")
            raw_value = (row.get(expression_column) or "").strip()
            try:
                value = float(raw_value)
            except ValueError as error:
                raise ValueError(
                    f"{source.name}:{line_number}: non-numeric expression for {gene_id}."
                ) from error
            if not np.isfinite(value):
                raise ValueError(
                    f"{source.name}:{line_number}: non-finite expression for {gene_id}."
                )
            values[gene_id] = value
    if not values:
        raise ValueError(f"{source.name} contains no gene-expression rows.")
    return values


def build_gdc_expression_matrix(
    metadata_csv: str | Path,
    data_directory: str | Path,
    gene_id_column: str = "gene_id",
    expression_column: str = "tpm_unstranded",
    output_npz: str | Path | None = None,
) -> ExpressionMatrix:
    """Build a strict sample-by-gene matrix from downloaded files in GDC metadata.

    Every selected metadata row must describe one open-access file linked to
    exactly one sample submitter ID. Files must exist locally and have the exact
    same gene set; columns are aligned to the first file's gene order.
    """
    metadata_path = Path(metadata_csv)
    data_root = Path(data_directory)
    if not metadata_path.exists():
        raise FileNotFoundError(metadata_path)
    if not data_root.is_dir():
        raise NotADirectoryError(data_root)
    required = {
        "file_name",
        "access",
        "sample_submitter_ids",
    }
    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Metadata CSV is missing columns: {', '.join(sorted(missing))}."
            )
        records = list(reader)
    if not records:
        raise ValueError("Metadata CSV contains no file records.")

    sample_ids: list[str] = []
    expression_by_sample: list[dict[str, float]] = []
    reference_genes: tuple[str, ...] | None = None
    for row in records:
        if (row.get("access") or "").strip().lower() != "open":
            raise PermissionError("Expression matrix assembly accepts open-access files only.")
        sample_ids_for_file = tuple(
            item.strip()
            for item in (row.get("sample_submitter_ids") or "").split("|")
            if item.strip()
        )
        if len(sample_ids_for_file) != 1:
            raise ValueError(
                "Each expression file must be linked to exactly one sample submitter ID."
            )
        sample_id = sample_ids_for_file[0]
        if sample_id in sample_ids:
            raise ValueError(f"Duplicate sample identifier in metadata: {sample_id}.")
        file_name = (row.get("file_name") or "").strip()
        if not file_name or Path(file_name).name != file_name:
            raise ValueError("GDC filename must be a plain filename without path components.")
        file_path = data_root / file_name
        if not file_path.is_file():
            raise FileNotFoundError(
                f"Downloaded expression file is missing for sample {sample_id}: {file_path}"
            )
        values = read_gdc_expression_table(
            file_path,
            gene_id_column=gene_id_column,
            expression_column=expression_column,
        )
        genes = tuple(values)
        if reference_genes is None:
            reference_genes = genes
        elif set(genes) != set(reference_genes):
            missing_genes = set(reference_genes) - set(genes)
            extra_genes = set(genes) - set(reference_genes)
            raise ValueError(
                "Gene sets differ between sample files "
                f"(missing={len(missing_genes)}, extra={len(extra_genes)})."
            )
        sample_ids.append(sample_id)
        expression_by_sample.append(values)

    assert reference_genes is not None
    matrix = np.asarray(
        [
            [sample[gene] for gene in reference_genes]
            for sample in expression_by_sample
        ],
        dtype=float,
    )
    result = ExpressionMatrix(
        sample_ids=tuple(sample_ids),
        gene_ids=reference_genes,
        values=matrix,
    )
    result.validate()
    if output_npz is not None:
        output_path = Path(output_npz)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            output_path,
            sample_ids=np.asarray(result.sample_ids, dtype=str),
            gene_ids=np.asarray(result.gene_ids, dtype=str),
            values=result.values,
            gene_id_column=np.asarray(gene_id_column),
            expression_column=np.asarray(expression_column),
        )
    return result
