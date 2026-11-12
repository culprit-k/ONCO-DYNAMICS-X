"""Privacy-preserving structural audit of the public I-SPY 2 MRI Excel workbooks.

The script reports workbook/sheet schemas and aggregate completeness only. It never
writes participant identifiers, cell values, row-level records, or source IDs.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


ROLE_PATTERNS: dict[str, tuple[str, ...]] = {
    "functional_tumour_volume": (r"\bftv\b", r"functional.?tumou?r.?volume"),
    "longest_diameter": (r"longest.?diam", r"\bld\b"),
    "sphericity": (r"spheric"),
    "background_parenchymal_enhancement": (r"background.*enhanc", r"\bbpe\b"),
    "visit_or_time": (r"visit", r"time.?point", r"scan.?date", r"acquisition.?date", r"day"),
    "treatment_or_regimen": (r"treatment", r"regimen", r"drug.?arm", r"therapy", r"arm"),
    "molecular_subtype_or_receptor": (r"subtype", r"\ber\b", r"\bpr\b", r"her2", r"hr.?status"),
    "pathologic_response": (r"pcr", r"pathologic.*response", r"response"),
    "patient_or_sample_key": (r"patient.?id", r"subject.?id", r"case.?id", r"participant.?id", r"study.?id"),
}


def classify_column(name: Any) -> list[str]:
    text = re.sub(r"[_-]+", " ", str(name)).strip().lower()
    return [
        role
        for role, patterns in ROLE_PATTERNS.items()
        if any(re.search(pattern, text, flags=re.IGNORECASE) for pattern in patterns)
    ]


def audit_workbook(path: str | Path) -> dict[str, Any]:
    """Return a value-free workbook audit suitable for sharing."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)

    book = pd.ExcelFile(source)
    report: dict[str, Any] = {
        "audit": "ispy2_excel_structure_and_completeness",
        "source_filename": source.name,
        "privacy_policy": (
            "Aggregate metadata only; no participant IDs, cell values, or row-level data."
        ),
        "sheets": [],
    }

    for sheet_name in book.sheet_names:
        frame = pd.read_excel(book, sheet_name=sheet_name)
        columns = []
        for column in frame.columns:
            series = frame[column]
            non_null = int(series.notna().sum())
            roles = classify_column(column)
            columns.append(
                {
                    "column_name": str(column),
                    "candidate_roles": roles,
                    "dtype": str(series.dtype),
                    "non_null_count": non_null,
                    "missing_count": int(series.isna().sum()),
                    "missing_fraction": (
                        float(series.isna().mean()) if len(frame) else None
                    ),
                    "unique_non_null_count": int(series.nunique(dropna=True)),
                }
            )

        report["sheets"].append(
            {
                "sheet_name": str(sheet_name),
                "row_count": int(len(frame)),
                "column_count": int(len(frame.columns)),
                "columns": columns,
                "candidate_role_counts": {
                    role: sum(role in col["candidate_roles"] for col in columns)
                    for role in ROLE_PATTERNS
                },
            }
        )

    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Local TCIA .xlsx workbook.")
    parser.add_argument(
        "--output",
        default="results/local_ispy2_workbook_audit.json",
        help="Aggregate JSON report; keep source workbook local.",
    )
    args = parser.parse_args()

    report = audit_workbook(args.input)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Sheets audited: {len(report['sheets'])}")
    for sheet in report["sheets"]:
        print(
            f"- {sheet['sheet_name']}: {sheet['row_count']} rows, "
            f"{sheet['column_count']} columns"
        )
    print(f"Aggregate audit written to: {output}")


if __name__ == "__main__":
    main()
