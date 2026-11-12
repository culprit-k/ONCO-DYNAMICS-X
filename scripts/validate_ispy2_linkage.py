"""Validate an explicit I-SPY 2 feature-to-clinical workbook join locally.

The default report contains counts and schema metadata only, never key values or
participant-level records. A merged file is written only when explicitly requested.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


def _read_table(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    suffix = source.suffix.casefold()
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(source)
    if suffix in {".csv", ".tsv"}:
        return pd.read_csv(source, sep="\t" if suffix == ".tsv" else ",", dtype="string")
    raise ValueError(f"Unsupported input extension: {source.suffix}")


def _normalise_key(series: pd.Series) -> pd.Series:
    key = series.astype("string").str.strip()
    return key.mask(key.eq(""))


def validate_linkage(
    feature: pd.DataFrame,
    clinical: pd.DataFrame,
    *,
    feature_key: str,
    clinical_key: str,
    cardinality: str = "one_to_one",
) -> tuple[dict[str, Any], pd.DataFrame]:
    if cardinality not in {"one_to_one", "many_to_one"}:
        raise ValueError(
            "cardinality must be one_to_one or many_to_one; many-to-many is prohibited"
        )
    if feature_key not in feature.columns:
        raise ValueError(f"Feature key column not found: {feature_key}")
    if clinical_key not in clinical.columns:
        raise ValueError(f"Clinical key column not found: {clinical_key}")
    if "__join_key" in feature.columns or "__join_key" in clinical.columns:
        raise ValueError("Input tables contain a reserved internal column name: __join_key")

    f = feature.copy()
    c = clinical.copy()
    f["__join_key"] = _normalise_key(f[feature_key])
    c["__join_key"] = _normalise_key(c[clinical_key])
    f_valid = f["__join_key"].notna()
    c_valid = c["__join_key"].notna()
    f_dup = int(f.loc[f_valid, "__join_key"].duplicated(keep=False).sum())
    c_dup = int(c.loc[c_valid, "__join_key"].duplicated(keep=False).sum())

    if c_dup:
        raise ValueError(
            "Clinical linkage key is not unique; resolve duplicate clinical records first."
        )
    if cardinality == "one_to_one" and f_dup:
        raise ValueError("Feature linkage key is not unique for one_to_one validation.")

    f_keys = set(f.loc[f_valid, "__join_key"].tolist())
    c_keys = set(c.loc[c_valid, "__join_key"].tolist())
    matched_keys = f_keys & c_keys
    matched_feature_rows = int((f_valid & f["__join_key"].isin(c_keys)).sum())
    matched_clinical_rows = int((c_valid & c["__join_key"].isin(f_keys)).sum())

    clinical_payload = c.drop(columns=[clinical_key])
    merged = f.merge(
        clinical_payload,
        on="__join_key",
        how="inner",
        validate="one_to_one" if cardinality == "one_to_one" else "many_to_one",
        suffixes=("", "__clinical"),
    )
    report = {
        "analysis": "ispy2_feature_clinical_linkage_validation",
        "cardinality_requested": cardinality,
        "feature_rows": int(len(feature)),
        "clinical_rows": int(len(clinical)),
        "feature_missing_key_rows": int((~f_valid).sum()),
        "clinical_missing_key_rows": int((~c_valid).sum()),
        "feature_duplicate_key_rows": f_dup,
        "clinical_duplicate_key_rows": c_dup,
        "unique_feature_keys": int(len(f_keys)),
        "unique_clinical_keys": int(len(c_keys)),
        "matched_unique_keys": int(len(matched_keys)),
        "matched_feature_rows": matched_feature_rows,
        "matched_clinical_rows": matched_clinical_rows,
        "unmatched_feature_rows": int(len(feature) - matched_feature_rows),
        "unmatched_clinical_rows": int(len(clinical) - matched_clinical_rows),
        "merged_rows_if_written": int(len(merged)),
        "feature_key_column": feature_key,
        "clinical_key_column": clinical_key,
        "feature_columns": [str(x) for x in feature.columns],
        "clinical_columns": [str(x) for x in clinical.columns],
        "privacy": "report contains no key values or row-level observations",
        "status": (
            "complete_overlap"
            if len(matched_keys) == len(f_keys) == len(c_keys)
            else "partial_overlap"
        ),
    }
    return report, merged.drop(columns=["__join_key"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, help="Local multi-feature workbook/CSV.")
    parser.add_argument("--clinical", required=True, help="Local clinical workbook/CSV.")
    parser.add_argument("--feature-key", required=True, help="Documented feature-table key column.")
    parser.add_argument(
        "--clinical-key", required=True, help="Documented clinical-table key column."
    )
    parser.add_argument(
        "--cardinality",
        choices=("one_to_one", "many_to_one"),
        default="one_to_one",
        help="Choose only after confirming source semantics; many-to-many is prohibited.",
    )
    parser.add_argument("--report", default="results/local/ispy2_linkage_audit.json")
    parser.add_argument(
        "--linked-output",
        help="Optional local CSV of matched rows. Never commit this file or its inputs.",
    )
    args = parser.parse_args()
    feature = _read_table(args.features)
    clinical = _read_table(args.clinical)
    report, merged = validate_linkage(
        feature,
        clinical,
        feature_key=args.feature_key,
        clinical_key=args.clinical_key,
        cardinality=args.cardinality,
    )
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.linked_output:
        linked_path = Path(args.linked_output)
        linked_path.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(linked_path, index=False)
    printable = {
        key: value
        for key, value in report.items()
        if key not in {"feature_columns", "clinical_columns"}
    }
    print(json.dumps(printable, indent=2))
    print(f"Aggregate linkage report: {report_path}")
    if args.linked_output:
        print(f"Local matched file written: {args.linked_output}")


if __name__ == "__main__":
    main()
