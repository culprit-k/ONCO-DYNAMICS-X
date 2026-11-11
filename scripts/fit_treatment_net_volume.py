"""Estimate treatment-exposed net tumour-volume rates from a local CSV."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any


from onco_dynamics_x.calibration.net_volume import (
    bootstrap_group_net_rates,
    fit_treatment_exposed_net_volume_rate,
)


def _key(identifier: str) -> str:
    return hashlib.sha256(identifier.encode("utf-8")).hexdigest()[:12]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Local CSV; not uploaded.")
    parser.add_argument("--id-column", default="subject_id")
    parser.add_argument("--time-column", default="time_days")
    parser.add_argument("--volume-column", default="ftv_mm3")
    parser.add_argument("--group-column", help="Optional regimen/subtype column.")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output",
        default="results/local_treatment_net_volume_rates.json",
    )
    args = parser.parse_args()
    if args.bootstrap < 100:
        parser.error("--bootstrap must be at least 100.")

    grouped: dict[str, dict[str, Any]] = {}
    with Path(args.input).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {args.id_column, args.time_column, args.volume_column}
        if args.group_column:
            required.add(args.group_column)
        missing = required - set(reader.fieldnames or [])
        if missing:
            parser.error(f"CSV is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            identifier = str(row.get(args.id_column, "")).strip()
            if not identifier:
                continue
            try:
                time = float(row[args.time_column])
                volume = float(row[args.volume_column])
            except (TypeError, ValueError):
                continue
            group = (
                str(row.get(args.group_column, "unknown")).strip() or "unknown"
                if args.group_column
                else "all"
            )
            item = grouped.setdefault(
                identifier, {"group": group, "observations": []}
            )
            if item["group"] != group:
                parser.error(
                    f"Subject {identifier} appears in multiple {args.group_column} groups."
                )
            item["observations"].append((time, volume))

    reports = []
    skipped: dict[str, int] = defaultdict(int)
    rates_by_group: dict[str, list[float]] = defaultdict(list)
    for identifier, item in grouped.items():
        try:
            observations = sorted(item["observations"], key=lambda pair: pair[0])
            fit = fit_treatment_exposed_net_volume_rate(
                [pair[0] for pair in observations],
                [pair[1] for pair in observations],
            )
        except ValueError as error:
            message = str(error)
            reason = (
                "censored_or_nonpositive_volume"
                if "left-censored" in message
                else "insufficient_or_invalid_visits"
            )
            skipped[reason] += 1
            continue
        rates_by_group[item["group"]].append(fit.net_log_rate_per_day)
        reports.append(
            {
                "subject_key": _key(identifier),
                "group": item["group"],
                "fit": fit.to_dict(),
            }
        )

    group_summary = {}
    for index, (group, rates) in enumerate(sorted(rates_by_group.items())):
        if len(rates) < 2:
            group_summary[group] = {
                "subject_series_count": len(rates),
                "status": "insufficient_for_group_bootstrap",
            }
        else:
            group_summary[group] = bootstrap_group_net_rates(
                rates,
                n_bootstrap=args.bootstrap,
                seed=args.seed + index,
            )

    output = {
        "analysis": "treatment_exposed_net_functional_tumour_volume_change",
        "input_policy": "local input only; no raw observations or source IDs exported",
        "model": "log(volume) = intercept + net_log_rate_per_day * elapsed_days",
        "units": {"time": "day", "volume": "mm^3"},
        "series_fitted": len(reports),
        "series_skipped": dict(skipped),
        "group_summary": group_summary,
        "series": reports,
        "limitations": [
            "The estimated slope is net volume change under the observed treatment context.",
            "It does not separate intrinsic proliferation from treatment-induced death.",
            "It is not a drug EC50, drug-killing rate, or PDE cancer_growth coefficient.",
            "Use actual acquisition days, not nominal visit labels, as time values.",
            "Zero/below-detection FTV is censored and requires a censored-data model.",
            "Patient-level outputs must remain in approved local/institutional storage.",
        ],
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"Series fitted: {len(reports)}")
    print(f"Series skipped: {dict(skipped)}")
    print(f"Aggregate output: {output_path}")


if __name__ == "__main__":
    main()
