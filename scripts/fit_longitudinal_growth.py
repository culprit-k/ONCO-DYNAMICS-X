"""Fit effective untreated tumour-growth kinetics from a longitudinal CSV."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from onco_dynamics_x.calibration.biological import fit_untreated_growth_curve


def _subject_key(value: str) -> str:
    """Create a stable pseudonymous output key; keep the source CSV private."""
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:12]


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Fit exponential/logistic/Gompertz curves to untreated tumour-volume "
            "series. Output is a local aggregate/pseudonymous summary."
        )
    )
    parser.add_argument("--input", required=True, help="Local CSV path; never uploaded.")
    parser.add_argument("--id-column", default="subject_id")
    parser.add_argument("--time-column", default="time_days")
    parser.add_argument("--volume-column", default="volume_mm3")
    parser.add_argument(
        "--treatment-column",
        help="Optional column identifying treatment status; filtered to untreated-value.",
    )
    parser.add_argument("--untreated-value", default="untreated")
    parser.add_argument("--minimum-visits", type=int, default=5)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", default="results/local_growth_calibration.json")
    args = parser.parse_args()

    if args.minimum_visits < 5:
        parser.error("--minimum-visits cannot be below five.")
    if args.bootstrap < 0:
        parser.error("--bootstrap cannot be negative.")

    grouped: dict[str, list[tuple[float, float]]] = defaultdict(list)
    with Path(args.input).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {args.id_column, args.time_column, args.volume_column}
        if args.treatment_column:
            required.add(args.treatment_column)
        missing = required - set(reader.fieldnames or [])
        if missing:
            parser.error(f"CSV is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            subject = str(row.get(args.id_column, "")).strip()
            if not subject:
                continue
            if args.treatment_column:
                observed_status = str(row.get(args.treatment_column, "")).strip().casefold()
                if observed_status != args.untreated_value.casefold():
                    continue
            try:
                time = float(row[args.time_column])
                volume = float(row[args.volume_column])
            except (TypeError, ValueError):
                continue
            grouped[subject].append((time, volume))

    reports = []
    skipped: dict[str, int] = defaultdict(int)
    for subject, observations in grouped.items():
        if len(observations) < args.minimum_visits:
            skipped["too_few_untreated_visits"] += 1
            continue
        observations.sort(key=lambda item: item[0])
        try:
            fit = fit_untreated_growth_curve(
                [item[0] for item in observations],
                [item[1] for item in observations],
                n_bootstrap=args.bootstrap,
                seed=args.seed,
            )
        except (ValueError, RuntimeError, FloatingPointError):
            skipped["invalid_or_non_growth_series"] += 1
            continue
        reports.append(
            {
                "subject_key": _subject_key(subject),
                "n_observations": fit.n_observations,
                "fit": fit.to_dict(),
            }
        )

    aggregate_by_model = {}
    for model in ("exponential", "logistic", "gompertz"):
        selected = [row["fit"] for row in reports if row["fit"]["model"] == model]
        rate_key = {
            "exponential": "growth_rate_per_day",
            "logistic": "logistic_rate_per_day",
            "gompertz": "gompertz_rate_per_day",
        }[model]
        rates = [
            row["parameters"][rate_key]
            for row in selected
            if rate_key in row["parameters"]
        ]
        if rates:
            import numpy as np

            aggregate_by_model[model] = {
                "series_count": len(rates),
                "median_rate_per_day": float(np.median(rates)),
                "rate_p25_p75_per_day": [
                    float(value) for value in np.quantile(rates, [0.25, 0.75])
                ],
            }

    output = {
        "analysis": "untreated_longitudinal_tumour_volume_growth_calibration",
        "model_scope": "scalar whole-tumour growth curves; not PDE calibration",
        "input_policy": "local input only; raw observations and source identifiers not exported",
        "assumption": (
            "all included observations belong to one untreated interval per subject; "
            "verify this from source metadata before interpretation"
        ),
        "time_unit": "day",
        "volume_unit": "mm^3",
        "series_fitted": len(reports),
        "series_skipped": dict(skipped),
        "aggregate_by_selected_model": aggregate_by_model,
        "series": reports,
        "limitations": [
            "AICc selects among descriptive curves; it does not prove a biological growth law.",
            "A fitted whole-tumour rate is not directly equal to the dimensionless PDE cancer_growth coefficient.",
            "Carrying capacity is weakly identified without observed saturation; inspect each fit warning and Jacobian condition number.",
            "Treatment-exposed series must be analysed with a treatment-effect model, not this untreated-growth estimator.",
            "Patient-level outputs remain local and should not be committed to the public repository.",
        ],
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Series fitted: {len(reports)}")
    print(f"Series skipped: {dict(skipped)}")
    print(f"Aggregate output: {output_path}")


if __name__ == "__main__":
    main()
