"""Run the public breast-cancer biological calibration target audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from onco_dynamics_x.data_engine.biological_calibration_audit import run_audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="results/biological_calibration/target_availability.json",
    )
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()
    report = run_audit(timeout=args.timeout, retries=args.retries)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for cohort in report["cohorts"]:
        print(f"{cohort['study_id']}: {cohort['matched_target_categories']}")
    print(f"Aggregate audit: {path}")


if __name__ == "__main__":
    main()
