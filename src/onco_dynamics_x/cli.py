"""Command-line entry points."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import simulate
from onco_dynamics_x.state import initial_state


def run_simulation(args: argparse.Namespace) -> Path:
    parameters = ModelParameters()
    state = initial_state(args.grid_size)
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="onco-dynamics",
        description="Run the ONCO-DYNAMICS X dimensionless 2D demonstration model.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    sim = subparsers.add_parser("simulate", help="run the reaction-diffusion demo")
    sim.add_argument("--grid-size", type=int, default=51)
    sim.add_argument("--steps", type=int, default=400)
    sim.add_argument("--dt", type=float, default=0.0003)
    sim.add_argument("--save-every", type=int, default=20)
    sim.add_argument("--output", default="results/simulation_summary.csv")
    sim.set_defaults(func=run_simulation)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.grid_size < 5:
        parser.error("--grid-size must be at least 5")
    if args.steps < 1 or args.dt <= 0:
        parser.error("--steps and --dt must be positive")
    args.func(args)


if __name__ == "__main__":
    main()
