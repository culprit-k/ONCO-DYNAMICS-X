"""Run reproducible synthetic numerical checks for the PDE prototype."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state
from onco_dynamics_x.validation.diagnostics import grid_refinement, timestep_refinement
from onco_dynamics_x.validation.reference_solver import (
    compare_explicit_euler_to_reference,
    simulate_adaptive_reference,
)


def run_numerical_validation() -> dict[str, Any]:
    parameters = ModelParameters()
    report: dict[str, Any] = {
        "analysis": "synthetic_pde_numerical_validation",
        "data": "deterministic synthetic initial conditions only; no patient data",
        "parameters": "default dimensionless ModelParameters; not clinically calibrated",
        "integrator": "SciPy adaptive DOP853 method of lines",
        "results": {},
        "interpretation": (
            "Checks numerical consistency and refinement behaviour on synthetic "
            "fields. It does not validate biological assumptions or clinical performance."
        ),
    }

    for dimensions, coarse_size, coarse_steps in ((2, 9, 8), (3, 5, 4)):
        state = initial_state(coarse_size, spatial_dimensions=dimensions)
        dx = 2.0 / (coarse_size - 1)
        dt = 0.001
        diffusion_limit = diffusion_timestep_limit(parameters, dx, dimensions)
        if dt >= diffusion_limit:
            dt = diffusion_limit * 0.1

        adaptive, adaptive_diagnostics = simulate_adaptive_reference(
            state,
            parameters,
            dx,
            final_time=dt * coarse_steps,
            method="DOP853",
            rtol=1e-8,
            atol=1e-10,
        )
        adaptive.validate()
        explicit_comparison = compare_explicit_euler_to_reference(
            state,
            parameters,
            dx=dx,
            dt=dt,
            steps=coarse_steps,
            method="DOP853",
            rtol=1e-8,
            atol=1e-10,
        )
        time_error = timestep_refinement(
            state, parameters, dx=dx, dt=dt, steps=coarse_steps
        )

        fine_size = (coarse_size - 1) * 2 + 1
        fine_state = initial_state(fine_size, spatial_dimensions=dimensions)
        space_error = grid_refinement(
            state,
            fine_state,
            parameters,
            coarse_dx=dx,
            coarse_dt=dt,
            coarse_steps=coarse_steps,
            refinement_ratio=2,
        )
        report["results"][f"{dimensions}d"] = {
            "grid_shape": list(state.shape),
            "coarse_dx": float(dx),
            "coarse_dt": float(dt),
            "coarse_steps": int(coarse_steps),
            "final_time": float(dt * coarse_steps),
            "diffusion_timestep_limit": float(diffusion_limit),
            "timestep_refinement_relative_l2_error": time_error,
            "grid_refinement_relative_l2_error": space_error,
            "explicit_vs_adaptive_reference": explicit_comparison,
            "adaptive_reference_solver": adaptive_diagnostics,
        }

    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="results/numerical_validation/synthetic_pde_validation.json",
    )
    args = parser.parse_args()
    report = run_numerical_validation()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for dimensions, result in report["results"].items():
        print(dimensions, "max explicit/reference error:",
              result["explicit_vs_adaptive_reference"]["maximum_relative_l2_error"])
        print(dimensions, "mean timestep refinement error:",
              sum(result["timestep_refinement_relative_l2_error"].values()) / 8)
        print(dimensions, "mean grid refinement error:",
              sum(result["grid_refinement_relative_l2_error"].values()) / 8)
    print("Aggregate synthetic validation report:", output)


if __name__ == "__main__":
    main()
