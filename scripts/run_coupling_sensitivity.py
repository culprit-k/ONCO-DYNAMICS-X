"""Run a synthetic parameter sweep for the bidirectional PDE-agent bridge.

This experiment quantifies software sensitivity to dimensionless coupling
controls. It is not parameter calibration and uses no patient data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from onco_dynamics_x.cellular_simulation.agent_based import (
    CANCER,
    FIBROBLAST,
    IMMUNE,
    initialize_agents,
)
from onco_dynamics_x.cellular_simulation.coupling import simulate_coupled
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state


COUPLING_STRENGTHS = (0.0, 0.25, 0.5, 0.75, 1.0)
AGENTS_PER_UNIT = (1.0, 2.0, 4.0, 8.0)


def _population_counts(population) -> dict[str, int]:
    return {
        "cancer": population.count(CANCER),
        "immune": population.count(IMMUNE),
        "fibroblast": population.count(FIBROBLAST),
        "total": len(population.x),
    }


def _field_summary(state, parameters: ModelParameters) -> dict[str, Any]:
    summaries = {}
    upper_bounds = {
        "cancer": parameters.cancer_capacity,
        "immune": parameters.maximum_field,
        "fibroblast": parameters.maximum_field,
        "oxygen": parameters.maximum_field,
        "nutrient": parameters.maximum_field,
        "drug": parameters.maximum_drug,
        "evasion": parameters.maximum_field,
        "vasculature": parameters.maximum_field,
    }
    for name, values in state.as_dict().items():
        array = np.asarray(values, dtype=float)
        upper = upper_bounds[name]
        saturated = (array <= 1e-10) | (array >= upper - 1e-10)
        summaries[name] = {
            "mean": float(np.mean(array)),
            "standard_deviation": float(np.std(array)),
            "minimum": float(np.min(array)),
            "maximum": float(np.max(array)),
            "fraction_at_or_near_bounds": float(np.mean(saturated)),
            "finite": bool(np.isfinite(array).all()),
        }
    return summaries


def run_coupling_sensitivity_sweep(
    *,
    steps: int = 20,
    seed: int = 42,
) -> dict[str, Any]:
    if steps < 1:
        raise ValueError("steps must be at least 1.")
    parameters = ModelParameters()
    scenarios = []
    baseline_by_dimension = {}

    for dimensions, grid_size in ((2, 9), (3, 5)):
        state = initial_state(grid_size, spatial_dimensions=dimensions)
        dx = 2.0 / (grid_size - 1)
        stable_limit = diffusion_timestep_limit(parameters, dx, dimensions)
        dt = min(0.01, stable_limit * 0.05)

        for coupling_strength in COUPLING_STRENGTHS:
            for agents_per_unit in AGENTS_PER_UNIT:
                population = initialize_agents(
                    cancer_count=40,
                    immune_count=25,
                    fibroblast_count=15,
                    seed=seed,
                    spatial_dimensions=dimensions,
                )
                initial_counts = _population_counts(population)
                times, states, populations = simulate_coupled(
                    state,
                    population,
                    parameters,
                    dx,
                    dt,
                    steps=steps,
                    coupling_strength=coupling_strength,
                    agents_per_unit=agents_per_unit,
                    save_every=max(1, steps // 4),
                    seed=seed,
                    maximum_agents=2000,
                )
                final_state = states[-1]
                final_population = populations[-1]
                final_fields = _field_summary(final_state, parameters)
                counts_by_snapshot = [
                    _population_counts(snapshot) for snapshot in populations
                ]
                baseline_key = f"{dimensions}d"
                if coupling_strength == 0.0:
                    baseline_by_dimension[baseline_key] = {
                        name: value["mean"]
                        for name, value in final_fields.items()
                    }
                baseline = baseline_by_dimension.get(baseline_key)
                delta_from_uncoupled = (
                    None
                    if baseline is None
                    else {
                        name: float(summary["mean"] - baseline[name])
                        for name, summary in final_fields.items()
                    }
                )
                scenarios.append(
                    {
                        "spatial_dimensions": dimensions,
                        "grid_shape": list(final_state.shape),
                        "grid_spacing": float(dx),
                        "time_step": float(dt),
                        "diffusion_stability_limit": float(stable_limit),
                        "steps": int(steps),
                        "final_time": float(times[-1]),
                        "coupling_strength": float(coupling_strength),
                        "agents_per_unit": float(agents_per_unit),
                        "initial_cell_counts": initial_counts,
                        "cell_counts_by_saved_snapshot": counts_by_snapshot,
                        "final_cell_counts": _population_counts(final_population),
                        "final_field_summary": final_fields,
                        "final_mean_field_delta_vs_uncoupled": delta_from_uncoupled,
                    }
                )

    return {
        "analysis": "synthetic_bidirectional_pde_agent_coupling_sensitivity",
        "data": "synthetic initial conditions and model-generated cell agents only",
        "seed": int(seed),
        "model_parameters": "default dimensionless ModelParameters; not biologically calibrated",
        "sweep": {
            "coupling_strengths": list(COUPLING_STRENGTHS),
            "agents_per_unit": list(AGENTS_PER_UNIT),
            "spatial_dimensions": [2, 3],
            "scenario_count": len(scenarios),
        },
        "scenarios": scenarios,
        "interpretation": (
            "This is a controlled software sensitivity sweep. It reports how the "
            "coupled numerical prototype changes as arbitrary dimensionless bridge "
            "parameters vary. It is not a mass-conservation proof, biological "
            "calibration, patient simulation, or clinical result."
        ),
        "limitations": [
            "The agent-to-grid occupancy scale and coupling strength are not "
            "calibrated from tissue observations.",
            "Nearest-grid deposition is used; no alternate deposition kernel is compared.",
            "The PDE-agent blend is not derived from a mass-conserving multiscale limit.",
            "The sweep does not estimate biological parameter uncertainty.",
            "Short synthetic runs do not validate cell kinetics or tumour progression.",
        ],
        "data_governance": {
            "patient_data_used": False,
            "patient_or_sample_identifiers_written": False,
            "aggregate_synthetic_metrics_only": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="results/coupling_sensitivity/synthetic_coupling_sweep.json",
    )
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    report = run_coupling_sensitivity_sweep(steps=args.steps, seed=args.seed)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("Scenarios:", report["sweep"]["scenario_count"])
    print("Dimensions:", report["sweep"]["spatial_dimensions"])
    print("Coupling strengths:", report["sweep"]["coupling_strengths"])
    print("Agents per unit:", report["sweep"]["agents_per_unit"])
    print("Aggregate synthetic report:", output)


if __name__ == "__main__":
    main()
