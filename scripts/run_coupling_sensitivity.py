"""Run a multi-seed deposition and coupling sensitivity experiment.

The experiment quantifies numerical sensitivity to dimensionless coupling
controls and two grid-deposition methods. It is not parameter calibration and
uses no patient data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

import onco_dynamics_x.cellular_simulation.agent_based as agent_based
import onco_dynamics_x.cellular_simulation.coupling as coupling
from onco_dynamics_x.parameters import ModelParameters
from onco_dynamics_x.solver import diffusion_timestep_limit
from onco_dynamics_x.state import initial_state


COUPLING_STRENGTHS = (0.0, 0.25, 0.5, 0.75, 1.0)
AGENTS_PER_UNIT = (1.0, 2.0, 4.0, 8.0)
DEPOSITION_METHODS = ("nearest_grid", "cloud_in_cell")
SEEDS = (42, 123, 2026)


def _population_counts(population) -> dict[str, int]:
    return {
        "cancer": population.count(agent_based.CANCER),
        "immune": population.count(agent_based.IMMUNE),
        "fibroblast": population.count(agent_based.FIBROBLAST),
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
        at_lower_bound = array <= 1e-10
        at_upper_bound = array >= upper - 1e-10
        at_either_bound = at_lower_bound | at_upper_bound
        summaries[name] = {
            "mean": float(np.mean(array)),
            "standard_deviation": float(np.std(array)),
            "minimum": float(np.min(array)),
            "maximum": float(np.max(array)),
            "fraction_at_lower_bound": float(np.mean(at_lower_bound)),
            "fraction_at_upper_bound": float(np.mean(at_upper_bound)),
            "fraction_at_or_near_bounds": float(np.mean(at_either_bound)),
            "finite": bool(np.isfinite(array).all()),
        }
    return summaries


def run_coupling_sensitivity_sweep(
    *,
    steps: int = 20,
    seeds: tuple[int, ...] = SEEDS,
) -> dict[str, Any]:
    if steps < 1:
        raise ValueError("steps must be at least 1.")
    if not seeds or len(set(seeds)) != len(seeds):
        raise ValueError("seeds must be non-empty and unique.")
    parameters = ModelParameters()
    scenarios = []
    baseline_by_condition = {}

    for dimensions, grid_size in ((2, 9), (3, 5)):
        state = initial_state(grid_size, spatial_dimensions=dimensions)
        dx = 2.0 / (grid_size - 1)
        stable_limit = diffusion_timestep_limit(parameters, dx, dimensions)
        dt = min(0.01, stable_limit * 0.05)

        for deposition_method in DEPOSITION_METHODS:
            for seed in seeds:
                for coupling_strength in COUPLING_STRENGTHS:
                    for agents_per_unit in AGENTS_PER_UNIT:
                        population = agent_based.initialize_agents(
                            cancer_count=40,
                            immune_count=25,
                            fibroblast_count=15,
                            seed=seed,
                            spatial_dimensions=dimensions,
                        )
                        initial_counts = _population_counts(population)
                        times, states, populations = coupling.simulate_coupled(
                            state,
                            population,
                            parameters,
                            dx,
                            dt,
                            steps=steps,
                            coupling_strength=coupling_strength,
                            agents_per_unit=agents_per_unit,
                            deposition_method=deposition_method,
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
                        baseline_key = (
                            dimensions,
                            deposition_method,
                            seed,
                            agents_per_unit,
                        )
                        if coupling_strength == 0.0:
                            baseline_by_condition[baseline_key] = {
                                name: value["mean"]
                                for name, value in final_fields.items()
                            }
                        baseline = baseline_by_condition.get(baseline_key)
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
                                "seed": int(seed),
                                "deposition_method": deposition_method,
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
        "analysis": "synthetic_bidirectional_pde_agent_deposition_and_coupling_sensitivity",
        "data": "synthetic initial conditions and model-generated cell agents only",
        "seeds": [int(seed) for seed in seeds],
        "model_parameters": "default dimensionless ModelParameters; not biologically calibrated",
        "sweep": {
            "coupling_strengths": list(COUPLING_STRENGTHS),
            "agents_per_unit": list(AGENTS_PER_UNIT),
            "deposition_methods": list(DEPOSITION_METHODS),
            "spatial_dimensions": [2, 3],
            "seed_count": len(seeds),
            "scenario_count": len(scenarios),
        },
        "scenarios": scenarios,
        "interpretation": (
            "This is a controlled multi-seed software sensitivity experiment. It "
            "compares nearest-grid and cloud-in-cell deposition as arbitrary "
            "dimensionless bridge parameters vary. It is not a mass-conservation "
            "proof for the full coupled dynamics, biological calibration, patient "
            "simulation, or clinical result."
        ),
        "limitations": [
            "The agent-to-grid occupancy scale and coupling strength are not calibrated from tissue observations.",
            "Cloud-in-cell conserves deposited agent count before field clipping, but the PDE-agent blend and agent kinetics are not a mass-conserving multiscale derivation.",
            "Only three random seeds and short synthetic trajectories are tested.",
            "Only small 2D and 3D grids are used.",
            "No histology or spatial-omics data were used to calibrate the bridge.",
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
    parser.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    args = parser.parse_args()
    report = run_coupling_sensitivity_sweep(
        steps=args.steps,
        seeds=tuple(args.seeds),
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("Scenarios:", report["sweep"]["scenario_count"])
    print("Dimensions:", report["sweep"]["spatial_dimensions"])
    print("Deposition methods:", report["sweep"]["deposition_methods"])
    print("Coupling strengths:", report["sweep"]["coupling_strengths"])
    print("Agents per unit:", report["sweep"]["agents_per_unit"])
    print("Seeds:", report["seeds"])
    print("Aggregate synthetic report:", output)


if __name__ == "__main__":
    main()
