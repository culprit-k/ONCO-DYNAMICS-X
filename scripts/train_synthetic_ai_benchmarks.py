"""Train PINN and FNO architectures on analytic synthetic heat-equation tasks.

This measures implementation/training behaviour only; no breast-cancer cohort,
patient, or clinical data are used.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch import Tensor

from onco_dynamics_x.ai_engine.neural_operator import FourierNeuralOperator2D
from onco_dynamics_x.ai_engine.pinn import FieldPINN, pde_residual
from onco_dynamics_x.ai_engine.spatial_gnn import SpatialCellGNN
from onco_dynamics_x.ai_engine.synthetic_benchmarks import (
    exact_neumann_heat_solution,
    periodic_heat_operator,
    synthetic_spatial_neighbour_graph,
)


def _relative_l2(prediction: Tensor, target: Tensor) -> float:
    numerator = torch.linalg.vector_norm(prediction - target)
    denominator = torch.linalg.vector_norm(target).clamp_min(1e-12)
    return float((numerator / denominator).detach().cpu())


def _coordinates(
    count: int,
    *,
    seed: int,
    device: torch.device,
) -> Tensor:
    generator = torch.Generator(device=device)
    generator.manual_seed(seed)
    values = torch.rand((count, 3), generator=generator, device=device)
    values[:, 0:2] = values[:, 0:2] * 2.0 - 1.0
    values[:, 2] = values[:, 2] * 0.25
    return values


def train_pinn(
    *,
    epochs: int,
    seed: int,
    device: torch.device,
) -> dict[str, Any]:
    diffusion = 0.02
    torch.manual_seed(seed)
    train_coordinates = _coordinates(1024, seed=seed + 1, device=device)
    train_targets = exact_neumann_heat_solution(
        train_coordinates, diffusion=diffusion
    )
    validation_coordinates = _coordinates(512, seed=seed + 2, device=device)
    validation_targets = exact_neumann_heat_solution(
        validation_coordinates, diffusion=diffusion
    )
    test_coordinates = _coordinates(1024, seed=seed + 3, device=device)
    test_targets = exact_neumann_heat_solution(test_coordinates, diffusion=diffusion)

    model = FieldPINN(output_fields=1, hidden_width=48, hidden_layers=3).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    def zero_reaction(prediction: Tensor) -> Tensor:
        return torch.zeros_like(prediction)
    best_validation = float("inf")
    best_state = None
    history = []

    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        data_prediction = model(train_coordinates)
        collocation = _coordinates(192, seed=seed + 10_000 + epoch, device=device)
        residual = pde_residual(
            model,
            collocation,
            diffusion_coefficients=(diffusion,),
            reaction_function=zero_reaction,
        )
        data_loss = torch.mean((data_prediction - train_targets) ** 2)
        physics_loss = torch.mean(residual**2)
        loss = data_loss + 0.1 * physics_loss
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        if epoch == 1 or epoch % 20 == 0 or epoch == epochs:
            model.eval()
            with torch.no_grad():
                validation_prediction = model(validation_coordinates)
                validation_loss = float(
                    torch.mean((validation_prediction - validation_targets) ** 2)
                    .detach()
                    .cpu()
                )
            history.append(
                {
                    "epoch": epoch,
                    "training_data_mse": float(data_loss.detach().cpu()),
                    "training_physics_mse": float(physics_loss.detach().cpu()),
                    "validation_mse": validation_loss,
                }
            )
            if validation_loss < best_validation:
                best_validation = validation_loss
                best_state = {
                    name: value.detach().cpu().clone()
                    for name, value in model.state_dict().items()
                }

    if best_state is None:
        raise RuntimeError("PINN training did not produce a checkpoint.")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        test_prediction = model(test_coordinates)
    test_error = _relative_l2(test_prediction, test_targets)

    diagnostic_coordinates = _coordinates(512, seed=seed + 4, device=device)
    diagnostic_coordinates.requires_grad_(True)
    diagnostic_residual = pde_residual(
        model,
        diagnostic_coordinates,
        diffusion_coefficients=(diffusion,),
        reaction_function=zero_reaction,
    )
    residual_rms = float(torch.sqrt(torch.mean(diagnostic_residual**2)).detach().cpu())
    initial_model = FieldPINN(
        output_fields=1, hidden_width=48, hidden_layers=3
    ).to(device)
    with torch.no_grad():
        initial_error = _relative_l2(initial_model(test_coordinates), test_targets)

    return {
        "task": "2D Neumann heat equation with known analytic solution",
        "architecture": "FieldPINN",
        "training_epochs_requested": int(epochs),
        "best_epoch": int(
            min(history, key=lambda row: row["validation_mse"])["epoch"]
        ),
        "train_samples": int(len(train_coordinates)),
        "validation_samples": int(len(validation_coordinates)),
        "test_samples": int(len(test_coordinates)),
        "diffusion_coefficient": diffusion,
        "time_domain": [0.0, 0.25],
        "relative_l2_test_error": test_error,
        "untrained_initialization_relative_l2_error": initial_error,
        "held_out_pde_residual_rms": residual_rms,
        "best_validation_mse": best_validation,
        "training_history_checkpoints": history,
        "model_configuration": {
            "hidden_width": 48,
            "hidden_layers": 3,
            "optimizer": "Adam",
            "learning_rate": 0.001,
            "loss": "data MSE + 0.1 * PDE residual MSE",
            "gradient_clip_norm": 5.0,
        },
    }


def _fno_dataset(
    count: int,
    *,
    size: int,
    seed: int,
    device: torch.device,
    diffusion: float,
    duration: float,
) -> tuple[Tensor, Tensor]:
    generator = torch.Generator(device=device)
    generator.manual_seed(seed)
    initial = torch.randn((count, size, size), generator=generator, device=device)
    initial = initial / torch.std(initial, dim=(1, 2), keepdim=True).clamp_min(1e-6)
    target = periodic_heat_operator(
        initial, diffusion=diffusion, duration=duration
    )
    return initial.unsqueeze(-1), target.unsqueeze(-1)


def train_fno(
    *,
    epochs: int,
    seed: int,
    device: torch.device,
) -> dict[str, Any]:
    size = 16
    diffusion = 0.05
    duration = 0.2
    train_x, train_y = _fno_dataset(
        128, size=size, seed=seed + 101, device=device,
        diffusion=diffusion, duration=duration,
    )
    validation_x, validation_y = _fno_dataset(
        32, size=size, seed=seed + 102, device=device,
        diffusion=diffusion, duration=duration,
    )
    test_x, test_y = _fno_dataset(
        64, size=size, seed=seed + 103, device=device,
        diffusion=diffusion, duration=duration,
    )

    model = FourierNeuralOperator2D(
        input_channels=1,
        output_channels=1,
        width=16,
        modes_x=6,
        modes_y=6,
        layers=2,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    best_validation = float("inf")
    best_state = None
    history = []
    batch_size = 16

    for epoch in range(1, epochs + 1):
        model.train()
        permutation = torch.randperm(len(train_x), device=device)
        losses = []
        for start in range(0, len(train_x), batch_size):
            indices = permutation[start : start + batch_size]
            optimizer.zero_grad(set_to_none=True)
            prediction = model(train_x[indices])
            loss = torch.mean((prediction - train_y[indices]) ** 2)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            losses.append(float(loss.detach().cpu()))

        if epoch == 1 or epoch % 10 == 0 or epoch == epochs:
            model.eval()
            with torch.no_grad():
                validation_prediction = model(validation_x)
                validation_loss = float(
                    torch.mean((validation_prediction - validation_y) ** 2)
                    .detach()
                    .cpu()
                )
            history.append(
                {
                    "epoch": epoch,
                    "training_mse": float(np.mean(losses)),
                    "validation_mse": validation_loss,
                }
            )
            if validation_loss < best_validation:
                best_validation = validation_loss
                best_state = {
                    name: value.detach().cpu().clone()
                    for name, value in model.state_dict().items()
                }

    if best_state is None:
        raise RuntimeError("FNO training did not produce a checkpoint.")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        test_prediction = model(test_x)
    test_error = _relative_l2(test_prediction, test_y)
    persistence_error = _relative_l2(test_x, test_y)

    return {
        "task": "2D periodic heat semigroup with exact Fourier solution",
        "architecture": "FourierNeuralOperator2D",
        "training_epochs_requested": int(epochs),
        "best_epoch": int(
            min(history, key=lambda row: row["validation_mse"])["epoch"]
        ),
        "train_samples": int(len(train_x)),
        "validation_samples": int(len(validation_x)),
        "test_samples": int(len(test_x)),
        "grid_size": size,
        "diffusion_coefficient": diffusion,
        "diffusion_duration": duration,
        "relative_l2_test_error": test_error,
        "persistence_baseline_relative_l2_error": persistence_error,
        "best_validation_mse": best_validation,
        "training_history_checkpoints": history,
        "model_configuration": {
            "width": 16,
            "modes_x": 6,
            "modes_y": 6,
            "layers": 2,
            "optimizer": "Adam",
            "learning_rate": 0.001,
            "batch_size": batch_size,
            "loss": "mean squared error",
            "gradient_clip_norm": 5.0,
        },
    }



def _batched_spatial_graphs(
    graph_count: int,
    *,
    seed: int,
    device: torch.device,
) -> tuple[Tensor, Tensor, Tensor]:
    node_features = []
    edge_indices = []
    targets = []
    offset = 0
    for graph_index in range(graph_count):
        features, edges, labels = synthetic_spatial_neighbour_graph(
            num_nodes=64,
            radius=0.4,
            seed=seed + graph_index,
        )
        node_features.append(features)
        edge_indices.append(edges + offset)
        targets.append(labels)
        offset += len(features)
    return (
        torch.cat(node_features, dim=0).to(device),
        torch.cat(edge_indices, dim=1).to(device),
        torch.cat(targets, dim=0).to(device),
    )


def train_gnn(
    *,
    epochs: int,
    seed: int,
    device: torch.device,
) -> dict[str, Any]:
    """Train node-level GNN to recover a synthetic local-neighbour statistic."""
    torch.manual_seed(seed + 201)
    train_x, train_edges, train_y = _batched_spatial_graphs(
        48, seed=seed + 210, device=device
    )
    validation_x, validation_edges, validation_y = _batched_spatial_graphs(
        12, seed=seed + 310, device=device
    )
    test_x, test_edges, test_y = _batched_spatial_graphs(
        12, seed=seed + 410, device=device
    )

    model = SpatialCellGNN(
        input_features=5,
        hidden_features=32,
        output_features=1,
        layers=2,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=2e-3)
    best_validation = float("inf")
    best_state = None
    history = []

    for epoch in range(1, epochs + 1):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        prediction = model(train_x, train_edges)
        loss = torch.mean((prediction - train_y) ** 2)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        if epoch == 1 or epoch % 10 == 0 or epoch == epochs:
            model.eval()
            with torch.no_grad():
                validation_prediction = model(validation_x, validation_edges)
                validation_mse = float(
                    torch.mean((validation_prediction - validation_y) ** 2)
                    .detach()
                    .cpu()
                )
            history.append(
                {
                    "epoch": epoch,
                    "training_mse": float(loss.detach().cpu()),
                    "validation_mse": validation_mse,
                }
            )
            if validation_mse < best_validation:
                best_validation = validation_mse
                best_state = {
                    name: value.detach().cpu().clone()
                    for name, value in model.state_dict().items()
                }

    if best_state is None:
        raise RuntimeError("GNN training did not produce a checkpoint.")
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        test_prediction = model(test_x, test_edges)
    test_mse = float(torch.mean((test_prediction - test_y) ** 2).cpu())
    test_relative_l2 = _relative_l2(test_prediction, test_y)
    constant_baseline = torch.full_like(test_y, float(train_y.mean().cpu()))
    baseline_mse = float(torch.mean((constant_baseline - test_y) ** 2).cpu())

    return {
        "task": "node-level prediction of local immune-neighbour fraction",
        "architecture": "SpatialCellGNN",
        "training_epochs_requested": int(epochs),
        "best_epoch": int(
            min(history, key=lambda row: row["validation_mse"])["epoch"]
        ),
        "train_graphs": 48,
        "validation_graphs": 12,
        "test_graphs": 12,
        "nodes_per_graph": 64,
        "train_nodes": int(len(train_x)),
        "validation_nodes": int(len(validation_x)),
        "test_nodes": int(len(test_x)),
        "relative_l2_test_error": test_relative_l2,
        "test_mse": test_mse,
        "constant_mean_baseline_mse": baseline_mse,
        "best_validation_mse": best_validation,
        "training_history_checkpoints": history,
        "model_configuration": {
            "input_features": ["x", "y", "cancer_one_hot", "immune_one_hot", "fibroblast_one_hot"],
            "hidden_features": 32,
            "message_passing_layers": 2,
            "output": "local immune-type fraction among radius-0.4 neighbours",
            "optimizer": "Adam",
            "learning_rate": 0.002,
            "loss": "mean squared error",
        },
    }

def run_benchmarks(
    *,
    pinn_epochs: int = 500,
    fno_epochs: int = 120,
    gnn_epochs: int = 80,
    seed: int = 42,
) -> dict[str, Any]:
    if min(pinn_epochs, fno_epochs, gnn_epochs) < 1:
        raise ValueError("Training epochs must be positive.")
    torch.set_num_threads(1)
    device = torch.device("cpu")
    pinn = train_pinn(epochs=pinn_epochs, seed=seed, device=device)
    fno = train_fno(epochs=fno_epochs, seed=seed, device=device)
    gnn = train_gnn(epochs=gnn_epochs, seed=seed, device=device)
    return {
        "analysis": "synthetic_ai_architecture_training_benchmark",
        "data": "analytic synthetic heat-equation and graph data only; no patients or cohorts",
        "device": str(device),
        "seed": int(seed),
        "pinn": pinn,
        "fno": fno,
        "gnn": gnn,
        "interpretation": (
            "Demonstrates that PINN, FNO and spatial GNN implementations can be "
            "trained and evaluated on separate analytic/synthetic mathematical tasks. "
            "Metrics are not breast-cancer accuracy and provide no clinical or biological "
            "validation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="results/ai_benchmarks/synthetic_pinn_fno_training.json",
    )
    parser.add_argument("--pinn-epochs", type=int, default=500)
    parser.add_argument("--fno-epochs", type=int, default=120)
    parser.add_argument("--gnn-epochs", type=int, default=80)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    report = run_benchmarks(
        pinn_epochs=args.pinn_epochs,
        fno_epochs=args.fno_epochs,
        gnn_epochs=args.gnn_epochs,
        seed=args.seed,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("PINN held-out relative L2:", report["pinn"]["relative_l2_test_error"])
    print("PINN held-out PDE residual RMS:", report["pinn"]["held_out_pde_residual_rms"])
    print("FNO held-out relative L2:", report["fno"]["relative_l2_test_error"])
    print("FNO persistence baseline:", report["fno"]["persistence_baseline_relative_l2_error"])
    print("GNN held-out relative L2:", report["gnn"]["relative_l2_test_error"])
    print("GNN test MSE:", report["gnn"]["test_mse"])
    print("GNN constant baseline MSE:", report["gnn"]["constant_mean_baseline_mse"])
    print("Aggregate synthetic AI benchmark:", output)


if __name__ == "__main__":
    main()
