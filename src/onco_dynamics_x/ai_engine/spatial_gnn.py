"""Optional PyTorch message-passing network for spatial cell graphs."""
from __future__ import annotations

import torch
from torch import Tensor, nn


class MessagePassingLayer(nn.Module):
    """Neighbour aggregation using directed edges and a residual update."""

    def __init__(self, input_features: int, output_features: int) -> None:
        super().__init__()
        if input_features < 1 or output_features < 1:
            raise ValueError("Feature dimensions must be positive.")
        self.message = nn.Sequential(
            nn.Linear(input_features * 2, output_features),
            nn.ReLU(),
            nn.Linear(output_features, output_features),
        )
        self.update = nn.Sequential(
            nn.Linear(input_features + output_features, output_features),
            nn.ReLU(),
            nn.Linear(output_features, output_features),
        )

    def forward(self, node_features: Tensor, edge_index: Tensor) -> Tensor:
        if node_features.ndim != 2:
            raise ValueError("node_features must have shape [nodes, features].")
        if edge_index.ndim != 2 or edge_index.shape[0] != 2:
            raise ValueError("edge_index must have shape [2, edges].")
        source, target = edge_index.long()
        if source.numel() and (
            int(source.max()) >= len(node_features) or int(target.max()) >= len(node_features)
        ):
            raise ValueError("edge_index contains a node index outside the feature array.")
        messages = self.message(
            torch.cat([node_features[source], node_features[target]], dim=1)
        )
        aggregated = torch.zeros(
            (len(node_features), messages.shape[1]),
            dtype=messages.dtype,
            device=messages.device,
        )
        aggregated.index_add_(0, target, messages)
        degree = torch.zeros(
            (len(node_features), 1), dtype=messages.dtype, device=messages.device
        )
        degree.index_add_(
            0, target, torch.ones((len(target), 1), dtype=messages.dtype, device=messages.device)
        )
        aggregated = aggregated / degree.clamp_min(1.0)
        return self.update(torch.cat([node_features, aggregated], dim=1))


class SpatialCellGNN(nn.Module):
    """Node-level model; nodes may represent segmented cells or cell regions."""

    def __init__(
        self,
        input_features: int,
        hidden_features: int = 64,
        output_features: int = 8,
        layers: int = 3,
    ) -> None:
        super().__init__()
        if layers < 1:
            raise ValueError("At least one message-passing layer is required.")
        modules = []
        current_features = input_features
        for _ in range(layers):
            modules.append(MessagePassingLayer(current_features, hidden_features))
            current_features = hidden_features
        self.layers = nn.ModuleList(modules)
        self.head = nn.Linear(hidden_features, output_features)

    def forward(self, node_features: Tensor, edge_index: Tensor) -> Tensor:
        hidden = node_features
        for layer in self.layers:
            hidden = layer(hidden, edge_index)
        return self.head(hidden)
