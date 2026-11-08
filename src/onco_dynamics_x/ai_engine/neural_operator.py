"""Optional Fourier neural operator for gridded multi-field simulations."""
from __future__ import annotations

import torch
from torch import Tensor, nn


class SpectralConv2d(nn.Module):
    """Fourier convolution over a limited set of low-frequency modes."""

    def __init__(self, in_channels: int, out_channels: int, modes_x: int, modes_y: int):
        super().__init__()
        if min(in_channels, out_channels, modes_x, modes_y) < 1:
            raise ValueError("Channel and Fourier mode counts must be positive.")
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes_x = modes_x
        self.modes_y = modes_y
        scale = 1.0 / (in_channels * out_channels)
        shape = (in_channels, out_channels, modes_x, modes_y)
        self.weight_positive = nn.Parameter(
            scale * torch.randn(*shape, dtype=torch.cfloat)
        )
        self.weight_negative = nn.Parameter(
            scale * torch.randn(*shape, dtype=torch.cfloat)
        )

    @staticmethod
    def _multiply(input_modes: Tensor, weights: Tensor) -> Tensor:
        return torch.einsum("bixy,ioxy->boxy", input_modes, weights)

    def forward(self, values: Tensor) -> Tensor:
        if values.ndim != 4:
            raise ValueError("SpectralConv2d input must have shape [batch, channels, x, y].")
        batch, _, size_x, size_y = values.shape
        modes_x = min(self.modes_x, size_x // 2)
        modes_y = min(self.modes_y, size_y // 2 + 1)
        transformed = torch.fft.rfft2(values)
        output = torch.zeros(
            batch,
            self.out_channels,
            size_x,
            size_y // 2 + 1,
            dtype=torch.cfloat,
            device=values.device,
        )
        output[:, :, :modes_x, :modes_y] = self._multiply(
            transformed[:, :, :modes_x, :modes_y],
            self.weight_positive[:, :, :modes_x, :modes_y],
        )
        output[:, :, -modes_x:, :modes_y] = self._multiply(
            transformed[:, :, -modes_x:, :modes_y],
            self.weight_negative[:, :, :modes_x, :modes_y],
        )
        return torch.fft.irfft2(output, s=(size_x, size_y))


class FourierNeuralOperator2D(nn.Module):
    """Map gridded input fields to gridded output fields."""

    def __init__(
        self,
        input_channels: int = 8,
        output_channels: int = 8,
        width: int = 32,
        modes_x: int = 12,
        modes_y: int = 12,
        layers: int = 4,
    ) -> None:
        super().__init__()
        if layers < 1:
            raise ValueError("At least one Fourier layer is required.")
        self.lift = nn.Conv2d(input_channels, width, kernel_size=1)
        self.spectral_layers = nn.ModuleList(
            [SpectralConv2d(width, width, modes_x, modes_y) for _ in range(layers)]
        )
        self.local_layers = nn.ModuleList(
            [nn.Conv2d(width, width, kernel_size=1) for _ in range(layers)]
        )
        self.project = nn.Sequential(
            nn.Conv2d(width, width, kernel_size=1),
            nn.GELU(),
            nn.Conv2d(width, output_channels, kernel_size=1),
        )

    def forward(self, values: Tensor) -> Tensor:
        if values.ndim != 4:
            raise ValueError("Input must have shape [batch, x, y, channels].")
        hidden = self.lift(values.permute(0, 3, 1, 2))
        for spectral, local in zip(
            self.spectral_layers, self.local_layers, strict=True
        ):
            hidden = torch.nn.functional.gelu(spectral(hidden) + local(hidden))
        output = self.project(hidden)
        return output.permute(0, 2, 3, 1)
