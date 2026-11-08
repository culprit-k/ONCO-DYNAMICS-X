"""Optional transformer encoder for consistently ordered bulk-expression vectors."""
from __future__ import annotations

import torch
from torch import Tensor, nn


class ExpressionTransformer(nn.Module):
    """Encode fixed-order gene-expression vectors for downstream research tasks.

    Input gene order must be identical between training and inference. This model
    is not a spatial transcriptomics model and cannot infer cell coordinates.
    """

    def __init__(
        self,
        number_of_genes: int,
        output_features: int,
        embedding_dim: int = 64,
        attention_heads: int = 4,
        transformer_layers: int = 2,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if min(number_of_genes, output_features, embedding_dim, attention_heads) < 1:
            raise ValueError("Model dimensions must be positive.")
        if embedding_dim % attention_heads:
            raise ValueError("embedding_dim must be divisible by attention_heads.")
        if transformer_layers < 1:
            raise ValueError("At least one transformer layer is required.")
        self.number_of_genes = number_of_genes
        self.expression_projection = nn.Linear(1, embedding_dim)
        self.gene_embedding = nn.Embedding(number_of_genes, embedding_dim)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embedding_dim,
            nhead=attention_heads,
            dim_feedforward=embedding_dim * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=transformer_layers)
        self.norm = nn.LayerNorm(embedding_dim)
        self.head = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(embedding_dim, output_features),
        )

    def forward(self, expression: Tensor) -> Tensor:
        if expression.ndim != 2 or expression.shape[1] != self.number_of_genes:
            raise ValueError(
                f"expression must have shape [batch, {self.number_of_genes}]."
            )
        values = expression.unsqueeze(-1)
        gene_ids = torch.arange(self.number_of_genes, device=expression.device)
        tokens = self.expression_projection(values) + self.gene_embedding(gene_ids)
        encoded = self.encoder(tokens)
        pooled = self.norm(encoded.mean(dim=1))
        return self.head(pooled)
