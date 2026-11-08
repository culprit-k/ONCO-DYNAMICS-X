# Spatial graph and leakage-safe cohort splits

## Cell graph construction
`cellular_simulation/graph.py` creates a radius-based directed cell graph from the 2D or 3D ABM coordinates. Each node feature vector contains spatial coordinates, one-hot cell phenotype, dimensionless agent age, and the eight local PDE fields. Edges represent geometric proximity under an analyst-chosen radius and neighbour cap. They are not direct evidence of molecular interaction or causality.

The graph builder returns NumPy arrays; convert them to PyTorch tensors before passing them to `SpatialCellGNN`. The current radius search is chunked but still quadratic in the number of cells. For large tissue images, replace it with a spatial index or radius-neighbour algorithm and validate the resulting graph.

## Patient-grouped splitting
`data_engine/splits.py` provides a reproducible group train/test split. All samples from a single patient remain in exactly one partition. It is not stratified; report patient counts and event balance, and use group-aware stratification when needed. Do not use sample-random splits for multi-sample patient datasets.
