# AI engine

Optional research models now include:
- `pinn.py`: coordinate-to-field MLP and automatic-differentiation PDE residual utilities.
- `spatial_gnn.py`: graph message passing for segmented cell or region nodes.
- `neural_operator.py`: Fourier neural operator for gridded multi-field input.
- `multiomics_transformer.py`: fixed-gene-order expression transformer for bulk expression vectors.

Install PyTorch separately with `pip install -e ".[pinn,spatial]"` or the relevant optional environment. The expression transformer and neural operator also use PyTorch. These are model components, not trained models. Supply compatible data and compare every model against transparent baselines. PINNs and neural operators must be benchmarked against numerical PDE solutions; a bulk expression transformer does not infer spatial cell locations.
