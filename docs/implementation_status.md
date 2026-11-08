# Implementation status

This document separates implemented software from work that still requires scientific validation or physical fabrication.

## Implemented and covered by base CI
- Dimensionless 2D and 3D reaction-diffusion solver with eight spatial fields, no-flux finite differences and dimension-aware explicit diffusion stability guard.
- Synthetic initial fields, state validation, 2D/3D CLI, timestep and nested-grid refinement diagnostics, and Streamlit 2D/3D visualization.
- NumPy 2D/3D agent-based model with cancer-agent division/death, local oxygen/nutrient dependence, immune recruitment and gradient-directed migration, fibroblast movement, reflecting boundaries and reproducible seeds.
- Spatial cell-graph builder using cell coordinates, phenotype, age and local PDE fields; patient-grouped train/test split utility.
- Public-only GDC metadata query, sample-manifest export, audit, selected open-file download with checksum verification, and strict TCGA expression-matrix assembly.
- Identifier-aware expression preprocessing, local parameter sensitivity and Monte Carlo parameter-range uncertainty.
- Validated JSON serial parser, normalized knob mapping, ESP32 firmware sketch and dashboard serial-input hook.
- GitHub Actions tests on Python 3.11 and 3.12.

## Implemented as optional model components
- PyTorch PINN with automatic-differentiation PDE residual and weighted observation/physics loss.
- Spatial graph message-passing network, Fourier neural operator and fixed-order expression transformer.
- PyMC Gaussian-likelihood posterior builder, DoWhy backdoor wrapper and CasADi mean-field control demonstration.

Optional modules need compatible data, training and validation; they are not trained models by default.

## Not yet completed
- Biological parameter calibration, identifiability and independent cohort validation.
- End-to-end PINN/FNO training and benchmarking against numerical baselines.
- Spatial GNN training with registered ground truth.
- PhysiCell integration and fully coupled PDE-ABM feedback.
- Real-observation Bayesian inference and defensible causal analysis.
- PDE-constrained optimal control and multi-omics external validation.
- Physical phantom printing/casting, ESP32 bench testing and integrated exhibit assembly.
- Hosted dashboard deployment; the repository owner must configure hosting.

## Evidence policy
Synthetic results demonstrate software behaviour only. No clinical utility, patient-specific prediction, treatment effect or therapeutic recommendation is established by this project.
