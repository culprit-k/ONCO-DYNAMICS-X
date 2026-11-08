# Implementation status

This document separates code that exists from work that remains planned or requires scientific validation.

## Implemented and covered by base CI
- Dimensionless 2D and 3D reaction-diffusion solver with eight spatial fields, no-flux finite differences and dimension-aware explicit diffusion stability guard.
- Synthetic initial fields, state validation, CLI simulation in 2D or 3D, timestep refinement, and Streamlit 2D/3D visualization.
- NumPy agent-based model with cancer-agent division/death, local oxygen/nutrient dependence, immune recruitment and gradient-directed migration, and fibroblast migration (2D environment).
- Public-only GDC API file-metadata query, sample-manifest export, manifest audit, and an individual open-file download method with checksum/size verification.
- Identifier-aware expression matrix validation and train-only variance filtering/scaling.
- Local finite-difference sensitivity and Monte Carlo propagation of user-specified parameter ranges.
- Validated JSON serial parser, normalized knob mapping, ESP32 firmware sketch and dashboard serial-input hook.
- GitHub Actions testing on Python 3.11 and 3.12.

## Implemented as optional model components
These architectures are present in source but need compatible data, training pipelines, benchmarking and validation before their outputs can be interpreted:
- Coordinate MLP PINN with automatic-differentiation PDE residual and weighted observation/physics loss.
- Spatial graph message-passing network.
- Fourier neural operator for gridded multi-field input.
- Transformer encoder for fixed-order bulk gene-expression vectors.
- PyMC Gaussian-likelihood posterior builder with caller-supplied forward model.
- DoWhy linear backdoor-estimation wrapper requiring analyst-declared common causes.
- CasADi mean-field optimal-control demonstration with bounded abstract input.

Optional modules are not installed by the core package. Install the appropriate extras in pyproject.toml.

## Not yet completed
- Calibration of parameters against breast-cancer datasets and biological measurements.
- Independent cohort validation and patient/sample linkage across modalities.
- End-to-end PINN/FNO training and numerical benchmark report.
- PhysiCell integration or externally validated agent rules.
- Spatial GNN training with registered spatial ground truth.
- Bayesian posterior estimation on real observations, including identifiability and convergence diagnostics.
- End-to-end multi-omics training and external cohort evaluation.
- Causal graph justification, negative-control analyses and unmeasured-confounding robustness.
- PDE-constrained optimal control; current CasADi example is a toy mean-field problem.
- Full 3D cell-agent model, clonal evolution, matrix mechanics and vascular topology.
- Physical breast phantom fabrication and hardware bench integration.
- Hosted dashboard deployment; instructions are provided, but repository owner must configure hosting.

## Evidence policy
Synthetic simulation behavior is software evidence only, not biological validation. No clinical utility, patient-specific prediction, treatment effect or therapeutic recommendation is established by this repository.
