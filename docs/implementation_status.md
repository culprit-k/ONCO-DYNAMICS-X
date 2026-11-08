# Implementation status

This document distinguishes code that exists from work that remains planned or needs scientific validation.

## Implemented and covered by the base CI

- Dimensionless 2D reaction-diffusion solver with eight spatial fields.
- Five-point no-flux finite-difference Laplacian and explicit Euler integration.
- Diffusion timestep guard, state validation, deterministic synthetic initial condition and CLI demo.
- Numerical timestep-refinement utility and scalar output comparison.
- NumPy stochastic agent-based model with cancer division/death, immune recruitment/chemotaxis, fibroblast motility and local environment coupling.
- Cohort manifest validation and gene-expression matrix validation with train-only variance filtering/scaling.
- Local finite-difference sensitivity and Monte Carlo parameter-range propagation.
- Validated JSON serial parser and normalized knob-to-simulation-parameter mapping.
- ESP32 firmware sketch and a dashboard serial-input hook.
- Streamlit dashboard with computed 2D maps and a 3D surface view of a 2D field.

## Implemented as optional model components

These architectures are present in source but need dedicated data, training pipelines, benchmarking and validation before their outputs can be interpreted:

- Coordinate MLP PINN with automatic-differentiation PDE residual and weighted data/physics loss.
- Spatial graph message-passing network.
- Fourier neural operator for gridded fields.
- Transformer encoder for fixed-order bulk gene-expression vectors.
- PyMC Gaussian-likelihood posterior builder with caller-supplied forward model.
- DoWhy linear backdoor estimation wrapper requiring analyst-declared common causes.
- CasADi mean-field optimal-control demonstration with bounded abstract input.

Optional modules are not automatically installed by the core package. Install the appropriate extras described in pyproject.toml.

## Not yet completed

- Calibration of parameters against breast-cancer datasets and biological measurements.
- Full independent cohort validation and patient/sample linkage across modalities.
- Spatial image/segmentation ingestion with assay-specific QC.
- Full coupled PDE-PINN training workflow and numerical comparison report.
- PhysiCell integration or fully validated agent rules.
- Spatial GNN training with registered spatial ground truth.
- Bayesian posterior estimation on real observations, including identifiability and convergence diagnostics.
- End-to-end multi-omics transformer training and external cohort evaluation.
- Causal graph justification, negative-control analysis and robustness/sensitivity to unmeasured confounding.
- PDE-constrained optimal control; current CasADi example is only a toy mean-field problem.
- A true 3D tissue solver and full interactive 3D tumour microenvironment.
- Physical breast phantom fabrication and hardware bench integration.
- A deployed hosted dashboard; deployment instructions are provided, but the owner must authorize and configure hosting.

## Evidence policy

Synthetic simulation behavior is software evidence only. It is not biological validation. No clinical utility, patient-specific prediction, treatment effect or therapeutic recommendation is established by the current repository.
