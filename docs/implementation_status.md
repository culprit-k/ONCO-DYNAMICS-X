# Implementation status

This status separates executable code from model components that still require training, scientific validation or physical integration.

## Implemented in the repository

### Mathematical simulation
- Eight dimensionless spatial fields: cancer, immune, fibroblast, oxygen, nutrient, drug, immune-evasion proxy and vasculature proxy.
- Finite-difference reaction-diffusion solver on 2D and 3D Cartesian grids, no-flux edges and explicit Euler time stepping.
- Dimension-aware explicit diffusion stability guard.
- Synthetic initial states, CLI 2D/3D simulations and Streamlit 2D maps / 3D voxel-volume visualization.
- Timestep refinement and nested-grid refinement diagnostics for 2D and 3D.

### Cellular and spatial AI inputs
- Stochastic 2D/3D agent-based simulation for cancer, immune and fibroblast cell proxies.
- Radius-based cell-neighbour graph construction and a GNN architecture.
- Patient-grouped train/test split utility.

### Data engineering and first real-cohort audit
- GDC public/open metadata query, sample-manifest export and audit.
- Sample-type-specific filtering, including Primary Tumor and Solid Tissue Normal.
- Selected open-file download with retry handling, byte-count verification and MD5 checksum verification.
- Strict expression matrix assembly with sample/gene alignment.
- Identifier-aware expression preprocessing.
- Paired TCGA-BRCA selection: one primary-tumour and one normal sample per shared patient.
- Reproducible tumour-versus-normal expression baseline using log1p(TPM), train-fold-only feature filtering/selection/scaling, patient-grouped holdout and StratifiedGroupKFold evaluation.
- Workflow artifacts contain selected metadata and metrics only; raw expression files and expression matrices are not uploaded.

### Statistical and model components
- Local finite-difference parameter sensitivity and Monte Carlo propagation of user-supplied parameter ranges.
- Optional PyTorch PINN, Fourier neural operator, spatial GNN and expression-transformer architectures.
- Optional PyMC posterior builder, DoWhy linear backdoor wrapper and CasADi mean-field control demonstration.

### Hardware and presentation
- ESP32 firmware, normalized serial protocol and dashboard control hook.
- Parametric OpenSCAD breast-phantom shell, removable lesion inserts, wiring map and physical build guide.

## Implemented as architecture only or optional components

PINN, GNN, neural operator, expression transformer, Bayesian, causal and optimal-control modules are code components. They are not trained or biologically validated models by simply being present in source.

## Still outstanding
- Confirm successful completion of the paired TCGA-BRCA download/matrix/model workflow and inspect the final metrics artifact.
- Increase cohort size and repeat the baseline with independent patient-level validation.
- Obtain and harmonize clinical outcome endpoints for survival/prognosis research; tumour-normal classification is not a survival or therapy-response task.
- Calibrate biological parameter values and assess identifiability and uncertainty using compatible observed data.
- Compare numerical solver with independent methods and document multi-level convergence.
- Fully couple PDE and agent-based dynamics.
- Validate spatial graphs against registered histology/spatial assay ground truth.
- Integrate additional compatible modalities/cohorts and perform independent validation.
- Implement the full PDE-constrained control problem.
- Fabricate the phantom, assemble ESP32 hardware and complete bench tests.

## Scientific and clinical boundary

All simulation parameters remain dimensionless teaching proxies. The TCGA baseline predicts known specimen labels in a small research subset only. No current output has clinical validation or utility for screening, diagnosis, prognosis, patient-specific treatment selection, dose calculation or treatment recommendations. Controlled, restricted or identifiable participant-level data must not be committed to this public repository.
