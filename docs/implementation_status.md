# Implementation status

This status separates executable code from model components that still require training, scientific validation or physical integration.

## Implemented in the repository

### Mathematical simulation
- Eight dimensionless spatial fields: cancer, immune, fibroblast, oxygen, nutrient, drug, immune-evasion proxy and vasculature proxy.
- Finite-difference reaction-diffusion solver on 2D and 3D Cartesian grids, no-flux edges and explicit Euler time stepping.
- Dimension-aware explicit diffusion stability guard.
- Synthetic initial states, CLI 2D/3D simulations and Streamlit 2D maps / 3D voxel-volume visualization.
- Timestep refinement and nested-grid refinement diagnostics for 2D and 3D. These measure numerical sensitivity, not biological validity.

### Cellular and spatial AI inputs
- Stochastic 2D/3D agent-based simulation for cancer, immune and fibroblast cell proxies, including movement, local environment coupling, recruitment and division/death rules.
- Radius-based cell-neighbour graph construction with coordinates, cell-type encoding, age and local PDE fields as node features.
- GNN architecture consuming these spatial graphs.
- Patient-grouped train/test split utility to keep a patient's samples in one partition.

### Data engineering
- GDC open-access metadata query with sample/case provenance.
- GDC sample-manifest export and required-column / duplicate-ID audit.
- Explicit selected-file open-access downloader with UUID/name checks, byte-count verification and MD5 checksum verification.
- TCGA expression-file parser and strict matrix assembler; requires exactly one sample link per file and an identical gene set across files.
- Identifier-aware expression matrix validation, training-only variance filtering and scaling.

### Statistical and model components
- Local finite-difference parameter sensitivity.
- Monte Carlo propagation of user-supplied parameter ranges (not a Bayesian posterior).
- Optional PyTorch PINN, Fourier neural operator and fixed-order expression transformer.
- Optional PyMC Gaussian-likelihood posterior builder with user-supplied forward model.
- Optional DoWhy linear backdoor wrapper requiring analyst-specified common causes.
- Optional CasADi dimensionless mean-field optimal-control example.

### Hardware and presentation
- ESP32 firmware sketch, normalized JSON serial protocol, parser and dashboard control hook.
- Parametric OpenSCAD breast-phantom shell, separate lesion inserts and stand ring.
- Wiring map, component list, physical build guide and dashboard deployment notes.

## Implemented as architecture only or optional components

The PINN, GNN, neural operator, expression transformer, Bayesian, causal and optimal-control modules are code components. They are not trained models and do not generate scientifically validated predictions by simply being imported. Optional dependencies are not part of the minimal base install.

## Still outstanding

- Execute and review an actual permitted TCGA-BRCA cohort query; inspect sample types, assay metadata, endpoints and sample linkage before selecting files.
- Train and benchmark AI models with suitable train/validation/test splits and an independent cohort.
- Calibrate biological parameter values and assess identifiability, posterior convergence and uncertainty using compatible observed data.
- Compare the numerical engine with an independent solver and document multi-level grid/time convergence.
- Fully couple PDE and agent-based cell dynamics through explicit source/sink feedback.
- Validate spatial graphs against registered histology/spatial assay ground truth.
- Integrate and validate additional compatible modalities/cohorts (METABRIC, GEO, proteomics or spatial data).
- Add clonal heterogeneity, detailed ECM mechanics, angiogenesis and more biologically resolved immune phenotypes.
- Implement the full PDE-constrained control problem; current CasADi module is a toy mean-field demonstration.
- Physically print/cast the phantom, wire the ESP32 console and complete bench integration tests.
- Configure and authorize a hosted Streamlit deployment; local deployment instructions are present.

## Scientific and clinical boundary

All current simulation parameters and initial fields are dimensionless teaching proxies. Synthetic tests establish software behaviour only. No current output has clinical validation or utility for screening, diagnosis, prognosis, patient-specific treatment selection, dose calculation or treatment recommendations. Controlled, restricted or identifiable participant-level data must not be committed to this public repository.
