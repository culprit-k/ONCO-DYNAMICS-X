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
- Bidirectional operator-split PDE–agent coupling: PDE fields drive agent updates; agent phenotypes are deposited back to the corresponding continuum density fields with explicit dimensionless coupling/occupancy controls.
- Radius-based cell-neighbour graph construction and a GNN architecture.
- Patient-grouped train/test split utility.

### Data engineering and real-cohort pipelines
- NCI GDC public/open metadata query, sample-manifest export and audit.
- Sample-type-specific filtering, including Primary Tumor and Solid Tissue Normal.
- Selected open-file download with retry handling, byte-count verification and MD5 checksum verification.
- Strict expression matrix assembly with sample/gene alignment and identifier-aware preprocessing.
- Paired TCGA-BRCA tumour/normal selection and successful 48-patient / 96-sample baseline workflow.
- Tumour-versus-normal baseline using log1p(TPM), train-fold-only feature filtering/selection/scaling and patient-grouped holdout/CV. Its perfect score is an engineering result on known specimen labels, not clinical diagnostic performance.
- Read-only cBioPortal API client with aggregate-only audits for METABRIC and TCGA-BRCA PanCancer Atlas.
- METABRIC audit: 2,509 patients/samples, 1,980 mRNA samples, 1,980 expression-linked patients with complete OS fields, 1,979 with positive OS follow-up.
- TCGA PanCancer Atlas audit: 1,084 patients/samples, 1,082 RNA-seq samples, 1,082 with complete OS fields, 1,069 with positive OS follow-up.
- Patient-level endpoint pair completeness, non-positive time diagnostics, sample-to-patient linkage, clinical-attribute coverage and candidate baseline covariate overlap.
- Endpoint-aware survival data assembler with explicit OS event/censor parsing and broad stage harmonization.
- Exploratory Cox proportional-hazards clinical baseline trained on TCGA (1,050 patients; 143 deaths) and evaluated once on METABRIC (1,454 patients; 821 deaths).
- External METABRIC Harrell C-index 0.616 (patient bootstrap 95% interval 0.592–0.641) for age + broad stage; age-only comparator C-index 0.591 (95% interval 0.570–0.617).
- Survival workflow uploads aggregate metrics only; no row-level clinical data, patient/sample identifiers or expression matrices are uploaded.

### Statistical and model components
- Local finite-difference parameter sensitivity and Monte Carlo propagation of user-supplied parameter ranges.
- Optional PyTorch PINN, Fourier neural operator, spatial GNN and expression-transformer architectures.
- Optional PyMC posterior builder, DoWhy linear backdoor wrapper and CasADi mean-field control demonstration.

### Hardware and presentation
- ESP32 firmware, normalized serial protocol and dashboard control hook.
- Parametric OpenSCAD breast-phantom shell, removable lesion inserts, wiring map and physical build guide.

## Implemented as architecture only or optional components

PINN, GNN, neural operator, expression transformer, Bayesian, causal and optimal-control modules are code components. They are not trained or biologically validated models by simply being present in source.

The survival Cox model is a clinical-only exploratory baseline. It does not use molecular expression and does not validate the digital twin.

## Still outstanding
- Review source-level endpoint definitions and raw TCGA/METABRIC stage-value distributions; aggregate category counts are now recorded, but TCGA patient-level AJCC pathologic stage and METABRIC sample-level generic tumour stage are not yet established as equivalent.
- Assess calibration and time-dependent survival performance, in addition to discrimination; resolve stage-field semantics and the below-0.5 stage-only transport warning before relying on stage-adjusted external estimates or advancing to molecular modelling.
- Build a fixed molecular feature/signature pipeline and solve cross-platform harmonization between TCGA RNA-seq RSEM and METABRIC Illumina HT-12 v3 microarray.
- Perform independent molecular evaluation without using METABRIC outcomes for feature selection or tuning.
- Calibrate biological parameter values and assess identifiability and uncertainty using compatible observed data.
- Compare numerical solver with independent methods and document multi-level convergence.
- Refine and biologically calibrate PDE–agent coupling: mass/number accounting, deposition-kernel comparison, coupling-strength sensitivity and operator-splitting convergence.
- Validate spatial graphs against registered histology/spatial assay ground truth.
- Integrate additional compatible modalities/cohorts.
- Implement the full PDE-constrained control problem.
- Fabricate the phantom, assemble ESP32 hardware and complete bench tests.

## Scientific and clinical boundary

The 48-patient TCGA tumour/normal classifier predicts known specimen labels only. The Cox model is an exploratory clinical-only cross-cohort evaluation, not a clinically validated prognostic tool. Neither result validates the spatial PDE/ABM, PINNs, GNNs, immune evasion, hypoxia fields or simulated drug response. No current output is suitable for screening, diagnosis, patient-specific prognosis, treatment selection, dose calculation or treatment recommendations. Controlled, restricted or identifiable participant-level data must not be committed to this public repository.
