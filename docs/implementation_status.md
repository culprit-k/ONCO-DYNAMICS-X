# Implementation status

This status separates executable code from model components that still require training, scientific validation or physical integration.

## Implemented in the repository

### Mathematical simulation
- Eight dimensionless spatial fields: cancer, immune, fibroblast, oxygen, nutrient, drug, immune-evasion proxy and vasculature proxy.
- Finite-difference reaction-diffusion solver on 2D and 3D Cartesian grids, no-flux edges and explicit Euler time stepping.
- Dimension-aware explicit diffusion stability guard.
- Synthetic initial states, CLI 2D/3D simulations and Streamlit 2D maps / 3D voxel-volume visualization.
- Timestep refinement and nested-grid refinement diagnostics for 2D and 3D.
- Independent adaptive SciPy DOP853 method-of-lines reference integrator and explicit-Euler comparison runner on identical spatial discretizations.
- Automated synthetic 2D/3D numerical validation workflow covering timestep refinement, grid refinement and adaptive-reference comparison.

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
- Endpoint-aware survival data assembler with explicit OS event/censor parsing and provisional broad-stage mapping; cross-cohort stage semantics remain unresolved and stage is excluded from the primary age-only model.
- Exploratory Cox proportional-hazards clinical baseline trained on TCGA (1,050 patients; 143 deaths) and evaluated once on METABRIC (1,454 patients; 821 deaths).
- External METABRIC Harrell C-index 0.616 (patient bootstrap 95% interval 0.592–0.641) for age + broad stage; age-only comparator C-index 0.591 (95% interval 0.570–0.617).
- Survival workflow uploads aggregate metrics only; no row-level clinical data, patient/sample identifiers or expression matrices are uploaded.
- Fixed 50-gene PAM50-panel research feature pipeline completed: cohort-wise mid-rank inverse-normal transformation, five-component PCA fitted in TCGA and applied frozen to METABRIC, and age + five-PC Cox external evaluation. All 50 genes were finite in both profiles; no imputation was required.
- Corrected cBioPortal DETAILED sample-object parsing and verified Primary Solid Tumor type for all 1,082 TCGA and 1,980 METABRIC profile samples; zero missing patient links and zero duplicate primary specimens. Full test suite and strict Ruff checks passed in the successful rerun.
- PAM50-panel external comparison: age-only C-index 0.5882; age + five PCs C-index 0.5992; paired delta +0.0111 with 95% bootstrap interval -0.0019 to +0.0241 (includes zero). This is exploratory and transductive; see docs/pam50_molecular_validation.md.
- Five-class intrinsic-subtype label transfer completed: train on 968 TCGA labelled patients; external test on 1,756 METABRIC patients with shared Basal/Her2/LumA/LumB/Normal labels. Accuracy 0.810 (95% CI 0.791–0.827), balanced accuracy 0.747 (0.724–0.768), macro-F1 0.778 (0.753–0.799), versus 0.399 majority-class accuracy. METABRIC claudin-low (217) and NC (6) labels excluded. This measures transfer of existing cohort annotations, not independent subtype discovery; see docs/subtype_transfer_benchmark.md and docs/subtype_transfer_results.md.
- A separate 40-scenario PDE-agent coupling sensitivity sweep completed for 2D/3D, five coupling strengths and four occupancy scales. At full coupling, mean lower-bound node occupancy for cell-density fields was 80–96%; the result indicates sparse-grid sensitivity, not biological calibration.
- Fixed 16-gene expression availability audit across TCGA-BRCA RSEM (1,082 samples) and METABRIC Illumina HT-12 v3 (1,980 samples): all marker-sample values finite and complete in the queried profiles; no missing/non-finite or duplicate records. No harmonization or model fitting performed.
- Censoring-aware age-only Cox evaluation trained in TCGA (1,069 patients; 151 deaths) and evaluated in METABRIC (1,979 patients; 1,143 deaths), with Harrell C-index 0.588 (bootstrap 95% interval 0.568–0.606).
- METABRIC time-dependent IPCW AUC at 12/36/60 months: 0.578 / 0.515 / 0.539; IPCW Brier scores: 0.0150 / 0.1014 / 0.1713; sparse-grid IBS 0.07935. Results and limitations are recorded in docs/age_only_external_results.md.

- Added two local longitudinal calibration tools: (i) untreated exponential/logistic/Gompertz volume fitting with AICc, bootstrap intervals and Jacobian/SVD identifiability diagnostics; (ii) treatment-exposed net log-volume slope fitting with a regression CI and patient-series bootstrap summaries. Added local-only CSV runners and aggregate-only clinical-target metadata audit. Synthetic recovery tests pass in CI; no real subject series has yet been fitted and no biological PDE coefficient is claimed. See docs/biological_calibration_protocol.md.
- Added an exploratory molecular-to-mechanistic feature interface: four prespecified expression modules (proliferation-associated, estrogen/luminal-associated, HER2-associated, basal/epithelial-associated) are computed as means of rank-normalised panel values with coverage/input checks. A machine-readable hypothesis registry documents candidate ABM phenotype/proliferation targets and explicitly prohibits direct PDE/ABM parameter updates pending matched-data calibration; see docs/molecular_mechanistic_bridge.md.

- Added an untreated longitudinal tumour-volume estimator comparing exponential, logistic and Gompertz curves using log-scale nonlinear least squares, AICc/BIC, residual bootstrap intervals and Jacobian/SVD diagnostics. Added a local CSV CLI and aggregate-only TCGA/METABRIC clinical-target metadata audit. These are calibration tools only: no real patient/animal series has yet been fitted and no biological PDE coefficient is claimed; see docs/biological_calibration_protocol.md.

- Added low-dimensional PDE parameter estimation in `src/onco_dynamics_x/calibration/pde_parameter_estimation.py`; test recovery of one synthetic `cancer_growth` coefficient using cancer/oxygen field snapshots. The fitter refuses more than three parameters per stage and withholds covariance when the Jacobian is rank-deficient/ill-conditioned. No fitted synthetic value changes the provisional model configuration.

### Statistical and model components
- CPU synthetic training/evaluation completed for FieldPINN, FourierNeuralOperator2D and SpatialCellGNN using analytic heat-equation and engineered graph targets. Latest boundary-aware PINN: held-out relative L2 0.02774, PDE residual RMS 0.02297 and Neumann boundary residual RMS 0.01188; FNO relative L2 0.01666; GNN relative L2 0.06995. These are software benchmarks, not biological training.
- Local finite-difference parameter sensitivity and Monte Carlo propagation of user-supplied parameter ranges.
- Optional PyTorch PINN, Fourier neural operator, spatial GNN and expression-transformer architectures.
- Optional PyMC posterior builder, DoWhy linear backdoor wrapper and CasADi mean-field control demonstration.

### Hardware and presentation
- ESP32 firmware, normalized serial protocol and dashboard control hook.
- Parametric OpenSCAD breast-phantom shell, removable lesion inserts, wiring map and physical build guide.

## Implemented as architecture only or optional components

PINN, GNN, neural operator, expression transformer, Bayesian, causal and optimal-control modules are code components. They are not trained or biologically validated models by simply being present in source.

The age-only and age-plus-stage Cox models are clinical-only exploratory baselines. The separate age-plus-PAM50-PC model is an exploratory molecular association analysis, not a validated clinical predictor. None of these analyses validates the spatial digital twin.

## Still outstanding
- Review source-level endpoint definitions and raw TCGA/METABRIC stage-value distributions; aggregate category counts are now recorded, but TCGA patient-level AJCC pathologic stage and METABRIC sample-level generic tumour stage are not yet established as equivalent.
- Extend calibration and time-dependent survival evaluation to future molecular models and alternative endpoints; review marginal independent-censoring assumptions. Keep the age-only model as the primary clinical reference until stage-field semantics are verified.
- Complete array-level METABRIC probe annotation and duplicate-probe summarisation; the current cBioPortal gene-level PAM50 profile did not rederive the underlying probe mapping.
- Replace transductive cohort-wise rank normalisation with a prespecified reference/frozen transform suitable for single-sample or prospective use; validate on a second independent cohort without using external outcomes for feature selection or tuning.
- Verify overall-survival time-origin and event-definition provenance against the original cohort guides. Primary Solid Tumor sample types are now verified using cBioPortal's DETAILED Sample projection; retain the full aggregate sample-type audit in the molecular result artifact.
- Verify METABRIC ER/PR/HER2/grade/size/cellularity field definitions and exact expression-sample overlap; the aggregate audit found candidate fields but no Ki-67 field. Metadata coverage is not a matched assay validation.
- Obtain approved serial tumour-volume or registered MRI cellularity data; fit only untreated intervals for the first effective growth-rate estimates, and reserve a later timepoint for prediction.
- Added a low-dimensional PDE parameter estimator for one to three positive reaction coefficients using registered field snapshots, weighted nonlinear least squares, Jacobian rank/conditioning and conditional local covariance. It has synthetic recovery tests only; no patient-derived parameter is estimated. PDE diffusion, oxygen, immune, stromal and drug parameters still require matched spatial/physical observations and a validated observation operator. Assess parameter confounding, profile likelihood/posterior correlation and uncertainty; no patient/animal-derived PDE parameter has yet been estimated.
- Expand the adaptive DOP853 comparison to a multi-resolution convergence ladder and report the resulting error trends; numerical consistency does not establish biological validity.
- Refine and biologically calibrate PDE–agent coupling: mass/number accounting, deposition-kernel comparison (nearest-grid versus cloud-in-cell/kernel), multi-seed coupling sensitivity and operator-splitting convergence. A single-seed 40-scenario sweep showed high lower-bound occupancy at strong coupling.
- Validate spatial graphs against registered histology/spatial assay ground truth.
- Integrate additional compatible modalities/cohorts.
- Implement the full PDE-constrained control problem.
- Fabricate the phantom, assemble ESP32 hardware and complete bench tests.

## Scientific and clinical boundary

The 48-patient TCGA tumour/normal classifier predicts known specimen labels only. The Cox model is an exploratory clinical-only cross-cohort evaluation, not a clinically validated prognostic tool. Neither result validates the spatial PDE/ABM, PINNs, GNNs, immune evasion, hypoxia fields or simulated drug response. No current output is suitable for screening, diagnosis, patient-specific prognosis, treatment selection, dose calculation or treatment recommendations. Controlled, restricted or identifiable participant-level data must not be committed to this public repository.
