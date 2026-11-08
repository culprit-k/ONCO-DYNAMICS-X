# Numerical and clinical validation plan

## Software verification already in CI
- Constant-field Laplacian is zero under the no-flux finite-difference operator.
- 2D and 3D state arrays remain finite and dimensionally consistent.
- Explicit diffusion timestep guard rejects a time increment above the dimension-dependent bound.
- The PDE solver is deterministic for a fixed initial state and parameter set.
- Timestep refinement compares dt and dt/2 at equal final time.
- Grid refinement compares nested grids at equal final time by sampling the fine-grid solution at coarse-grid nodes; this is tested for 2D and 3D.
- Agent simulation tests cover deterministic seeds, boundary constraints, cell phenotypes and 2D/3D environment matching.
- GDC/cBioPortal tests cover open-access enforcement, sample linkage, metadata parsing, checksums, expression-matrix alignment, survival endpoint parsing, stage normalization and aggregate-only reporting.

## Required numerical studies
1. Run a grid-refinement ladder (at least three resolutions) and report field-wise error norms.
2. Compare the current explicit finite-difference solver against an independent method-of-lines stiff integrator and, later, an independent spatial discretization.
3. Test positivity, clipping frequency, mass-balance effects and reaction stiffness.
4. Check parameter identifiability and confidence/credible intervals under synthetic known-parameter recovery experiments.
5. Benchmark PINNs and neural operators against the finite-difference reference on held-out initial conditions and parameter settings.

## Clinical survival analysis: current baseline and required checks

An exploratory clinical-only Cox model has been fitted on the TCGA-BRCA PanCancer Atlas cohort and evaluated once on METABRIC. It uses age and broad stage groups, not gene expression. The METABRIC Harrell C-index is 0.616 (patient-bootstrap 95% interval 0.592–0.641); this is modest discrimination and does not establish clinical utility.

Required before stronger interpretation:
1. Verify OS time origin, event/censor coding and follow-up units against each source cohort's documentation.
2. Verify broad-stage mapping across TCGA AJCC pathologic stage and METABRIC TUMOR_STAGE; report ambiguous categories.
3. Perform sensitivity analyses for complete-case stage selection and missingness.
4. Assess proportional-hazards assumptions, calibration at prespecified time horizons and time-dependent AUC/Brier score.
5. Compare with simple age-only and clinically appropriate reference models.
6. Keep all preprocessing and model selection in TCGA development data; do not tune on METABRIC outcomes.
7. Treat the result as an exploratory retrospective cohort analysis, not clinical validation.

## Biological validation

Numerical convergence is not biological validity. Biological validation requires permitted, compatible, independent datasets with correctly linked sample identifiers, endpoint definitions and modality-specific ground truth. Spatial claims require spatially resolved observations. Bulk RNA expression and a clinical survival score cannot validate cell positions, cell contacts, oxygen gradients, drug penetration or treatment benefit. Report limitations and failure cases.
