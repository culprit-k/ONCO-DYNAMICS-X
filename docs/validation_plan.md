# Numerical validation plan

## Software verification already in CI
- Constant-field Laplacian is zero under the no-flux finite-difference operator.
- 2D and 3D state arrays remain finite and dimensionally consistent.
- Explicit diffusion timestep guard rejects a time increment above the dimension-dependent bound.
- The PDE solver is deterministic for a fixed initial state and parameter set.
- Timestep refinement compares dt and dt/2 at equal final time.
- Grid refinement compares nested grids at equal final time by sampling the fine-grid solution at coarse-grid nodes; this is tested for 2D and 3D.
- Agent simulation tests cover deterministic seeds, boundary constraints, cell phenotypes and 2D/3D environment matching.
- Cohort/GDC tests cover open-access enforcement, sample linkage, metadata parsing, checksums and expression-matrix alignment.

## Required next numerical studies
1. Run a grid-refinement ladder (at least three resolutions) and report field-wise error norms.
2. Compare the current explicit finite-difference solver against an independent method-of-lines stiff integrator and, later, an independent spatial discretization.
3. Test positivity, clipping frequency, mass-balance effects and reaction stiffness.
4. Check parameter identifiability and confidence/credible intervals under synthetic known-parameter recovery experiments.
5. Benchmark PINNs and neural operators against the finite-difference reference on held-out initial conditions and parameter settings.

## Biological validation
Numerical convergence is not biological validity. Biological validation requires permitted, compatible, independent datasets with correctly linked sample identifiers, endpoint definitions and modality-specific ground truth. Spatial claims require spatially resolved observations. Report limitations and failure cases; do not infer treatment benefit from prognostic associations.
