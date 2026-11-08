# Implementation roadmap

All originally requested workstreams remain in scope. Phases describe implementation dependency, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI and tests | Reproducible install; automated checks | Implemented |
| 1 | Coupled 2D/3D PDE/ODE engine | Positivity/stability checks, timestep and grid refinement, independent solver | Reference prototype and diagnostics implemented; full convergence report pending |
| 2 | PINNs, Bayesian inference and neural operators | Trained models benchmarked against numerical reference; posterior diagnostics | Optional architectures/builders implemented; training and calibration pending |
| 3 | Agent-based model and spatial GNN | Reproducible cell-level experiments and spatial validation | Initial 2D/3D NumPy ABM, radius graph builder and GNN architecture implemented; PDE feedback and biological calibration pending |
| 4 | Multi-omics and cohort integration | Audited metadata, sample-level linkage, preprocessing provenance and external cohort | Public GDC metadata/file workflow, TCGA matrix assembly, expression preprocessing and transformer architecture implemented; additional cohorts/modalities pending |
| 5 | Optimal control, causal analysis and uncertainty | Explicit assumptions, identifiable estimands, sensitivity and uncertainty reporting | Sensitivity, range-based Monte Carlo, optional DoWhy and toy CasADi OCP implemented; real-data analysis pending |
| 6 | 3D digital twin and dashboard | Computed volumetric fields, interaction, experiment logs and uncertainty | True 3D PDE field and volume rendering implemented; complete integrated twin pending |
| 7 | Breast phantom and hardware console | Physical fabrication, safe wiring, serial and bench tests | ESP32 sketch/parser, dashboard serial hook, OpenSCAD source and build guide implemented; physical construction/test pending |
| 8 | University demonstration and report | End-to-end reproducible demonstration and limitations report | Pending |

## Immediate next tasks
1. Review a live open-access TCGA-BRCA file query and inspect the exported sample manifest.
2. Download a small, deliberately selected set of open expression files and assemble the gene matrix.
3. Run patient-grouped cohort partitions and establish a leakage-safe baseline analysis.
4. Run a three-resolution grid-refinement ladder and an independent integrator comparison.
5. Train the PINN/FNO/GNN components first on synthetic ground-truth tasks, then evaluate their suitability for biological observations.
6. Build the fully coupled PDE-ABM interaction loop and validate it on controlled synthetic experiments.
7. Connect spatially resolved histology or spatial-omics data only after registration and sample identity are confirmed.
8. Bench-test the ESP32 and physically fabricate the phantom before the live university demonstration.
