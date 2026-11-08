# Implementation roadmap

All originally requested workstreams remain in scope. Phases describe implementation dependency, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI and tests | Reproducible install; automated checks | Implemented |
| 1 | Coupled 2D/3D PDE/ODE engine | Positivity/stability checks, timestep and grid refinement, independent solver | Reference prototype and diagnostics implemented; full convergence report pending |
| 2 | PINNs, Bayesian inference and neural operators | Trained models benchmarked against numerical reference; posterior diagnostics | Optional architectures/builders implemented; training and calibration pending |
| 3 | Agent-based model and spatial GNN | Reproducible cell-level experiments and spatial validation | Initial 2D/3D NumPy ABM, radius graph builder and GNN architecture implemented; PDE feedback and biological calibration pending |
| 4 | Multi-omics and cohort integration | Audited metadata, sample-level linkage, endpoint-specific analysis and external evaluation | GDC expression workflow; cBioPortal TCGA/METABRIC audits; exploratory TCGA-trained/METABRIC-tested clinical Cox baseline implemented. Molecular harmonization and external molecular validation pending |
| 5 | Optimal control, causal analysis and uncertainty | Explicit assumptions, identifiable estimands, sensitivity and uncertainty reporting | Sensitivity, range-based Monte Carlo, optional DoWhy and toy CasADi OCP implemented; real-data causal/control analysis pending |
| 6 | 3D digital twin and dashboard | Computed volumetric fields, interaction, experiment logs and uncertainty | True 3D PDE field and volume rendering implemented; complete integrated twin pending |
| 7 | Breast phantom and hardware console | Physical fabrication, safe wiring, serial and bench tests | ESP32 sketch/parser, dashboard serial hook, OpenSCAD source and build guide implemented; physical construction/test pending |
| 8 | University demonstration and report | End-to-end reproducible demonstration and limitations report | Pending |

## Immediate next tasks

1. Verify the broad-stage mapping against both cohort data dictionaries and perform sensitivity analyses for missing/unrecognized stage.
2. Add survival calibration and time-dependent discrimination metrics, and review proportional-hazards assumptions.
3. Build a fixed, biologically justified molecular gene/signature feature set with a common gene identifier mapping.
4. Compare platform-harmonization approaches for TCGA RNA-seq RSEM and METABRIC Illumina HT-12 v3 microarray without using external outcomes for tuning.
5. Evaluate the frozen molecular score in METABRIC and report uncertainty, calibration and failure cases.
6. Run a three-resolution PDE grid-refinement ladder and independent integrator comparison.
7. Train PINN/FNO/GNN components first on synthetic ground-truth tasks, then evaluate suitability for biological observations.
8. Build the fully coupled PDE-ABM interaction loop and validate it on controlled synthetic experiments.
9. Connect spatially resolved histology or spatial-omics data only after registration and sample identity are confirmed.
10. Bench-test the ESP32 and physically fabricate the phantom before the live university demonstration.
