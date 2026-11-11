# Implementation roadmap

All originally requested workstreams remain in scope. Phases describe implementation dependency, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI and tests | Reproducible install; automated checks | Implemented |
| 1 | Coupled 2D/3D PDE/ODE engine | Positivity/stability checks, timestep and grid refinement, independent solver | Explicit solver, timestep/grid diagnostics and adaptive DOP853 method-of-lines reference implemented; automated 2D/3D numerical validation running; multi-resolution convergence report pending |
| 2 | PINNs, Bayesian inference and neural operators | Trained models benchmarked against numerical reference; posterior diagnostics | Optional architectures/builders implemented; training and calibration pending |
| 3 | Agent-based model and spatial GNN | Reproducible cell-level experiments and spatial validation | 2D/3D NumPy ABM, radius graph builder, GNN architecture and operator-split bidirectional PDE-agent coupling implemented; mass-consistent coupling and biological calibration pending |
| 4 | Multi-omics and cohort integration | Audited metadata, sample-level linkage, endpoint-specific analysis and external evaluation | GDC/cBioPortal audits; age-only TCGA-trained/METABRIC-tested Cox model with IPCW AUC, Brier and calibration summaries implemented. Molecular harmonization and external molecular validation pending |
| 5 | Optimal control, causal analysis and uncertainty | Explicit assumptions, identifiable estimands, sensitivity and uncertainty reporting | Sensitivity, range-based Monte Carlo, optional DoWhy and toy CasADi OCP implemented; real-data causal/control analysis pending |
| 6 | 3D digital twin and dashboard | Computed volumetric fields, interaction, experiment logs and uncertainty | True 3D PDE field and volume rendering implemented; complete integrated twin pending |
| 7 | Breast phantom and hardware console | Physical fabrication, safe wiring, serial and bench tests | ESP32 sketch/parser, dashboard serial hook, OpenSCAD source and build guide implemented; physical construction/test pending |
| 8 | University demonstration and report | End-to-end reproducible demonstration and limitations report | Pending |

## Immediate next tasks

1. Verify the broad-stage mapping against both cohort data dictionaries and perform sensitivity analyses for missing/unrecognized stage.
2. Extend survival calibration/time-dependent metrics to molecular models and additional endpoints; assess censoring assumptions and horizon stability.
3. Build a fixed, biologically justified molecular gene/signature feature set with a common gene identifier mapping.
4. Compare platform-harmonization approaches for TCGA RNA-seq RSEM and METABRIC Illumina HT-12 v3 microarray without using external outcomes for tuning.
5. Evaluate the frozen molecular score in METABRIC and report uncertainty, calibration and failure cases.
6. Complete a three-resolution PDE grid-refinement ladder and summarize adaptive-integrator agreement.
7. Train PINN/FNO/GNN components first on synthetic ground-truth tasks, then evaluate suitability for biological observations.
8. Validate PDE-ABM mass/number consistency, deposition kernels, coupling-strength sensitivity and splitting-order effects on controlled synthetic experiments.
9. Connect spatially resolved histology or spatial-omics data only after registration and sample identity are confirmed.
10. Bench-test the ESP32 and physically fabricate the phantom before the live university demonstration.
