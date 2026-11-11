# Implementation roadmap

All originally requested workstreams remain in scope. Phases describe implementation dependency, not feature removal.

| Phase | Workstream | Acceptance evidence | Status |
|---|---|---|---|
| 0 | Repository, docs, package, CI and tests | Reproducible install; automated checks | Implemented |
| 1 | Coupled 2D/3D PDE/ODE engine | Positivity/stability checks, timestep and grid refinement, independent solver | Explicit solver, timestep/grid diagnostics and adaptive DOP853 method-of-lines reference implemented; automated 2D/3D numerical validation running; multi-resolution convergence report pending |
| 2 | PINNs, Bayesian inference and neural operators | Trained models benchmarked against numerical reference; posterior diagnostics | PINN/FNO/GNN synthetic training benchmarks completed; PINN now includes initial/Neumann boundary losses. Biological training, Bayesian calibration and independent validation pending |
| 3 | Agent-based model and spatial GNN | Reproducible cell-level experiments and spatial validation | 2D/3D NumPy ABM, radius graph builder, GNN architecture and operator-split bidirectional PDE-agent coupling implemented; 40-scenario 2D/3D parameter sweep completed. Mass/number consistency, deposition-kernel comparisons and biological calibration pending |
| 4 | Multi-omics and cohort integration | Audited metadata, sample-level linkage, endpoint-specific analysis and external evaluation | GDC/cBioPortal audits; age-only Cox evaluation; fixed 50-gene PAM50-panel-PC model evaluated in METABRIC (C=0.599 vs age-only 0.588; paired delta CI includes zero). Probe-level mapping, deployable transform and independent-cohort confirmation remain pending |
| 5 | Optimal control, causal analysis and uncertainty | Explicit assumptions, identifiable estimands, sensitivity and uncertainty reporting | Sensitivity, range-based Monte Carlo, optional DoWhy and toy CasADi OCP implemented; real-data causal/control analysis pending |
| 6 | 3D digital twin and dashboard | Computed volumetric fields, interaction, experiment logs and uncertainty | True 3D PDE field and volume rendering implemented; complete integrated twin pending |
| 7 | Breast phantom and hardware console | Physical fabrication, safe wiring, serial and bench tests | ESP32 sketch/parser, dashboard serial hook, OpenSCAD source and build guide implemented; physical construction/test pending |
| 8 | University demonstration and report | End-to-end reproducible demonstration and limitations report | Pending |

## Immediate next tasks

1. Verify the broad-stage mapping against both cohort data dictionaries and perform sensitivity analyses for missing/unrecognized stage.
2. Extend survival calibration/time-dependent metrics to molecular models and additional endpoints; assess censoring assumptions and horizon stability.
3. Reconstruct and verify the METABRIC array probe-to-gene mapping and duplicate-probe summarisation; the present gene-level analysis does not replace this step.
4. Define a deployable frozen reference transform for the fixed molecular feature panel; the current cohort-wise rank transform is transductive.
5. Compare platform-harmonization approaches for TCGA RNA-seq RSEM and METABRIC Illumina HT-12 v3 microarray without using external outcomes for tuning.
6. Confirm the molecular result in a second independent cohort with frozen processing and report uncertainty, calibration and failure cases.
7. Complete a three-resolution PDE grid-refinement ladder and summarize adaptive-integrator agreement.
8. Extend PINN/FNO/GNN synthetic tests to cross-resolution and out-of-distribution settings before evaluating suitability for biological observations.
9. Use the completed 40-scenario coupling sweep as an initial sensitivity record; next test PDE-ABM mass/number consistency, multiple deposition kernels, multi-seed robustness and splitting-order effects.
10. Connect spatially resolved histology or spatial-omics data only after registration and sample identity are confirmed.
11. Bench-test the ESP32 and physically fabricate the phantom before the live university demonstration.
