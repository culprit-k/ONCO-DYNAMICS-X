# Synthetic AI architecture training benchmark

**Status:** completed successfully  
**Date:** 11 November 2026  
**Workflow:** [Synthetic PINN, FNO and GNN training run 42982733597](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42982733597)  
**Source commit evaluated:** 93f4941faf14d7fba32ca3e92519edf3c3a338f6  
**Aggregate result:** results/ai_benchmarks/synthetic_pinn_fno_gnn_summary.json

## Scope

Three optional architecture implementations were trained and evaluated on distinct, deterministic synthetic tasks with known mathematical targets:

1. FieldPINN: analytic two-dimensional no-flux heat-equation solution, trained with data loss plus a PDE residual penalty.
2. FourierNeuralOperator2D: exact periodic heat-semigroup mapping generated through Fourier multipliers.
3. SpatialCellGNN: generated spatial graphs with a target equal to the fraction of immune-labelled neighbours within a fixed radius.

All runs used CPU, Python 3.11 and seed 42. Training, validation and test samples/graphs were generated separately. No breast-cancer patient, clinical, molecular, histology, imaging or treatment-response data were used.

## Results

| Architecture | Synthetic held-out task | Held-out metric | Baseline / diagnostic |
|---|---|---:|---:|
| FieldPINN | 2D Neumann heat solution | Relative L2 error **0.0909** | Untrained initialization relative L2: 1.3800; held-out PDE residual RMS: 0.05325 |
| FNO | 2D periodic heat operator | Relative L2 error **0.01363** | Persistence baseline relative L2: 8.1083 |
| Spatial GNN | Local immune-neighbour fraction on generated graphs | Test MSE **0.000687**; relative L2: 0.06995 | Constant-mean baseline MSE: 0.04054 |

### FieldPINN

- Training / validation / test: 1,024 / 512 / 1,024 coordinate samples.
- 500 epochs; best validation MSE 0.0009141.
- Held-out relative L2 error: 9.09%.
- Held-out PDE residual RMS: 0.05325.

The network improves substantially over its untrained initialization, but the remaining approximation error and PDE residual are not negligible. This is a working training pipeline, not a claim of high-accuracy PINN performance.

### Fourier neural operator

- Training / validation / test: 128 / 32 / 64 fields on 16 × 16 grids.
- 120 epochs; best validation MSE 3.00 × 10⁻⁶.
- Held-out relative L2 error: 1.36%.

The FNO learns the selected periodic heat operator and outperforms the no-evolution persistence baseline on this synthetic task. This does not establish transfer to tumour dynamics or irregular tissue geometry.

### Spatial GNN

- Training / validation / test: 48 / 12 / 12 independently generated graphs.
- Each graph has 64 nodes; 3,072 / 768 / 768 nodes across the respective splits.
- 80 epochs; best validation MSE 0.000702.
- Test MSE: 0.000687, compared with constant-mean baseline MSE 0.04054.

The target is computed directly from the synthetic graph's node labels and radius-based adjacency. This is a software benchmark of message passing and generalization to new generated graphs, not a measurement or prediction of real tumour immune infiltration.

## CI result

The dedicated GitHub Actions workflow completed successfully:

- CPU PyTorch and development dependencies installed.
- Full test suite passed.
- Ruff lint passed.
- All three training jobs completed.
- Aggregate JSON artifact uploaded.

## Interpretation and limitations

These runs demonstrate that the PINN, FNO and GNN implementations can be trained, checkpointed using validation performance and evaluated on held-out synthetic examples. They do not rank the architectures because each architecture solves a different task.

No biological model was trained here. There is no evidence from these metrics about breast-cancer subtype, immune escape, tumour progression, drug response, recurrence or patient prognosis. The PINN still has measurable residual error; its training objective and architecture should be improved before it is used as a numerical surrogate. The GNN target is intentionally synthetic and should not be described as an inferred tumour-microenvironment label.

## Next AI validation gate

1. Add several analytic solutions and coefficient regimes, including heterogeneous diffusion and nontrivial reaction terms.
2. Compare against a simple numerical baseline on the same held-out tasks.
3. Report error by time, spatial frequency and boundary region, not only aggregate relative L2.
4. Test stability across seeds and grid resolutions.
5. For the GNN, vary graph density, radius and node count, and include a non-message-passing baseline.
6. Only after those software checks, design a separate cohort-specific molecular or spatial task with documented data provenance, patient-level splits and external validation.
