# Synthetic PDE numerical validation — results

**Status:** completed successfully  
**Date:** 11 November 2026  
**Workflow:** [PDE numerical validation run 42980986063](https://github.com/culprit-k/ONCO-DYNAMICS-X/actions/runs/42980986063)  
**Source commit evaluated:** 4b381346e54c876bbac2d7836dd8f5d12f53386a  
**Aggregate summary:** results/numerical_validation/synthetic_pde_validation_summary.json

## What was tested

The existing explicit finite-difference reaction-diffusion solver was checked in small 2D and 3D synthetic systems using three comparisons:

1. **Adaptive time-integration reference:** SciPy DOP853 integrated the same semi-discrete finite-difference equations with adaptive time stepping (relative tolerance 1e-8; absolute tolerance 1e-10). The adaptive solver does not apply the explicit solver's post-step clipping.
2. **Time-step refinement:** explicit Euler runs at dt and dt/2 were compared at the same final time.
3. **Nested-grid refinement:** coarse and fine grids were compared at common coarse-grid nodes with dx and dt divided by two and step count multiplied by four.

All experiments used default dimensionless parameters and deterministic synthetic initial conditions. There are no patient observations in this validation.

## Results

### Explicit Euler versus adaptive DOP853

| Dimension | Grid | Final time | Maximum per-field relative L2 error | Mean per-field relative L2 error |
|---|---:|---:|---:|---:|
| 2D | 9 × 9 | 0.008 | 2.10 × 10⁻⁶ | 5.68 × 10⁻⁷ |
| 3D | 5 × 5 × 5 | 0.004 | 2.35 × 10⁻⁷ | 9.54 × 10⁻⁸ |

The two time integrators closely agree for these short, small-grid runs. The same spatial discretization and reaction equations are used in both implementations, so this primarily checks time-integration consistency.

### Time-step refinement

The maximum relative L2 differences across the eight fields were approximately:

- 2D: 1.05 × 10⁻⁶ (immune field).
- 3D: 1.18 × 10⁻⁷ (immune field).

The drug field had zero difference because no drug input was applied and the initial drug field was zero.

### Nested-grid refinement

The maximum relative L2 differences after sampling the fine solution at coarse-grid nodes were:

- 2D: 1.69 × 10⁻³ (immune field).
- 3D: 1.16 × 10⁻³ (cancer field).

These are single coarse/fine comparisons. They do not establish a formal spatial convergence order.

## Interpretation

The explicit solver is internally consistent with an adaptive time integrator on these demonstration cases. The small errors support the correctness of the implemented time stepping for these chosen synthetic conditions, but do not prove stability for all parameter regimes, coupling strengths, infusion inputs, or grid sizes.

The nested-grid comparison is encouraging but limited: only one refinement ratio was tested in each dimension, the same finite-difference stencil was used, and the grids were deliberately small. A three-resolution sequence and an independent spatial discretization are still required before making stronger convergence claims.

## Limitations

- Initial conditions and parameters are synthetic and dimensionless.
- No patient-level or tissue-level data were used.
- No biological parameter calibration or identifiability analysis was performed.
- The adaptive reference shares the same spatial stencil and reaction terms; it is not an independent biological model.
- The explicit solver clips fields after each step; the adaptive reference does not. This difference was not activated in a way that made the tested results diverge, but should be monitored in more stressed cases.
- The coupled PDE–agent system has separate numerical risks; this workflow validates the uncoupled PDE solver only.
- The tests do not validate PINN, FNO, GNN, immune escape, hypoxia mechanisms, drug transport, or clinical predictions.

## Next numerical milestone

Run at least three nested grid resolutions in each dimension, compare more than one spatial discretization or boundary treatment, quantify observed convergence order where the solution is in the asymptotic regime, and repeat the tests over a predeclared parameter/stiffness range. For the coupled PDE–agent model, separately test population-number accounting, deposition kernel sensitivity, coupling-strength sensitivity and operator-splitting order.
