# Bidirectional PDE–agent-based coupling

**Implementation:** added in src/onco_dynamics_x/cellular_simulation/coupling.py  
**Validation status:** software coupling tests added; biological calibration is not established.

## Coupling scheme

The original reaction-diffusion PDE and cell-agent model can now exchange information at each macro time step using an explicit operator-splitting scheme:

1. Advance the eight continuum fields by one PDE solver step.
2. Advance discrete cancer, immune and fibroblast agents in the updated PDE environment. Agents sample local cancer, oxygen, nutrient, immune, drug, fibroblast and immune-evasion fields.
3. Deposit the updated agents onto the PDE grid by nearest-grid mapping, separately for cancer, immune and fibroblast phenotypes.
4. Convert per-voxel counts into dimensionless density proxies using a user-specified agents-per-unit scale.
5. Blend the agent-derived density proxies with the corresponding PDE cell-density fields using a user-specified coupling strength.
6. Carry the coupled fields and updated agents into the next macro step.

The other five fields (oxygen, nutrient, drug, immune-evasion and vasculature) remain continuum fields in the current implementation.

## Coupling equations

For phenotype k in cancer, immune or fibroblast:

```
rho_agent,k(x_i) = count_k(x_i) / agents_per_unit

u_k,new(x_i) =
    clip((1-alpha)*u_k,PDE(x_i) + alpha*rho_agent,k(x_i), 0, u_k,max)
```

where alpha is the coupling strength in [0,1]. At alpha=0 the cell fields remain purely continuum; at alpha=1 they are replaced by the grid-aggregated agent density proxies.

The occupancy scale has units of agents per grid node per dimensionless density unit. It is intentionally an explicit input because the current project has no registered histology or spatial assay from which to calibrate an agent-to-continuum conversion.

## Software interface

- agent_density_fields(population, grid_shape, agents_per_unit=...)
- couple_agents_to_pde(state, population, parameters, coupling_strength=..., agents_per_unit=...)
- coupled_step(state, population, parameters, dx, dt, ...)
- simulate_coupled(initial_state, initial_population, parameters, dx, dt, steps, ...)

The simulator returns time points, continuum-state snapshots and agent-population snapshots. Randomness is reproducible through the seed argument, which advances deterministically by macro-step index.

## Numerical and biological limitations

- Nearest-grid deposition is simple and reproducible but may introduce voxel-scale aliasing. A future version should compare nearest-grid, cloud-in-cell and kernel deposition.
- The blend is an operator-splitting bridge, not a mass-conserving derivation of the continuum equations from the agent model.
- The PDE and ABM use overlapping phenomenological processes; without a carefully derived conservation/consistency framework, blending can double-count births, deaths or recruitment.
- Cell-agent motility, division, death and recruitment rules are illustrative, not calibrated to patient measurements.
- The occupancy scale and coupling strength are not clinically meaningful parameter estimates.
- The PDE explicit diffusion timestep guard does not establish stability for the coupled ABM/PDE system. Timestep, grid and coupling-strength sensitivity tests are required.
- This module is a research/software prototype, not a patient simulator or tool for treatment decisions.

## Acceptance tests currently included

- 2D and 3D phenotype-specific grid deposition.
- Zero-coupling identity.
- Full-coupling replacement of cell-density fields while preserving other PDE fields.
- Deterministic one-step coupling under a fixed seed.
- Validation of invalid occupancy scales and coupling strengths.

## Next calibration and validation tasks

1. Add mass/number accounting across births, deaths, migration and deposition.
2. Compare grid deposition kernels and spatial convergence.
3. Run coupling-strength and occupancy-scale sensitivity analyses.
4. Test operator-splitting order and time-step refinement.
5. Validate cell counts and spatial distributions against registered tissue imaging before assigning biological units.
6. Compare results with an independent cell-based simulator after the biological rules and units are defined.
