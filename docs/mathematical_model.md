# Mathematical model: dimensionless 2D/3D prototype

## State vector

For spatial position `(x,y)` or `(x,y,z)` and dimensionless time `t`:

```
u = [C, I, F, O, N, D, E, V]
```

- `C`: cancer-cell density proxy
- `I`: immune-cell density proxy
- `F`: fibroblast / stromal-cell density proxy
- `O`: oxygen availability proxy
- `N`: nutrient availability proxy
- `D`: simulated drug-field concentration proxy
- `E`: immune-evasion activity proxy
- `V`: vascular-density proxy

The model has the coupled reaction-diffusion form:

```
du_j/dt = D_j * Laplacian(u_j) + R_j(u; theta)
```

A dimensionless input may be applied to the drug-field equation. Diffusion is constant and diagonal in the current prototype.

## Local reaction terms

Let `gO = O/(K_O+O)`, `gN = N/(K_N+N)`, and `gD = D/(EC50+D)`.

```
dC/dt = rC*C*(1-C/KC)*gO*gN - [kI/(1+rhoE*E)]*I*C - etaD*gD*C - deltaC*C

dI/dt = sI + alphaI*C/(hI+C) - muI*I - gammaF*F*I - gammaE*E*I

dF/dt = rF*C/(hF+C) - muF*F

dO/dt = -qO*C*O/(KO+O) + lambdaO*(Oref*V-O)

dN/dt = -qN*C*N/(KN+N) + lambdaN*(Nref-N)

dD/dt = qD(t) - kD*D

dE/dt = aE*C/(hE+C) - bE*E

dV/dt = aV*max(Otarget-O,0) - bV*(V-1)
```

These are simplified phenomenological terms, not a complete description of breast-cancer biology. Fibroblast and immune variables aggregate cell states and mechanisms. The current model does not explicitly represent ECM fibre mechanics, pharmacokinetics, immune exhaustion subtypes, or clonal evolution.

## Spatial discretization

The solver supports 2D slices and 3D Cartesian voxel grids. It uses:
- A uniform grid with spacing `dx`.
- A nearest-neighbour finite-difference Laplacian.
- Edge-value padding to approximate zero normal derivative (no-flux) boundaries.
- Explicit forward Euler time stepping.
- Post-step clipping to declared numerical bounds.

For `d` spatial dimensions and the d-dimensional nearest-neighbour Laplacian, the conservative explicit diffusion restriction is:

```
dt <= dx^2 / (2 * d * max_j(D_j))
```

This restriction controls the diffusion operator only. Reaction stiffness may impose a stricter limit. Every substantive experiment must include timestep refinement, grid refinement, and comparison with an independent numerical method. The current timestep-refinement utility compares dt and dt/2 on the same grid.

## Inputs and interventions

`qD(t)` is an abstract dimensionless input to the simulated drug field. It is not a pharmacological dose, exposure recommendation, or treatment schedule. The optimal-control module is a separate toy mean-field simulation. It must not be used to guide patient care.

## Calibration and validation requirements

Define units and scaling; document initial and boundary conditions; assess parameter identifiability; estimate uncertainty; compare grid and timestep refinements; compare an independent solver; evaluate on independent compatible cohorts; report failure cases. A survival association does not prove treatment benefit. Bulk molecular data cannot validate a spatial oxygen or cell-contact field.
