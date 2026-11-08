# Mathematical model: dimensionless prototype

## State vector
For spatial position `(x,y)` and dimensionless time `t`:
```
u = [C, I, F, O, N, D, E, V]
```
`C`: cancer-cell density proxy; `I`: immune-cell density; `F`: fibroblast/stromal density; `O`: oxygen availability; `N`: nutrient availability; `D`: simulated drug-field concentration; `E`: immune-evasion activity; `V`: vascular density.

The current engine implements `du_j/dt = D_j * Laplacian(u_j) + R_j(u; theta)`, with a simulated drug input added to the drug equation. Diffusion is constant and diagonal.

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

## Numerical method
The baseline uses a 2D uniform Cartesian grid, five-point finite-difference Laplacian, edge-value padding to approximate zero normal gradient, and explicit forward Euler updates. Post-step clipping enforces declared display bounds.

For the 2D five-point diffusion operator, the conservative explicit diffusion restriction used is `dt <= dx^2 / (4 * max_j(D_j))`. This controls diffusion stability only; reaction stiffness may impose a stricter limit. Run timestep refinement, grid refinement and an independent-solver comparison before interpreting results.

## Inputs and interventions
`qD(t)` is an abstract dimensionless input to the simulated drug field. It is not a pharmacological dose, exposure recommendation or treatment schedule. Future optimal-control experiments must remain mathematical simulations.

## Calibration and validation requirements
Define units and scaling; document initial and boundary conditions; assess parameter identifiability; estimate uncertainty; compare grid/timestep refinements; compare an independent solver; evaluate on independent compatible cohorts; report failure cases. A survival association does not prove treatment benefit. Bulk molecular data cannot validate a spatial oxygen or cell-contact field.
