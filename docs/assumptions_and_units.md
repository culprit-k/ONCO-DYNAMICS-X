# Assumptions, units and parameter policy

## Current prototype
All fields, parameters, space, time and inputs are **dimensionless**. Values are bounded numerical proxies. No coefficient has been estimated from patient data, a clinical trial or a validated biological experiment. A value of 0.5 does not mean 50% tumour cellularity, 50 mmHg oxygen, a specific cell count, or a drug concentration.

The mathematical solver supports 2D and 3D Cartesian grids. The current agent-based model is 2D only. The dashboard's 3D mode uses the 3D PDE solver; the 2D mode can additionally show a surface plot of a 2D field.

## Assumptions
- Spatial coefficients are uniform and diffusion is diagonal and constant.
- Boundaries are no-flux.
- Initial fields are synthetic and deterministic.
- Growth is reduced by oxygen and nutrient limitation.
- Immune activity is reduced by fibroblast and evasion proxies.
- Hypoxia drives a simple vascular-adaptation proxy.
- Drug exposure is a synthetic field driven by a user-defined abstract source.
- The agent model uses stochastic division/death and simplified migration/recruitment rules.
- Post-step clipping prevents negative values but can alter mass balance and must be assessed in solver comparisons.

## Parameter policy
For each future calibrated parameter, record its mathematical definition, biological interpretation, dimensional unit and scaling; source paper/dataset/experiment; estimation method and prior assumptions; identifiability assessment and uncertainty; validated range and extrapolation limits. Do not copy prototype values into biological or clinical claims.

## Future dimensionalization
Specify length, time, density and concentration scales and transform each parameter before interpreting time or distance as days, hours, millimetres or micrometres.
