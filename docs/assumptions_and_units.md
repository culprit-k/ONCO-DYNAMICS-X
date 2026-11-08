# Assumptions, units and parameter policy

## Current prototype
All fields, parameters, space, time and inputs are **dimensionless**. Values are bounded numerical proxies. No coefficient has been estimated from patient data, a clinical trial or a validated biological experiment. A value of 0.5 does not mean 50% tumour cellularity, 50 mmHg oxygen, a specific cell count, or a drug concentration.

## Assumptions
- The domain is two-dimensional with uniform spatial coefficients.
- Diffusion uses a diagonal constant coefficient matrix.
- Boundaries are no-flux.
- Initial fields are synthetic and deterministic.
- Growth is reduced by oxygen and nutrient limitation.
- Immune activity is reduced by fibroblast and evasion proxies.
- Hypoxia drives a simple vascular-adaptation proxy.
- Drug exposure is a synthetic field driven by a user-defined abstract source.
- Post-step clipping prevents negative values but can alter mass balance and must be assessed in solver comparisons.

## Parameter policy
For each future calibrated parameter, record its mathematical definition, biological interpretation, dimensional unit and scaling; source paper/dataset/experiment; estimation method and prior assumptions; identifiability assessment and uncertainty; validated range and extrapolation limits. Do not copy prototype values into biological or clinical claims.

## Future dimensionalization
Specify length, time, density and concentration scales and transform each parameter before interpreting time or distance as days, hours, millimetres or micrometres.
