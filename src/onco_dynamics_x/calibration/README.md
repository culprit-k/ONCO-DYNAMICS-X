# Calibration and uncertainty

Implemented:
- Local central finite-difference parameter sensitivity for the final mean cancer-field summary.
- Monte Carlo propagation of user-specified parameter ranges through the deterministic PDE solver.
- Optional PyMC Gaussian-likelihood posterior construction for analyst-supplied, PyTensor-compatible forward models.

These tools currently operate on the illustrative dimensionless model. Monte Carlo parameter-range uncertainty is not Bayesian posterior uncertainty. Neither can be interpreted biologically until parameter ranges/priors are justified and observations support identifiability.
