# Cellular simulation

The initial NumPy agent-based implementation is in agent_based.py. It includes cancer-agent division/death, local oxygen/nutrient dependence, immune-agent recruitment and chemotaxis proxy, fibroblast movement affected by local stromal density, and reproducible random seeds. Its rules are illustrative and not calibrated biological kinetics.

Future work: phenotype-specific states, cell-cell contact rules, matrix mechanics, explicit migration constraints, coupling to the PDE at every substep, and validation against spatially resolved observations. PhysiCell remains a candidate integration engine.
