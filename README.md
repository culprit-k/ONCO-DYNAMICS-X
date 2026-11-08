# Cellular simulation

The NumPy agent-based implementation supports both 2D and 3D spatial domains. It includes cancer-agent division/death, local oxygen/nutrient dependence, immune recruitment and gradient-directed migration, fibroblast migration influenced by local stromal density, reproducible random seeds, reflecting boundaries and a population cap.

The agent rules are illustrative and are not calibrated biological kinetics. The PDE fields act as the local environment, but the ABM is not yet integrated into a fully coupled multiscale co-solver with cell secretion/consumption feeding back to every PDE time step.

Future work: phenotype-specific immune states, contact rules, matrix mechanics, explicit migration constraints, clonal evolution, full PDE-ABM coupling and validation against spatially resolved observations. PhysiCell remains a candidate external integration engine.
