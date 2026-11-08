# Causal analysis

An optional DoWhy wrapper is implemented in observational.py. It requires declared common causes, checks required variables and missingness, and estimates a linear backdoor effect. It does not discover or validate the causal graph. Analysts must defend temporal ordering, confounding control, selection assumptions, consistency and positivity before interpreting any result causally.
