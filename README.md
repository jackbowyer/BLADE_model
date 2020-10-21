### BLADE_model ###
Repository for code relating to the BLADE platform developed by Weinberg et al. 2017

BLADE_model.SBML is the SBML code describing the BLADE model
- All model variables (species) are listed first, followed by the model parameters.
- Each reaction that takes place within the BLADE system is listed individually.
- The naming convention for all species and parameters is consistent with the accompanying published paper.
- This SBML model was NOT used to produce any of the results in the paper, but it can be used in conjunction with the published parameter values to reproduce the results.

gillespie_algorithm_BLADE.m is the MATLAB code used to simulate Gillespie algorithm simulations of the BLADE model.
- See Gillespie 1977 for a full description of the Gillespie algorithm.
- Each of the 17 'cases' listed in the script corresponds to the 17 reactions that constitute the BLADE model. 
- When run with the inputs provided, the script will produce Gillespie simulations of the two circuits used as testing data in the accompanying paper.
- The adapted angular metric used to score circuit performance is calculated with respect to the corresponding ideal truth table outputs.
- This score and the error compared to experimental performance is printed.
