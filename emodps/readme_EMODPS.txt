README EMODPS

EMODPS searches the optimal operating policy within a pre-specified family of parameterized functions. The policy parameters are determined via simulation-based optimization, with function evaluate_objective.m providing the interface between each alternative (parameters vector) and its simulated performance.

DEFAULT SETTINGS:

- A piece-wise linear function (i.e., standard operating policy, defined in std_operating_policy.m) is adopted as the parametrized class for the control policy. The standard operating policy is a three-piecewise linear function of the level with a constant, fixed central piece corresponding to the delivery target, in this case taken equal to the downstream water demand (sys_param.simulation.w).

- NSGAII evolutionary algorithm is used to optimize the policy parameters. The adopted setting includes 50 generations with a population of 40 individuals.

SUGGESTIONS FOR MORE COMPLEX PROBLEMS:

To solve more complex problems we recommend to increase the number of function evaluations and therefore the combination of generations and individuals. Alternative MOEAs can be used and we suggest using the Borg MOEA, which has been demonstrated to outperform other algorithms in solving complex EMODPS problems (see Zatarain et al. (2016), A diagnostic assessment of evolutionary algorithms for multi-objective surface water reservoir control, Advances in Water Resources).
The parameterized class for the control policy is case-study dependent. In case no prior knowledge on the policy shape is available, a universal approximator can be adopted, and for this purpose we recommend using Gaussian Radial Basis Functions.
