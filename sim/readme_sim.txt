README SIM

This is the folder which holds functions describing the characteristics of the target problem, and therefore it should be re-adapted to users' own case study.

You should find a list of following files, namely:

- `construct_rel_matrices.m`: this function builds the matrices of minimum and maximum release to be used in dynamic programming optimization scheme.
- `extractor_ref.m`: function that retrieves a single optimal release decision value and the corresponding index from a set of potential good decisions, by choosing the one that is closest to the delivery target defined by sys_param.simulation.w.
- `immediate_costs.m`: this function calculates the step costs for flooding and irrigation objectives.
- `interp_lin_scalar.m`: given a set of points (X,Y) and a new independent variable x, this function interpolates the values of y in x.
- `levelToStorage.m`: function that converts the reservoir level into storage. Cylinder shape is assumed for the reservoir in this case.
- `storageToLevel.m`: function that converts the reservoir storage into the corresponding level. Cylinder shape is assumed for the reservoir in this case.
- `massBalance.m`: function used to integrate the reservoir dynamics over 24 hours interval (the simulation step).
- `max_release.m`: function that calculates the maximum daily release from the reservoir. The maximum release value is used to constrain the daily decision variable as the upper bound.
- `min_release.m`: function that calculates the minimum daily release from the reservoir. The minimum release value is used to constrain the daily decision variable as the lower bound.
- `simLake.m`: the main function that simulates the performance and system dynamics given the policy optimized by a specific method of M3O toolbox.


HOW TO CREATE A NEW CASE STUDY

In order to connect the M3O toolbox with your own case study, it is crucial that the user correctly describes his/her problem regarding to the dynamics of the system, the state variable, the objectives (and step cost for DP family) as well as constraints.
Specifically:
- to change the objectives formulation, modify immediate_costs.m and extractor_ref.m
- for a different reservoir shape, modify levelToStorage.m and storageToLevel.m
- to change the integration step for the reservoir dynamics, modify massBalance.m
