README run_iso.m

run_iso.m performs an Implicit Stochastic Optimization (ISO) to design the reservoir operating policy by regressing the state variables (levels) towards the decision variables (releases). To collect the samples for building the regression models, run_iso.m relies on a set of deterministic optimizations obtained running opt_ddp.m before hand.

DEFAULT SETTINGS:

- in order to explore the Pareto front, different weights combinations are used.
Specifically: lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1]

- the implemented regressor is a piecewise linear interpolation performed with slmengine.m function setting 'degree' as 1.

SUGGESTIONS FOR MORE COMPLEX PROBLEMS:

To adapt the code to a different case-study we recommend to tune the weights combinations in order to efficiently explore the Pareto front during the deterministic optimization.
Different regressors can be used, provided that they are defined in run_iso.m.
For more complex problems user should increase the number of samples over which the regression model is built and/or apply a multiple regressor.
Also the weights combinations used for the aggregation of the objectives can be changed, in order to explore the Pareto front.