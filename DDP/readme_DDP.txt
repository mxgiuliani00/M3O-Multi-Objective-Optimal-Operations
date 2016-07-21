README DDP

Deterministic Dynamic Programming (DDP) designs the optimal operating policy assuming a perfect knowledge on the entire trajectory of reservoir inflow.

DEFAULT SETTING:

- the discretized domains of state (relative storage wrt a reference level), control (release decision), and disturbance (inflow) are defined as follows: 
	- 55 values of storage (-7295000, 1118600000)
	- 500 values of release decision (1, 500)
	- 500 values of inflow (1, 500)

- in order to explore the Pareto front, different weights combinations are used.
Specifically: lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1]


SUGGESTIONS FOR MORE COMPLEX PROBLEMS:

To solve more complex applications we recommend adapting the discretization of state, decision, and disturbance vectors. Finer discretization allows better performance but implies higher computational costs. Also the weights combinations used for the aggregation of the objectives can be changed, in order to explore the Pareto front.

