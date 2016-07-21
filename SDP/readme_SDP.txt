README SDP

Stochastic Dynamic Programming (SDP) designs the optimal operating policy describing the reservoir inflow as a stochastic disturbance characterized by a log-normal pdf.

DEFAULT SETTING:

- the discretized domains of state (relative storage wrt a reference level), control (release decision), and disturbance (inflow) are defined as follows: 
	- 55 values of storage (-7295000, 1118600000)
	- 500 values of release decision (1, 500)
	- 500 values of inflow (1, 500)

- the lognormal distribution is estimated on the trajectory of inflow available, which is also used for the simulation-based evaluation of the policy.

- in order to explore the Pareto front, different weights combinations are used.
Specifically: lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1]


SUGGESTIONS FOR MORE COMPLEX PROBLEMS:

To solve more complex applications we recommend adapting the discretization of state, decision, and disturbance vectors. Finer discretization allows better performance but implies higher computational costs. Also the weights combinations used for the aggregation of the objectives can be changed, in order to explore the Pareto front. Finally, the longer the time-series of inflow used for the estimation of the associated pdf, the better the characterization of their dynamics and, consequently, the quality of the designed policies. It is worth noting that if the inflows are temporally correlated, it is suggested to estimate a dynamic model of the inflow (e.g., autoregressive model) and condition the policy on the new enlarged state, i.e. storage of the reservoir and inflow, describing the residuals of the inflow model as stochastic disturbances with an associated pdf. 