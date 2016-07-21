README SSDP

Sampling Stochastic Dynamic Programming (SSDP) designs the optimal operating policy by using an ensemble for reservoir inflow trajectories, either resampled from historical observations or using the modern streamflow forecast, to characterize the stochasticity of the inflow process. This replaces the SDP characterization of the inflow using pdf.

DEFAULT SETTING:

- the discretized domains of state (relative storage wrt a reference level), control (release decision) are defined as follows: 
	- 55 values of storage (-7295000, 1118600000)
	- 500 values of release decision (1, 500)

- an ensemble of 25 members is used for describing the inflow

- in order to explore the Pareto front, different weights combinations are used.
Specifically: lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1]


SUGGESTIONS FOR MORE COMPLEX PROBLEMS:

To solve more complex applications we recommend adapting the discretization of state, and decision vectors. Finer discretization allows better performance but implies higher computational costs. Also the weights combinations used for the aggregation of the objectives can be changed, in order to explore the Pareto front. Finally, the larger and the more accurate the ensemble of inflow is, the better the characterization of their dynamics and, consequently, the quality of the designed policies will be.