README MPC 
MPC performs deterministic Model Predictive Control over the finite receding horizon defined by the variable sys_param.algorithm.P. The prediction accuracy for disturbance in MPC is affected by a random noise proportional to the input percentage value errorLevel.


DEFAULT SETTINGS:
- MPC is performed via Matlab GlobalSearch fmincon algorithm. Default settings for the algorithm are ‘interior-point’ algorithm, with maximum number of iterations equal to 2000 and maximum number of function evaluations equal to 1000. Fmincon option for the ‘Hessian’ is set to ‘lbfgs’ and the sub-problem algorithm is set to ‘ldl-factorisation’. Tolerances on objective function improvement and controls of fmincon algorithm are set equal to ‘eps’. These default settings can be customised in file mpc_Test.m.

- the random noise level affecting the accuracy of MPC prediction is set to 0 by default (i.e. perfect prediction).- in order to explore the Pareto front, different weights combinations are used.
Specifically: lambda = [1 0; .75 .25; .5 .5 ; .35 .65; .2 .8; .1 .9;  0 1]

SUGGESTIONS FOR MORE COMPLEX PROBLEMS:
To solve more complex applications we recommend increasing the number of function evaluations or iterations in the optimization algorithm. 
If real predictions of the disturbance are available (e.g., from a model or from agencies) those should be injected as MPC input. In order to test the effect of noise on disturbance prediction accuracy we recommend setting the errorLevel to values larger than 0.
Also the weights combinations used for the aggregation of the objectives can be changed, in order to explore the Pareto front.
