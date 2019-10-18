### M3O-Multi-Objective-Optimal-Operations
M3O is a Matlab toolbox for designing the optimal operations of multipurpose water reservoir systems. M3O allows users to design Pareto optimal (or approximate) operating policies for managing water reservoir systems through several alternative state-of-the-art methods. Version 1.0 of M3O includes Deterministic and Stochastic Dynamic Programming, Implicit Stochastic Optimization, Sampling Stochastic Dynamic Programming, fitted Q-iteration, Evolutionary Multi-Objective Direct Policy Search, and Model Predictive Control. The toolbox is designed to be accessible to practitioners, researchers, and students, and to provide a fully commented and customizable code for more experienced users.


### List of methods available
	- Deterministic Dynamic Programming (DDP);
	- Stochastic Dynamic Programming (SDP);
	- Implicit Stochastic Optimization (ISO);
	- Sampling Stochastic Dynamic Programming (SSDP);
	- Evolutionary Multi-Objective Direct Policy Search (EMODPS);
	- Fitted Q-Iteration (FQI);
	- Model Predictive Control (MPC);


### How to handle `<filename>.mex` files in MATLAB
M3O uses some `<filename>.mex` files built from the associated `<filename>.c` files (i.e., `nakeinterp1mx.c` and `rtenslearn_c.c`) which are used to speed up the computation efficiency.

To build up the `.mex` files from `<filename>.c`, run `mex -v <filename>.c` in MATLAB. If MATLAB cannot find a valid compiler, then you may need to install one on your machine by following this [link](http://it.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html).


### `Simlake` test case and `sys_param` 
M3O comes with a test-case using synthesized data to demonstrate the value of each algorithms. To run the test case, simple run `main_test.m` in MATLAB. A global structure variable `sys_param` is shared by the inner functions of the toolbox via the `global` command. Under the `sys_param` there are two subfields, namely the '.simulation' which holds generic settings for running the dynamic simulation of supplied case study, and `.algorithm` which holds the specific algorithmic settings for each algorithm. For more details see the [HELP FILE](http://www.nrm.deib.polimi.it/wp-content/uploads/2016/07/main_script_demo.html).

### Libraries
M3O relies on external code in the `lib` directory, including the [rtree-c code](http://www.montefiore.ulg.ac.be/~geurts/Software.html) by P. Geurts for the implementation of the Extremely Randomized Trees (Extra-Trees), SLMtools by John D'Errico, and NSGAII implementation by Aravind Seshadri. 

### References
Giuliani, M., Y. Li, A. Cominola, S. Denaro, E. Mason, and A. Castelletti (2016), [A Matlab toolbox for designing Multi-Objective Optimal Operations of water reservoir systems](https://www.sciencedirect.com/science/article/pii/S1364815216305345), Environmental Modelling & Software, 85, 293–298

----
### Copyright:

Copyright 2016 NRM group - Politecnico di Milano

Developers: Matteo Giuliani, Yu Li, Andrea Cominola, Simona Denaro, Emanuele Mason, Andrea Castelletti, Francesca Pianosi, Daniela Anghileri, Stefano Galelli.

M3O is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

The code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with M3O.  If not, see <http://www.gnu.org/licenses/licenses.en.html>.
