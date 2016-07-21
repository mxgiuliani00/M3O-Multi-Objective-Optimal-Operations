function [JJ_mpc, Ompc] = run_mpc(errorLevel)
%RUN_MPC evaluates an approximation of the Pareto optimal set using the
%        model predictive control scheme.
%
%   [JJ_mpc, Ompc] = RUN_MPC(errorLevel) returns the sequence of controls,
%   releases, storage state, maximum and minimum release, and step costs,
%   for the considered time horizon, as fields of the structure Ompc. The
%   approximated Pareto set is returned in matrix N-by-2 matrix JJ_mpc, for
%   N weight combinations and 2 objectives. Step costs in the optimization
%   problem are combined as a weighted sum defined by the weights in the
%   1-by-2 vector lambda. The prediction accuracy for disturbance in MPC is
%   affected by a random noise proportional to the input percentage value
%   errorLevel. Perfect prediction is assumed if errorLevel is set to 0.
%
% Output:
%     JJ_mpc - performance of multi-objective. 
%       Ompc - Structure with trajectory of optimal controls, releases,
%              storage states, maximum and minimum releases, step costs for
%              the whole time horizon.
%
% Input:
% 		errorLevel - random noise level affecting MPC prediction.
%
% See also SIM_MPC_TEST MPC_TEST 

% Copyright 2016 NRM group - Politecnico di Milano
%
% This file is part of Multi-Objective Optimal Operation (M3O) toolbox
%
% Multi-Objective Optimal Operation toolbox is free software: you can
% redistribute it and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
% 
% Multi-Objective Optimal Operation (M3O) toolbox  is distributed in the
% hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with Multi-Objective Optimal Operation. If not, see
% http://www.gnu.org/licenses/


global sys_param;

h_in  = sys_param.simulation.h_in;
e = sys_param.simulation.q;

% --- System Parameters
S0 = levelToStorage(h_in); % Initial storage

% --- Simulation and Optimization
% Run MPC
[ g_flo, g_irr, g, s, u, r ,v, V] = sim_mpc_Test( e, S0, errorLevel);

% Collect outputs
Ompc.u = u;
Ompc.r = r;
Ompc.s = s;
Ompc.v = v;
Ompc.V = V;
Ompc.g_flo = g_flo;
Ompc.g_irr = g_irr;

% Pareto front objectives
J_flo = mean(g_flo);
J_irr = mean(g_irr);

% Pareto front objectives
JJ_mpc(1) = J_flo;
JJ_mpc(2) = J_irr;

end