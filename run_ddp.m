function [JJ, HH] = run_ddp()
%RUN_DDP  Run Deterministic Dynamic programming algorithm using the
%         supplied test case
%
% [JJ, HH] = RUN_DDP()
% Deterministic dynamic programming (DDP) formulate the operating policy
% design problem as a sequential decision-making process. Specifically, the
% long-term cost of an operating policy is computed for each state by means
% of the value function:
%
%     H(x) = min G(x, u, e) + H(x1)
%
% where H() is the optimal cost-to-go function defined over a discrete grid
% of states and G() is the immediate (time separable) cost function
% associated to the transition from state x to state x1 under the decision
% u and the known disturbance value e.
%
% Output:
%       JJ - performance of multi-objective 
%       HH - optimized control policy
%
%
% See also OPT_DDP, SIMLAKE

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

%-- Initialization --
disturbance   = sys_param.simulation.q;
initial_state = sys_param.simulation.h_in;

%-- Run optimization --
policy.H = opt_ddp( disturbance );
HH = policy.H;

%-- Run simulation --
[Jflo, Jirr] = simLake( disturbance, initial_state, policy );
JJ(1) = Jflo;
JJ(2) = Jirr;

end
