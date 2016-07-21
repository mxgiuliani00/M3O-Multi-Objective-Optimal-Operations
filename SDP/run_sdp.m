function [JJ, HH] = run_sdp()
%RUN_SDP  Run Stochastic Dynamic Programming algorithm using the supplied
% test case
%
% [JJ, HH] = RUN_SDP()
% Stochastic dynamic programming (DDP) formulate the operating policy design
% problem as a sequential decision-making process. Specifically, the expected 
% long-term cost of an operating policy is computed for each state by means of
% the value function:
%
%     H(x) = min E[G(x, u, e) + gamma*H(x1)]
%
% where H() is the optimal cost-to-go function defined over a discrete grid of
% states and G() is the immediate (time separable) cost function associated to
% the transition from state x to state x1 under the decision u and E[] 
% represents the expected value over the stochastic disturbances e. 
%
% Output:
%       JJ - performance of multi-objective 
%       HH - optimized control policy
%
%
% See also OPT_SDP, SIMLAKE

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
q    = sys_param.simulation.q;
h_in = sys_param.simulation.h_in;

tol = -1;    % accuracy level for termination 
max_it = 10; % maximum iteration for termination 

%-- Run optimization --
policy.H = opt_sdp( tol, max_it );
HH = policy.H;

%-- Run simulation --
[Jflo, Jirr] = simLake( q, h_in, policy );
JJ(1) = Jflo;
JJ(2) = Jirr;

end