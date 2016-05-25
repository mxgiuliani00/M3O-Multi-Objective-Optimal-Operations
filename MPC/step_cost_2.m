function G = step_cost_2( s, r, v, V)
%STEP_COST_2 calculates the aggregated step cost for a finite future 
%            horizon.
% 
%     G = STEP_COST_2( s, r, v, V) evaluates the aggregate step cost for
%     the finite MPC horizon P, and returns it in the 1-by-P vector G.
%     State, release, and maximum and minimum release trajectories for the
%     finite horizon P are given, respectively, by 1-by-P input vectors
%     s,r,v,V. Step costs aggregation is evaluated as a weighted sum
%     defined by the weights in the 1-by-2 vector lambda. Step costs are
%     normalized before aggregation. Penalties for the final state achieved
%     in time horizon P, and physical constraint violations are added to
%     the step-cost to be minimized.
%
% Output:
%     G - aggregated (weighted) step cost trajectory with penalties and
% 		    constraints violations.
%
% Input:
%       s - state (storage) trajectory.
%       r - release trajectory.
%       v - minimum release allowed trajectory.
%       V - maximum release allowed trajectory.
%
% See also RUN_MPC SIM_MPC_TEST MPC_TEXT

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

w    = sys_param.simulation.w;
hFLO = sys_param.simulation.hFLO;

weights = sys_param.algorithm.weights;

% compute step-costs:
H = length(s)-1;

[g_flo, g_irr] = deal(zeros(H,1)); % pre-allocate the memory

for i = 1: H
  ht = storageToLevel(s(i));
  [g_flo(i), g_irr(i)] = immediate_costs(ht, r(i));
end

% Cost for violation of constraint
g_MaR = max( r - V ,0);
g_MiR = max( v - r ,0);

% Penalties
g_floP = 100*max(storageToLevel(s(end)) - hFLO, 0);
g_irrP = max(w - max_release(s(end)), 0);

% aggregate step-costs:
G = g_flo.*weights(1) + g_irr.*weights(2) + ...
    g_floP.*weights(1) + g_irrP.*weights(2) + ...
    g_MaR(:) + g_MiR(:);
