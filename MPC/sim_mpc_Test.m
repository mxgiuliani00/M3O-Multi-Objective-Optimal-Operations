function [ g_flo, g_irr, g, s, u, r, v, V] = sim_mpc_Test(e, S0, errorLevel)
%SIM_MPC_TEST run MPC simulation along the full time horizone.
% 
%     SIM_MPC_TEST(e, S0, errorLevel) simulates the dynamic system given as
%     input along the full time horizon, subject to optimal controls
%     provided by MPC optimization. Step cost trajectories for flooding and
%     irrigation, and their weighted stepcost aggregation trajectory are
%     returned, respectively, in 1-by-H vectors g_flo, g_irr, and g, being
%     H the length of the simulation horizon. System state trajectory,
%     controls trajectory, release trajectory, and minimum and maximum
%     release trajectories are returned in 1-by-H vectors s, u, r, v, V.
%     The initial condition for the system state is given by S0. Step costs
%     aggregation is evaluated as a weighted sum defined by the weights in
%     the 1-by-2 vector lambda. The prediction accuracy for disturbance e
%     in MPC is affected by a random noise proportional to the input
%     percentage value errorLevel. Perfect prediction is assumed if
%     errorLevel is set to 0.
%
% Output:
%     g_flo - flood step-cost.
%     g_irr - irrigation deficit step cost.
%         g - aggregated (weighted) step cost.
%         s - state (storage) trajectory.
%         u - control trajectory.
%         r - release trajectory.
%         v - minimum release allowed trajectory.
%         V - maximum release allowed trajectory.
%
% Input:
%           e - disturbance (inflow) trajectory to the reservoir.
%          S0 - initial condition (storage).
%  errorLevel - random noise level affecting MPC prediction.
%
% See also RUN_MPC MPC_TEST MIN_RELEASE MAX_RELEASE

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

% -- Initialization --
delta   = sys_param.simulation.delta;

p       = sys_param.algorithm.P;
weights = sys_param.algorithm.weights;

N = length(e);  % Length of the simulation

[u, r, s] = deal( nan(N+1,1) );
[v, V, g_flo, g_irr, g] = deal( zeros(N-p+1,1) );

% define initial conditions and the integration time-step
s(1) = S0;

% mpc input trajectory for simulation
e_sim = [nan; e];

% -- Simulation --
for t = 1: N-p+1
  disp(num2str(t))
  
  % prediction (b-a).*rand(1000,1) + a;
  if (t+p) <= N
    e_pred = e(t:t+(p-1))+errorLevel*(max(e)-min(e))*(2.*rand(p,1)-1);
  else
    e_pred = e(t:end)+errorLevel*(max(e)-min(e))*(2.*rand(length(e(t:end)),1)-1);
  end
  
  % Define the initial condition (i.e. current state s(t))
  s1_init = s(t);
  
  % Determine the trajectory of the optimal controls
  [x]  = mpc_Test(s1_init, e_pred, weights);
  u(t) = x(1);
  
  % Compute release and mass balance equation
  v(t) = min_release(s(t));
  V(t) = max_release(s(t));
  
  r(t+1) = u(t);
  s(t+1) = s(t) + delta * ( e_sim(t+1) - r(t+1) );
  ht = storageToLevel(s(t));
  
  [g_flo(t), g_irr(t)] = immediate_costs(ht, r(t+1));
  g(t) = g_flo(t)*weights(1) + g_irr(t)*weights(2);
end