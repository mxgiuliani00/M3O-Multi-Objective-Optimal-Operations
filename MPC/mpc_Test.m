function [x] = mpc_Test(s1_init, e_pred, lambda)
% MPC_TEST solves the optimization problem for the finite MPC optimization
% orizon P under disturbance prediction e_pred, and returns the 1-by-P
% optimal control vector x. The constrained optimization problem is solved
% through Matlab GlobalSearch and fmincon.
% The initial state of the system is set equal to s1_init. Step costs
% aggregation in the optimization problem is evaluated as a weighted sum
% defined by the weights in the 1-by-2 vector lambda.
%
% Output:
%       x - optimal control provided by MPC optimization.
%
% Input:
%      s1_init - state of the system at the beginning of the finite
%                prediction horizon.
%       e_pred - disturbance (inflow) to be predicted by MPC.
%   errorLevel - random noise level affecting MPC prediction.
%
% See also RUN_MPC SIM_MPC_TEST

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
w     = sys_param.simulation.w;
delta = sys_param.simulation.delta;

% define the length of the prediction horizon, over which the optimal
% decision vector x must be optimized
e_pred = [nan ; e_pred];
H      = length(e_pred);

% initialize vectors
[r, s, v, V] = deal(nan(H,1));

% define initial conditions and the integration time-step
s(1) = s1_init;

% minimum and maximum value of the control variables
min_u = sys_param.simulation.r_min;
max_u = sys_param.simulation.r_max;

% -- Optimization of the decision(s) vector --
% define constraints
A = [eye(H-1); -eye(H-1)];
b = [max_u*ones(H-1,1); -min_u*ones(H-1,1)];

% Determine the initialization vector
xo = ones(H-1,1).*w;

% Define Options
options = optimset;
options.Algorithm = 'interior-point';
options.GradObj = 'off';
options.Hessian = 'lbfgs';
options.Display = 'off';
options.MaxIter = 2000;
options.MaxFunEvals = 1000;
options.TolFun = eps;
options.TolX = eps;
options.SubproblemAlgorithm = 'ldl-factorization';

% Optimization
% Multi-start version
% problem = createOptimProblem('fmincon','objective',@optfun,'x0',xo,'Aineq',A, 'bineq',b,'options',options);
% ms = MultiStart('UseParallel',false);
% x = run(ms,problem,5);

problem = createOptimProblem('fmincon','objective',@optfun,'x0',xo,'Aineq',A, 'bineq',b,'options',options);
ms = GlobalSearch;
x  = run(ms,problem);

%[x] = fmincon(@optfun,xo,A,b,[],[],[],[],[],options);

% Nested function that computes the objective function
  function f = optfun(x)
    
    % simulation
    for t = 1: H-1      
      % compute release and mass balance equation
      v(t) = min_release(s(t));
      V(t) = max_release(s(t));
      r(t+1) = x(t);
      s(t+1) = s(t) + delta * ( e_pred(t+1) - r(t+1) ) ;
    end
    
    % compute the step-costs over the simulation horizon and the aggregated
    % cost, which correspond to the objective function
    g = step_cost_2( s(1:end),r(2:end), v(1:end-1), V(1:end-1));
    f = mean(g);
  end

end
