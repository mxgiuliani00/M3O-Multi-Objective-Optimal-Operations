function [JJ, policy, err_perc] = run_iso(regressor)
% RUN_ISO evaluates an approximation of the Pareto optimal set resulting
% from Implicit Stochastic Optimization (ISO) (Labadie 2004)
%
% RUN_ISO first employs a Deterministic Dynamic Programming (DDP) to find
% optimal controls (i.e. release decisions) under several possible inflow
% realizations. Then it performs a regression based on the optimized data
% to retrieve a release policy. The regression method can be selected to
% best suit the given problem. The regressor already implemented is a
% piecewise linear interpolation.
%
%   [JJ, policy, err_perc] = run_iso(regressor)
%
% Outputs :
%         JJ - the approximated Pareto set: N-by-2 matrix for N weight
%              combinations and 2 objectives.
%     policy - the policy struct variable obtained from regression
%   err_perc - the regression percentage error on the performance
% 
% Inputs :
%   regressor - string containing the regression method:
%               'linear_spline' -> linear piecewise interpolator
% 
%
% See also OPT_DDP, SLMENGINE
%
% REFERENCES:
%     Labadie, J. W. (2004). Optimal operation of multireservoir systems:
%     state-of-the-art review. Journal of water resources planning and
%     management, 130(2), 93-111.

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

% -- Error check --
if strcmp(regressor,'linear_spline')
  sys_param.regressor = regressor ;
else
  error(['Regressor not defined. Please check or modify',...
    ' this function to use a different regression method']);
end

sys_param.algorithm.name = 'ddp';

% -- Initialization --
q    = sys_param.simulation.q;
h_in = sys_param.simulation.h_in;

% -- Run optimization to collect samples -- 
policy.H = opt_ddp( q );
[JJ_DDP(1), JJ_DDP(2), h, u, ~] = simLake( q, h_in, policy );

% -- Run regression to find the approximate decision rules --
X = h(2:end); % regressors
Y = u(2:end); % regressand

if strcmp(regressor,'linear_spline') 
  % Fit a piecewise linear function to data
  policy = slmengine(X,Y,'degree',1,'plot','off') ;
else
  error(['Regressor not defined. Please check or modify',...
    ' this function to use a different regression method']);
end

% -- Run simulation to collect performance --
sys_param.algorithm.name = 'iso';
[JJ(1), JJ(2)] = simLake( q, h_in, policy );

% compute error
err_perc = 100.*(abs(JJ_DDP - JJ)./JJ_DDP);
