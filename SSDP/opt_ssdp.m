function H = opt_ssdp( ESP )
%OPT_SSDP Optimization scheme of sampling stochastic dynamic programming
%        (SSDP)
%
%  H = OPT_SSDP( ESP ) calculate the optimal Bellman value function
%  H_Bellman, given the sampled ensembles or historical ensemble streamflow
%  prediction (ESP). In SSDP method, the ESP trajectories are used to
%  estiamte transition probability of state variable due to the external
%  random disturbances, e.g., inflow. In this case, the performance of SSDP
%  highly depends on the quality of ESP.
%
%
% Output:
%       H - optimal value function. H is a three-dimensional matrix, with
%           first and second dimension correspond to the discretization
%           level for lake state and hydrological one (i.e., inflow). The
%           third dimension is total optimization horizon.
%
% Input:
%       ESP - ensemble streamflow prediction scenarios. ESP is a M-by-N
%             matrix, where M is the number of simulation time steps, and N
%             is the number of ensembles.
%
% See also REOPT_SSDP, BELLMAN_SSDP

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

%% Initilization
discr_q = sys_param.algorithm.discr_q;
discr_s = sys_param.algorithm.discr_s;
min_rel = sys_param.algorithm.min_rel;
max_rel = sys_param.algorithm.max_rel;

% Error check and pre-configuration

try
  nakeinterp1(1,1,1);
  sys_param.algorithm.interp_foo = @nakeinterp1;
catch
  warning(['Cannot use ''nakeinterp1'' function for fast interpolation.',...
    ' To enable the ''nakeinterp1'' please run ''mex -v nakeinterp1mx.c'' to compile it first.',...
    ' For now, I will use the MATLAB version for interpolation.']);
  sys_param.algorithm.interp_foo = @interp1qr;
end
interp_foo = sys_param.algorithm.interp_foo;

% Setup parameters

horizon_T = sys_param.algorithm.T;         % optimization time horizon
Pr_mode   = sys_param.algorithm.Pr_mode;   % mode for computing transition of ESP from one trace to the other

no_ensemble = size(ESP, 2);

n_s = length(discr_s);

% Pre-compute the state transition matrix

switch Pr_mode
  
  case  1 % deterministic case, i.e., no transition allowed across different scenarios
    p_ji = eye(no_ensemble);
    
  case  2 % equal transition probability across different scenarios
    p_ji = ones(no_ensemble) / no_ensemble;
    
  case  3 % applying the method in <J.R. Stedinger 2001>
    cycle_T    = sys_param.algorithm.cycle_T;
    forecast_T = sys_param.algorithm.forecast_T;
    
    p_ji = zeros(no_ensemble);
    [reg_param, error_std] = regreEstimate(cycle_T, forecast_T, ESP);
    
  otherwise
    error(['Undefined ''Pr_mode'' value for computing transition probability.',...
      ' Valid value is ''1'', ''2'',or ''3''.'])
    
end


% Initialize Bellman value matrix

H = zeros(n_s , no_ensemble, horizon_T + 1);
H(:, :, end) = sys_param.algorithm.Hend;

%% Backward recursive optimization

for t = horizon_T: -1: 1
  
  H_ = H(:, :, t+1) ;  % [ n_s x no_ensemble ]
  
  % Compute the p_ji before hand
  
  if Pr_mode == 3
    for l = 1: no_ensemble
      e_t = ESP(t, l);
      p_ji(:, l) = transitionPrOpt(reg_param, error_std, e_t);
    end
  end
  sys_param.algorithm.p_ji = p_ji;
  
  % Compute Bellman value function for each state variable
  
  for i = 1: n_s
    
    storage_curr = discr_s(i); % scalar
    
    for l = 1: no_ensemble
      
      e_series = ESP(t, l);
      
      sys_param.simulation.vv = interp_foo(discr_q(:), min_rel(i, :)', e_series);
      
      sys_param.simulation.VV = interp_foo(discr_q(:), max_rel(i, :)', e_series);
      
      H(i,l,t) = Bellman_ssdp(H_, storage_curr, ESP(t,:), l);
      
    end
  end
end


end

