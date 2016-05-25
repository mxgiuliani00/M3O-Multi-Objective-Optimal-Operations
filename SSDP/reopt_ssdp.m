function u_idx = reopt_ssdp( H_Bellman_futur, storage_curr, ESP )
% REOPT_SSDP Reoptimization of Bellman value using real-time streamflow
%            forecast to determine the actual optimal decisions at each
%            stage.
%
% u_idx = REOPT_SSDP( H_Bellman_futur, storage_curr, ESP ) update the
% Bellman value function at each decision stage, given the current forecast
% information which can be different from the one used in calculation of
% future value function
%
% Output:
%       u_idx - the indice of optimal decision from discritized feasbile 
%               decisions
%
% Input:
%    H_Bellman_futur - future value function
%       storage_curr - current storage value
%                ESP - ensemble streamflow prediction scenarios
%
% See also BELLMAN_SSDP, OPT_SSDP

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

Pr_mode    = sys_param.algorithm.Pr_mode;
gamma      = sys_param.algorithm.gamma;
weights    = sys_param.algorithm.weights;
interp_foo = sys_param.algorithm.interp_foo;
discr_s    = sys_param.algorithm.discr_s;
discr_u    = sys_param.algorithm.discr_u;

delta = sys_param.simulation.delta;
vv    = sys_param.simulation.vv;
VV    = sys_param.simulation.VV;

n_u = length(discr_u);

no_ensemble = size(ESP, 2);


% Compute the transitional probability

switch Pr_mode
  
  case  2 % Equal probability across scenarios
    p_ji = ones(no_ensemble, 1) / no_ensemble;
    
    %TODO improve the regreEstimate to allow different use of ESP/Hist data
    %   case  3 % applying the method in <J.R. Stedinger 2001>
    %     p_ji = zeros(1, no_ensemble) ;
    %     [reg_param, error_std] = regreEstimate(cycle_T, forecast_T, ESP) ;
    %
    %     for j = 1: no_ensemble
    %       e_t = ESP(1, j) ;
    %       p_ji(:, j) = transitionPrOpt(reg_param, error_std, e_t) ;
    %     end
    
  otherwise
    error('undefined mode for computing transition probability in reoptimization')
    
end

[H_curr, H_futur] = deal( zeros(n_u, no_ensemble) );

for j = 1: no_ensemble
  
  e_series = ESP(1, j); % forecast initiated at decision stage
  
  %==========================================================================
  % Calculate the state transition; TO BE ADAPTAED ACCORDING TO
  % YOUR OWN CASE STUDY
  
  release_curr = min( VV, max( vv, discr_u) ); 
  storage_next = storage_curr + delta * (e_series - release_curr);
  
  %==========================================================================
  % Compute immediate costs and aggregated one; TO BE ADAPTAED ACCORDING TO
  % YOUR OWN CASE STUDY
  
  [g1, g2] = immediate_costs(...
    storageToLevel(storage_curr) .* ones(size(release_curr)), release_curr);
  H_curr(:, j) = ( g1*weights(1) + g2*weights(2) )';
  
  % Compute future return 
  H_futur(:,j) = interp_foo(discr_s(:) , H_Bellman_futur(:,j), storage_next(:)); % ( n_u x 1)
  
end


% Find target decision

H_Bellman = ( H_curr + gamma * H_futur ) * p_ji; % n_u x 1
H_min = min(H_Bellman);

u_idx = find( H_Bellman <=  H_min + eps );
