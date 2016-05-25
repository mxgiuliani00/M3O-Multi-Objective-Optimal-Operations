function [H, u_opt_idx] = Bellman_ssdp( H_, storage_curr, ESP, i )
%BELLMAN_SSDP Solver for Bellman value function.
%
% [H, u_opt_idx] = BELLMAN_SSDP(H_, storage_curr, ESP, i) return the
% optimal Bellman value H, and the indice of the optimal decision  from
% discritized feasible decisions.
%
% Output:
%           H - optimal value function
%   u_opt_idx - the indice of optimal decision from discritized feasbile decisions
%
% Input:
%            H_ - future value function
%  storage_curr - current storage value
%           ESP - ensemble streamflow prediction scenarios
%             i - indice for current scenario trace
%
% See also REOPT_SSDP, OPT_SSDP
 
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

interp_foo = sys_param.algorithm.interp_foo;

% Initialize the parameters

weights = sys_param.algorithm.weights;
p_ji    = sys_param.algorithm.p_ji;
gamma   = sys_param.algorithm.gamma;
discr_u = sys_param.algorithm.discr_u;
discr_s = sys_param.algorithm.discr_s;

delta = sys_param.simulation.delta;
vv    = sys_param.simulation.vv;
VV    = sys_param.simulation.VV;

n_u = length(discr_u);


no_ensemble = length(ESP);
e_series    = ESP(i);

%==========================================================================
% Calculate the state transition; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUD

release_curr = min( VV, max( vv, discr_u) );
storage_next = storage_curr + delta * (e_series - release_curr);

%==========================================================================
% Compute immediate costs and aggregated one; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY

[g1, g2] = immediate_costs(...
  storageToLevel(storage_curr) .* ones(size(release_curr)), release_curr);
H_curr = ( g1*weights(1) + g2*weights(2) )';

% Compute future return; Linear interpolation is used

H_futur = zeros(n_u, no_ensemble); % (n_u, no_ensemble)
for j = 1: no_ensemble    
  H_futur(:,j) = interp_foo(discr_s(:) , H_(:,j), storage_next(:)); % ( n_u x 1)
end


% Find target decision under current state

H_sum = H_curr +  gamma * H_futur * p_ji(:, i); % [ n_u x 1 ]
[~, u_opt_idx] = min(H_sum);


% Update the Bellman value function

H = H_curr(u_opt_idx) + H_futur(u_opt_idx, i);