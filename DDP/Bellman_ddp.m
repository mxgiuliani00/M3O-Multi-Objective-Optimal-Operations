function [ H , idx_u ] = Bellman_ddp( H_ , s_curr , q_curr )
%BELLMAN_DDP Update future value function and optimal decisions for
%            Deterministic Dynamic Programming
%
%  [ H , idx_u ] = BELLMAN_DDP( H_ , s_curr , q_curr ) update the estimated
%  future value function given the previous one H_, the current state s_curr
%  and the current value of disturbance.
%
% Output:
%       H - update future value function.
%   idx_u - index of optimal decision
%
% Input:
%          H_ - future value function from previous estimation
%      s_curr - value of state variable at current stage
%      q_curr - observed trajectory of the disturbance at current stage
%
% See also OPT_SDP, BELLMAN_SDP, BELLMAN_DDP

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
discr_s = sys_param.algorithm.discr_s ;
discr_u = sys_param.algorithm.discr_u ;
weights = sys_param.algorithm.weights ;

VV = sys_param.simulation.VV; % max release
vv = sys_param.simulation.vv; % min release

delta = sys_param.simulation.delta;

%-- Calculate actual release contrained by min/max release rate --
R = min( VV , max( vv , discr_u ) ) ;

%==========================================================================
% Calculate the state transition; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY
s_next = s_curr + delta* ( q_curr - R );

%==========================================================================
% Compute immediate costs and aggregated one; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY
[g1, g2] = immediate_costs( storageToLevel(s_next), R ) ;
G = (g1*weights(1) + g2*weights(2))';

%-- Compute cost-to-go given by Bellman function --
% apply linear interpolation to update the Bellman value H_
H_ = interp1qr( discr_s , H_ , s_next );

%-- Compute resolution of Bellman value function --
Q     = G + H_;
H     = min(Q);
sens  = eps;
idx_u = find( Q <= H + sens );



end
