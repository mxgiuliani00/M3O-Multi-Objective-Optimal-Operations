function [ H , idx_u ] = Bellman_sdp( H_ , s_curr )
%BELLMAN_SDP Update future value function and optimal decisions for
%            Stochastic Dynamic Programming
%
%  [ H , idx_u ] = BELLMAN_SDP( H_ , s_curr ) update the estimated future
%  value function given the previous one H_ and the current state s_curr.
%
% Output:
%       H - update future value function.
%   idx_u - index of optimal decision
%
% Input:
%       H_ - future value function from previous estimation
%   s_curr - value of state variable at current step
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
discr_s = sys_param.algorithm.discr_s;
discr_u = sys_param.algorithm.discr_u;
discr_q = sys_param.algorithm.discr_q;

n_u = length(sys_param.algorithm.discr_u);
n_q = length(sys_param.algorithm.discr_q);

weights = sys_param.algorithm.weights;
gamma   = sys_param.algorithm.gamma;


VV = repmat(sys_param.simulation.VV', 1 , n_u);
vv = repmat(sys_param.simulation.vv', 1 , n_u);

delta = sys_param.simulation.delta;

mi_q    = sys_param.algorithm.q_stat(1); % mean of log(disturbance)
sigma_q = sys_param.algorithm.q_stat(2); % std of log(disturbance)

%-- Calculate actual release contrained by min/max release rate --
discr_u = repmat(discr_u, n_q, 1);
R = min( VV , max( vv , discr_u ) );

%==========================================================================
% Calculate the state transition; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY

qq = repmat( discr_q', 1, n_u );
s_next = s_curr + delta*( qq - R );

%==========================================================================
% Compute immediate costs and aggregated one; TO BE ADAPTAED ACCORDING TO
% YOUR OWN CASE STUDY
[g1, g2] = immediate_costs( storageToLevel(s_next), R ) ;
G = g1*weights(1) + g2*weights(2);

%-- Compute cost-to-go given by Bellman function --
H_ = interp1qr( discr_s , H_ , s_next(:) ) ;
H_ = reshape( H_, n_q, n_u )         ;

%-- Compute resolution of Bellman value function --
% compute the probability of occourence of inflow that falls within the
% each bin of descritized inflow level
cdf_q      = logncdf( discr_q , mi_q , sigma_q );  
p_q        = diff(cdf_q);                          
p_diff_ini = 1-sum(p_q);                           
p_diff     = [ p_diff_ini ; p_q'];                 

Q     = (G + gamma.*H_)'*p_diff ;
H     = min(Q)                  ;
sens  = eps                     ;
idx_u = find( Q <= H + sens )   ;

end