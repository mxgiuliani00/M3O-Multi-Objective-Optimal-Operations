function [Jflo, Jirr, h, u, r, g_flo, g_irr] = simLake( q, h_in, policy )
% SIMLAKE Simulation routine that return the final performance of flooding
%         and irrigation, as well as the a number of time-series variables
%         related to the lake dynamics. 
%
% [Jflo, Jirr, h, u, r, g_flo, g_irr] = SIMLAKE( q, h_in, policy )
% 
% Output:
%       Jflo - final performance for flooding objective. ([cm])
%       Jirr - final performance for irrigation objective. ([m3/s])
%          h - time-series of lake level. ([cm])
%          u - time-series of release decision. ([m3/s])
%          r - time-series of actual release. ([m3/s])
%      g_flo - time-series of step cost of flooding. ([cm])
%      g_irr - time-series of step cost of irrigation. ([cm3/s])
% 
% Input:
%          q - inflow time series. ([m3/s])
%       h_in - initial lake level. ([m])
%     policy - structure variable containing the policy related parameters.
%
% See also LEVELTOSTORAGE, EXTRACTOR_REF, INTERP_LIN_SCALAR, HOURLY_INTEG,
%     STORAGETOLEVEL, IMMEDIATE_COSTS, STD_OPERATING_POLICY, BELLMAN_DDP,
%     BELLMAN_SDP, PREDICTWITHANENSEMBLE_R, READQ, REOPT_SSDP

%  Copyright 2016 NRM group - Politecnico di Milano
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

% Simulation setting
q_sim = [ nan; q ];
H = length(q_sim) - 1;

% Initialization
[h,s,r,u] = deal(nan(size(q_sim)));

% Start simulation
h(1) = h_in;
s(1) = levelToStorage(h(1));

for t = 1: H
  
  % Compute release decision
  
  switch (sys_param.algorithm.name)
    
    case 'rand'      
      r_min = sys_param.simulation.r_min;
      r_max = sys_param.simulation.r_max;
      
      uu = r_min + rand*(r_max - r_min);
      
    case 'doe'
      uf = rand;
      r_min = sys_param.simulation.r_min;
      r_max = sys_param.simulation.r_max;
      w = sys_param.simulation.w;
      
      if uf < 0.2
        uu = r_min + rand*(r_max*0.66 - r_min);
      elseif uf < 0.40
        uu = r_max*0.66 + rand*(r_max - r_max*0.66);
      else
        uu = w;
      end
      
    case 'ddp'      
      discr_s = sys_param.algorithm.discr_s;
      discr_q = sys_param.algorithm.discr_q;
      discr_u = sys_param.algorithm.discr_u;
      
      min_rel = sys_param.algorithm.min_rel;
      max_rel = sys_param.algorithm.max_rel;
      
      w = sys_param.simulation.w;
     
      [ ~ , idx_q ] = min( abs( discr_q - q_sim(t+1) ) );
      
      % Minimum and maximum release for current storage and inflow:
      sys_param.simulation.vv = interp1qr( discr_s , min_rel(: , idx_q) , s(t) );
      sys_param.simulation.VV = interp1qr( discr_s , max_rel(: , idx_q) , s(t) );
      
      [ ~ , idx_u ] = Bellman_ddp( policy.H(:,t+1) , s(t) , q_sim(t+1) );
      
      % Choose one decision value (idx_u can return multiple equivalent
      % decisions)
      uu = extractor_ref( idx_u , discr_u , w );
      
    case 'sdp'
      discr_s = sys_param.algorithm.discr_s;
      discr_q = sys_param.algorithm.discr_q;
      discr_u = sys_param.algorithm.discr_u;
      
      min_rel = sys_param.algorithm.min_rel;
      max_rel = sys_param.algorithm.max_rel;
      
      w = sys_param.simulation.w;
      
      % Minimum and maximum release for current storage and inflow:
      [ ~ , idx_q ] = min( abs( discr_q - q_sim(t+1) ) );
      
      v =interp1qr( discr_s , min_rel( : , idx_q ) , s(t) );
      sys_param.simulation.vv = repmat( v, 1, length(discr_q) );
      
      V = interp1qr( discr_s , max_rel( : , idx_q ) , s(t) );
      sys_param.simulation.VV = repmat( V, 1, length(discr_q) );
      [ ~, idx_u ] = Bellman_sdp( policy.H , s(t) );

      % Choose one decision value (idx_u can return multiple equivalent
      % decisions)
      uu = extractor_ref( idx_u , discr_u , w );
      
    case 'emodps'
      policy_class = sys_param.algorithm.policy_class;
      
      if strcmp(policy_class, 'stdOP')
        uu = std_operating_policy(h, policy);
      else
        error(['Policy class not defined.',...
          ' Please check or modify this function to use a different',...
          ' class of parameterized functions']);
      end
      
    case 'fqi'
      discr_s = sys_param.algorithm.discr_s;
      discr_u = sys_param.algorithm.discr_u;
      
      w = sys_param.simulation.w;
      
      [ ~, idx_u] = readQ(s(t), discr_s, discr_u, policy.Q);
      uu = extractor_ref( idx_u , discr_u , w );
      
    case 'ssdp'    
      interp_foo = sys_param.algorithm.interp_foo;
      T          = sys_param.algorithm.T;
      discr_u    = sys_param.algorithm.discr_u;
      discr_s    = sys_param.algorithm.discr_s;
      discr_q    = sys_param.algorithm.discr_q;      
      min_rel    = sys_param.algorithm.min_rel;
      max_rel    = sys_param.algorithm.max_rel;
      esp_sample = sys_param.algorithm.esp_sample;
      
      t_idx = mod(t-1, T) + 1;
      
      w = sys_param.simulation.w;
      
      [ ~, idx_q ] = min(abs( discr_q - q_sim(t+1) ));
      
      sys_param.simulation.vv = interp_foo( discr_s(:), min_rel( : , idx_q ) , s(t) );
      sys_param.simulation.VV = interp_foo( discr_s(:), max_rel( : , idx_q ) , s(t) );
      
      idx_u = reopt_ssdp( policy.H(:,:,t_idx+1), s(t), esp_sample(t_idx,:) );
      uu = extractor_ref( idx_u , discr_u, w );
      
    case 'iso'
      regressorName = sys_param.algorithm.regressorName;
      
      if strcmp(regressorName,'linear_spline')
        uu = slmeval(h(t), policy);
      else
        error(['Regressor not defined. Please check or modify',...
          ' this function to use a different regression method']);
      end
      
    otherwise
      uu = nan;
      
  end
  u(t) = uu;
  
  % Hourly integration of mass-balance equation
  [s(t+1), r(t+1)] = massBalance( s(t), u(t), q_sim(t+1) );
  h(t+1) = storageToLevel(s(t+1));
  
end

% Calculate objectives (daily average of immediate costs)
[g_flo, g_irr] = immediate_costs(h(2:end), r(2:end));

Jflo = mean(g_flo);
Jirr = mean(g_irr);
