function H = opt_ddp( disturbance )
%OPT_DDP Run deterministic dynamic programming on the given trajectory of
%        disturbances
%
%  H = OPT_DDP( disturbance ) calculate the optimal Bellman value function
%  H, given observed trajectory of inflow.
%
% Output:
%       H - optimal value function.
%
% Input:
%       disturbance - observed trajector of disturbance of any length.
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
discr_q = sys_param.algorithm.discr_q;

min_rel = sys_param.algorithm.min_rel;
max_rel = sys_param.algorithm.max_rel;
Hend = sys_param.algorithm.Hend;

N   = length(disturbance);
n_s = length(discr_s);


%-- Backward recursive optimization -- 
% note that 2 loops are run, with the first used for computing a penalty
% over the final state used in the second loop
H = zeros(n_s, N+1); % create Bellman
H(:,end) = Hend;     % initialize Bellman to penalty function
for t = N: -1: 1
  for i = 1: n_s; % for each discretization level of state 
    % calculate the min/max storage-discharge relationship. See README.md
    % file for the description of 'min_rel' and 'max_rel' ;
    sys_param.simulation.vv = interp1qr( discr_q , min_rel( i, : ) , disturbance(t) );
    sys_param.simulation.VV = interp1qr( discr_q , max_rel( i, : ) , disturbance(t) );
    
    H(i,t) = Bellman_ddp( H(:,t+1), discr_s(i), disturbance(t) );
  end
end
Hend = H(:,end);
H = zeros(n_s, N+1); % create Bellman
H(:,end) = Hend;     % initialize Bellman to penalty function
for t = N: -1: 1
  for i = 1: n_s; % for each discretization level of state 
    % calculate the min/max storage-discharge relationship. See README.md
    % file for the description of 'min_rel' and 'max_rel' ;
    sys_param.simulation.vv = interp1qr( discr_q , min_rel( i, : ) , disturbance(t) );
    sys_param.simulation.VV = interp1qr( discr_q , max_rel( i, : ) , disturbance(t) );
    
    H(i,t) = Bellman_ddp( H(:,t+1), discr_s(i), disturbance(t) );
  end
end

end
