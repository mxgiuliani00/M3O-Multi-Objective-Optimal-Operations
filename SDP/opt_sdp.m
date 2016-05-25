function H = opt_sdp( tol, max_iter )
%OPT_SDP Run stochastic dynamic programming.
%
%  H = OPT_SDP( tol, max_iter  ) calculate the optimal Bellman value
%  function H usign stochatic dynamic programming (SDP). In SDP method, the
%  external disturbance is stochastic and can be described by a certain
%  probability density function (pdf). The Bellman value function is
%  calculated as : 
% 
%     H(x) = min E[G(x, u, e) + gamma*H(x1)]
%
% where the new value function H(x) is updated as the optimal expected sum
% of immediate cost G(x,u,e) and discounted future value gamma*H(x1). x, u,
% e, and x1 correspond to the current state, decision variable, external
% stochastic disturbance and the next state, respectively. H(x) will be
% updated repeatly until maximum iteration 'max_iter' or accuracy level
% 'tol' is reached.
%
% Output:
%       H - optimal value function.
%
% Input:
%        tol - tolerence value for termination test.
%   max_iter - maximal iteration of Bellman value evaluations.
%
% See also OPT_DDP, BELLMAN_SDP, BELLMAN_DDP

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
T       = sys_param.algorithm.T;
n_s     = length(discr_s);

min_rel = sys_param.algorithm.min_rel;
max_rel = sys_param.algorithm.max_rel;

H = zeros( n_s , T );

%-- Backward recursive optimization --
diff_H = Inf;
count  = 1;

while diff_H >= tol
  H_ = H;
  
  for i = 1: n_s
    sys_param.simulation.vv = min_rel( i, : );
    sys_param.simulation.VV = max_rel( i, : );
    
    H(i) = Bellman_sdp( H_, discr_s(i) );
  end
  
  diff_H = max( abs( H_ - H ));
  count  = count+1;
  
  if count > max_iter;
    return;
  end
  
end


end