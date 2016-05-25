function [Qh, idx_u] = readQ(x1, grid_x, grid_u, Q)
%READQ interrogate the action-value function
%
%  [Qh, idx_u] = readQ(x1, grid_x, grid_u, Q)
%
% Output:
%      Qh - optimal value for a given state x1
%   idx_u - index of optimal decision
%
% Input:
%       x1 - value of state variable
%   grid_x - discretized domain of state variable
%   grid_u - discretized domain of decision variable
%        Q - action-value function

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

% Read action-value functions for state x1
n_u = length( grid_u );
n_s  = size(x1,1);

Q_ = nan(n_s, n_u);

for i = 1: n_u
  % Interpolate the value for any state that falls in between the grid of
  % descritized domain
  Q_(:,i) = interp1qr( grid_x , Q(:,i), x1 );
  
  % Extrapolation for any state value that is outside the descritized
  % domain
  idx = x1 >= grid_x(end);
  Q_(idx,i) = Q(end, i) ;
end

% Extract optimal value (min cost) and corresponding optimal decision
Qh = min(Q_, [], 2);

if (length(x1)>1)
  idx_u = nan;
else
  idx_u = find( Q_ <= Qh + eps ) ;
end

end

% This code has been written by Matteo Giuliani (matteo.giuliani@polimi.it)
