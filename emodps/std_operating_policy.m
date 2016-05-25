function r = std_operating_policy(h, policy)
%STD_OPERATING_POLICY evaluates the release decision for the current level.
%
% function r = STD_OPERATING_POLICY(h, policy)
%
% Evaluates the release decision for the current level according to the
% Standard Operating Policy (Maass et al., 1962) and value contained in the
% parameter vector. 
% 
% For the demonstration, the Standard Operating Policy
% here is a three-piecewise linear function of the level with a constant,
% fixed central piece corresponding to the delivery target, in this case
% taken equal to ''sys_param.simulation.w''. User shall define their own
% policy class, e.g., non-linear function such as radial basis functions, 
% by following the same structure of this function.
%
% Output:
%       r - release decision. Same dimension as h.
%
% Input:
%        h - scalar or vector of reservoir levels.
%   policy - a struct with a field theta that is a vector of four parameter 
%            values. Specifically:
% 
%            policy.theta(1): level at which water supply edging is stopped
%            policy.theta(2): level at which release > demand for flood
%                             protection
%            policy.theta(3): angle of the water supply edging segment
%            policy.theta(4): angle of the flood protection segment
%
% See also RUN_EMODPS

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

w = sys_param.simulation.w;

% -- Get policy parameters --
h1 = policy.theta(1);
h2 = policy.theta(2);
m1 = policy.theta(3);
m2 = policy.theta(4);

% -- Construct the policy using piecewise linear functions --
L1 = w + tan(m1) * ( h - h1 );
L2 = w + tan(m2) * ( h - h2 );
r  = max( [ min( L1 , w ) ; L2 ] );

r( r < 0 ) = 0;

end