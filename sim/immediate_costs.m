function [g_flo, g_irr] = immediate_costs(ht, rt)
% IMMEDIATE_COSTS calculate the step costs for flooding and irrigation 
%                 objectives, respectively.
%
% [g_flo, g_irr] = IMMEDIATE_COSTS(ht, rt)
%
% Output:
%       g_flo - step cost for flood. ([cm])
%       g_irr - step cost for irrigation. ([m3/s])
% 
% Input:
%       ht - lake level at time step t. ([m])
%       rt - water release at time step t. ([m3/s])

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

hFLO = sys_param.simulation.hFLO;
w = sys_param.simulation.w;

g_flo = max( ht - hFLO, 0 )*100; % exceedance of water level above the flooding threshold
g_irr = max( w-rt, 0 ); % deficit between the demand and supply