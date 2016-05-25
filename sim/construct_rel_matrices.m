function [vv, VV] = construct_rel_matrices(discr)
% CONSTRUCT_REL_MATRICES constructs the matrices of minimum and maximum
%                        release to be used in dynamic programming
%                        optimization scheme.
%
% [vv, VV] = CONSTRUCT_REL_MATRICES( discr )
%
% Decision and inflow are considered constant over all the simulation steps,
% although actual release is calculated dynamically following the lake 
% evolution.
%
% Output:
%       vv - minimal release for all combination of storage and inflow.
%       VV - maximal release for all combination of storage and inflow.
% 
% Input:
%    discr - structure variable containing the pre-defined descritized
%            storage, inflow and decision space.
%
% See also HOURLY_INTEG

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



discr_s = discr.discr_s;
discr_q = discr.discr_q;

% min release for each value of storage and inflow, assuming the release
% decision is equal to 0
vv = nan(length(discr_s), length(discr_q));
for i = 1:length(discr_s)
  for j = 1:length(discr_q)
    [~, r1] = massBalance( discr_s(i), 0, discr_q(j) );
    vv(i,j) = r1;
  end
end

% min release for each value of storage and inflow, assuming the release
% decision is equal to the maximum
VV = nan(length(discr_s), length(discr_q));
for i = 1: length(discr_s)
  for j = 1: length(discr_q)
    [~, r1] = massBalance( discr_s(i), discr.discr_u(end), discr_q(j) );
    VV(i,j) = r1;
  end
end