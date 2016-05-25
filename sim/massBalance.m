function [s1,r1] = massBalance( s, u, q )
% MASSBALANCE simulate the reservoir dynamics over 24 hours interval.
%
% v = MASSBALANCE( s, u, q )
%
% Decision and inflow are considered constant over all the simulation steps,
% although actual release is calculated dynamically following the lake
% evolution.
%
% Output:
%      s1 - final storage. ([m3])
%      r1 - average release over the 24 hours. ([m3/s])
%
% Input:
%       s - initial storage. ([m3])
%       u - release decision ([m3/s])
%       q - inflow rate ([m3/s])
%
% See also MIN_RELEASE, MAX_RELEASE

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



HH = 24;
delta = 3600;
s_ = nan(HH+1,1);
r_ = nan(HH+1,1);

s_(1) = s;
for i=1:HH
  qm = min_release(s_(i));
  qM = max_release(s_(i));
  r_(i+1) = min( qM , max( qm , u ) );
  s_(i+1) = s_(i) + delta*( q - r_(i+1) );
end

s1 = s_(HH);
r1 = mean(r_(2:end));
