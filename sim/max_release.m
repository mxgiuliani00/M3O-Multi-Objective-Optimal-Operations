function V = max_release(s)
% MAX_RELEASE calculate the maximum release rate [m3/s]
%
% v = MAX_RELEASE(s)
%
% Output:
%       v - maximum release rate. ([m3/s])
% 
% Input:
%       s - value of storage. ([m3])
%
% See also MIN_RELEASE, STORAGETOLEVEL

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



h = storageToLevel(s);

if (h <= -0.5)
  q = 0.0;
elseif (h <= -0.40)
  q = 1488.1*h + 744.05;
else
  q = 33.37*(h + 2.5).^2.015;
end
V = q;