function v = min_release(s)
% MIN_RELEASE calculate the minimum release rate [m3/s]
%
% v = MIN_RELEASE(s)
%
% Output:
%       v - minimum release rate. ([m3/s])
% 
% Input:
%       s - value of storage. ([m3])
%
% See also MAX_RELEASE, STORAGETOLEVEL

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

if (h <= 1.25)
  q = 0.0;
else
  q = 33.37*(h + 2.5).^2.015;
end

v = q;