function Xsf = tuples_shuffling(X)
%TUPLES_SHUFFLING shuffles a dataset of tuples.
%
%  Xsf = tuples_shuffling(X)
%
% Output:
%      Xsf - the input dataset (X) after random shuffling performed
%            on the rows of the matrix
%
% Input:
%       X - a dataset of tuples (observations on the rows and
%           variables on the columns)

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


% Get dimensions of the input matrix
r = size(X, 1);

% Random shuffling of the tuples (i.e., by row index)
random_idx = randperm(r);
Xsf = X(random_idx, :);

end

% This code has been written by Matteo Giuliani (matteo.giuliani@polimi.it)