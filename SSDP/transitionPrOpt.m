function p_ji = transitionPrOpt(reg_param, error_std, Q_curr)
%TRANSITIONPROPT  calculate the transitional probability from one trace to
%                 another.
%
% p_ji = TRANSITIONPROPT(reg_param, error_std, Q_curr) calculate the
% transitional probability from trace i to trace j, given current flowrate
% Q_curr and the estimated regression parameters reg_param and error_std.
% The calculation is based on eq. (12) in <J.R. Stedinger 2001>.
%
% Output:
%        p_ji - transitional probability from trace i to trace j.
%
% Input:
%      reg_param - the regression parameters [slope and intercept]
%      error_std - the standard error from regression prediction
%         Q_curr - current flowrate
%
%

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


no_ensembles = length(error_std) ;

p_j  = 1 / no_ensembles ;

fcast_mean = [1, Q_curr] * reg_param ;

%TODO: check if the use of pdf() may underestimate the probability
pr = p_j * pdf( 'norm', Q_curr(ones(no_ensembles,1), :),...
  fcast_mean(:), error_std(:) );

p_ji = pr ./ sum(pr(:, ones(no_ensembles,1)))' ;