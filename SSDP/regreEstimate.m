function [reg_param, error_std] = regreEstimate(T, T_forecast, esp_new)
%REGREESTIMATE Estimate the parameters from the given ensembles for
%               computing the conditional probability.
%
% [reg_param, error_std] = REGREESTIMATE(T, T_forecast, esp_new) calculates
% the regression parameters and the standard error given the historical
% forecast ensembles. In this case, the regressor is the flowrate forecast
% at time t, while the regressand is the total streamflow volume from t+1
% to t + T_forecast.It is assumed that the regression prediction and the
% associated error are used as the mean and standard deviation for
% computing the conditional probability. See eq. (13) in <J.R. Stedinger
% 2001>.
%
%
% Output:
%      reg_param - the regression parameters [slope and intercept]
%      error_std - the standard error from regression prediction
% 
% Input:
%            T - the length of each forecast trace, e.g., 365 days 
%   T_forecast - the forecast horizon used to estimate the total flow
%                volume
%      esp_new - the ensemble forecast used to estimate the regression
%                parameters
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


[horizon_tot, no_ensembles] = size(esp_new);

if ( T + T_forecast >= horizon_tot )
  error('T must be less than the total horizon of single trace minus T_forecast')
end

reg_param = zeros(2, no_ensembles);
error_std = zeros(no_ensembles, 1);

for i = 1: no_ensembles
  esp_trace = esp_new(:, i);
  
  X = esp_trace(1: T);  % X_t --> fQ_t
  Y = filter(ones(T_forecast, 1), 1, esp_trace); % rolling sum of flow rate over T_forecast
  
  regressand = Y(2:T+1); % Y_t+1
  regressor  = [ones(size(X)), X];
  
  reg_param(:, i) = regressor \ regressand; % least square estimation
  error_std(i) = std(regressand - regressor * reg_param(:,i));
end