function [JJ, Q] = run_fqi(F, G, reg_param)
% RUN_FQI run Fitted Q-Iteration algorithm over a given sample dataset
%
% function [JJ, Q] = RUN_FQI(F, G, reg_param)
%
% The function designs the optimal operating policy for a pre-defined
% tradeoff, via Fitted-Q Iteration algorithm (Ernst et al. 2005;
% Castelletti et al. 2010), where the action-value function is approximated
% through Extremely Randomized Trees.
%
% Output:
%       JJ - objective function vector
%        Q - estimated action-value function
%
% Input:
%           F - sample dataset in the form < x(t), u(t), x(t+1) > obtained
%               from design of experiments (see RUN_DOE)
%           G - sample dataset of the immediate costs < g1, g2 > for each
%               tuple of F
%   reg_param - type of regressor used for approximating the
%               action-value function and associated parameters. Currently,
%               possible values are: 'ET': Extremely Randomized Trees, see
%               lib/RT for additional details.
%
% See also: TUPLES_SHUFFLING, READQ, RUN_DOE.
% 
% REFERENCES:
%     D. Ernst, P. Geurts, L. Wehenkel (2005). "Tree-based batch mode
%     reinforcement learning". Journal of Machine Learning Research 6 (1),
%     503-556. 
%     A. Castelletti, S. Galelli, M. Restelli, R. Soncini-Sessa. (2010).
%     "Tree-based reinforcement learning for optimal water reservoir
%     operation". Water Resources Research 46 (W09507).

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

% -- Initialization --
q    = sys_param.simulation.q;
h_in = sys_param.simulation.h_in;

weights = sys_param.algorithm.weights;
discr_s = sys_param.algorithm.discr_s;
discr_u = sys_param.algorithm.discr_u;
gamma   = sys_param.algorithm.gamma;

n_s = length(discr_s);
n_u = length(discr_u);

% Sample dataset - aggregated into single objective
GG = G * weights';
FG = [F, GG] ;

% Get configuration for selected regressor
if strcmp(reg_param.name, 'ET')
  
  % Shuffle the data (it improves ET regression)
  FF = tuples_shuffling(FG);
  
  % Extra-Trees parameters
  maxIter                        = reg_param.maxIter;         
  rtensparam                     = init_extra_trees();
  rtensparam.nbterms             = reg_param.M;
  rtensparam.rtparam.nmin        = reg_param.nmin;
  rtensparam.rtparam.extratreesk = size(F,2); % Number of random cuts
  
  % RUN FQI with ET
  % initialization of Q = 0 for each state-action pair
  Q_curr = zeros(n_s, n_u);
  ls = int32(1:size(FF,1));
  
  for i = 1: maxIter
    
    % Construction of training set: input = < x(t), u(t) >, output =
    % r(t+1) + gamma*Q(x(t+1),u(t))
    x  = FF(:,1);
    u  = FF(:,2);
    x1 = FF(:,3);
    
    H  = readQ(x1, discr_s, discr_u, Q_curr);
    Qh = FF(:,end) + gamma*H;
    
    tuple_in  = single( [x u] );
    tuple_out = single( Qh );
    
    % Regression to estimate new state-action function Q(x,u)
    [SS, UU] = meshgrid( discr_s, discr_u );
    SS = SS'; 
    UU = UU';
    XU = single( [ SS(:) UU(:) ] );
    
    Qout   = rtenslearn_c(tuple_in,tuple_out,ls,[],rtensparam,XU,0);
    Q_curr = double(reshape(Qout, n_s, n_u));
  end
  
  % Approximate Q function
  Q = Q_curr;
  policy.Q = Q;
  
else
  error(['Regressor not defined. Please check or modify'...
    ' this function to use a different regressor.']);
end

% -- RUN SIMULATION --
[Jflo, Jirr] = simLake( q, h_in, policy );
JJ(1) = Jflo;
JJ(2) = Jirr;

end

% This code has been written by Matteo Giuliani (matteo.giuliani@polimi.it)