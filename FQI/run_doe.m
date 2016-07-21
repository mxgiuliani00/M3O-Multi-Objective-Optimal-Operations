function [F, G] = run_doe(No_exp)
% RUN_DOE run design of experiments to construct sample dataset for batch-mode Reinforcement Learning algorithms (i.e. FQI)
%
% function [F, G] = run_doe(No_exp)
%
% The function runs a design of experiments phase for constructing two
% sample datasets of tuples of the form F = < x(t), u(t), x(t+1) > and G =
% < g1(t+1), g2(t+1) > which will be used to learn the optimal operating
% policy via bacth mode Reinforcement Learning algorithms (i.e. FQI). DOE will be
% performed by running 'Nexp' simulations with random decisions.
%
% Output:
%       F - sample dataset in the form < x(t), u(t), x(t+1) >
%       G - sample dataset of the immediate costs < g1(t+1), g2(t+1) >
%
% Input:
%       No_exp - number of experiments (simulations) to perform
%
% Copyright 2016 NRM group - Politecnico di Milano
%
% This file is part of Multi-Objective Optimal Operation (M3O)
%
%     Multi-Objective Optimal Operation is free software: you can redistribute
%     it and/or modify it under the terms of the GNU General Public License as
%     published by the Free Software Foundation, either version 3 of the
%     License, or (at your option) any later version.
%
%     Multi-Objective Optimal Operation (M3O) is distributed in the hope that
%     it will be useful, but WITHOUT ANY WARRANTY; without even the implied
%     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License along
%     with Multi-Objective Optimal Operation.  If not, see http://www.gnu.org/licenses/


% Run simulations with pseudo-random decisions (irregular sampling of
% decision space to ensure good coverage of the state-action-value space)

global sys_param;

q    = sys_param.simulation.q;
h_in = sys_param.simulation.h_in;

n_q = length(q);

[h, u] = deal( nan(n_q+1, No_exp) );
[g_flo, g_irr] = deal( nan(n_q, No_exp) );

for i = 1: No_exp
  [~, ~, h(:,i), u(:,i), ~, g_flo(:,i), g_irr(:,i)] = simLake( q, h_in, [] );
end

% Construct sample datasets
xt  = levelToStorage(h(1:end-1,:));
ut  = u(1:end-1,:);
xt1 = levelToStorage(h(2:end,:));
F   = [ xt(:), ut(:), xt1(:) ];
G   = [ g_flo(:), g_irr(:) ];

end

% This code has been written by Matteo Giuliani (matteo.giuliani@polimi.it)