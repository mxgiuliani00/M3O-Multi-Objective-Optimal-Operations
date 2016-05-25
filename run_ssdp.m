function [JJ, HH] = run_ssdp()
%RUN_SSDP Run SSDP algorithm using the supplied test case 
%
% [JJ, HH] = RUN_SSDP() performs the sampling stochastic dynamic
% programming, and return the array of objective values ''JJ'' and the
% corresponding value function ''HH''.
%
% Output:
%       JJ - performance of multi-objective 
%       HH - optimized control policy
%
%
% see also REOPT_SSDP, BELLMAN_SSDP, REOPT_SSDP

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

%-- Initialization --
q    = sys_param.simulation.q;
h_in = sys_param.simulation.h_in;

esp_sample = sys_param.algorithm.esp_sample;

%-- Run optimization --
policy.H = opt_ssdp( esp_sample );
HH = policy.H;

%-- Run simulation --
[Jflo, Jirr] = simLake( q, h_in, policy );
JJ(1) = Jflo;
JJ(2) = Jirr;

end
