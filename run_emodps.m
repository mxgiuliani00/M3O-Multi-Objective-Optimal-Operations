function [JJ, PS] = run_emodps(policy_class, moea_param)
%RUN_EMODPS calculate an approximate pareto front using Evolutionary
%           Multi-Objective Direct Policy Search.
%
% function [JJ, PS] = RUN_EMODPS(policy_class, moea)
%
% Designs a Pareto approximate set of operating policies via Evolutionary
% Multi-Objective Direct Policy Search (Giuliani et al., 2016), where the
% policy is assumed to be of a certain function faimly 'policy_class' (e.g,
% piecewise linear) and shall be parameterized accordingly, where the
% state-of the art multi-objective evolutionary algorithm (i.e.,
% NSGAII) can be used to find the optimal parameter sets that best describe 
% the Pareto Optimal solutions. In principle, the performance of EMODPS is
% subjected to the 'shape' of policy class and the MOEA algorithm. See 
% < M. Giuliani,et al 2016 > for more information.
%
% Output:
%       JJ - values of the objective in the best Pareto front approximation
%       PS - values of the parameter of the policies in the Pareto front
%
% Input:
%  policy_class - type of policy. Currently, possible values are:
%                'stdOP': Standard Operating Policy, see
%                 std_operating_policy.m for references and other 
%                 information.
%          moea - configuration of evolutionary algorithm. Currently, 
%                 possible values are:
% 
%             		'NSGAII': Non Sorting Genetic Algorithm II, see
%                           lib/NSGA2/nsga_2.m for additional details.
%
% See also EVALUATE_OBJECTIVE
%
% REFERENCES:
% M. Giuliani, A. Castelletti, F. Pianosi, E. Mason, and P. M. Reed (2016). Curses,
% tradeoffs, and scalable management: advancing evolutionary multi-objective
% direct policy search to improve water reservoir operations. Journal of water
% resources planning and management, vol. 142, iss. 2.

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

% -- POLICY SETTING --
if strcmp(policy_class,'stdOP')
  sys_param.algorithm.policy_class = policy_class ;
else
  error(['Policy class not defined. Please check or modify',...
    ' this function to use a different class of parameterized functions']);
end

% -- OPTIMIZATION SETTING --
M = 2 ;     % objective number
V = 4 ;     % decision variables number (policy parameters)
min_range = [ -.5 -.5 0 0 ] ;
max_range = [ 1.25 1.25 pi/2 pi/2 ] ;

% -- RUN OPTIMIZATION --
if strcmp(moea_param.name, 'NSGAII')
  % The algorithm returns the initial (ch0) and the final population (chF)
  [ch0, chF] = nsga_2(moea_param.pop,moea_param.gen,M,V,min_range,max_range);
  
  % Save Pareto approximate set (PS) and front (JJ)
  PS = chF(:, 1:V);
  JJ = chF(:, V+1:V+M);
else
  error(['MOEA not defined. Please check or modify this function',...
    'to run a different optimization algorithm.']);
end

end

% This code has been written by Matteo Giuliani (matteo.giuliani@polimi.it)
