function f = evaluate_objective(x, M, V)
%EVALUATE_OBJECTIVE evaluate the objective functions values
%
% function f = EVALUATE_OBJECTIVE(x, M, V)
%
% Function to evaluate the objective functions for the given input vector
% x. x is an array of decision variables and f(1), f(2), etc are the
% objective functions. The algorithm always minimizes the objective
% function hence if you would like to maximize the function then multiply
% the function by negative one. M is the numebr of objective functions and
% V is the number of decision variables.
%
% This function has basically to be written by the user who defines his/her
% own objective function. Make sure that the M and V matches your initial
% user input.
%
% Outputs:
%     f - array containing the performance of each objective ;
%
% Inputs:
%     x - parameters to be optimzied. The length of x shall equal to V.
%     M - the number of objective
%     V - the number of input variable
%
% See also RUN_EMODPS

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

x = x(1:V);
x = x(:);
if ( length(x) ~= V )
  error(['The number of decision variables does not match' ,...
    ' you previous input. Kindly check your objective function'])
end

% --------------------------------------
% insert here your function f = f(x):
global sys_param;

q    = sys_param.simulation.q;
h_in = sys_param.simulation.h_in;

% -- Get policy parameters --
policy.theta = x;

% -- Run simulation to collect the results --
[ J1, J2 ] = simLake( q, h_in, policy );
f = [J1, J2];
% --------------------------------------

% Check for error
if ( length(f) ~= M )
  error(['Incorrect number of outptu objectives. Expecting to solve' ,...
    ' %d objectives formulation. Please check your objective function again'],...
    M );
end