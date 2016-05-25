function f = initialize_variables(N, M, V, min_range, max_range)
%
% function f = initialize_variables(N, M, V, min_range, max_range) 
%
% This function initializes the chromosomes. Each chromosome has the
% following at this stage
%       * set of decision variables
%       * objective function values
% 
% where,
% N - Population size
% M - Number of objective functions
% V - Number of decision variables
% min_range - A vector of decimal values which indicate the minimum value
% for each decision variable.
% max_range - Vector of maximum possible values for decision variables.

%  Copyright (c) 2009, Aravind Seshadri
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions are 
%  met:
%
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%      
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%  POSSIBILITY OF SUCH DAMAGE.

K = M + V        ;
% K is the total number of array elements. For ease of computation decision
% variables and objective functions are concatenated to form a single
% array. For crossover and mutation only the decision variables are used
% while for selection, only the objective variable are utilized.

% Initialize the decision variables based on the minimum and maximum
% possible values. V is the number of decision variable. A random
% number is picked between the minimum and maximum possible values for
% the each decision variable.
f = nan*ones(N,K);
min_range = min_range(:)' ;
max_range = max_range(:)' ;
f(:,1:V) = repmat( min_range, N, 1 ) + repmat( max_range - min_range, N, 1 ).*rand( N, V );

% For ease of computation and handling data the chromosome also has the
% value of the objective function concatenated at the end. The elements
% V + 1 to K has the objective function valued. 
% The function evaluate_objective takes one chromosome at a time,
% infact only the decision variables are passed to the function along
% with information about the number of objective functions which are
% processed and returns the value for the objective functions. These
% values are now stored at the end of the chromosome itself.
for i = 1 : N
    f(i,V+1:K) = evaluate_objective(f(i,1:V),M,V);
end
