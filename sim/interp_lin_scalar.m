function y = interp_lin_scalar( X , Y , x )
% INTERP_LIN_SCALAR linearly interpolate points X and Y to estimate y in x.
%
% y = INTERP_LIN_SCALAR( X , Y , x )
%
% Performs:
%
%            Y(k+1) + Y(k)
% y = Y(k) + ------------- ( x - X(k) ) 
%            X(k+1) - X(k)
%
% with 'k' such that X(k) <= x < X(k+1)
%
% Output:
%       y - the interpolated value.
% 
% Input:
%       X - vector of independent variables.
%       Y - vector of dependent variables.
%       x - independent variable value to interpolate on.
%
% See also INTERP1

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



% extreme cases
if x <= X( 1 )
	y = Y( 1 ) ; 
	return
end
if x >= X(end)
	y = Y(end) ;
	return
end

% otherwise

% Find index 'k' of subinterval [ X(k) , X(k+1) ] s.t. X(k) <= x < X(k+1)
[ ignore , i ] = min( abs( X - x ) ) ;

% If X( i ) = x     then   y = Y( i ) :
if X( i ) == x
	y = Y( i ) ;
	return
end

% Else :
% if X( i ) < x     then   k = i  
% if X( i ) > x     then   k = i - 1
k = i - ( X( i ) > x ) ;       
% Line joining points ( X(k) , Y(k) ) and ( X(k+1) , Y(k+1) ) 
Dy = Y( k + 1 ) - Y( k ) ;
Dx = X( k + 1 ) - X( k ) ;
m  = Dy / Dx             ; % slope
% Interpolate :
y = Y( k ) +  m * ( x - X( k ) ) ;