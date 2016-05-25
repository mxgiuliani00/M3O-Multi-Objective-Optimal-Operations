function slm = hermite2slm(harray)
% hermite2slm: converts an array of coefficients in Hermite form to a slm
% usage: slm = hermite2slm(harray)
%
% arguments: (input)
%  harray - nx2 or nx3 array of coefficients
%           harray(:,1) == knots
%           harray(:,2) == function values
%
%           and if they are present:
%           harray(:,3) == first derivatives
%
%           If additional columns, they define the
%           higher order derivatives of a (2*k-3)
%           degree hermite interpolant.
%
%           The knots (or breaks) must be an increasing
%           sequence of numbers, or the spline will be
%           a problem to evaluate in an SLM form.
%
% arguments: (output)
%  slm - slm struct, usable by slmeval and plotslm
%
%
% See also: slmset, slmengine, slmeval, ppval, slmfit, plotslm
%
%
% Author: John D'Errico
% E-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/6/07

[n,k] = size(harray);

% check that the knots are sorted in increasing order
% if not, then the spline will fail eventually.
if any(diff(harray(:,1)) <= 0)
  error('HERMITE2SLM:nonmonotonic','Knots were non-monotonic')
end

% create the basic struct
slm.form = 'slm';
slm.knots = harray(:,1);

% we don't know how it was created
slm.prescription = [];
% nor what it was derived from
slm.x = [];
slm.y = [];

if k == 2
  % linear Hermite
  slm.degree = 1;
  slm.coef = harray(:,2);
  
elseif k == 3
  % cubic Hermite
  slm.degree = 3;
  slm.coef = harray(:,2:3);
  
elseif k > 3
  % quintic Hermite or higher
  slm.degree = 2*k-3;
  slm.coef = harray(:,2:end);

else
  error 'harray is not in my Hermite standard form'
end

