function y = slmeval(x,slm,evalmode)
% slmeval: evaluates a Hermite function, its derivatives, or its inverse
% usage 1: y = slmeval(x,slm)           % evaluates y = f(x)
% usage 2: y = slmeval(x,slm,evalmode)  % general form
%
% As opposed to extrapolation as I am, slmeval will not extrapolate
% unless explicitly told to do so. If you wanted to extrapolate,
% then you really should have built the spline differently, with knots
% that extend to the limits as far as you will expect to predict.
% This way, the user will be able to control tha behavior of the curve
% in the extrapolation region.
%
% Extrapolation MAY be done however by SLMEVAL, if the extrapolation
% parameter was set properly by either slmset or slmengine. There
% are several choices. See slmset for more description.
% 
% arguments: (input)
%  x       - array (any shape) of data to be evaluated through model
%
%            All points are assumed to lie inside the range of
%            the knots. I simply don't ever recommend extrapolating
%            a spline beyond the limits over which it was built.
%
%            However, IF you insist on extrapolating the spline,
%            the default behavior of SLMEVAL is to use the value
%            of the spline at the first or last know as is appropriate.
%            This is an implicit constant extrapolation.
%
%            Other extrapolation types are (grudgingly) allowed, but
%            are controlled by setting the Extrapolation field in
%            the model. It is best to simply build the spline with a
%            desired shape using knots that go out as far as you wish
%            using SLMENGINE in the first place. Then you have complete
%            control of the shape.
%
%  slm     - a shape language model structure, normally constructed
%            by slmfit or slmengine.
%
%            The fields in this struct are:
%            slm.type  = either 'cubic', 'linear', or 'constant'
%
%            slm.knots = a vector of knots (distinct & increasing)
%               There must be at least two knots.
%
%            slm.coef  = an array of coefficients for the hermite
%               function.
%
%            If the function is linear or constant Hermite, then
%            slm.coef will have only one column, composed of the
%            value of the function at each corresponding knot.
%            The 'constant' option uses a function value at x(i)
%            to apply for x(i) <= x < x(i+1).
%
%            If the function is cubic Hermite, then slm.coef
%            will have two columns. The first column will be the
%            value of the function at the corresponding knot,
%            the second column will be the corresponding first
%            derivative.
%
% evalmode - (OPTIONAL) numeric flag - specifies what evaluation
%            is to be done.
%
%            DEFAULT VALUE: 0
%
%            == 0 --> evaluate the function at each point in x.
%            == 1 --> the first derivative at each point in x.
%            == 2 --> the second derivative at each point in x.
%            == 3 --> the third derivative at each point in x.
%            == -1 --> evaluate the inverse of the function at
%             each point in x, thus y is returned such that x=f(y)
%
%            Note 1: Piecewise constant functions will return zero
%            for all order derivatives, since I ignore the delta
%            functions at each knot.
%            The inverse operation is also disabled for constant
%            functions.
%
%            Note 2: Linear hermite functions will return zero
%            for second and higher order derivatives. At a knot
%            point, while technically the derivative is undefined
%            at that location, the slope of the segment to the
%            right of that knot is returned.
%
%            Note 3: Inverse computations will return the
%            LEFTMOST zero (closest to -inf) in the event that
%            more than one solution exists.
%
%            Note 4: Inverse of points which fall above the
%            maximum or below the minimum value of the function
%            will be returned as a NaN. NO extrapolation will
%            ever be done for the inverse.
%
% Arguments: (output)
%  y       - Evaluated result, the predicted value, i.e., f(x)
%            or f'(x), f''(x), f'''(x), or the functional inverse
%            such that x = f(y). y will have the same shape and
%            size as x.
%
% Example:
%




% See also: slmpar, ppval, slmengine, slmset, slm2sym
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 6/10/09

% get the shape of x, so we can restore it at the end
sizex = size(x);
x=x(:);
nx = length(x);

% check for the minimum required information
if ~isstruct(slm) || ~isfield(slm,'knots') || ...
    ~isfield(slm,'coef') || ~isfield(slm,'degree')
  error('slm must be a struct, with "knots" "coef" and "degree" fields')
end

% if an older spline is evaluated, then the Extrapolation field may not
% have been defined. Use the default here, of 'constant'
if ~isfield(slm,'Extrapolation')
  slm.Extrapolation = 'constant';
end

extraps = {'error','warning','constant','linear','cubic','nan'};
extraptype = find(strncmpi(slm.Extrapolation,extraps,length(slm.Extrapolation)));
if numel(extraptype) < 1
  error('SLMEVAL:invalidextraptype', ...
    'Extrapolation field is not one of the valid options. see slmset')
elseif numel(extraptype) > 1
  error('SLMEVAL:invalidextraptype', ...
    'Extrapolation field is ambiguous. see slmset')
end

knots = slm.knots(:);
coef = slm.coef;
nk = length(knots);
dx = diff(knots);

% default evalmode is 0, i.e., evaluate f(x)
if (nargin<3) || isempty(evalmode)
  evalmode = 0;
end

% what extrapolation class was indicated, and is any necessary?
if (evalmode >= 0) 
  extrapx = (x < knots(1)) | (x > knots(end));
  if any(extrapx(:))
    switch extraptype
      case 1 % error
        % just kick the user out with an error
        error('SLMEVAL:extrapolationdisallowed', ...
          'Extrapolation has been disallowed by the user, and one or more points fell outside of the knots')
        
      case 2 % warning
        % issue a warning as requested
        warning('SLMEVAL:extrapolation', ...
          'One or more points fell outside of the knots. constant extrapolation performed')
        % then perform constant extrpolation. this is most simply
        % accomplished by shifting the x values to the corresponding
        % end knots of the spline.
        x = min(max(x,knots(1)),knots(end));
        
      case 3 % constant
        % just shift x to the min or max knot. no warning issued
        x = min(max(x,knots(1)),knots(end));
        
      case 6 % NaN
        % stuff the prediction with Nan where extrapolation would have been
        % attempted.
        x(extrapx) = NaN;
    end
    % extrapolation types 4 and 5 will be accomplished by the code
  end
end

% which knot interval does each point fall in? This
% question only applies to forward evaluation, thus
% for evalmode >= 0.
% Any needed interval clipping has already been done
% above, so if a point falls outside the knot interval,
% it will be extrapolated.
if evalmode >= 0
  % use histc to bin the points.
  [junk,xbin] = histc(x,knots);

  % any point which falls at the top end, is said to
  % be in the last bin.
  xbin(xbin==nk)=nk-1;
  
  % catch any points which fell outside the bins. Note
  % that we have already clipped the points with min/max
  % unless extraptype was 4 or 5.
  % by considering a point to lie in the first/last bin,
  % EVEN if it falls below/above the first/last knot, we
  % can thus easily allow the function to be extrapolated.
  xbin(x < knots(1)) = 1;
  xbin(x > knots(end)) = nk - 1;
 
else
  % we need to treat the inverse carefully, and it will
  % be specific to the spline order. So just create xbin
  % for now as zeros.
  xbin = zeros(nx,1);
end

% evaluate the appropriate class of piecewise hermite function
switch slm.degree
  % we don't do anything special for extrapolation mode here.
  % everything has already been fixed in advance.
  case 1
    % verify that coef is the right size
    sc=size(coef);
    if (sc(1) ~= nk)
      error('Improper size of coef field for these knots')
    elseif (sc(2)~=1)
      error('Improper size of coef field for a linear Hermite')
    end
    
    % Evaluation mode:
    switch evalmode
      case 0
        % f(x)
        t = (x-knots(xbin))./dx(xbin);
        y = coef(xbin) + (coef(xbin+1)-coef(xbin)).*t;

      case 1
        % first derivative
        y = (coef(xbin+1)-coef(xbin))./dx(xbin);

      case {2 3}
        % higher derivatives of a piecewise linear function
        % are all zero
        y = zeros(sizex);

      case -1
        % functional inverse
        % here the problem of which bin we fall into is not
        % trivial, since the function need not be monotone.
        % we might also have constant segments.
        % Identify the leftmost interval which contains the
        % point in question.

        % first, clip the data in terms of function value
        miny = min(coef);
        maxy = max(coef);
        y = repmat(NaN,nx,1);
        k = find((x<=maxy) & (x>=miny));

        % determine which knot interval to look in for
        % each point.
        ybin = nmbs(x(k),coef);

        % rule out the divide by zero cases for intervals
        % where the function was constant
        L = (coef(ybin+1) == coef(ybin));
        if any(L)
          y(k(L)) = knots(ybin(L));
          ybin(L)=[];
          k(L)=[];
        end

        % inverse interpolation
        y(k) = knots(ybin) + (x(k) - coef(ybin)).* ...
          (knots(ybin+1)-knots(ybin))./(coef(ybin+1)-coef(ybin));

      otherwise
        % anything else
        error('Evalmode must be one of [0 1 2 3 -1]')
    end

  case 3
    % verify shape of coef
    sc=size(coef);
    if (sc(1) ~= nk)
      error('Improper size of coef field for these knots')
    elseif (sc(2)~=2)
      error('Improper size of coef field for a cubic Hermite')
    end

    % extrapolation mode need only be dealt with here if
    % it was set at 4 or 5. All other extrapolation cases
    % have already been treated. And for mode 5, nothing
    % special needs to be done, as in that case, we just evaluate
    % the corresponding polynomial, even for a point that falls
    % below/above the first/last bin.
    
    % Evaluation mode:
    switch evalmode
      case 0
        % f(x)
        t = (x-knots(xbin))./dx(xbin);
        t2 = t.^2;
        t3 = t.^3;
        s2 = (1-t).^2;
        s3 = (1-t).^3;
        y = (-coef(xbin,2).*(s3-s2) + ...
          coef(xbin+1,2).*(t3-t2)).*dx(xbin) + ...
          coef(xbin,1).*(3*s2-2*s3) + ...
          coef(xbin+1,1).*(3*t2-2*t3);
        
        % catch the points that fell outside, IF extraptype
        % indicates a linear extrapolation in those areas.
        if extraptype == 4
          k = (x < knots(1));
          y(k) = coef(1,1) + coef(1,2)*(x(k) - knots(1));
          
          k = (x > knots(end));
          y(k) = coef(end,1) + coef(end,2)*(x(k) - knots(end));
        end
        
      case 1
        % first derivative for the cubic case
        t = (x-knots(xbin))./dx(xbin);
        t2 = t.^2;
        s = 1-t;
        s2 = (1-t).^2;
        y = -coef(xbin,2).*(-3*s2+2*s) + ...
          coef(xbin+1,2).*(3*t2-2*t) + ...
          (coef(xbin,1).*(-6*s+6*s2) + ...
          coef(xbin+1,1).*(6*t-6*t2))./dx(xbin);
        
        % catch the points that fell outside, IF extraptype
        % indicates a linear extrapolation in those areas.
        if extraptype == 4
          k = (x < knots(1));
          y(k) = coef(1,2);
          
          k = (x > knots(end));
          y(k) = coef(end,2);
        end

      case 2
        % second derivative of a cubic
        t = (x-knots(xbin))./dx(xbin);
        s = 1-t;
        y = (-coef(xbin,2).*(6*s - 2) + ...
          coef(xbin+1,2).*(6*t - 2))./dx(xbin) + ...
          (coef(xbin,1).*(6 - 12*s) + ...
          coef(xbin+1,1).*(6 - 12*t))./(dx(xbin).^2);
        
        % catch the points that fell outside, IF extraptype
        % indicates a linear extrapolation in those areas.
        if extraptype == 4
          k = (x < knots(1)) | (x > knots(end));
          y(k) = 0;
        end

      case 3
        % third derivative
        y = 6*(coef(xbin,2) + coef(xbin+1,2))./(dx(xbin).^2) + ...
          12*(coef(xbin,1) - coef(xbin+1,1))./(dx(xbin).^3);
        
        % catch the points that fell outside, IF extraptype
        % indicates a linear extrapolation in those areas.
        if extraptype == 4
          k = (x < knots(1)) | (x > knots(end));
          y(k) = 0;
        end
        
      case -1
        % functional inverse
        % here the problem of which bin we fall into is not
        % trivial, since the function need not be monotone.
        % we might also have constant segments.
        
        % first, convert the spline into a pp form.
        pp = slm2pp(slm);
        coefs = pp.coefs;
        
        % scale the cubic polys so they live on [0,1].
        % this will stabilize things a bit, as well as make
        % the tests easier later on. I could do this in one
        % line with a variety of tricks, but feeling lazy...
        coefs(:,1) = coefs(:,1).*(dx.^3);
        coefs(:,2) = coefs(:,2).*(dx.^2);
        coefs(:,3) = coefs(:,3).*dx;
        
        % Identify the leftmost interval which contains the
        % point in question. The problem is the cubic segments
        % might not be monotone. nmbs will try though.
        binedges = slm.coef(:,1);
        ybin = nmbs(x,binedges);
        
        % a simple solution is to test every interval up to
        % and including the interval found in ybin, but no
        % further. We can stop there since we know a solution
        % must exist in that interval. If no interval was found,
        % then ybin will be a NaN. In that event, we must search
        % every interval, since the curve need not be monotone.
        ybin(isnan(ybin) | (ybin == nk)) = nk - 1;
        
        % just a loop over roots here. not terribly efficient,
        % but I don't terribly want to code a vectorized cubic
        % solver for this problem.
        y = NaN(size(x));
        tol = 1000*eps;
        for j = 1:nx
          flag = true;
          i = 1;
          while flag && (i <= ybin(j))
            Ci = coefs(i,:);
            
            % offset the constant term. a root of this
            % cubic is what we want.
            Ci(4) = Ci(4) - x(j);
            
            % scale the poly coefficients to improve things
            % yet some more.
            Ci = Ci./max(abs(Ci(1:3)));
            
            % get the roots. at least one must be real.
            Ri = roots(Ci);
            
            % of the real roots, is one of them strictly
            % inside [0,1]? If more than one is, take the
            % smallest root.
            k = (abs(imag(Ri)) == 0) & (real(Ri) >= 0) & (real(Ri) <= 1);
            if any(k)
              k = find(k,1,'first');
              y(j) = knots(i) + dx(i)*real(Ri(k));
              flag = false;
              continue
            end
            
            % if we did not succeed in the last test,
            % of the real roots, is one of them within 
            % inside [-tol,1+tol]? If more than one is,
            % take the smallest root.
            k = (abs(imag(Ri)) == 0) & (real(Ri) >= -tol) & (real(Ri) <= (1 + tol));
            if any(k)
              k = find(k,1,'first');
              y(j) = knots(i) + dx(i)*real(Ri(k));
              flag = false;
              continue
            end
            
            % increment i, check the next interval.
            i = i + 1;
          end
        end
        
      otherwise
        % anything else
        error('Evalmode must be one of [0 1 2 3 -1]')
    end
  case 0
    % piecewise constant function (discontinuous at the knots)
    % verify that coef is the right size
    sc=size(coef);
    if (sc(1) ~= (nk-1))
      error('Improper size of coef field for these knots')
    elseif (sc(2)~=1)
      error('Improper size of coef field for a piecewise constant function')
    end
    
    % note that the extrapolation mode is essentially ignored
    % for piecewise constant functions. xbin is already fixed correctly.
    
    % Evaluation mode:
    switch evalmode
      case 0
        % f(x)
        y = coef(xbin);
        
      case {1 2 3}
        % all derivatives of a piecewise constant function will
        % be zero if we ignore the delta function singularities
        % at each knot.
        y = zeros(sizex);
        
      case -1
        % functional inverse - no inverse will be offered for a
        % discontinuous function.
        y = reshape(NaN,sizex);
        
      otherwise
        % anything else
        error('Evalmode must be one of [0 1 2 3 -1]')
    end
    
  otherwise
    error('slmeval only handles the constant, linear, or cubic Hermite case')
end

% when all done, make sure the size of the output is the
% same as the input was.
y = reshape(y,sizex);


% ===================================================
% ================ begin subfunction ================
% ===================================================
function bind = nmbs(x,binedges)
% nmbs: non-monotone bin search
%
% Finds the leftmost bin that contains each x. Assumes
% that x has already been clipped so it must lie in
% some bin.
nb = length(binedges);
nx = length(x);

% if the bins are actually monotone, then just use histc
db = diff(binedges);
if all(db>0)
  % increasing bins
  [junk,bind] = histc(x,binedges);
  bind(bind==nb)=nb-1;

elseif all(db<0)
  % decreasing sequence of edges
  [junk,bind] = histc(-x,-binedges);
  bind(bind==nb)=nb-1;

else
  % non-monotone sequence of edges. Do this one the
  % hard way. Find the first bin that fits.
  i = 1;
  j = 1:nx;
  bind = ones(nx,1);
  while ~isempty(j) && (i<nb)
    k = ((binedges(i)>=x(j)) & (binedges(i+1)<=x(j))) | ...
      ((binedges(i)<=x(j)) & (binedges(i+1)>=x(j)));

    if any(k)
      bind(j(k)) = i;
      j(k)=[];
    end

    % increment bin we will look in
    i=i+1;

  end

end






