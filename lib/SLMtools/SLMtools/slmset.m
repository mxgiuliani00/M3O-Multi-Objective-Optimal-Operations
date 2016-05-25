function prescription=slmset(varargin)
% slmset: defines the shape prescription for a model
% usage 1: prescription=slmset(prop1,val1,pro2,val2,...)
% usage 2: prescription=slmset(prescription,prop1,val1,pro2,val2,...)
%
%  A set of property/value pairs are parsed into a structure
%  to use as the prescription for slmengine and slmfit.
% 
% arguments:
%  prop, val - property/value pairs (see below for the complete list)
%          There is no limit on the upper number of property/value
%          pairs allowed.
%          Property names may be shortened, as long as they are
%          unambiguous. Case is ignored.
%
%  prescription - shape prescription structure, defining that
%          which the user is willing to say about the model.
%          The fields of the prescription structure reflect
%          the defaults, plus any modifications specified.
%          This structure will control the behaviours built
%          a curve fit by slmengine or the defaults for slmfit.
%
%
% Property/value pairs:
%
% Properties are character strings, chosen to be mnemonic of their
% purpose. Any property name may be shortened, as long as the
% shortened string is unambiguous. Thus since no other property
% starts with the letter k, any of these alternatives are acceptable
% shortenings for the 'knots' property:  'knot', 'kno', 'kn' or 'k'.
% In the event that a given property is assigned more than once in the
% list of property/value pairs, only the last value in the list is
% assigned.
%
% Property names & admissable values:
%
% 'C2     - Most users want their cubic spline curves to be as
%         smooth and differentiable as possible. For a cubic
%         spline, this is called a twice continuously differentiable
%         function, here simply 'C2'. Of course, if your spline is not
%         cubic, this specification is blithely ignored.
%
%         = 'on' --> Causes the result to be twice differentiable
%           across the knots.
%
%         = 'off' --> Allows a cubic spline to be only C1. The second
%           derivative of the spline MAY have discontinuities across
%           each knot of the spline.
%
%         DEFAULT VALUE:  'on'
%
%         Comment: I would rarely recommend changing this option
%         from the default.
%
% 'concavedown' - controls curvature of the function
%         = 'off'  --> No part of the spline is constrained to be a
%           concave down function (i.e., a negative second derivative.)
%         = 'on' --> f''(x) >= 0 over the entire domain of the spline.
%         = vector of length 2 denoting the start and end points of a
%           region of the spline over which the second derivative is
%           negative.
%         = array of size nx2, each row of which denotes the start and
%           end points of a region of the spline over which the second
%           derivative is negative.
%
%         DEFAULT VALUE:  'off'
%
%         Comment: curvature properties do not apply to piecewise
%         constant functions.
%
%         Comment: This constraint is equivalent to a monotone
%         decreasing slope over the specified region.
%
% 'concaveup' - controls curvature of the function
%         = 'off'  --> No part of the spline is constrained to be a
%           concave up function (i.e., a positive second derivative.)
%         = 'on' --> f''(x) >= 0 over the entire domain of the spline.
%         = vector of length 2 denoting the start and end points of a
%           region of the spline over which the second derivative is
%           positive.
%         = array of size nx2, each row of which denotes the start and
%           end points of a region of the spline over which the second
%           derivative is positive.
%
%         DEFAULT VALUE:  'off'
%
%         Comment: curvature properties do not apply to piecewise
%         constant functions.
%
%         Comment: This constraint is equivalent to a monotone
%         increasing slope over the specified region.
%
% 'constantregion' - defines a region over which the curve is forced
%         be constant, although the exact level is not defined.
%         = [] --> No region of the spline is forced to be a constant
%           function.
%         = vector of length 2 denoting the start and end points of a
%           region of the spline over which it is a constant function.
%         = array of size nx2, each row of which denotes the start and
%           end points of a region of the spline over which it is a
%           constant function.
%
%         DEFAULT VALUE:  []
%
%         Comments: A segment which is forced to be constant over
%         only part of a knot interval must necessarily be constant
%         over that entire interval, since the curve is composed of
%         a single polynomial segment in a knot interval.
%
% 'decreasing' - controls monotonicity of the function
%         = 'off'  --> No part of the spline is constrained to be a
%           decreasing function.
%         = 'on' --> the function will be decreasing over its entire domain.
%         = vector of length 2 --> denotes the start and end points of a
%           region of the curve over which it is monotone decreasing.
%         = array of size nx2 --> each row of which denotes the start
%           and end points of a region of the curve over which it is
%           monotone decreasing.
%
%         DEFAULT VALUE:  'off'
%
%         Comments: in actuality this property should be named
%         'non-increasing', since a constant function is admissible.
%         In addition, it is a sufficient constraint for monotonicity.
%         It is not a necessary constraint. There may exist another
%         spline which has a slightly lower sum of squares and is also
%         monotone.
%
% 'degree' - controls the degree of the piecewise Hermite function
%         = 'constant' --> Use a piecewise constant "Hermite" 
%         = 'linear' --> Use a piecewise linear Hermite 
%         = 'cubic' --> Use a piecewise cubic Hermite 
%         
%         As a concession to memory, valid synonyms for 'constant'
%         'linear', and 'cubic' are respectively the integers 0, 1, & 3
%         
%         DEFAULT: 'cubic'
%
%         Comment: Some properties are inappropriate for all degrees
%         of function. E.g., it would be silly to specify a specific
%         value for the left hand end point slope of a piecewise
%         constant function. All information supplied will be used
%         to whatever extent possible.
%
%         The "order" of a form, as used by the spline toolbox, will
%         be one more than the degree. 'order' is also a valid property
%         in these tools, but it results only in setting the degree
%         field, where degree = order - 1.
%
% 'endconditions' - controls the end conditions applied to the spline
%         = 'natural' --> The "natural" spline conditions will be applied.
%           I.e., f''(x) = 0 at end end of the spline.
%         = 'notaknot' --> Not-a-knot end conditions applied.
%         = 'periodic' --> Periodic end conditions applied.
%         = 'estimate' --> end conditions are estimated from the data.
%
%         DEFAULT VALUE:  'estimate
%
%         Comment: Except for periodicity, end conditions are not
%         relevant to any degree model below cubic.
%
%         Periodic end conditions mean that the function values are
%         the same at each end of the curve for linear and cubic fits.
%
%         For cubic fits, periodicity means that the function has
%         also first and second derivative continuity across the
%         wrapped boundary.
%
%         For piecewise constant fits, end conditions do not apply,
%         and are ignored.
%
% 'envelope' - allows the user to solve for an envelope of the data
%         = 'off' --> the curve will be a simple least squares spline
%         = 'supremum' --> comute a model such that all residuals
%           (yhat - y) are positive. In effect, the curve will be a
%           "supremum" (least upper bound) function.
%         = 'infimum' --> comute a model such that all residuals
%           (yhat - y) are negative. In effect, the curve will be a
%           "infimum" (greatest lower bound) function.
%
%         DEFAULT VALUE:  'off'
%
%         Comment: this option should be rarely used, but its a cute
%         idea when it does come up.
%
% 'errorbar' - defines a set of lower and upper bounds for the function
%         value of the spline at each data point. If there are n data
%         points, then if the corresponding value is:
% 
%         = [] --> No errorbar constraints will be imposed
%         = a scalar E --> error bars will be set at [Y-E,Y+E]
%         = a vector E --> error bars will be set at [Y-E,Y+E]
%         = an nx2 array --> error bars will be set at [Y-E(:,1),Y+E(:,2)]
%
%         DEFAULT VALUE:  []
%
%         Comment: It is possible that depending on the choice of
%         knots, it will be impossible to satisfy some sets of
%         error bars. It may be best to simply use every single data
%         point as a knot. Of course, replicate data points with
%         non-overlappng error bars will always cause a failure.
%
% 'extrapolation' - tells SLMEVAL how to extrapolate the function outside
%         the defined knot points of the spline. This parameter has
%         NOTHING to do with how the curve will be fit, but only
%         with how it will be evaluated after the fit is done.
%
%         = 'error' --> No extrapolation will be allowed. An error
%            will be returned if the user tries to evaluate the
%            spline at any location outside of the knots.
%
%         = 'warning' --> A warning message will be generated when
%            extrapolation is attempted, but extrapolation will still
%            done as a constant value outside of the end point knots.
%
%         =  'constant' --> Constant extrapolation will be performed,
%            using the value of the spline at the associated end point
%            knot. This is the default, to be consistent with the
%            behavior of this tool in previous releases.
%
%         =  'linear' --> Linear extrapolation will be performed,
%            using the slope of the spline at the associated end point
%            knot.
%
%         =  'cubic' --> Cubic extrapolation will be performed,
%            using the shape of the spline at the associated end point
%            knot.
%
%         = 'NaN' --> NaNs will be inserted for any attempted extrapolation.
%            although no error or warning message will be generated.
%
%         DEFAULT value: 'constant'
%
% 'increasing' - controls monotonicity of the function
%         = 'off'  --> No part of the spline is constrained to be an
%           increasing function.
%         = 'on' --> the function will be increasing over its entire domain.
%         = vector of length 2 --> denotes the start and end points of a
%           region of the curve over which it is monotone increasing.
%         = array of size nx2 --> each row of which denotes the start
%           and end points of a region of the curve over which it is
%           monotone increasing.
%
%         DEFAULT VALUE:  'off'
%
%         Comments: in actuality this property should be named
%         'non-decreasing', since a constant function is admissible.
%         In addition, it is a sufficient constraint for monotonicity.
%         It is not a necessary constraint. There may exist another
%         spline which has a slightly lower sum of squares and is also
%         monotone.
%
% 'integral' - known aim value for the integral of the curve over
%         its domain.
%         
%         DEFAULT VALUE:  []
%
% 'interiorknots' - allows for free knot placement (of the interior knots)
%         = 'fixed' 
%           spaced knots. In this case, the first data point will
%           be the first knot.
%         = 'free' --> uses fmincon to optimize the interior knot
%           placement to minimize the overall rmse.
%
%         DEFAULT VALUE:  'fixed'
%
%         Comment: The initial values for the knot placement are
%         taken from the knots property.
%
%         Comment: Since the free knot placement is done by an
%         optimizer (fmincon), it is not allowed to both choose
%         a set of free knots and set the rmse of the fit.
%
%         Comment: The first and last knots are not adjusted by
%         the optimization, so there must be at least 3 knots.
%
%         Comment: Piecewise constant functions sometimes have
%         difficulty with free knots.
%
% 'jerk' - Controls the sign of the jerk function over the support
%         of the spline, to be either 'positive' or 'negative' (or no
%         constraint at all applied.) This a global sign
%         constraint on the third derivative of the function.
%
%         = 'positive' --> causes the third derivative to be
%           everywhere positive over the support of the function
%
%         = 'negative' --> causes the third derivative to be
%           everywhere negative over the support of the function
%
%         = '' --> No constraint applied
% 
%         DEFAULT VALUE: ''
%
%         Comment: This property can be used to constrain the
%         curvature of the function to be increasing (or decreasing)
%         over the support, in combination with the concaveup or
%         concavedown properties as appropriate.
%
%         Comment: The third derivative is sometimes known as the
%         'jolt', but from my experience, 'jerk' seems to be the
%         most common name. I will add that the third derivative
%         quite difficult to estimate for a spline model, especially
%         if the data is at all noisy.t
%
% 'knots' - controls the number of knots used, or the number of
%           equally spaced knots.
%
%         = A scalar integer which denotes the number of equally
%           spaced knots. In this case, the first and last data
%           points will be the first and last knots.
%
%         = A vector containing the list of knots themselves.
%           The knots must be distinct and MUST wholly contain
%           the data or an error will be generated. No replicate
%           knots are allowed.
%
%         = A negative value K (no larger in absolute value than
%           the number of data points - 1) causes every K'th data
%           point to be used as a knot.
%
%         DEFAULT VALUE:  6
%
% 'leftmaxslope' - controls the maximum slope allowed at the left hand
%           end of the curve.
%         = [] --> No explicit value provided for the maximum slope
%           of the spline at its left hand end point.
%         = A numeric scalar --> the slope of the function will be
%           constrained to not rise above this value at its left
%           hand end point (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
% 
% 'leftmaxvalue' - controls the maximum valued allowed at the left hand
%           end of the curve.
%         = [] --> No explicit value provided for the maximum value
%           of the spline at its left hand end point.
%         = A numeric scalar --> the function will be constrained
%           to not rise above this value at its left hand end point
%           (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
%
% 'leftminslope' - controls the minimum slope allowed at the left hand
%           end of the curve.
%         = [] --> No explicit value provided for the minimum slope
%           of the spline at its left hand end point.
%         = A numeric scalar --> the slope of the function will be
%           constrained to not fall below this value at its left
%           hand end point (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
%
% 'leftminvalue' - controls the minimum valued allowed at the left hand
%           end of the curve.
%         = [] --> No explicit value provided for the minimum value
%           of the spline at its left hand end point.
%         = A numeric scalar --> the function will be constrained
%           to not fall below this value at its left hand end point
%           (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
%
% 'leftslope' - controls the function slope at the left hand endpoint.
%         = [] --> No explicit value provided for the slope of the
%           curve at its left hand end point.
%         = A numeric scalar --> the function will be assigned this
%           slope at its left hand end point (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
%
% 'leftvalue' - controls the function value at its left hand endpoint.
%         = [] --> No explicit value provided for the value of the
%           curve at its left hand end point.
%         = A numeric scalar --> the function will be assigned this
%           value at its left hand end point (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
%
% 'linearregion' - defines a region over which the curve is forced
%         be linear, although the exact level is not defined.
%         = [] --> No region of the spline is forced to be a linear
%           function.
%         = vector of length 2 denoting the start and end points of a
%           region of the spline over which it is a linear function.
%         = array of size nx2, each row of which denotes the start and
%           end points of a region of the spline over which it is a
%           linear function.
%
%         DEFAULT VALUE:  []
%
%         Comment 1: A segment which is forced to be linear over
%         only part of a knot interval must necessarily be constant
%         over that entire interval, since the curve is composed of
%         a single polynomial segment in a knot interval.
%         Comment 2: A linear region may extend across knots, in
%         which case the slope will take on the same value across
%         knot boundaries.
%         
% 'maxslope' - controls the globally maximum slope of the function
%         = [] --> No explicit value provided for the globally maximum
%           slope of the curve.
%         = A numeric scalar --> the globally maximum slope of the spline
%
%         DEFAULT VALUE:  []
%
%         Comment: This is a sufficient constraint for the maximum
%         slope of the spline. It is not a necessary constraint.
%         There may exist another spline which has a slightly lower
%         sum of squares and also has the same maximum slope.
%
% 'maxvalue' - controls the globally maximum value of the curve
%         = [] --> No explicit value provided for the globally
%           maximum value of the spline.
%         = A numeric scalar --> the function must lie no higher
%           than this maximum value.
%
%         DEFAULT VALUE:  []
%
%         Comment 1: This constraint is only a necessary constraint.
%         It is not sufficient. In some circumstances the spline may
%         pass slightly above this maximum value
%         Comment 2: The location of the global minimizer is unspecified.
%
% 'minslope' - controls the globally minimum slope of the function
%         = [] --> No explicit value provided for the globally minimum
%           slope of the curve.
%         = A numeric scalar --> the globally minimum slope of the spline
%
%         DEFAULT VALUE:  []
%
%         Comment: This is a sufficient constraint for the minimum
%         slope of the spline. It is not a necessary constraint.
%         There may exist another spline which has a slightly lower
%         sum of squares and also has the same minimum slope.
%
% 'minvalue' - controls the globally minimum value of the curve
%         = [] --> No explicit value provided for the globally
%           minimum value of the spline.
%         = A numeric scalar --> the function will lie no lower
%           than this minimum value.
%
%         DEFAULT VALUE:  []
%
%         Comment 1: This constraint is only a necessary constraint.
%         It is not sufficient. In some circumstances the spline may
%         pass slightly below this minimum value
%         Comment 2: The location of the global minimizer is unspecified.
%
% 'negativeinflection' - controls the existence and placement of a
%           inflection point of the final curve. The second derivative
%           of the function will pass through zero at this point.
%
%         = [] --> no inflection point constraint employed
%         = a numeric scalar --> the function will have a point of
%           inflection at the supplied x location
%
%         DEFAULT VALUE:  []
%
%         Comment: NegativeInflection is really just a composite property
%         of a curve. A negative inflection point at x == a is equivalent
%         to a ConcaveUP function for x<=a, and a ConcaveDown function
%         for x>=a. As such, NegativeInflection will override any other
%         curvature properties one has previously set.
%
%         Comment: Inflection points only apply to linear or cubic models
%
% 'order' - An implicit synonym for 'degree', controls the degree of
%         the piecewise Hermite function. Order MUST be an numeric
%         integer, from the set [1, 2, 4].
%
%         = 1 --> Use a piecewise constant "Hermite" 
%         = 2 --> Use a piecewise linear Hermite 
%         = 4 --> Use a piecewise cubic Hermite 
%         
%         DEFAULT: 4
%
%         Setting the order to some value has the effect of setting
%         the degree of the Hermite spline fit to one less than the
%         order.
%
%         Comment: Some properties are inappropriate for all degrees
%         of function. E.g., it would be silly to specify a specific
%         value for the left hand end point slope of a piecewise
%         constant function. All information supplied will be used
%         to whatever extent possible.
%
% 'plot' - controls whether a final plot is generated of the curve
%         = 'off' --> No plot
%         = 'on' --> plot the curve and data using slmplot 
%         
%         DEFAULT VALUE: 'off'
%
% 'positiveinflection' - controls the existence and placement of a
%           inflection point of the final curve. The second derivative
%           of the function will pass through zero at this point.
%
%         = [] --> no inflection point constraint employed
%         = a numeric scalar --> the function will have a point of
%           inflection at the supplied x location
%
%         DEFAULT VALUE:  []
%
%         Comment: PositiveInflection is really just a composite property
%         of a curve. An inflection point at x == a is equivalent to a
%         ConcaveDown function for x<=a, and a ConcaveUp function for
%         x>=a. As such, PositiveInflection will override any other
%         curvature properties one has previously set.
%
%         Comment: Inflection points only apply to linear or cubic models
%
% 'predictions' - The number of points to evaluate the curve at. Assumed
%         to be equally spaced. The supplied value must be a positive
%         scalar, integer value, >= 2, or empty. If empty, no predictions
%         will be generated.
%         
%         DEFAULT VALUE: []
%
%         Comment: plotslm generates 1001 points along the curve.
%
% 'regularization' - 
%         = [] --> Uses the default regularization parameter of 0.0001.
%
%         = A Non-negative scalar value --> explicitly defines the weight
%           given to smoothness of the resulting curve.
%
%         = 'crossvalidation' --> Use cross validation method to choose
%           the regularization parameter.
%
%         = A NEGATIVE scalar value --> attempts to choose a
%           regularization parameter which has as its rmse the absolute
%           value of the supplied value.
%
%         = 'smoothest' --> Finds the smoothest curve that satisfies
%           the supplied prescription. This is not really a least
%           squares model, since the goal is purely to maximize the
%           smoothness. This option would very often be used in
%           conjunction with errorbars.
%
%         DEFAULT VALUE:  0.0001
%
%         Comment: Smaller values will yield less smoothing, larger
%         values more smoothing. In most cases this parameter should
%         be left alone. It is used to prevent numerical singularities
%         in the linear algebra, as well as help in the case of
%         extrapolation and intrapolation. Smoothness of the resulting
%         spline can be far better controlled by changing the number
%         of knots and their placement. Specifically, the regularization
%         parameter is a scale factor applied to the integral of the
%         squared second derivative of the spline.
%
%         Comment: It is possible that no value for the regularization
%         parameter will yield the given rmse. In this case the sky will
%         fall down.
%
%         Comment: Since the cross validation option and matching a
%         given rmse are both optimizations, it is not allowed to use
%         these options together with the interiorknots estimation.
%
%         Comment: Both the cross validation and rmse options can
%         both be quite slow. Remember that they are optimizations.
%
% 'result' - controls the output structure style
%         = 'pp'  --> Returns a pp struct, use ppval to evaluate
%         = 'slm' --> Returns a slm struct in a Hermite form. Evaluate
%         using slmeval.
%
%         DEFAULT VALUE:  'slm'
%
% 'rightmaxslope' - controls the maximum slope allowed at the right
%           hand end of the curve.
%         = [] --> No explicit value provided for the maximum slope
%           of the spline at its right hand end point.
%         = A numeric scalar --> the slope of the function will be
%           constrained to not rise above this value at its right
%           hand end point (i.e., the last knot.)
%
%         DEFAULT VALUE:  []
% 
% 'rightmaxvalue' - controls the maximum valued allowed at the right
%           hand end of the curve.
%         = [] --> No explicit value provided for the maximum value
%           of the spline at its right hand end point.
%         = A numeric scalar --> the function will be constrained
%           to not rise above this value at its right hand end point
%           (i.e., the first knot.)
%
%         DEFAULT VALUE:  []
%
% 'rightminslope' - controls the minimum slope allowed at the right
%           hand end of the curve.
%         = [] --> No explicit value provided for the minimum slope
%           of the spline at its right hand end point.
%         = A numeric scalar --> the slope of the function will be
%           constrained to not fall below this value at its right
%           hand end point (i.e., the last knot.)
%
%         DEFAULT VALUE:  []
%
% 'rightminvalue' - controls the minimum valued allowed at the right
%           hand end of the curve.
%         = [] --> No explicit value provided for the minimum value
%           of the spline at its right hand end point.
%         = A numeric scalar --> the function will be constrained
%           to not fall below this value at its right hand end point
%           (i.e., the last knot.)
%
%         DEFAULT VALUE:  []
%
% 'rightslope' - controls the function slope at the right hand endpoint.
%         = [] --> No explicit value provided for the slope of the
%           curve at its right hand end point.
%         = A numeric scalar --> the function will be assigned this
%           slope at its right hand end point (i.e., the last knot.)
%
%         DEFAULT VALUE:  []
%
% 'rightvalue' - controls the function value at its right hand endpoint.
%         = [] --> No explicit value provided for the value of the
%           curve at its right hand end point.
%         = A numeric scalar --> the function will be assigned this
%           value at its right hand end point (i.e., the last knot.)
%
%         DEFAULT VALUE:  []
%
% 'robust' - Controls the use of a simple robust solver, employing
%           iteratively reweighted least squares.
%         = 'off' --> Simple ordinary least squares is employed
%         = 'on' --> Performs iteratively re-weighted least squares
%
%         DEFAULT VALUE:  'off'
%
%         Comment: when weights are also supplied, iterative
%         re-weighting is applied on top of those weights.
%
%         Comment: Beware use of the robust fitting option combined
%         with other options such as:
%
%            free interior knots
%            an rmse goal
%            cross-validation as a choice of a regularization parameter
%
%         In combination with the robust fitting option, any of those
%         options may introduce instabilities in the result and will at
%         best be slow.
%
% 'scaling' - controls data scaling to avoid numerical problems
%         = 'on' --> data is shifted and scaled so as to minimize
%           any numerical problems that may result in the solution.
%         = 'off' --> No scaling is done.
%
%         DEFAULT VALUE:  'on'
%
%         Comments: There is no transformation that will positively
%         eliminate all problems, but this transformation will
%         somtimes make the solution a bit kess prone to problems.
%         All scalings are undone in the final spline coefficients,
%         so these parameters are completely transparent to the user.
%         In most cases you would see no difference in the solution,
%         with or without scaling.
%
%         Any scaling is automatically applied to both the data as
%         well as all provided prescription parameters.
%
%         The actual shift and scale parameters applied will be reported
%         as part of the stats.
%
% 'segmentconstant' - Allows the user to specify the index of a constant
%         segment or segments. This is for the rare case where the knots
%         may be free, yet there is a known internal constant region, so
%         perhaps one always needs the second polynomial segment to be
%         exactly constant even though the knots vary in their placement.
%
%         The index of the knot interval which is to be held fixed is
%         supplied. More than one interval may be so held fixed. If there
%         are N knots, then the supplied index must be an integer between 1
%         and N-1.
%
%         Comment: This property is meaningless for piecewise constant
%         splines.
%         
%         Comment: Cubic splines with constant segments will still have
%         the overall degree of continuity. This may result in strange
%         looking curves if there are several segments held constant.
%         
%         DEFAULT Value:  []
%
% 'segmentlinear' - Allows the user to specify the index of a linear
%         segment or segments. This is for the rare case where the knots
%         may be free, yet there is a known internal linear region.
%
%         The index of the knot interval which is to be held fixed is
%         supplied. More than one interval may be so held fixed. If there
%         are N knots, then the supplied index must be an integer between 1
%         and N-1.
%
%         This property is meaningless for piecewise constant or linear
%         splines.
%
%         DEFAULT Value:  []
%
% 'simplepeak' - controls the existence and placement of a single
%           maximizer of the final curve
%
%         = [] --> no peak placement constraint employed
%         = a numeric scalar --> the function will attain its maximum
%           at the supplied x location
%
%         DEFAULT VALUE:  []
%
%         Comment: SimplePeak is really just a composite property of
%         a curve. A peak at x == a is equivalent to a monotone increasing
%         function for x<=a, and a monotone decreasing function for
%         x>=a. As such, simplepeak will override any other monotonicity
%         properties one has previously set.
%
% 'simplevalley' - controls the existence and placement of a single
%           minimizer of the final curve
%
%         = [] --> no valley placement constraint employed
%         = a numeric scalar --> the function will attain its minimum
%           at the supplied x location
%
%         DEFAULT VALUE:  []
%
%         Comment: SimpleValley is really just a composite property
%         of a curve. A valley at x == a is equivalent to a monotone
%         decreasing function for x<=a, and a monotone increasing
%         function for x>=a. As such, simplevalley will override any
%         other monotonicity properties one has previously set.
%
% 'sumresiduals' - Allows the user to specify an explicit goal for the
%         sum of the residuals. This is a property that virtually
%         NOBODY should ever have a need for, since normally you get
%         a zero sum implicitly, as a freebie. This is due to the
%         presence of an effective constant term in the regression
%         model. Under some circumstances however, that property is
%         circumvented, for example by forcing the model to explicitly
%         pass through a given point. Use of non-unit weights in the
%         model could also cause a non-zero sum of residuals.
%
%         In practice, any real scalar value for the sum of the
%         residuals can be specified, but a value of 0 is the only
%         value that makes any sense that I can see. Weights on the
%         data points are disregarded when the residual sum is
%         computed.
%
%         The computation of residual here is defined as [yhat - y].
%         
%         DEFAULT VALUE: []
%
% 'verbosity' - controls commandline feedback to the user
%         = 0 --> No response at the commandline
%         = 1 --> Basic fit statistics reported at the command line
%         = 2 --> Debug level output
%         
%         DEFAULT VALUE: 0
%
% 'weights' - defines a weight vector 
%         = [] --> all data points are assigned equal (unit) weight.
%         = vector of the same length as length(x), denotes relative
%           weights for each data point. If supplied, the length of
%           this vector of weights must be the same as the number of
%           data points.
%
%         DEFAULT VALUE:  []
%
% 'xy'  - Forces the curve through an individual point or set of points
%         = [] --> no points are forced
%         = a 1x2 vector --> an x-y pair that the curve must pass
%           through with no error.
%         = an nx2 array --> each row corresponds to a single point
%           that the curve passes through.
%         
%         DEFAULT VALUE: []
%
%         Comment 1: The curve will pass through the desired point
%         to within computational error, IF it is possible to do so.
%         Comment 2: Multiple points that are inconsistent with each
%         other, or inconsistent with other parameters that are set
%         will cause failure of the least squares.
%
% 'xyp'  - Forces the curve to have a specified slope at an individual
%         point or set of points
%         = [] --> no points have their slope enforced
%         = a 1x2 vector --> an x-yprime pair that the curve must
%           satisfy through with no error.
%         = an nx2 array --> each row corresponds to an x-yprime
%           pair that the curve must satisfy.
%         
%         DEFAULT VALUE: []
%
%         Comment 1: The curve will satisfy the desired slope to
%         within computational error, IF it is possible to do so.
%         Comment 2: Multiple points that are inconsistent with each
%         other, or inconsistent with other parameters that are set
%         will cause failure of the least squares.
%         Comment 3: Setting the slope at some point of a zero'th
%         degree function will be ignored.
%
% 'xypp' - Forces the curve to have a specified second derivative
%         at an individual point or set of points
%         = [] --> no points have their 2nd derivative enforced
%         = a 1x2 vector --> an x-y'' pair that the curve must
%           satisfy through with no error.
%         = an nx2 array --> each row corresponds to an x-y''
%           pair that the curve must satisfy.
%         
%         DEFAULT VALUE: []
%
%         Comment: The curve will satisfy the desired 2nd derivative
%         to within computational error, IF it is possible to do so.
%         Comment: Multiple points that are inconsistent with each
%         other, or inconsistent with other parameters that are set
%         will cause failure of the least squares.
%         Comment: Setting the 2nd derivative at some point of a
%         zero'th degree or linear function will be ignored.
%         Comment: This property only applies to cubic models.
%
% 'xyppp' - Forces the curve to have a specified third derivative
%         at an individual point or set of points.
%
%         Really, the only meaningful use of this property that I
%         see under normal circumstances is to force a segment
%         to be only quadratic, rather than cubic. Thus, by forcing
%         the third derivative to zero at some point in the segment,
%         since the third derivative of a cubic polynomial is a
%         constant, will reduce the cubic to a quadratic.
%
%         = [] --> no points have their 3rd derivative enforced
%         = a 1x2 vector --> an x-y''' pair that the curve must
%           satisfy through with no error.
%         = an nx2 array --> each row corresponds to an x-y'''
%           pair that the curve must satisfy.
%         
%         DEFAULT VALUE: []
%
%         Comment: Setting the 3rd derivative at some point of a
%         zero'th degree or linear function will be ignored.
%         Comment: This property only applies to cubic models.


% was an initial prescription struct provided?
if (nargin==0) || ~isstruct(varargin{1})
  % defaults for all properties
  prescription.C2 = 'on';
  prescription.ConcaveDown = 'off';
  prescription.ConcaveUp = 'off';
  prescription.ConstantRegion = [];
  prescription.Decreasing = 'off';
  prescription.Degree = 3;
  prescription.EndConditions = 'estimate';
  prescription.Envelope = 'off';
  prescription.ErrorBar = [];
  prescription.Extrapolation = 'constant';
  prescription.Increasing = 'off';
  prescription.Integral = [];
  prescription.InteriorKnots = 'fixed';
  prescription.Jerk = '';
  prescription.Knots = 6;
  prescription.LeftMaxSlope = [];
  prescription.LeftMaxValue = [];
  prescription.LeftMinSlope = [];
  prescription.LeftMinValue = [];
  prescription.LeftSlope = [];
  prescription.LeftValue = [];
  prescription.LinearRegion = [];
  prescription.MaxSlope = [];
  prescription.MaxValue = [];
  prescription.MinSlope = [];
  prescription.MinValue = [];
  prescription.NegativeInflection = [];
  prescription.Order = [];
  prescription.Plot = 'off';
  prescription.PositiveInflection = [];
  prescription.Predictions = 1001;
  prescription.Regularization = 0.0001;
  prescription.Result = 'slm';
  prescription.RightMaxSlope = [];
  prescription.RightMaxValue = [];
  prescription.RightMinSlope = [];
  prescription.RightMinValue = [];
  prescription.RightSlope = [];
  prescription.RightValue = [];
  prescription.Robust = 'off';
  prescription.Scaling = 'on';
  prescription.SegmentConstant = [];
  prescription.SegmentLinear = [];
  prescription.SimplePeak = [];
  prescription.SimpleValley = [];
  prescription.SumResiduals = [];
  prescription.Verbosity = 0;
  prescription.Weights = [];
  prescription.XY = [];
  prescription.XYP = [];
  prescription.XYPP = [];
  prescription.XYPPP = [];
elseif isstruct(varargin{1})
  % there is a struct provided. use it for defaults
  prescription = varargin{1};
  varargin(1)=[];
end

% The YShift and YScale parameters cannot be passed in.
% They are created on the fly whenever 'scaling' is set as 'on'.
prescription.YScale = [];
prescription.YShift = [];

% begin processing property/value pairs
if isempty(varargin)
  % just use the defaults already present in prescription
else
  prescription = parse_pv_pairs(prescription,varargin);
end

% check that all properties were set, also check
% for some simple errors

% C2 must be: 'off', 'on', or ''
prescription = value_check(prescription,'C2', ...
   {'off' 'on'},[],0);

% ConcaveDown must be: 'off', 'on', an nx2 array
prescription = value_check(prescription,'ConcaveDown', ...
   {'off' 'on'},{[NaN,2]},0);
if (ischar(prescription.ConcaveDown) && ...
    strcmp(prescription.ConcaveDown,'on')) || ...
   (isnumeric(prescription.ConcaveDown) && ...
   ~isempty(prescription.ConcaveDown))
  % 
  % override any prior inflection point settings
  prescription.PositiveInflection = [];
  prescription.NegativeInflection = [];
  
end

% ConcaveUp must be: 'off', 'on', an nx2 array
prescription = value_check(prescription,'ConcaveUp', ...
   {'off' 'on'},{[NaN,2]},0);
if (ischar(prescription.ConcaveUp) && ...
    strcmp(prescription.ConcaveUp,'on')) || ...
   (isnumeric(prescription.ConcaveUp) && ...
   ~isempty(prescription.ConcaveUp))
  % 
  % override any prior inflection point settings
  prescription.PositiveInflection = [];
  prescription.NegativeInflection = [];
  
end

% ConstantRegion must be: [], or an nx2 array
prescription = value_check(prescription,'ConstantRegion', ...
   {},{[NaN,2]},1);

% Decreasing must be: 'off', 'on', or an nx2 array
prescription = value_check(prescription,'Decreasing', ...
   {'off' 'on'},{[NaN,2]},0);
if (ischar(prescription.Decreasing) && ...
    strcmp(prescription.Decreasing,'on')) || ...
   (isnumeric(prescription.Decreasing) && ...
   ~isempty(prescription.Decreasing))
  % 
  % override any prior peak or valley settings
  prescription.SimplePeak = [];
  prescription.SimpleValley = [];
  
end

% EndConditions must be: 'natural', 'notaknot', 'periodic' 'estimate'
prescription = value_check(prescription,'EndConditions', ...
   {'natural', 'notaknot', 'periodic' 'estimate'},{},0);

% Envelope must be: 'off', 'supremum', 'infinmum'
prescription = value_check(prescription,'Envelope', ...
   {'off' 'supremum' 'infimum'},{},0);
 
% ErrorBar must be: [], scalar, vector, or an nx2 array
prescription = value_check(prescription,'ErrorBar', ...
     {},{[1 1] [NaN 1] [1 NaN] [NaN 2]},1);

% Extrapolation must be: 'error', 'warning', 'constant', 'linear', 'cubic', 'nan'
prescription = value_check(prescription,'Extrapolation', ...
   {'error', 'warning', 'constant', 'linear', 'cubic', 'nan'},[],0);

% InteriorKnots must be: 'fized', 'free'
prescription = value_check(prescription,'InteriorKnots', ...
   {'fixed', 'free'},{},0);

% Increasing must be: 'off', 'on', or an nx2 array
prescription = value_check(prescription,'Increasing', ...
   {'off' 'on'},{[NaN,2]},0);
if (ischar(prescription.Increasing) && ...
    strcmp(prescription.Increasing,'on')) || ...
   (isnumeric(prescription.Increasing) && ...
   ~isempty(prescription.Increasing))
  % 
  % override any prior peak or valley settings
  prescription.SimplePeak = [];
  prescription.SimpleValley = [];
  
end

% Integral must be a numeric scalar or empty
prescription = value_check(prescription,'Integral', ...
   {},{[1 1]},1);

% Jerk must be: '', 'positive', or 'negative'
prescription = value_check(prescription,'Jerk', ...
   {'positive' 'negative'},{},1);

% Knots must be a numeric scalar or a vector
prescription = value_check(prescription,'Knots', ...
   {},{[1 1], [NaN,1], [1,NaN]},0);

% LeftMaxSlope must be a numeric scalar or empty
prescription = value_check(prescription,'LeftMaxSlope', ...
   {},{[1 1]},1);

% LeftMaxValue must be a numeric scalar or empty
prescription = value_check(prescription,'LeftMaxValue', ...
   {},{[1 1]},1);

% LeftMinSlope must be a numeric scalar or empty
prescription = value_check(prescription,'LeftMinSlope', ...
   {},{[1 1]},1);

% LeftMinValue must be a numeric scalar or empty
prescription = value_check(prescription,'LeftMinValue', ...
   {},{[1 1]},1);

% LeftSlope must be a numeric scalar or empty
prescription = value_check(prescription,'LeftSlope', ...
   {},{[1 1]},1);

% LeftValue must be a numeric scalar or empty
prescription = value_check(prescription,'LeftValue', ...
   {},{[1 1]},1);

% LinearRegion must be: [], or an nx2 array
prescription = value_check(prescription,'LinearRegion', ...
   {},{[NaN,2]},1);

% MaxSlope must be a numeric scalar or empty
prescription = value_check(prescription,'MaxSlope', ...
   {},{[1 1]},1);

% MaxValue must be a numeric scalar or empty
prescription = value_check(prescription,'MaxValue', ...
   {},{[1 1]},1);

% MinSlope must be a numeric scalar or empty
prescription = value_check(prescription,'MinSlope', ...
   {},{[1 1]},1);

% MinValue must be a numeric scalar or empty
prescription = value_check(prescription,'MinValue', ...
   {},{[1 1]},1);

% NegativeInflection must be a numeric scalar or empty
prescription = value_check(prescription,'NegativeInflection', ...
   {},{[1 1]},1);
if ~isempty(prescription.NegativeInflection)
  % override any prior curvature settings
  prescription.ConcaveUp = 'off';
  prescription.ConcaveDown = 'off';
  prescription.PositiveInflection = [];
end

% Degree must be: 'constant', 'linear', 'cubic', or 0, 1, 3
prescription = value_check(prescription,'Degree', ...
   {'constant' 'linear', 'cubic' '0', '1', '3'},{[1 1]},0);
switch prescription.Degree
case {'constant' '0' 0}
  prescription.Degree = 0;
case {'linear' '1' 1}
  prescription.Degree = 1;
case {'cubic' '3' 3}
  prescription.Degree = 3;
otherwise
  error 'degree may be any of: {0,1,3, ''constant'', ''linear'', ''cubic''}'
end

% Order must be: 1, 2, 4
prescription = value_check(prescription,'Order', ...
   {'1', '2', '4'},{1 2 4},1);
if isempty(prescription.Order)
  % Do not override Degree, already defined.
else
  switch prescription.Order
    case {'1' 1}
      prescription.Degree = 0;
    case {'2' 2}
      prescription.Degree = 1;
    case {'4' 4}
      prescription.Degree = 3;
    otherwise
      error('Order may be any of: {[],1,2,4}')
  end
end

% Plot must be: 'off', 'on' 
prescription = value_check(prescription,'Plot', ...
   {'off' 'on'},{},0);

% PositiveInflection must be a numeric scalar or empty
prescription = value_check(prescription,'PositiveInflection', ...
   {},{[1 1]},1);
if ~isempty(prescription.PositiveInflection)
  % override any prior curvature settings
  prescription.ConcaveUp = 'off';
  prescription.ConcaveDown = 'off';
  prescription.NegativeInflection = [];
end

% Predictions must be a numeric scalar or empty
prescription = value_check(prescription,'Predictions', ...
   {},{[1 1]},1);
% make sure its an integer scalar, >= 2, or empty
if ~isempty(prescription.Predictions) && ...
    ((numel(prescription.Predictions) > 1) || ...
    (prescription.Predictions < 2) || ...
    (rem(prescription.Predictions,1) ~= 0))
  error('Predictions must be scalar, >= 2 if supplied.')
end
 
% Regularization must be: 'crossvalidation', 'smoothest', or a numeric scalar
prescription = value_check(prescription,'Regularization', ...
   {'crossvalidation' 'smoothest'},{[1 1]},0);
 
% Result must be: 'pp', 'slm', or []
prescription = value_check(prescription,'Result', ...
   {'pp' 'slm'},{},0);

% RightMaxSlope must be a numeric scalar or empty
prescription = value_check(prescription,'RightMaxSlope', ...
   {},{[1 1]},1);

% RightMaxValue must be a numeric scalar or empty
prescription = value_check(prescription,'RightMaxValue', ...
   {},{[1 1]},1);

% RightMinSlope must be a numeric scalar or empty
prescription = value_check(prescription,'RightMinSlope', ...
   {},{[1 1]},1);

% RightMinValue must be a numeric scalar or empty
prescription = value_check(prescription,'RightMinValue', ...
   {},{[1 1]},1);

% RightSlope must be a numeric scalar or empty
prescription = value_check(prescription,'RightSlope', ...
   {},{[1 1]},1);

% RightValue must be a numeric scalar or empty
prescription = value_check(prescription,'RightValue', ...
   {},{[1 1]},1);

% Robust must be: [], 'off', 'on', or a positive scalar value
prescription = value_check(prescription,'Robust', ...
   {'off' 'on'},{},1);
if isempty(prescription.Robust)
  prescription.Robust = 'off';
end

% Scaling must be: 'off', 'on'
prescription = value_check(prescription,'Scaling', ...
   {'off' 'on'},{},0);

% SegmentConstant must be an integer numeric vector, scalar, or empty
prescription = value_check(prescription,'SegmentConstant', ...
   {},{[NaN 1] [1 NaN]},1);
if ~isempty(prescription.SegmentConstant) && ...
    (any(rem(prescription.SegmentConstant,1) ~= 0) || ...
        any(prescription.SegmentConstant < 1))
      
  error('Segment Intervals held constant must be referred to by an integer segment index.')
end

% SegmentLinear must be an integer numeric vector, scalar, or empty
prescription = value_check(prescription,'SegmentLinear', ...
   {},{[NaN 1] [1 NaN]},1);
if ~isempty(prescription.SegmentLinear) && ...
    (any(rem(prescription.SegmentLinear,1) ~= 0) || ...
        any(prescription.SegmentLinear < 1))
      
  error('Segment Intervals held linear must be referred to by an integer segment index.')
end

% Segmentlinear must be an integer numeric vector, scalar, or empty
prescription = value_check(prescription,'SegmentLinear', ...
   {},{[NaN 1] [1 NaN]},1);

% SimplePeak must be a numeric scalar or empty
prescription = value_check(prescription,'SimplePeak', ...
   {},{[1 1]},1);
if ~isempty(prescription.SimplePeak)
  % override any prior monotonicity settings
  prescription.Decreasing = 'off';
  prescription.Increasing = 'off';
  prescription.SimpleValley = [];
end

% SimpleValley must be a numeric scalar or empty
prescription = value_check(prescription,'SimpleValley', ...
   {},{[1 1]},1);
if ~isempty(prescription.SimpleValley)
  % override any prior monotonicity settings
  prescription.Decreasing = 'off';
  prescription.Increasing = 'off';
  prescription.SimplePeak = [];
end

% SumResiduals must be a numeric scalar or empty
prescription = value_check(prescription,'SumResiduals', ...
   {},{[1 1]},1);

% Verbosity must be 0 or 1 (small talk) or 2(garrulous)
prescription = value_check(prescription,'Verbosity', ...
   {},{[1 1]},0);

% Weights must be a vector or empty
prescription = value_check(prescription,'Weights', ...
   {},{[NaN 1] [1 NaN]},1);

% XY must be an nx2 array or empty
prescription = value_check(prescription,'XY', ...
   {},{[NaN 2]},1);

% XYP must be an nx2 array or empty
prescription = value_check(prescription,'XYP', ...
   {},{[NaN 2]},1);

% XYPP must be an nx2 array or empty
prescription = value_check(prescription,'XYPP', ...
   {},{[NaN 2]},1);

% XYPPP must be an nx2 array or empty
prescription = value_check(prescription,'XYPPP', ...
   {},{[NaN 2]},1);

if strcmpi(prescription.Robust,'on') && strcmpi(prescription.Regularization,'smoothest')
  error('A robust fit combined with the smoothest possible fit would seem to be incompatible goals')
end
  

% ======================================================
% =========== begin value_check subfunction ============
% ======================================================
function prescription=value_check(prescription,fieldname,legalchar,legalnumeric,legalempty)
% checks that:
% 1. the field exists
% 2. it contains an acceptable character string or a
%    numeric array of valid shape for this field

% first test that the field exists
if ~isfield(prescription,fieldname)
  error(['Field not found: ',fieldname])
end

% grab that field
Field = getfield(prescription,fieldname);

% is an empty array a legal option?
if isempty(Field)
  if ~legalempty
    error(['Empty array is not legal for this field ',fieldname])
  else
    % empty is acceptable, so just return
    return
  end
end

% if its a character string, is it a legal option?
if ischar(Field)
  if ~isempty(legalchar)
    Field = lower(Field);
    ind = strmatch(Field,legalchar);
    if isempty(ind)
      error(['Illegal value for ',fieldname,', see help slmset'])
    elseif length(ind)>1
      error(['Ambiguous value for ',fieldname,', see help slmset'])
    else
      prescription = setfield(prescription,fieldname,legalchar{ind});
    end
  else
    error(['Char is not legal for: ',fieldname,', see help slmset'])
  end
end

% is a numeric array a legal value?
if isnumeric(Field)
  if isempty(legalnumeric)&&~isempty(Field)
    error(['Numeric is not legal for this field ',fieldname])
  else
    % we've already checked about empty fields
    s = size(Field);
    for i=1:length(legalnumeric)
      % verify that the size of Field fits one of the
      % legal size templates in legalnumeric
      
      
      
      
      
      
      
      
      
      
    end
  end
end



