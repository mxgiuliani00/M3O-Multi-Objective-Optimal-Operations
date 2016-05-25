%%
% The following script is a tutorial on the methodology of what I call 
% "Shape Prescriptive Modeling". Its a tool for modeling a function of
% a single variable
%
%    y = f(x) + error
%
% The actual model will be in the form of a least squares spline. Shape
% prescriptive modeling is really a way of thinking about a model - it
% has a lot of very Bayesian underpinnings. I'll get to that in a bit.
% First, lets talk about modeling, who builds models, and why we choose
% a specific model form.
%
% Why do we, as scientists and engineers, fit a curve to data? There are
% really two basic reasons why, if we ignore the common "My boss told me
% to do this." The first such reason is for purely predictive value. We may
% wish to reduce some set of data to a simply accessible functional relationship
% for future use. The process of modeling may allow us to smooth the data,
% producing a smooth curve where only distinct data points fell before.
% As well, this smooth curve is now a platform for interpolation. Note that
% if we have smoothed the data, removing unwanted noise that was present
% in the data, then the interpolant is more an approximant. I.e., one
% will not interpolate the data, but a smoother form of it.
% 
% The second common goal of curve fitting is often to predict/estimate some
% parameter of the relationship, such as a minimum or maximum value where
% the curve approaches an asymptote, or perhaps a maximum slope or a rate
% parameter. 
%
% If we are to fit a model to data, what kind of models do we use, and
% how does the model class reflect our goals in the modeling process?
% To answer this question, I'll first look at the common modeling tools
% we might see. First and most primeval is polyfit. Polyfit allows us to
% fit a polynomial model to our data. Such a model tells us little about
% our process beyond perhaps a simple slope. These models are really of
% value only for very simple relationships, where we are willing to
% invest little in the modeling process. Many of the individuals who
% use polynomial modeling tools do so for the wrong reasons. They use
% higher order polynomial fits to get a better approximation to their
% data, not realizing the problems inherent with those higher order
% polynomials.
% 
% At the other end of the spectrum, one sees individuals using nonlinear
% regression to fit a large variety of curves to data. Exponential terms,
% gaussian-like modes, sigmoid shapes, as well as many more. Did these
% model forms result from valid mechanistic arguments? I.e., has the
% analyst developed a real model of the process under study, where all
% that is left is to estimate the parameters of that model from data?
% 
% Sometimes it is so, but far more often the analyst individual knows
% something about the basic underlying process, and the model chosen
% has the correct basic functional shape. I'll argue that truly valid
% mechanistic models are nearly as rare as hen's teeth. 
% 
% There is a subtle variation on the mechanistic model. I call it the
% metaphoric or metaphorical model. Here we use a mathematical model
% of one process that we do understand, used as a model for a process
% that we do not really understand. My favorite examples of such
% metaphorical models are cubic splines, where a mathematical model
% of a thin flexible beam is used to predict many other processes. A
% second good example is the use of epidemiological models as used to
% predict sales of some product. There are many other examples of course.
% 
% An advantage of a metaphorical model is it may help us to understand/
% predict/even extrapolate future behavior based on our knowledge of the
% metaphor. Many of the nonlinear regression models that we see built in
% Matlab are of the shape variety. I.e., we see a bell shaped curve, so
% we fit a gaussian to it. We see a function with a lower and upper
% constant asymptote, and we fit some variety of logistic function, or
% an erf variant to it. If the fit is inadequate for our purposes, we
% try another functional form. I'll argue that this group of individuals
% is often using nonlinear regression for the wrong reason. They simply
% want some curve that behaves the way they know it should. Examples of
% this might be a photographic film modelor, who knows the film darkens
% as more photons land upon it, or a biologist, who knows that the count
% of bacteria in a petri dish will grow with time. (All of these examples
% have limits, but a good scientist/engineer will know those limits.)
% 
% The fact is, we often understand something about a process we are trying
% to model. What we may not have is the ability to build that knowledge and
% understanding into our model. Even if we know the underlying functional
% relationship is monotone increasing/decreasing, how do we estimate a
% model that has that behavior built into it? This is where shape
% prescriptive modeling shines. It is really a very Bayesian concept.
% It is the process of quantifying our basic knowledge of a process in
% mathematical terms, and then estimating a predictive model of the
% process that has that knowledge built into it. Our knowledge becomes
% a prescription for the final shape of the process. Sometimes we know
% much about a process, in which case we have a detailed prescription.
% Other times we have a very sketchy prescription.
% 
% The slm tools in this toolkit are designed to help the user to embed
% their knowledge into a model using a langauage of shape. One builds
% a prescription for the desired shape, then combines it with their data
% to build a model. These steps can be done distinctly, where an explicit
% prescription is created, then applied with a modeling engine to data,
% or the prescription can be built in a fluid manner using a gui. Here
% a point and click interface can provide the information necessary.
% 
% Some final notes before we begin actual experimentation: This toolkit
% assumes possession of the optimization toolbox. 
% 
% Also, all examples will use the commandline engine for modeling,
% instead of the gui interface. One good reason is that at the time
% of writing this tutorial, I had not finished writing the gui. A
% better reason is that anything that can be done via the gui is
% also possible to do using the computational engine, and cell mode
% favors this approach.
% 
% Finally, ...

%% Some data
% Lets start out by creating some data. A simple functional form that
% we understand the shape of.
n = 100;
x1 = sort(rand(n,1));
y1 = exp(x1) + randn(n,1)/20;
plot(x1,y1,'ro')

%% Polyfit?
% First, lets build a simple model using polyfit.
P = polyfit(x1,y1,6);
xev = 0:.001:1;
yev = polyval(P,xev);

%%
% and plot the results on top of the data
plot(xev,yev,'g-',x1,y1,'ro')

%%
% Are you happy with the resulting curve? Perhaps. Much depends upon
% our goals in building a model for this process. How good of a fit
% do we really need, how much noise is there in the data? Does the
% model behave as we believe it should? Does this model have lack of
% fit? Actually, it looks pretty good to me, at least given the amount
% of noise in the data. But see what happens when I use a high order
% polynomial model.

%% 
% This one was just too much, with condition problems.
P = polyfit(x1,y1,15);
xev = 0:.001:1;
yev = polyval(P,xev);
%%
% and plot the results on top of the data
plot(xev,yev,'g-',x1,y1,'ro')

%% Less data, more noise, with polyfit
% Lets try again, only this time we'll use less data, with a much
% higher level of noise.
close
n = 30;
x2 = sort(rand(n,1));
y2 = exp(x2) + randn(n,1)/3;
plot(x2,y2,'ro')

P = polyfit(x2,y2,6);
xev = linspace(min(x2),max(x2),100);
yhat = polyval(P,xev);

hold on
plot(xev,yhat,'g-')
hold off
%%
% Its a clear example of overfitting this noisy data.

%% A piecewise linear SLM fit
% Now lets try using the slm tools. After all, its why we are
% working through this tutorial.
slm = slmengine(x2,y2,'degree',1,'plot','on')
%%
% note that the result from slmengine is a structure. It contains
% 7 fields, and has much the same information as the 'pp' form does.
%
% 'form' is (in this case) 'slm'. This means that the function
% this structure contains is a slm model. We could have had it
% return a pp form instead.
%
% 'degree' is the polynomial degree of the piecewise segments.
% 1 indicates a piecewise linear function, continuous across the
% breakpoints (or knots.) Note that the splines toolbox tends to
% use the term 'order', one more than the degree.
%
% 'knots' contains the set of knots or breakpoints. The default
% is to use 6 equally spaced knots, from min to max of the data.
%
% 'coef' is the set of coefficients of the spline model. I've used
% a Hermite form here. Its an easy one to understand, as well as
% to look at the coefficients and visualize the shape of the curve.
% The field coef contains the value of the function at the
% corresponding knot.
%
% 'prescription' is the default prescription. A bayesian might call
% it an uninformative prior. Its the list of all settings for every
% optional parameter one can set for a model. Since we specified
% nothing but the data in our call to slmengine, we get all defaults.
%
% Finally, there are fields for x and y. The model thus contains
% the data used to build it. This serves a documentation purpose
% but it also allows plotslm to work well later.


%% Shape prescriptions
% We can look at the parameters in slm.prescription, but now is as
% good a tiem as any to see what slmset does for us. Here are the
% default set of properties.
prescription = slmset

%%
% slmset is the tool that allows us to specify a variety of pieces
% of information to the engine.


%% Use of slmset
% Help slmset will list all the options. It works on a property/value
% pair interface. There are a wide variety of properties that one can set,
% enough to satisfy almost any type of information a user can bring to
% a modeling problem. Here are the options one can set:
help slmset

%% Evaluating a model
% One more thing before we return to fitting data. We can use the 
% models that we build in a predictive mode. 
slmeval(.5,slm)

%%
% slmeval can also differentiate the function (although piecewise
% constant models would have a zero derivative everywhere, with a
% singularity at every knot.) The third argument controls the order
% ofthe derivative to be generated.
slmeval(.5,slm,1)

%%
% and even compute an inverse at any point. 
slmeval(1.75,slm,-1)

%%
% We could plot the function as fit uing slmeval,
yhat = slmeval(xev,slm,0);
plot(x2,y2,'ro',xev,yhat,'b-')

%% Plotting a model using plotslm
% Or we could have plotted it using a nice utility that I've provided.
% Plotslm plots the curve, the data that it was built from, and the
% knots/breaks of the spline model. Plotslm also allows you to plot
% the residuals of the model to the data, or the derivatives of the
% curve. Finally, you can control many features of this plot using the
% plot menu.
plotslm(slm)

%% slmengine can plot too
% Or we could have just told our engine to plot after it finished
% estimating its curve fit. This time we will fit it as a cubic spline.
close all
slm = slmengine(x2,y2,'degree',3,'knots',4,'plot','on');

%% Computations on a spline model, min max, integral, etc.
% plotslm can plot the function derivatives or integral, and
% slmeval can evaluate the function, its derivatives or the inverse
% function at specific points. However, it is occasionally useful to be
% able to find the minimum (or maximum) value of a spline over some interval, or
% find the minimum or maximum value of the slope, or compute the integral of the
% curve between a pair of points.
%
% The slmpar function does these things, both for a slm form as well as a pp form.
x = -5:.1:5;
y = exp(-x.^2/2)/sqrt(2*pi) + randn(size(x))/100;
slm = slmengine(x,y,'plot','on','knots',15,'integral',1,'deg',3,'minval',0)
%%
% Of course, here the integral was forced to 1, so the integral should be correct.
res = slmpar(slm,'in')
%%
% The maximum value should occur very near zero.
[res,loc] = slmpar(slm,'maxfun')

% The next few examples show slmpar at work on a pp form, as created by spline:
x = linspace(-pi,pi);
y = sin(x);
pp = spline(x,y);
%%
% The overall integral should be zero, or as close as floating point
% computations can give us for the integral of a spline.
res = slmpar(pp,'integral')
%%
% The integral from 0 to pi should be 2.
res = slmpar(pp,'integral',[0,pi])
%%
% The maximum of a sine wave over its period is 1, and should occur at
% x = pi/2.
[res,loc] = slmpar(pp,'maxfun')
%%
% The minimum of a sine wave over [-1,1] should occur at the lower end point.
[res,loc] = slmpar(pp,'minfun',[-1,1])

%%
% A nice thing about plotslm is it overlays the knot positions
% on top of the curve. The knots are at the vertical green dotted
% lines.

%% Curvature constraints
% Having seen the fit above, one might wonder why one would use
% a shape prescriptive model at all. It appears to be little better
% than that 6'th order polynomial model we fit before using polyfit.
%
% What do we know about the underlying functional relationship?
% We want to think in terms of fundamental shape primitives.
%
% I'd suggest that the function is known to be positively curved.
% Its was a simple exponential after all. The second derivatve
% should never be negative.
close all
slm = slmengine(x2,y2,'plot','on','concaveup','on');

%% Monotonicity constraints
% Having seen the fit above, one might wonder why one would use
% a shape prescriptive model at all. It appears to be little better
% than that 6'th order polynomial model we fit before using polyfit.
%
% What do we know about the underlying functional relationship?
% We want to think in terms of fundamental shape primitives.
%
% I'd suggest that the function is known to be positively curved.
% Its was a simple exponential after all. The second derivatve
% should never be negative.
close all
slm = slmengine(x2,y2,'plot','on','increasing','on');

%% Both curvature and monotonicity constraints
% Having seen the fit above, one might wonder why one would use
% a shape prescriptive model at all. It appears to be little better
% than that 6'th order polynomial model we fit before using polyfit.
%
% What do we know about the underlying functional relationship?
% We want to think in terms of fundamental shape primitives.
%
% I'd suggest that the function is known to be positively curved.
% It was a simple exponential after all. The second derivatve
% should never be negative.
close all
slm = slmengine(x2,y2,'plot','on','concaveup','on','increasing','on');

%% Minimum slope or maximum slope constraints
% Try an erf function, with just a bit of noise. But then, show that the
% curve can be constrained to have an overall minimum and maximum slope.
n = 100;
x = 4*rand(n,1) - 2;
y = erf(x) + randn(size(x))/100;
%%
% First, with no slope constraints at all
slm = slmengine(x,y,'plot','on','knots',15);
%%
% Plot the first derivative. See that it should strictly be bounded by
% in the interval [0,1].
xev = linspace(-2,2,1000);
d = slmeval(xev,slm,1);
plot(xev,d,'-')
grid on

%%
% Redo the fit, but with just a touch of a constraint on the slopes.
slm = slmengine(x,y,'plot','on','knots',15,'minslope',0.1,'maxslope',0.9);
d = slmeval(xev,slm,1);
figure
plot(xev,d,'-')
grid on


%% Extrapolation in the fit, non-uniform knot spacings
% We might also decide to fit the curve over a wider knot range,
% extrapolating here beyond the data. Note that it is BETTER to build a
% spline to have the shape you want over a desired region then to try
% to extrapolate a spline that was built over a smaller region.
%
% In fact, I try to force the user to do no extrpolation at all,
% unless they insist on doing so. The extrapolation options will
% be shown in some detail later on.
close all
slm = slmengine(x2,y2,'plot','on','concaveup','on','increasing','on','knots',[-1,0:.2:1,2]);

%% Piecewise constant models
% Or even to use a lower order model. This version of my toolkit only
% allows piecewise constant or piecewise linear besides the cubic models.
%
% Also, for brevity, we can always shorten any property names,
% as long as the name remains unambiguous. Capitalization is ignored.
slm = slmengine(x1,y1,'plot','on','kn',0:.1:1,'deg',0);

%%
% Note that I also used more knots in this last fit. I also returned
% to our original set of data with 100 points.

%% Fixed specific values
% There are other shape parameters we could have set in our prescription.
% Kowing that the underlying relationship is an exponential, we might
% remember that at x == 0, that y == 1. Specific values like this are
% often known. If our data was from an MTF curve, we might choose to
% assume a modulation of 100% at a frequency of 0. Likewise, population
% growth data might have a known population at time == 0.
slm = slmengine(x1,y1,'plot','on','concaveup','on','knots',0:.2:1,'leftvalue',1);

%% Flawed information, "Garbage in, garbage out"
% What if our information is flawed? The old saying always applies.
slm = slmengine(x1,y1,'plot','on','concaveup','on','knots',0:.2:1,'leftvalue',2);

%%
% Or suppose we specify a bit too small a value for the maximum slope. The
% result is lack of fit.
slm = slmengine(x1,y1,'plot','on','maxslope',1.5,'knots',0:.1:1);

%%
% Garbage in, garbage out. If we specify a monotone decreasing function
% on increasing data, then we have wasted a lot of machine cycles to
% compute the mean of our data.
slm = slmengine(x1,y1,'plot','on','decreasing','on','knots',0:.2:1,'degree',1);

%%
% We could have done the same with somewhat less work...
slm = slmengine(x1,y1,'plot','on','degree',0,'knots',2);

%%
% The really lazy would just use mean(y1). Perhaps efficient is a
% better word.

%% 
% lets look at a function with slightly more interesting behaviour. 
% still keep it simple enough that we know what it should look like.
close all
n = 500;
x = linspace(0,1,n)';
y = sin(x*pi) + randn(n,1)/5;

slm = slmengine(x,y,'plot','on','knots',10);

%%
% If I remember correctly, sin(x) will have a negative second derivative
% over that interval. If I choose to apply this information to our
% curve, the fit gets quite good.
slm = slmengine(x,y,'plot','on','knots',10,'concavedown','on');

%%
% A simple way of specifying a curve shape with a single maximum
% value at some point is just this:
slm = slmengine(x,y,'plot','on','knots',0:.1:1,'simplepeak',.5);

%%
% all it does is to automatically specify increasing and decreasing
% regions of the curve.

%%
% lets change our data again now. Make it a fairly noisy sine wave.
clear
close all
n = 250;
x = linspace(0,1,n)';
y = sin(2*x*pi) + randn(n,1)/2;

%%
% this time we might choose to assume periodicity.
slm = slmengine(x,y,'plot','on','knots',0:.1:1,'endconditions','periodic');

%%
% That last curve had a lot of noise in it, I might not be surpised
% if it was a little bumpy. We also know that somewhere in the middle
% of the curve, this function goes through an inflection point. I'll
% just pick somehere to put that inflection point. Note that the
% curvature would be negative to the left, and positive to the right,
% so set a positive inflection point.
slm = slmengine(x,y,'plot','on','kn',0:.1:1,'endc','per','positiveinflection',.5);

%%
% There are certainly many ways to define a prescription for the 
% shape of a curve. We have only covered a few of them here.
%
% Lets return to slmset though. Note that the two lines below are
% the same as the direct call to slmengine.
prescription = slmset('plot','on','kn',0:.1:1,'endc','per','positiveinflection',.5);
slm = slmengine(x,y,prescription);

%% 
% The model as it is returned includes a field with the shape prescription
% that defined that model. It also includes the original data in named fields.
% Each of these fields is useful for documentation purposes.
slm

%%
% We can also modify an existing prescription, adding new pieces of information.
slm = slmengine(x,y,prescription,'maxvalue',0.5);

%% Too many knots?
% Sometimes we have little data, but use too many knots. A
% traditional spline model will often fail in this case, with
% a singular matrix warning. Slmengine will not fail, although
% you may not always be happy with the result. Remember, when
% there is little information brought to a problem by the data,
% a reasonable fit may still be obtained if you can add enough
% information from your own kowledge of the process. In this
% particular example, I'm tryin only to show that the use of too
% many kots did not result in a singularity.
clear
close all
n = 10;
x = rand(n,1);
y = sin(2*x*pi) + randn(n,1)/10;
slm = slmengine(x,y,'plot','on','knots',0:.05:1);

%% 
% Fewer knots may be better, even though we did survive having
% way too many knots.
slm = slmengine(x,y,'plot','on','knots',[0 1]);

%%
% And while it may be a little slow, one can use cross
% validation to help smooth your curve. Don't try it with too
% many knots or too many points.
slm = slmengine(x,y,'plot','on','knots',0:.05:1,'regularization','cross');

%%
% Another option, if you had a decent estimate of the error
% standard deviation, is to supply that value. slmengine will
% build a fit with the desired residual error.
slm = slmengine(x,y,'plot','on','knots',0:.05:1,'regularization',-0.10);

%% Use in practice
% for some more realistic data, try this one out.
% The function is known to be increasing, and everywhere non-negative.
% There is a little curvature at the bottom that is probably real.
% At the top end, we expect the curve to be moderately well behaved
% so extrapolation to x = 255 is linear.
close all
load scann
slm = slmengine(R,L,'plot','on','knots', ...
  [0 5 10 15 20 30 50 80 120 160 200 255],'leftminvalue',0,'increasing','on','deg',3,'linearregion',[200,255])

%%
% I often see people using splines wrongly. For example, given a
% relationship with a derivative singularity, how best to build a
% spline model?
x = linspace(0,1,20);
truefun = @(x) nthroot(x,5);
y = truefun(x);
xfine = linspace(0,1,1000);
yfine = truefun(xfine);
plot(x,y,'-o')
title 'A singularity at x = 0'
xlabel X
ylabel Y
%%
% Obviously, fitting a spline to this data, in the form y(x), we
% should expect a problem. Even with no noise present in the system.
slm = slmengine(x,y,'plot','on','knots',10)
hold on
plot(xfine,yfine,'b--')
hold off

%%
% We might force the spline to be monotone, or even choose a non-uniform
% knot spacing, but the derivative singularity at x == 0 is too much for
% polynomial segments to handle well.
slm = slmengine(x,y,'plot','on','knots',[0 .01 .02 .05 .1 .3 .6 1],'incr','on','concavedown','on')
hold on
plot(xfine,yfine,'b--')
hold off

%%
% A better solution is to recognize the essential singularity in the problem.
% Here, swap the dependent and independent variables. What was once a singularity
% is no longer so. Se that the spline fit and the original function now overlay
% each other almost perfectly here.
slm = slmengine(y,x,'plot','on','knots',10,'incr','on','concaveup','on','deg',3)
hold on
plot(yfine,xfine,'b--')
hold off

%%
% Of course, prediction will be more difficult, since the resulting model
% must be interpolated in an inverse form. slmevel can do so. The result
% here would be 0.7, IF the spline were an exact representation of the
% fifth root of its input. Since we live in a floating point world, the
% result is as close as I might expect.
format long g
slmeval(0.7.^5,slm,-1)
format short g

%% Returning a set of points, evaluated through the fitted curve.
% SLM has the ability to generate a list of points, then evaluate
% them through the resulting model. While this is easily enough done
% using slmeval, this is a friendly option that will save a set for
% the user. By default, 1001 points are generated along the curve,
% however, any number of equally spaced points can be generated. For
% the user who wishes for an unequal (user specified) spacing, there
% is always slmeval to fall back upon. In this example, only 11 equally
% spaced points are generated along the curve to be returned. The user
% can then plot these points as they choose. Of course, plotslm has
% significantly more functionality than a simple plot, since you can
% view more than just the simple function values.
close all
format long g
[slm,xp,yp] = slmengine(y,x,'knots',10,'incr','on','concaveup','on','rightvalue',1,'predictions',11)

%% Combining various constraints to produce a better model
% On this curve, we know the actual form. It should be a hyperbolic shape,
% that approaches a linear asymptote.
x = rand(40,1);
y = x - 1./(10*(x+.1)) + randn(size(x))/25;
plot(x,y,'o')

%%
% We can model the curve as simply a monotone increasing function.
% I'll use more knots than I really need to exacerbate the problem.
% (There really are way too many knots in this spline, possibly causing
% numerical problems.)
slm = slmengine(x,y,'incr','on','knots',35,'plot','on')

%%
% Now, try enforcing both monotonicity, and a concave down function, 
% (a negative second derivative.) 
slm = slmengine(x,y,'incr','on','concavedown','on','knots',35,'plot','on')

%%
% Our second try is better, but still not well enough behaved. See that
% by requiring that the jerk be positive everywhere, we get a nice, smooth
% curve, even with far too many knots for this simple curve. In fact, it
% even appears to be approaching a linear asymptote.
slm = slmengine(x,y,'incr','on','jerk','positive','concavedown','on','knots',35,'plot','on')

%% Extrapolation of a spline model
% Normally, I try to push a user to NOT extrapolate a spline model.
% Cubic polynomial segments do squirrely things away from where they
% were used in a fit. So if you seriously tried to extrapolate a cubic
% spline, then don't be surprised to see it do silly things. For this
% reason the default for SLMEVAL is to only extrapolate as a constant
% function. Thus below the first knot or above the top knot, the spline
% prediction will take on the values of the function at the corresponding
% knot, but no more. Essentially only constant extrapolation is done.
%
% SLMSET does provide for other options in extrapolation, but only if you
% ask for it. The 'extrapolation' prescription may take on the values:
%
%  {'error', 'warning', 'constant', 'linear', 'cubic', 'NaN'}
%
% The default is 'constant'. For this case, constant extrapolation is
% done with with no warning generated for all points outside of the knots.
% The other options seem pretty stratight forward.
%
%  'error' generates an error if any point falls outside the knots.
%
%  'warning' generates a warning message, then does constant extrapolation.
%
%  'constant' does constant extrapolation. (The default)
%
%  'linear' does linear extrapolation.
%
%  'cubic' does cubic extrapolation.
%
%  'nan' inserts NaN values for all points outside the knots.
close all
x = 0:.1:1;
y = sin(2*pi*x) + randn(size(x))/100;

%%
% See that the Extrapolation field has been set to 'constant' by default.
slm = slmengine(x,y,'knots',x,'plot','on')

%%
% See what happens when various extrapolation modes are specified though.
% Note that the extrapolation type is indicated when you create the spline
% itself. I won't show what happens when 'error' is chosen, since that will
% screw up this published document. You can guess though.
xev = -1:.01:2;

%%
% Warning throws a warning if extrapolation is tried, then does constant
% extrapolation, like the default.
slm = slmengine(x,y,'knots',x,'extrapolation','warning');
yev = slmeval(xev,slm,0);
plot(x,y,'o',xev,yev,'r-')
title 'Constant extrapolation (warning generated)'

%%
% Linear simply does linear extrapolation, with no warnings or errors.
slm = slmengine(x,y,'knots',x,'extrapolation','linear');
yev = slmeval(xev,slm,0);
plot(x,y,'o',xev,yev,'r-')
title 'Linear extrapolation'

%%
% Cubic simply does cubic extrapolation, with no warnings or errors. I
% would rarely ever recommend 'cubic' extrpolation, certainly not if you
% extraapolate too far. See what happens here.
slm = slmengine(x,y,'knots',x,'extrapolation','cubic');
yev = slmeval(xev,slm,0);
plot(x,y,'o',xev,yev,'r-')
title 'Cubic extrapolation'

%% You can get information about the fit itself
% slmengine always returns this information in the stats field of the returned model.
% However, the verbosity propert has three levels {0,1,2}. 0 tells
% slmengine to be quiet, with no output to the command line. 1 tells it to
% report fit statistics at the command line,
close all
slm = slmengine(x,y,'verbosity',1);

%%
% and verbosity of 2 tells slmengine to be more expressive yet.
slm = slmengine(x,y,'verbosity',1);

%% Knot placement optimization
% You can let an optimizer choose the knots. It may give you better
% results than a blind choice of knots. But, the fact is all optimizers can
% get lost at times.
close all
x = -5:.1:5;
y = erf(x) + randn(size(x))/100;
slm = slmengine(x,y,'increasing','on','plot','on','degree',1, ...
  'interiorknots','free','knots',10)

%%
% If however, you have some idea where the knots might belong based on the
% shape of your function, then it is a good idea to provide some input.
% Here for example, there is little happening with this function on both
% ends of the curve. So why waste the knots out there?
slm = slmengine(x,y,'increasing','on','plot','on','degree',1,...
  'interiorknots','free','knots',[-5 -1.25 -1 -.75 -.5 .5 .75 1 1.25 5])

%%
% slmpar has some other options, that are not terribly useful for
% computation, but are nice to understand how splines work, as well as
% how the coefficients are stored in these spline containers. These last
% options in slmpar require the symbolic toolbox.
%
% First, I'll create a piecewise linear spline with 4 segments.
x = linspace(-pi,pi,100);
y = sin(x);
slm = slmengine(x,y,'knots',5,'result','slm','degree',1)

%%
% The SLM storage form can be called a Hermite form. It contains a list of
% knots, and the function values of the curve at those knots. So we see
% below the first column isthe list of knots, while the second column the
% values of the sine function approximation at those knots.
[slm.knots,slm.coef]

%%
% However, we can see the actual linear polynoial segments using slmpar.
% This option tells slmpar to create a cell array that lists the knot
% intervals over which each segment is defined, as well as the actual
% linear polynomial used over that interval, as a symbolic polynomial
% in the variable x.
p = slmpar(slm,'symabs')
p{:}

%%
% We can use slmpar to understand the pp form of a spline, and how the
% spline coefficients are stored. They are stored as the coefficients of a
% polynomial in standard matlab poly form, so highest order coefficient
% first. However, the pp form uses what I'll call a relative polynomial.
pp = slm2pp(slm)

%%
% We can list out the breaks (or knots if you choose) of the spline.
% between each pair of breaks, there is here a linear polynomial segment
% defined. However, the i'th such polynomial will be evaluated as if it
% lives in the interval [0,pp.breaks(i+1) - pp.breaks(i)]
pp.breaks
pp.coefs

%%
% So the first linear polynomial segment actually lives on the interval
[pp.breaks(1),pp.breaks(2)]

%%
% However, it is evaluated on the translated interval
[0,pp.breaks(2)-pp.breaks(1)]

%%
% This translation of the interval is not terribly important when one is
% working with linear spline segments, but when cubic (or even higher order)
% segments are employed, there may well be numerical issues if MATLAB is
% forced to take powers of very large numbers. Suppose for example, our
% spline was forced to be evaluated for numbers on the order of 1e10? Now
% forming the cube of such large numbers will cause massive loss of
% accuracy in the floating point numbers.

%%
% Similarly, we might create a cubic spline. Here I'll use only 4 knots,
% thus three cubic segments, but I'll generate data that lies in a
% different domain.
x = linspace(1e6,1e6 + 3,100);
y = sin(x);
slm = slmengine(x,y,'knots',4,'result','slm','degree',3,'plot','on')

%%
% First, lets understand the Hermite form here. The first column below is
% the list of knots. The second column the function values at those knots,
% and the third column contains the first derivatives of the spline function
% at those knots. Personally, I've always had a preference for this form
% because one can easily visualize much about the shape of the cubic
% segments directly from those coefficients.
[slm.knots,slm.coef]

%%
% Converting the spline to its pp form, we get:
pp = slm2pp(slm)

%%
% There are 4 breaks, so 3 segments. These coefficients presume the relative
% form for evaluation
pp.coefs
crel = slmpar(pp,'symrel')
crel{:}

%%
% See however, what happens to the coefficients of the polynomials in the
% absolute form. Those coefficients are now a bit nasty for a cubic
% polynomial. If x is on the order of 1e6, I would expect to see a serious
% loss of accuracy.
cabs = slmpar(pp,'symabs')
cabs{:}

%%
% To compare the accuracy problems one would expect, what is the sin of
% a moderately large number?
sin(sym('1e6 + 0.5'))
sin(1e6 + 0.5)
slmeval(1e6 + 0.5,slm)

%%
% So the 4 knot spline fit itself was only accurate to about 3 significant
% digits. Lets now use polyval to evaluate those polynomials at that same
% effective location. The relative form returns a value that matches
% slmeval exactly, as well as being a viable value for the sine function.
% However, the absolute form just returns numerical junk.
polyval(fliplr(double(coeffs(crel{2,1}))),0.5)
polyval(fliplr(double(coeffs(cabs{2,1}))),1e6 + 0.5)

close all





