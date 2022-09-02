function [ k1, k2, x1, x2, f1, f2, expF1, expF2 ] = rcmaesVScmaes( fname, lb, ub, i, lambda, maxIter )
% rcmaesVScmaes( fname, lb, ub, i, lambda, maxIter )
% compares R-CMA-ES result with CMA-ES result
%
% INPUT
% fname		objective function name
% lb		lower bounds on variables
% ub		upper bounds on variables
% i		index of stochastic component with fix std.-dev. (ub(i)-lb(i))/2
% lambda	number of samples per iteration
% maxIter	number of optimization iterations
%
% OUTPUT
% k1, k2	number of iterations
% x1, x2	solutions found by cmaes resp. rcmaes
% f1, f2	f(x1) resp. min(f(x2))
% expF1, expF2  expectation for f(x) when x(i) is N(m,s)-distributed
%		with m=(lb(i)+ub(i))/2 and s=(ub(i)-lb(i))/2
%		and the other components of x agree with x1 resp. x2

opts.LBounds = lb';
opts.UBounds = ub';
opts.DispModulo = 0;
opts.DispFinal = 'off';
opts.LogModulo = 0;
opts.PopSize = lambda;
opts.maxIter = maxIter;
opts.TolX = 1e-4;
m = ( lb + ub ) / 2;
s = 0.3 * ( ub - lb )';
[ x1, f1, k1 ] = cmaes( fname, m, s, opts );
expF1 = expectedFitness( fname, x1, i, lb( i ), ub( i ) );

[ x2, k2 ] = rcmaes( fname, lb, ub, i, lambda, maxIter );
xLine = linspace( lb( i ), ub( i ), 10000 );
x = x2;
for k = 1 : 10000
	x( i ) = xLine( k );
	fLine( k ) = feval( 'testf', fname, x );
end
f2 = min( fLine );	
expF2 = expectedFitness( fname, x2, i, lb( i ), ub( i ) );
