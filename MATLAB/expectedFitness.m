function expF = expectedFitness( fname, x, iStoch, lb, ub )

% expectedFitness calculate "expected" objective value
%
% expF = expectedFitness( fname, x, iStoch, lb, ub )
%
% OUTPUT
% expF   "expected" objective value
%
% INPUT
% fname  name of objective function, e.g., 'linear' or 'rosenbrock'
% x      some fix argument vector
% iStoch component of x to be varied
% lb, ub lower and upper limit of variation
% 
% The expectation is calculated w.r.t. a modification of N(m,s)
% where m=(lb+ub)/2 and s=(ub-lb)/2,
% see our paper draft for the exact definition of the distribution

expF = 0;
x( iStoch ) = ( lb + ub ) / 2; % mean of component interval
for k = 1 : 100000
	w = x;
	r = randn;
	while r < -1 || r > 1
		r = randn;
	end;
	w( iStoch ) = x( iStoch ) + r * ( ub - lb ) / 2;
	expF = expF + feval( 'testf', fname, w );
end
expF = expF / 100000;
