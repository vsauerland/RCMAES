function [ xmean, countIter ] = rcmaes( fname, lb, ub, iStoch, lambda, maxIter, illustrate, colorLimits )
% RCMAES: Evolution Strategy with Covariance Matrix Adaptation allowing for a ramdom parameter
%
% [ xmean, countIter ] = rcmaes( fname, lb, ub, iStoch, lambda, maxIter, illustrate, colorLimits )
%
% OUTPUT
% xmean		mean of final normal distribution = optimized x value
% countIter	number of iterations
%
% INPUT
% fname         objective function name, default 'rosenbrock'
% lb            lower bounds on variables, default [0 0]
% ub            upper bounds on variables, default [1 1]
% iStoch        index of stochastic component with fix std.-dev. (ub-lb)/2, default 2
% lambda        number of samples per iteration, default 10
% maxIter	number of optimization iterations, default 50
% illustrate    indicates if an illustration of the optimization
%               procedure is desired, default 0
% colorLimits   minimum and maximum color scale limits for objective function, default [0 100]
%
% EXAMPLE with illustration
% xmean = rcmaes( 'rosenbrock', [0 0], [1 1], 2, 10, 50, 1, [0 100] );
%
% USES "EXTERNAL": testf and hybridRank
% USES "INTERNAL": xIntoUnitCube, scale

%rng( 'default' ); % initialize random number generator
if ( nargin < 1 )
	fname = 'rosenbrock';
end
if ( nargin < 3 )
	lb = [ 0 0 ];
	ub = [ 1 1 ]
end
if ( nargin < 4 )
	iStoch = 2;
end
if ( nargin < 5 )
	lambda = 10;
end
if ( nargin < 6 )
	maxIter = 50;
end
if ( nargin < 7 )
	illustrate = 0;
end
if ( nargin < 8 )
	colorLimits = [ 0 0 ];
end

% Strategy parameter setting: Selection
N = length( lb );
N1 = N - 1;
isStochastic = zeros( N, 1 );
isStochastic( iStoch ) = 1; % indicator vector for stochastic index
% lambda = 4 + floor( 3 * log( N ) ); % population size, offspring number
mu = lambda / 2;
% VSA: fix even population size and good half portion mu (lambda = 2*mu)
weights = log( mu + 1 / 2 ) - log( 1 : mu )'; % muXone recombination weights
%mu = floor( mu ); % number of parents/points for recombination
weights = weights / sum( weights ); % normalize recombination weights array
%weights = 1 / mu * ones( mu, 1 ); % VSA: equal weights alternative should suffice
mueff = sum( weights ) ^ 2 / sum( weights .^ 2 ); % variance-effective size of mu
s = 0.3; % coordinate wise standard deviation (step-size)
chi = N ^ 0.5 * ( 1 - 1 / ( 4 * N ) + 1 / ( 21 * N ^ 2 ) );
% chi is expectation of || N( 0, I ) ||

% Strategy parameter setting: Adaptation
cs = ( mueff + 2 ) / ( N + mueff + 5 ); % t-const for cumulation for s control
%ds = 1 + 2 * max( 0, sqrt( ( mueff - 1 ) / ( N + 1 ) ) - 1 ) + cs; % damping for s
ds = 1 + cs; % (the same as above with, e.g., mu=5, equal weights, and N>=3)
cc = ( 4 + mueff / N ) / ( N + 4 + 2 * mueff / N ); % time constant for cumulation for C
c1 = 2 / ( ( N + 1.3 ) ^ 2 + mueff ); % learning rate for rank-one update of C
cmu = 2 * ( mueff + 1 / mueff - 2 ) / ( ( N + 2 ) ^ 2 + mueff ); % learning rate for rank-mu update
cmu = min( 1 - c1, cmu );

xmean = 0.5 * ones( N, 1 ); % initial mean 

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros( N1, 1 ); % evolution path for C
ps = zeros( N1, 1 ); % evolution path for s
B = eye( N1 ); % B defines the coordinate system
D = eye( N1 ); % diagonal matrix D defines the scaling
C = eye( N1 ); % covariance matrix

% -------------------- Generation Loop --------------------------------
countIter = 0;
maxStdDevX = 1;
xmeanOld = 2 * xmean;
while countIter < maxIter && ( maxStdDevX > 1e-4 || norm( xmean - xmeanOld ) > 1e-4 )
	xmeanOld = xmean;
	countIter = countIter + 1;
	% Generate and evaluate lambda offspring
	z = zeros( N, lambda );
	y = zeros( N, lambda );
	x = zeros( N, lambda );
	arfitness = zeros( lambda, 1 );
	for i = 1 : N1
		if B( 1, i ) < 0
	        	B( :, i ) = -B( :, i );
		end
	end
	% Sample probability distribution
	% (pseudo code lines 6-10)
	for k = 1 : mu
		z( : , k  ) = randn( N , 1 ); % standard normally distributed vector
		% z( :, mu + k ) = randn( N, 1 );
		% keep the stochastic component within its bounds:
		while z( iStoch, k ) < -1 || z( iStoch, k ) > 1
			z( iStoch, k ) = randn;
		end;
		% let z( :, mu + k ) be the mirrored (at stochastic axis) version of z( :, mu + k ):
		z( :, mu + k ) = -z( :, k );
		z( iStoch, mu + k ) = z( iStoch, k );
	end
	% (pseudo code lines 11-15)
	for k = 1 : lambda
                y( ~isStochastic, k ) = B * D * z( ~isStochastic, k );
		x( ~isStochastic, k ) = xmean( ~isStochastic ) + s * ( B * D * z( ~isStochastic, k ) );
		y( iStoch, k ) = z( iStoch, k ); % actually not used (in calculation of Cmu below)
		x( iStoch, k ) = 0.5 + 0.5 * z( iStoch, k );
		[ tx, ti ]  = xIntoUnitCube( x( :, k ) ); % closest feasible point
		xdiff = tx - x( :, k );
		arfitness( k ) = feval( 'testf', fname, scale( x( :, k ), lb, ub ) ) ... % objective function call
				+ 1e6 * norm( xdiff( ~isStochastic ) ) ^ 2; % penalty term (simple boundary handling)
	end
	maxStdDevX = max( std( x( ~isStochastic, : )' ) );
	arindex = hybridRank( arfitness, x, zeros( N, 1 ), ones( N, 1 ), isStochastic, mu );
	arfitness = arfitness( arindex );
	z_mu = z( :, arindex( 1 : mu ) ); % the "good z_k-s"
	y_mu = y( :, arindex( 1 : mu ) ); % the "good y_k-s"
	x_mu = x( :, arindex( 1 : mu ) ); % the "good x_k-s"

	% Update mean (pseudo code lines 17-19)
	xmean = x_mu * weights; % recombination
	zmean = z_mu * weights; % = D^ -1 * B * ( xmean - "former xmean" ) / s	
	xmean( iStoch ) = 0.5;
	zmean( iStoch ) = 0;   
 
	% Update evolution paths (pseudo code lines 20-22) 
	ps = ( 1 - cs ) * ps + ( sqrt( cs * ( 2 - cs ) * mueff ) ) * ( B * zmean( ~isStochastic ) );
	pc = ( 1 - cc ) * pc + sqrt( cc * ( 2 - cc ) * mueff ) * ( B * D * zmean( ~isStochastic ) );
    
	% Update covariances and scaling (pseudo code lines 23-26)
	s = s * exp( ( cs / ds ) * ( norm( ps ) / chi - 1 ) );
	Cmu = y_mu( ~isStochastic, : ) * diag( weights ) * y_mu( ~isStochastic, : )';
	C1 = pc * pc';
	C = ( 1 - c1 - cmu ) * C ... % regard old matrix
	    + c1 * C1 ... % plus rank one update
	    + cmu * Cmu'; % plus rank mu update
	C = triu( C ) + triu( C, 1 )'; % enforce symmetry
	% Derive B and D from C (pseudo code line 27)
	[ B, D ] = eig( C ); % eigen decomposition, B="normalized eigenvectors", D="eigen values"
	% obtained D corresponds to D^2 in our paper and Hansens tutorial, 2016
	D = diag( sqrt( diag( D ) ) ); % D contains roots of eigenvalues, now
%	disp( [ num2str( countIter ) ': ' num2str( sum( arfitness ) / N ) ] );
end % while, end generation loop
xmean = scale( xmean, lb, ub );

% ------------------ Auxiliary Functions -----------------------
function x_scaled = scale( x, lb, ub )
% scale sample from [0,1]^N to original feasible box to apply objective function
if size( lb, 1 ) ~= size( x, 1 )
	lb = lb';
	ub = ub';
end
x_scaled = lb + ( ub - lb ) .* x;
% ---------------------------------------------------------------
function [ x, ix ] = xIntoUnitCube( x ) % map x to closest point in [0,1]^N
N = length( x );
lb = zeros( N, 1 );
ub = ones( N, 1 );
ix = x < lb;
x( ix ) = lb( ix );
ix2 = x > ub;
x( ix2 ) = ub( ix2 );
ix = ix2 - ix;
