function val = testf( fname, x, y )
% TESTF different two-dimensional test functions
%
% val = testf( fname, x, y )
%
% fname  name of fitness function
% x      first scalar argument or two-dimensional argument vector
% y      optional second scalar argument
%--------------------------------------------------------
if nargin == 3
	val = feval( fname, x, y );
else
	val = feval( fname, x );
end

function val = ackley( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = -20 * exp( -0.2 * sqrt( 0.5 * ( x ^ 2 + y ^ 2 ) ) ) - exp( 0.5 * ( cos( 2 * pi * x ) + cos( 2 * pi * y ) ) ) + exp( 1 ) + 20;
elseif nargin == 1 && length( x ) == 2
	val = -20 * exp( -0.2 * sqrt( 0.5 * ( x( 1 ) ^ 2 + x( 2 ) ^ 2 ) ) ) - exp( 0.5 * ( cos( 2 * pi * x( 1 ) ) + cos( 2 * pi * x( 2 ) ) ) ) + exp( 1 ) + 20;
else
	error( 'Ackley function requires 2 scalars or a 2-vector as argument' );
end

function val = griewank( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = 1 + ( x ^ 2 + y ^ 2 ) / 4000 - cos( x ) * cos( 1 / sqrt( 2 ) * y );
else
	val = 1 + sum( x .^ 2 ) / 4000;
	product = 1;
	for i = 1 : length( x )
		product = product * cos( 1 / sqrt( i ) * x( i ) );
	end
	val = val - product;
end

function val = griewank2drot( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	% first rotate by 90 degree around first axis
	xr = x * cos( pi / 4 ) - y * sin( pi / 4 );
	yr = x * sin( pi / 4 ) + y * cos( pi / 4 );
	% then apply griewank	
	val = 1 + ( xr ^ 2 + yr ^ 2 ) / 4000 - cos( xr ) * cos( 1 / sqrt( 2 ) * yr );
elseif nargin == 1 && length( x ) == 2
	% first rotate by 90 degree around first axis
	xr = x( 1 ) * cos( pi / 4 ) - x( 2 ) * sin( pi / 4 );
	yr = x( 1 ) * sin( pi / 4 ) + x( 2 ) * cos( pi / 4 );
	% then apply griewank	
	val = 1 + ( xr ^ 2 + yr ^ 2 ) / 4000 - cos( xr ) * cos( 1 / sqrt( 2 ) * yr );
else
	error( 'griewank2dot function requires 2 scalars or a 2-vector as argument' );
end

function val = himmelblau( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = ( x ^ 2 + y - 11 ) ^ 2 + ( x + y ^ 2 - 7 ) ^ 2;
elseif nargin == 1 && length( x ) == 2
	val = ( x( 1 ) ^ 2 + x( 2 ) - 11 ) ^ 2 + ( x( 1 ) + x( 2 ) ^ 2 - 7 ) ^ 2;
else
	error( 'Himmelblau function requires 2 scalars or a 2-vector as argument' );
end

function val = linear1( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = -x - y;
else
	val = -sum( x );
end

function val = rastrigin( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = 20 + x ^ 2 + y ^ 2 - 10 * cos( 2 * pi * x ) - 10 * cos( 2 * pi * y );
else
	val = 10 * length( x );
	for i = 1 : length( x )
		val = val + x( i ) ^ 2 - 10 * cos( 2 * pi * x( i ) );
	end
end

function val = rosenbrock( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = ( 1 - x ) ^ 2 + 100 * ( y - x ^ 2 ) ^ 2;
else
	val = 0;
	for i = 1 : length( x ) - 1
		val = val + ( 1 - x( i ) ) ^ 2 + 100 * ( x( i + 1 ) - x( i ) ^ 2 ) ^2;
	end
end

function val = sphere( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = x ^ 2 + y ^ 2;
else
	val = sum( x .^ 2 );
end
