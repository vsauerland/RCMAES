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
