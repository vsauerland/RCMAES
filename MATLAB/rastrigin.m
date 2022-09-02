function val = rastrigin( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = 20 + x ^ 2 + y ^ 2 - 10 * cos( 2 * pi * x ) - 10 * cos( 2 * pi * y );
else
	val = 10 * length( x );
	for i = 1 : length( x )
		val = val + x( i ) ^ 2 - 10 * cos( 2 * pi * x( i ) );
	end
end
