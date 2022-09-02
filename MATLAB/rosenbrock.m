function val = rosenbrock( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = ( 1 - x ) ^ 2 + 100 * ( y - x ^ 2 ) ^ 2;
else
	val = 0;
	for i = 1 : length( x ) - 1
		val = val + ( 1 - x( i ) ) ^ 2 + 100 * ( x( i + 1 ) - x( i ) ^ 2 ) ^2;
	end
end
