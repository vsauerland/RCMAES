function val = linear1( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = -x - y;
else
	val = -sum( x );
end
