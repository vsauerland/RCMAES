function val = sphere( x, y )
if nargin == 2 && length( x ) == 1 && length( y ) == 1
	val = x ^ 2 + y ^ 2;
else
	val = sum( x .^ 2 );
end
