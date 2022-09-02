function ix = hybridRank( arfitness, arx, lb, ub, isStochastic, mu )
%
% function ix = hybridRank( fitnesses, vectors, isStochastic )
% 
% ranking of vectors w.r.t.
% 1) a certain component
% 2) objective function values
%
% (i)   calculate ranks acc. to 1)
% (ii)  partition consecutively ranked pairs w.r.t. 2) into good A and bad B
% (iii) rank A and B by 2) again and increase B ranks by |A|
%
% ix (output)   index vector; column ix( i ) of arx has rank i
% arfitness:    objective function values of the vectors
% arx:          matrix of the vectors as columns
% isStochastic: indicator vector for the ONE ranking component
% mu:           if at leat mu columns of arx are within bounds,
%               also the mu best ranked colums of arx should be within bounds

n = length( arfitness );
m = length( isStochastic );
s = size( arx );
if s( 1 ) ~= m || s( 2 ) ~= n
  error( "size( arx ) should be [ length( arfitness ) length( isStochastic )" );
end
sIndex = find( isStochastic == 1 );
if length( sIndex ) > 1
  error( "isStochastic must only indicate a single stochastic index" );
end
%arx
v = arx( sIndex, : );
[ v1, ix1 ] = sort( v );
%v1
%ix1
A = [];
B = [];
for i = 1 : n / 2
  i1 = ix1( 2 * i - 1 );
  i2 = ix1( 2 * i );
  if arfitness( i1 ) < arfitness( i2 )
    A = [ A i1 ];
    B = [ B i2 ];
  else
    A = [ A i2 ];
    B = [ B i1 ];
  end
end
if mod( n, 2 ) == 1
  B = [ B, ix1( n ) ];
end
[ va, ixa ] = sort( arfitness( A ) );
[ vb, ixb ] = sort( arfitness( B ) );
A = A( ixa );
B = B( ixb );
% column vector ix( i ) has rank i
ix = [ A B ];
