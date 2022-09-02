#include "auxiliaries.hpp"
#include <iostream>
#include "stdio.h"
#include <queue>

double normaldistribution( double m, double s )
{
	double r1 = ( double )rand() / RAND_MAX;
	double r2 = ( double )rand() / RAND_MAX;
	double pi = 3.141592653589793;
	return ( s * sqrt( -2 * log( r1 ) ) * sin( 2 * pi * r2 ) + m );
}

double normalCDF( double value )
// cumulative density function of normal distribution
// VS: added 02.02.2022
{
   return 0.5 * erfc( -value * M_SQRT1_2 );
}

void merge( double *v, double *x, int a, int m, int b )
// merge two consecutive segments of a vector 
// and arrange a second vector correspondingly
// 
// v, x     vectors
// a, m, b  component indizes ( counting from 0 )
//
// both subvectors of consecutive components from a up to m 
// and from m+1 up to b of v are supposed to be in ascending order 
// merging brings all components from a up to b in ascending order
// and arranges the components of x correspondingly  
{
	double	*w, *y;
	int	i, j, k;
	w = ( double* ) calloc( m + 1 - a, sizeof( double ) );
	y = ( double* ) calloc( m + 1 - a, sizeof( double ) );
	// copy left segments of v and x to auxially vectors w and y:
	for( i = 0; i < m + 1 - a; i++ )
	{
		*( w + i ) = *( v + a + i );
		*( y + i ) = *( x + a + i );
	}
	i = 0;
	j = m + 1;
	k = a;
	// copy back the next greater elements:
	while( k < j && j <= b )
	{
		if( *( w + i ) <= *( v + j ) )
		{
			*( v + k ) = *( w + i );
			*( x + k ) = *( y + i );
			k++;
			i++;
		} else
		{
			*( v + k ) = *( v + j );
			*( x + k ) = *( x + j );
			k++;
			j++;	
		} 
	}
	// copy back remaining elements of w:
	while ( k < j )
	{
		*( v + k ) = *( w + i );
		*( x + k ) = *( y + i );
		k++;
		i++;
	}
	free( w );
	free( y );
}

void mergeSort( double *v, double *x, int a, int b )
// sort one double array and permute a second double array accordingly
//
// v: array to be sorted in ascending order
// x: array to be permuted like v
// a,b: components a to b (0 <= a <= b <= length(v) - 1 ) are sorted
{
	int m;
	if( a < b )
	{
		m = ( a + b ) / 2;
		mergeSort( v, x, a, m );
		mergeSort( v, x, m + 1, b );
		merge( v, x, a, m, b );
	}	
}

void sortVectors( VectorXld &v, VectorXld &x )
// sort a vector and permute a second vector accordingly
{
	int n = v.size();
	// copy vectors to some double arrays
	double *va, *xa;
	va = ( double* )calloc( n, sizeof( double ) );
	xa = ( double* )calloc( n, sizeof( double ) );
	for ( int i = 0; i < n; i++ )
	{
		va[ i ] = v( i );
		xa[ i ] = x( i );
	}
	// sort double arrays
	mergeSort( va, xa, 0, n - 1 );
	// copy back double arrays to vectors
	for ( int i = 0; i < n; i++ )
	{
		v( i ) = va[ i ];
		x( i ) = xa[ i ];
	}
	free( va );
	free( xa );
}

double median( VectorXld v ) // Computes median of a vector v
{
	int n = v.size();
	VectorXld dummy( n );
	sortVectors( v, dummy );
	if ( n % 2 == 1 ) return( v( n / 2 ) );
	else return( 0.5 * ( v( n / 2 - 1 ) + v( n / 2 ) ) );
}

double percentile( VectorXld v, double p ) // Computes the percentile p from vector v
{
	int n = v.size();
	double r;	
	VectorXld dummy( n );
	sortVectors( v, dummy );
	if ( p <= 100 * ( 0.5 / n ) ) r = v( 0 );
	else if ( p >= 100 * ( ( n - 0.5 ) / n ) ) r = v( n - 1 );
	else
	{
		// find largest index smaller than required percentile
		VectorXld indices = VectorXld::LinSpaced( n, 1, n );
		VectorXld a = 100 * ( indices.array() - 0.5 ) / n;
		int i = ( a.array() >= p ).select( 0, indices ).lpNorm<Eigen::Infinity>();
		// vsa debugging 18.09.2017: introduced "the non interpolation" cases
		// to cope with small lambda
		if ( i == 0 )
		{
			r = v( 0 );
		}
		else if ( i == n - 1 )
		{
			r = v( n - 1 );
		}
		else
		{
			// interpolate linearly
			assert( a( i + 1 ) > a( i ) );
			r = v( i ) + ( v( i + 1 ) - v( i ) ) * ( p - a( i ) ) / ( a( i + 1 ) - a( i ) );
		}
	}
	return( r );
}

void xintobounds( VectorXld x, VectorXld lb, VectorXld ub, VectorXld &bx, VectorXld &ix )
// brings a vector x between lower and upper bound vectors lb and ub
// bx is the result, ix indicates the out of bounds components of x  
{
	ix = VectorXld::Zero( ix.size() );
	VectorXld ix2 = VectorXld::Zero( ix.size() );
	ix = ( x.array() < lb.array() ).select( 1, ix );
	bx = ( x.array() < lb.array() ).select( lb, x );
	ix2 = ( bx.array() > ub.array() ).select( 1, ix2 );
	bx = ( bx.array() > ub.array() ).select( ub, bx );
	ix = ix2 - ix;
}

void xIntoUnitCube( VectorXld x, VectorXld &bx, VectorXld &ix )
// map a vector x to unit cube (maping outside components onto closest bound)
// bx is the result, ix indicates the out of bounds components of x  
{
	ix = VectorXld::Zero( ix.size() );
	VectorXld ix2 = VectorXld::Zero( ix.size() );
	ix = ( x.array() < 0 ).select( 1, ix );
	bx = ( x.array() < 0 ).select( 0, x );
	ix2 = ( bx.array() > 1 ).select( 1, ix2 );
	bx = ( bx.array() > 1 ).select( 1, bx );
	ix = ix2 - ix;
}

VectorXld scale( VectorXld x, VectorXld lb, VectorXld ub )
// scale a vector x from unit cube
// to rectangle described by lower bounds lb and upper bounds ub
{
	return( lb.array() + ( ub.array() - lb.array() ) * x.array() );
}

VectorXld calcWeights( MatrixXld Z, VectorXld lb, VectorXld ub, VectorXi isStochastic )
// distribute weights acc. to the "probabilities" of sample vectors
// stochastic component (sampled from N(0,1))
//
// X             matrix of sampled vectors
// lb, ub        vectors of lower and upper bounds
// isStochastic  indicator vector for the stochastic component (only 1 allowed, here)
// 
// VS: added 03.02.2022
{
	int N = Z.col( 0 ).size();
	int mu = Z.row( 0 ).size();
	int row = 0;
	while ( isStochastic( row ) != 1 && row < N - 1 ) row++;
	long double stochasticRange = ub( row ) - lb( row );
	long double stochasticCenter = ( lb( row ) + ub( row ) ) / 2.0;
	VectorXld x = Z.row( row ); // stochastic components of the sampled vectors
	x = x.array() - stochasticCenter;
	x = ( 2 / stochasticRange ) * x; // scale to [-1,1]
	VectorXld ix = VectorXld::LinSpaced( mu, 0, mu - 1 );
	sortVectors( x, ix ); // sort by value, ix is the applied permutation
	VectorXld weights = VectorXld::Zero( mu );
	VectorXld wa = VectorXld::Zero( mu ); // first weights term
	VectorXld wb = VectorXld::Zero( mu ); // second weights term
	for ( int i = 0; i < mu; i++ )
	{
		if ( i == 0 || x( i - 1 ) < -1 )
		{
			wa( i ) = normalCDF( x( i ) ) - normalCDF( -1 );
			if ( wa( i ) < 0 ) wa( i ) = 0;
		}
		else if ( x( i ) <= 1 )
		{
			wa( i ) = ( normalCDF( x( i ) ) - normalCDF( x( i - 1 ) ) ) / 2;
		}
		else wa( i ) = 0;
		if ( i == mu - 1 || x( i + 1 ) > 1 )
		{
			wb( i ) = normalCDF( 1 ) - normalCDF( x( i ) );
			if ( wb( i ) < 0 ) wb( i ) = 0;
		}
		else if ( x( i ) >= -1 )
		{
			wb( i ) = ( normalCDF( x( i + 1 ) ) - normalCDF( x( i ) ) ) / 2;
		}
		else wb( i ) = 0;
		weights( i ) = wa( i ) + wb( i );
	}
	sortVectors( ix, weights );
	cout << "stochastic index = " << row << endl << endl;
	cout << "x:" << endl << x << endl << endl;
	cout << "wa:" << endl << wa << endl << endl;
	cout << "wb:" << endl << wb << endl << endl;
	weights = weights / weights.sum();
	return( weights );
}

VectorXld hybridSort( VectorXld arfitness, MatrixXld arx, VectorXld lb, VectorXld ub, VectorXi isStochastic, int mu )
// ranking of vectors w.r.t.
// 1) a certain ("stochastic") component
// 2) objective function values and bound violations
//
// (i)   calculate ranks acc. to 1)
// (ii)  partition consecutively ranked pairs w.r.t. 2) into good A and bad B
// (iii) swap "off bound ranks <= mu" with "in-bound ranks > mu
// (iii) rank A and B by 2) again and increase B ranks by |A|
//
// returns index vector ix such that column ix( i ) of arx has rank i
// arfitness:    objective function values of the vectors
// arx:          matrix of the vectors as columns
// lb:           lower bounds for feasible vectors
// ub:           upper bounds for feasilbe vectors
// isStochastic: indicator vector for the UNIQUE component used by 1)
// mu:           off-bound columns of arx only then get ranks <= mu
//               if less then mu columns of arx are within bounds
{
	int n = arfitness.size();
	int m = isStochastic.size();
	assert( arx.col( 0 ).size() == m );
	assert( arx.row( 0 ).size() == n );
	int sIndex = 0;
	while ( isStochastic( sIndex ) != 1 && sIndex < m - 1 ) sIndex++;
	VectorXld v = arx.row( sIndex );
	VectorXld ix1 = VectorXld::LinSpaced( n, 0, n - 1 );
	sortVectors( v, ix1 );
	queue<long double> qa, qb;
	for ( int i = 0; i < n / 2; i++ )
	{
		int i1 = ( int )ix1( 2 * i );
		int i2 = ( int )ix1( 2 * i + 1 );
		if ( arfitness( i1 ) < arfitness( i2 ) )
		{
			qa.push( i1 );
			qb.push( i2 );
		}
		else
		{
			qa.push( i2 );
			qb.push( i1 );
		}
	}
	int na = n / 2;
	int nb = n / 2;
	if ( 2 * ( n / 2 ) != n )
	{
		nb = nb + 1;
		qb.push( ix1( n - 1 ) );
	}
	assert( qa.size() == na );
	assert( qb.size() == nb );
	VectorXld a = VectorXld::Zero( na );
	VectorXld b = VectorXld::Zero( nb );
	VectorXld va = VectorXld::Zero( na );
	VectorXld vb = VectorXld::Zero( nb );
	for ( int i = 0; i < na; i++ )
	{
		a( i ) = qa.front();
		qa.pop();
		va( i ) = arfitness( int( a( i ) ) );
	}
	for ( int i = 0; i < nb; i++ )
	{
		b( i ) = qb.front();
		qb.pop();
		vb( i ) = arfitness( int( b( i ) ) );
	}
	VectorXld ixa = VectorXld::LinSpaced( na, 0, na - 1 );
	VectorXld ixb = VectorXld::LinSpaced( nb, 0, nb - 1 );
	sortVectors( va, ixa );
	sortVectors( vb, ixb );
//	cout << "ixa" << endl << ixa << endl << endl;
//	cout << "ixb" << endl << ixb << endl << endl;
	VectorXld ix = VectorXld::Zero( n );
	for ( int i = 0; i < na; i++ ) ix( i ) = a( int( ixa( i ) ) );
	for ( int i = 0; i < nb; i++ ) ix( na + i ) = b( int( ixb( i ) ) );
	return( ix );
}

void correctColumn( MatrixXld &B, MatrixXld &D, int iStoch )
{
	int N = ( int )B.rows();
	VectorXld eStoch = VectorXld::Zero( N );
	eStoch( iStoch ) = 1;
	// find column which is (nearly) +/- eStoch
	int iSwap = 1;
	double bestErrorNorm = min( ( B.col( 1 ) - eStoch ).norm(), ( B.col( 1 ) + eStoch ).norm() );
	for ( int i = 1; i < N; i++ )
	{
		double errorNorm = min( ( B.col( i ) - eStoch ).norm(), ( B.col( i ) + eStoch ).norm() );
		if( errorNorm < bestErrorNorm )
		{
			bestErrorNorm = errorNorm;
			iSwap = i;
		}
	}
	if ( iSwap != iStoch )
	{
		// swap colums of B and rows and columns of D
		MatrixXld P = MatrixXld::Identity( N, N );
		P( iSwap, iSwap ) = 0;
		P( iStoch, iStoch ) = 0;
		P( iSwap, iStoch ) = 1;
		P( iStoch, iSwap ) = 1;
		B = B * P;
		D = P * D * P;
	}
}
