// Implementation of the Covariance Matrix Adaption Evolution Strategy (CMAES)
// This implements the basic strategy with box constraints using C++ and the
// Eigen package for matrix calculations
// We follow [1] for the pure CMAES extended by a simple "boundary handling"

#include <iostream>
#include <fstream>
#include "eigenreq.hpp"
#include "auxiliaries.hpp"

using namespace std;

void writeMuDB( int nI, int N, VectorXld mu, long double sigma, MatrixXld D, MatrixXld B )
{
	stringstream fileName1;
	stringstream fileName2;
	stringstream fileName3;
	fileName1 << "mu" << ( nI + 1 ) << ".txt";
	fileName2 << "ew" << ( nI + 1 ) << ".txt";
	fileName3 << "ev" << ( nI + 1 ) << ".txt";
	ofstream f1( fileName1.str().c_str() );
	ofstream f2( fileName2.str().c_str() );
	ofstream f3( fileName3.str().c_str() );
	f1.precision( 16 );
	f2.precision( 16 );
	f3.precision( 16 );
	for ( int i = 0; i < N; i++ )
	{
		f1 << mu( i ) << endl;
		f2 << sigma * D( i, i ) * sigma * D( i, i ) << endl;
		for ( int j = 0; j < N; j++ ) f3 << " " << B( i, j );
		f3 << endl;
	}
	f1.close();
	f2.close();
	f3.close();
}

int writeScaledSamples( int nI, int lambda, int N, VectorXld lb, VectorXld ub, MatrixXld arx )
{
	for ( int j = 0; j < lambda; j++ )
	{
		stringstream fileName;
		fileName << "parameters_" << ( nI + 1 ) << "_" << ( j + 1 ) << ".txt";
		ofstream f( fileName.str().c_str() );
		f.precision( 43 );
		for ( int i = 0; i < N; i++ )
		{
			f << lb( i ) + ( ub( i ) - lb( i ) ) * arx( i, j ) << "\n";
		}
		f.close();
	}
	return( 0 );
}

int readResults( int nI, int lambda, VectorXld penalty, VectorXld &arfitness )
{
	for ( int i = 0; i < lambda; i++ )
	{
		string line;
		stringstream ss;
		stringstream fileName;
		fileName << "fitness_" << nI << "_" << ( i + 1 ) << ".txt";
		printf( "read results from file %s\n", fileName.str().c_str() );
		ifstream f( fileName.str().c_str() );
		getline( f, line );
		if ( line.find( "." ) == string::npos ) arfitness( i ) = 1e12;
		else
		{
			ss << line;
			assert( ss >> arfitness( i ) );
		}
		f.close();
		arfitness( i ) = arfitness( i ) + penalty( i );
	}
	return( 0 );
}

int writeAlgVars( int nI, int lambda, int N, int counteval, long double sigma, VectorXld pc, VectorXld ps, MatrixXld C, MatrixXld B, MatrixXld D, MatrixXld arz, MatrixXld ary, MatrixXld arx, VectorXld penalty, VectorXld xmean )
{
	stringstream fileName;
	fileName << "algVars_" << ( nI + 1 ) << ".txt";
	ofstream f( fileName.str().c_str() );
	f.precision( 43 );
	// 1.) write counteval
	f << "counteval" << nI << "\n"; 
	f << counteval << "\n";
	// 2.) write sigma
	f << "sigma" << nI << "\n"; 
	f << sigma << "\n";
	// 3.) write pc
	f << "pc" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << pc( i ) << "\n";
	}
	// 4.) write ps
	f << "ps" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << ps( i ) << "\n";
	}
	// 5.) write C
	f << "C" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << C( i, j ) << " ";
		}
		f << "\n";
	}
	// 6.) write B
	f << "B" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << B( i, j ) << " ";
		}
		f << "\n";
	}
	// 7.) write D
	f << "D" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < N; j++ )
		{
			f << D( i, j ) << " ";
		}
		f << "\n";
	}
	// 8.) write arz
	f << "arz" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << arz( i, j ) << " ";
		}
		f << "\n";
	}
	// 9.) write ary
	f << "ary" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << ary( i, j ) << " ";
		}
		f << "\n";
	}
	// 10.) write arx
	f << "arx" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		for ( int j = 0; j < lambda; j++ )
		{
			f << arx( i, j ) << " ";
		}
		f << "\n";
	}
	// 11.) write penalty
	f << "penatly" << nI << "\n"; 
	for ( int i = 0; i < lambda; i++ )
	{
		f << penalty( i ) << "\n";
	}
	// 12.) write xmean
	f << "xmean" << nI << "\n"; 
	for ( int i = 0; i < N; i++ )
	{
		f << xmean( i ) << "\n";
	}
	f.close();
	return( 0 );
}

int readAlgVars( int nI, int lambda, int N, int *counteval, long double *sigma, VectorXld &pc, VectorXld &ps, MatrixXld &C, MatrixXld &B, MatrixXld &D, MatrixXld &arz, MatrixXld &ary, MatrixXld &arx, VectorXld &penalty )
{
	string line;
	stringstream ss;
	stringstream fileName;
	fileName << "algVars_" << nI << ".txt";
	ifstream f( fileName.str().c_str() );
	// 1.) read counteval
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *counteval );
	// 2.) read sigma
	getline( f, line );
	getline( f, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> *sigma );
	// 3.) read pc
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> pc( i ) );
	}
	// 4.) read ps
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> ps( i ) );
	}
	// 5.) read C
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line.c_str();
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> C( i, j ) );
		}
	}
	// 6.) read B
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> B( i, j ) );
		}
	}
	// 7.) read D
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < N; j++ )
		{
			assert( ss >> D( i, j ) );
		}
	}
	// 8.) read arz
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> arz( i, j ) );
		}
	}
	// 9.) read ary
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> ary( i, j ) );
		}
	}
	// 10.) read arx
	getline( f, line );
	for ( int i = 0; i < N; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		for ( int j = 0; j < lambda; j++ )
		{
			assert( ss >> arx( i, j ) );
		}
	}
	// 11.) read penalty 
	getline( f, line );
	for ( int i = 0; i < lambda; i++ )
	{
		getline( f, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss >> penalty( i ) );
	}
	f.close(); 
	return( 0 );
}



int main( int argc, char* argv[] )
{
	int N, N1, lambda, nI; // problem size, population size and iteration number
	int seed; // random number seed
	string functionName; // default fitness function
	if ( argc < 3 ) 
	{
		printf( "program requires 2 command line arguments:\ninstance filename (string) and current iteration (int)\nthe associated file is supposed to contain:\n - problem size,\n - fitness function name,\n - lower and upper bounds\n\n" );
	}
	assert( argc == 3 );
//	assert( argc == 2 );
	nI = atoi( argv[ 1 ] );
	string inFileName = argv[ 2 ];
	string line;
	stringstream ss;
	ifstream inFile( inFileName.c_str() );
//	inFile >> N;
//	inFile >> lambda;
//	inFile >> functionName;
	getline( inFile, line );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> functionName );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> N );
	N1 = N;
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> lambda );
	if ( lambda < 2 )
	{
		printf( "the population size must be at least 2\n" );
	}
	assert( lambda >= 2 );
	getline( inFile, line );
	getline( inFile, line );
	ss.str( "" ); ss.clear(); ss << line;
	assert( ss >> seed );
	getline( inFile, line );
	VectorXld lb( N );
	VectorXld ub( N );
	for ( int i = 0; i < N; i++ )
	{
		getline( inFile, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss  >> lb( i ) >> ub( i ) );
	}
	// VSA (2.3.2018) introduced indicator for stochastic components as parameter
	getline( inFile, line );
	VectorXi isStochastic = VectorXi::Zero( N );
	for ( int i = 0; i < N; i++ )
	{
		getline( inFile, line );
		ss.str( "" ); ss.clear(); ss << line;
		assert( ss  >> isStochastic( i ) );
	}
	inFile.close();

	// derived operational parameters: Selection
//	int lambda = 4 + floor( 3 * log( N ) ); // population size, offspring number (vsa: original, now set in nIter.txt)
	// find stochastic index iStoch, define accordant unity vector eStoch
	int iStoch = 0;
	for ( int i = 0; i < N; i++ ) if ( isStochastic( i ) ) iStoch = i;
	VectorXld eStoch = VectorXld::Zero( N );
	eStoch( iStoch ) = 1;

	int mu = floor( ( long double )lambda / 2 ); // number of parents/points for recombination
	VectorXld weights = VectorXld::LinSpaced( mu, 1, mu );
	weights = log( mu + 0.5 ) - weights.array().log();// muXone recombination weights
	weights = weights / weights.sum(); // normalize recombination weights array
	long double mueff = weights.sum() * weights.sum() / ( weights.array() * weights.array() ).sum(); // variance-effective size of mu
	long double sigma; // coordinate wise standard deviation (step-size)
	VectorXld sVec( N );
	long double chiN = sqrt( N ) * ( 1 - ( long double )1 / ( 4 * N ) + ( long double )1 / ( 21 * N * N ) ); // expectation of || N( 0, I ) || == norm( randn( N, 1 ) )

	// Strategy parameter setting: Adaptation
	long double cs = ( mueff + 2 ) / ( N1 + mueff + 5 ); // t-const for cumulation for sigma control
	long double damps = 1 + cs;
	long double cc = ( 4 + mueff / N1 ) / ( N1 + 4 + 2 * mueff / N1 ); // time constant for cumulation for C
	long double c1 = 2 / ( ( N1 + 1.3 ) * ( N1 + 1.3 ) + mueff ); // learning rate for rank-one update of C
	long double cmu = 2 * ( mueff - 2 + 1 / mueff ) / ( ( N1 + 2 ) * ( N1 + 2 ) + mueff ); // learning rate for rank-mu update
	cmu = min( 1 - c1, cmu );

	// declare dynamic (internal) strategy parameters
	VectorXld pc( N ); // evolution path for C
	VectorXld ps( N ); // evolution path for sigma
	MatrixXld B( N, N );  // B defines coordinate system
	MatrixXld D( N, N );  // diagonal matrix D defines scaling
	MatrixXld C( N, N );  // covariance matrix
	int counteval; // count total number of function evaluations
	MatrixXld arz( N, lambda );
	MatrixXld ary( N, lambda );
	MatrixXld arx( N, lambda );
	MatrixXld arz_sub( N, mu );
	MatrixXld ary_sub( N, mu );
	MatrixXld arx_sub( N, mu );
	VectorXld arfitness( lambda );
	VectorXld penalty( lambda ); // penalties for samples
	VectorXld arindex( lambda );
	VectorXld xmean( N );
	VectorXld zmean( N );

	if ( nI == 0 )
	{
		xmean = VectorXld::LinSpaced( N, 0.5, 0.5 );
		sigma = 0.3;
		sVec = VectorXld::LinSpaced( N, sigma, sigma );
		sVec( iStoch ) = 0.5;
		pc = VectorXld::Zero( N );
		ps = VectorXld::Zero( N );
		B = MatrixXld::Identity( N, N );
		D = MatrixXld::Identity( N, N );
		C = MatrixXld::Identity( N, N );
		cout << "xmean0" << endl << xmean << endl << endl;
		cout << "sigma0" << endl << sigma << endl << endl;
		cout << "C0" << endl << C << endl << endl;
		cout << "D0" << endl << D << endl << endl;
		cout << "B0" << endl << B << endl << endl;
		
		// vsa (17.10.2021): need to initialize penatly, too
		// and doing it also for other "eigen" type variables
		penalty = VectorXld::Zero( lambda );
		arz = MatrixXld::Zero( N, lambda );
		ary = MatrixXld::Zero( N, lambda );
		arx = MatrixXld::Zero( N, lambda );
		arz_sub = MatrixXld::Zero( N, mu );
		ary_sub = MatrixXld::Zero( N, mu );
		arx_sub = MatrixXld::Zero( N, mu );

		counteval = 0;
		printf( "initialized\n" );
	}
	else // ( nI > 0 )
	{
		// read current dynamic strategy parameters and fitness values of
		// previous iteration samples
		readAlgVars( nI, lambda, N, &counteval, &sigma, pc, ps, C, B, D, arz, ary, arx, penalty );
		readResults( nI, lambda, penalty, arfitness );

		// Sort and compute weighted mean into xmean
		arindex = hybridSort( arfitness, arx, lb, ub, isStochastic, mu );
		cout << "arfitness" << nI << endl << arfitness << endl << endl;
		cout << "arindex" << nI << endl << arindex << endl << endl;
		for ( int k = 0; k < mu; k++ )
		{
			arx_sub.col( k ) = arx.col( arindex( k ) );
			ary_sub.col( k ) = ary.col( arindex( k ) );
			arz_sub.col( k ) = arz.col( arindex( k ) );
		}
		xmean = arx_sub * weights; // recombination [1] (38),(39)
		zmean = arz_sub * weights; // = D ^ -1 * B * ( xmean - "xmean_old" ) / sigma
		// VSA (2.3.2018) we keep stochastic components of the mean
		// in the center of their box
		for ( int i = 0; i < N; i++ ) if ( isStochastic( i ) == 1 )
		{
			xmean( i ) = 0.5; // we work in the unit cube!  * ( lb( i ) + ub( i ) ); 
			zmean( i ) = 0;
		}
		cout << "scaled xmean" << nI << endl << scale( xmean, lb, ub ) << endl << endl;
		cout << "zmean" << nI << endl << zmean << endl << endl;
		cout << "cs =" << cs << endl;
		cout << "cc =" << cc << endl << endl;
		
		// Cumulation: Update evolution paths
		ps = ( 1 - cs ) * ps + ( sqrt( cs * ( 2 - cs ) * mueff ) ) * ( B * zmean ); // [1] (40)
		pc = ( 1 - cc ) * pc + sqrt( cc * ( 2 - cc ) * mueff ) * ( B * D * zmean ); // [1] (42)
		cout << "ps" << nI << endl << ps << endl << endl;
		cout << "pc" << nI << endl << pc << endl << endl;

		// Adapt covariance matrix C and step size sigma
		sigma = sigma * exp( ( cs / damps ) * ( ps.norm() / chiN - 1 ) ); // [1] (41)
		cout << "sigma" << nI << endl << sigma << endl << endl;
		MatrixXld C1 = pc * pc.transpose();
		MatrixXld Cmu = ary_sub * weights.asDiagonal() * ary_sub.transpose();
		C =	  ( 1 - c1 - cmu ) * C // [1] (43)
			+ c1 * C1 // plus rank one update
			+ cmu * Cmu; // plus rank mu update
		// Enforce symmetry
		for ( int i = 0; i < N - 1; i++ ) for ( int j = i + 1; j < N; j++ ) C( j, i ) = C( i, j );
		// Manipulate C to yield std.-dev. ellipse with a "stochastic axis"
		// that is parallel to the stochastic axis of the coordinate system
		C.col( iStoch ) = eStoch;
		C.row( iStoch ) = eStoch.transpose();
		// Update B and D from C
		SelfAdjointEigenSolver<MatrixXld> es( C );
		D = es.eigenvalues().asDiagonal();
		B = es.eigenvectors(); // B==normalized eigenvectors
		D = D.array().sqrt(); // D contains standard deviations now
		// eigen decomposition sorts B and D w.r.t. eigenvalues
		// => rearrange system such that column iStoch of B is eStoch
		correctColumn( B, D, iStoch );
		cout << "C" << nI << endl << C << endl << endl;
		cout << "B" << nI << endl << B << endl << endl;
		cout << "diag( D )'" << nI << endl << D.diagonal().transpose() << endl << endl;
		cout << counteval << ": " << arfitness( 0 ) << endl;
	} // else ( nI > 0 )

	srand( seed + counteval );
	printf( "calling srand( %i ) in iteration %i\n", seed + counteval, nI );
	sVec = VectorXld::LinSpaced( N, sigma, sigma );
	sVec( iStoch ) = 0.5;
	// Generate and evaluate lambda offspring by mirrored sampling
	for ( int k = 0; k < mu; k++ )
	{
		for ( int i = 0; i < N; i++ )
		{
			arz( i, k ) = normaldistribution( 0.0, 1.0 );
			arz( i, mu + k ) = -arz( i, k );
		}
		// keep the stochastic component within its bounds
		while ( arz( iStoch, k ) < -1 || arz( iStoch, k ) > 1 ) arz( iStoch, k ) = normaldistribution( 0.0, 1.0 );
		arz( iStoch, mu + k ) = arz( iStoch, k );
	}
	for ( int k = 0; k < lambda; k++ )
	{
		ary.col( k ) = B * D * arz.col( k );
		arx.col( k ) = xmean.array() + sVec.array() * ( B * D * arz.col( k ) ).array();
		VectorXld dummy( N );
		VectorXld v( N );
		xIntoUnitCube( arx.col( k ), v, dummy );
		// penalize non-stochastic out of bound components:
		penalty( k ) = 0;
		for ( int i = 0; i < N; i++ ) if ( i != iStoch )
		{
			long double q = v( i ) - arx.col( k )( i );
			q = q * q;
			penalty( k ) = penalty( k ) + q;
		}
		penalty( k ) = 1e6 * penalty( k );
		counteval = counteval + 1;
	}

	writeScaledSamples( nI, lambda, N, lb, ub, arx );
	writeAlgVars( nI, lambda, N, counteval, sigma, pc, ps, C, B, D, arz, ary, arx, penalty, xmean );
	writeMuDB( nI, N, xmean, sigma, D, B ); // data for illustrations with standard deviation ellipses
	// update nIter.txt
	ofstream f;
	f.open( "nIter.txt", ofstream::out | ofstream::app );
	f << nI + 1 << " " << arfitness( 0 ) << " 0" << "\n";
	f.close();
	return( 0 );
}

/* -------------------------------- Literature --------------------------------

[1] N. Hansen (2011). The CMA Evolution Strategy: A Tutorial

-----------------------------------------------------------------------------*/
