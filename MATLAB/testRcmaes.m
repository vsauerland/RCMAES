function testRCMAES

% testRCMAES: a set of test instances to compare R-CMA-ES and CMA-ES
% (we apply rcmaesVScmaes to each instance)

fnames={ 'linear1', 'linear1', 'sphere', 'sphere', 'rosenbrock', 'rosenbrock', 'griewank', 'griewank', 'rastrigin', 'rastrigin' };
lbs={ zeros( 1, 6 ), zeros( 1, 10 ), -ones( 1, 6 ), -ones( 1, 10 ), [ -1 -1 0.1 -1 -1 -1 ], [ 0.5 -ones( 1, 9 ) ], -10 * ones( 1, 6 ), -10 * ones( 1, 10 ), -5.12 * ones( 1, 6 ), -5.12 * ones( 1, 10 ) };
ubs={ ones( 1, 6 ), ones( 1, 10 ), ones( 1, 6 ), ones( 1, 10 ), [ 3 3 0.3 3 3 3 ], [ 0.8 3 * ones( 1, 9 ) ], 10 * ones( 1, 6 ), 10 * ones( 1, 10 ), 5.12 * ones( 1, 6 ), 5.12 * ones( 1, 10 ) };
iStochs={ 3 1 3 1 3 1 3 1 3 1};
lambdas={ 10, 12, 10, 12, 10, 12, 10, 12, 10, 12 };
maxIters={ 200, 300, 200, 300, 200, 300, 200, 300, 200, 300 };

numInstances = length( fnames );
numTrials = 10;
k1 = zeros( numInstances, numTrials );
k2 = zeros( numInstances, numTrials );
k3 = zeros( numInstances, numTrials );
f1 = zeros( numInstances, numTrials );
f2 = zeros( numInstances, numTrials );
f3 = zeros( numInstances, numTrials );
expF1 = zeros( numInstances, numTrials );
expF2 = zeros( numInstances, numTrials );
expF3 = zeros( numInstances, numTrials );
meanK1 = zeros( numInstances );
bestK1 = zeros( numInstances );
meanK2 = zeros( numInstances );
bestK2 = zeros( numInstances );
meanK3 = zeros( numInstances );
bestK3 = zeros( numInstances );
meanF1 = zeros( numInstances );
bestF1 = zeros( numInstances );
meanF2 = zeros( numInstances );
bestF2 = zeros( numInstances );
meanF3 = zeros( numInstances );
bestF3 = zeros( numInstances );
meanExpF1 = zeros( numInstances );
bestExpF1 = zeros( numInstances );
meanExpF2 = zeros( numInstances );
bestExpF2 = zeros( numInstances );
meanExpF3 = zeros( numInstances );
bestExpF3 = zeros( numInstances );
%rng( 1, 'twister' );
rng( 'default' );
for i = 1 : numInstances
	fname = fnames{ i }
	lb = lbs{ i }
	ub = ubs{ i }
	iStoch = iStochs{ i };
	lambda = lambdas{ i };
	maxIter = maxIters{ i };
	for k = 1 : numTrials
		[ k1( i, k ), k2( i, k ), x1, x2, f1( i, k ), f2( i, k ), expF1( i, k ), expF2( i, k ) ] = rcmaesVScmaes( fname, lb, ub, iStoch, lambda, maxIter );
	end
	meanK1( i ) = sum( k1( i, : ) ) / numTrials;
	bestK1( i ) = min( k1( i, : ) );
	meanK2( i ) = sum( k2( i, : )  ) / numTrials;
	bestK2( i ) = min( k2( i, : )  );
	meanF1( i ) = sum( f1( i, : ) ) / numTrials;
	bestF1( i ) = min( f1( i, : ) );
	meanF2( i ) = sum( f2( i, : )  ) / numTrials;
	bestF2( i ) = min( f2( i, : )  );
	meanExpF1( i ) = sum( expF1( i, : )  ) / numTrials;
	bestExpF1( i ) = min( expF1( i, : )  );
	meanExpF2( i ) = sum( expF2( i, : )  ) / numTrials;
	bestExpF2( i ) = min( expF2( i, : )  );
	display( "Results: best and mean values out of trials" );
	display( "columns 1,2 for iterations, columns 3,4 for minimum value, columns 5,6 for expected value" );
	display( "first row for CMAES, second row for SCMAES2, third row for SCMAES3" );
	[ bestK1( i ) meanK1( i ) bestF1( i ) meanF1( i ) bestExpF1( i ) meanExpF1( i ) ]
	[ bestK2( i ) meanK2( i ) bestF2( i ) meanF2( i ) bestExpF2( i ) meanExpF2( i ) ]
end
