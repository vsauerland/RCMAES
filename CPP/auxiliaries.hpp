#include "eigenreq.hpp"

using namespace Eigen;

double normalCDF( double value );

double normaldistribution( double m, double s );

void merge( double *v, double *x, int a, int m, int b );

void mergeSort( double *v, double *x, int a, int b );

void sortVectors( VectorXld &v, VectorXld &x );

double median( VectorXld v );

double percentile( VectorXld v, double p );

void xintobounds( VectorXld x, VectorXld lb, VectorXld ub, VectorXld &bx, VectorXld &ix );

void xIntoUnitCube( VectorXld x, VectorXld &bx, VectorXld &ix );

VectorXld scale( VectorXld x, VectorXld lb, VectorXld ub );

VectorXld calcWeights( MatrixXld Z, VectorXld lb, VectorXld ub, VectorXi isStochastic );

VectorXld hybridSort( VectorXld arfitness, MatrixXld arx, VectorXld lb, VectorXld ub, VectorXi isStochastic, int mu );

void correctColumn( MatrixXld &B, MatrixXld &D, int iStoch );
