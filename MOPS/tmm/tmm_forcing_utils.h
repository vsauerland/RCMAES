#define MAX_MATRIX_PERIODS 20
#define MAX_FORCING_PERIODS 2000
typedef struct {
  Vec *up;
  PetscBool firstTime;
  PetscInt numPeriods;
} PeriodicVec;

typedef struct {
  Mat Ap[MAX_MATRIX_PERIODS];
  PetscBool firstTime;
  PetscInt numPeriods;
} PeriodicMat;

typedef struct {
  PetscScalar *up[MAX_FORCING_PERIODS];
  PetscInt arrayLength;
  PetscBool firstTime;
  PetscInt numPeriods;
} PeriodicArray;

extern PetscErrorCode calcInterpFactor(PetscInt nmax,PetscScalar t,PetscScalar tarray[],PetscInt *itf,PetscScalar *alpha);
extern PetscErrorCode calcPeriodicInterpFactor(PetscInt n,PetscScalar t,PetscScalar tparr[],PetscInt *itf1,PetscInt *itf2,PetscScalar *al1,PetscScalar *al2);
extern PetscInt findindex(PetscScalar tarr[],PetscInt nmax,PetscScalar t);
extern PetscErrorCode interpPeriodicMatrix(PetscScalar tc, Mat *A, PetscScalar cyclePeriod, PetscInt numPeriods, 
                                           PetscScalar *tdp, PeriodicMat *user, char *filename);
extern PetscErrorCode interpPeriodicVector(PetscScalar tc, Vec *u, PetscScalar cyclePeriod, PetscInt numPeriods, 
                                           PetscScalar *tdp, PeriodicVec *user, char *filename);

extern PetscErrorCode destroyPeriodicVec(PeriodicVec *user);
extern PetscErrorCode destroyPeriodicMat(PeriodicMat *user);
extern PetscErrorCode destroyPeriodicArray(PeriodicArray *user);

extern PetscErrorCode interpTimeDependentVector(PetscScalar tc, Vec *u, PetscInt numTracers, PetscInt nt, PetscScalar *t, Vec **ut);

extern PetscErrorCode writeBinaryScalarData(char *fileName, PetscScalar *arr, PetscInt N, PetscBool appendToFile);
