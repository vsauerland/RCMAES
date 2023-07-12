[![DOI](https://zenodo.org/badge/532043505.svg)](https://zenodo.org/badge/latestdoi/532043505)


The repository provides C++ and Matlab implementations of R-CMA-ES,
a variant of the Covariance Matrix Adaption Evolution Strategy.
R-CMA-ES allows to declare random components.
The algorithm is described in [1].

# PURPOSE

R-CMA-ES is supposed to optimize the parameters of models, e.g. BGC models.
It supports the declaration of a random parameter, e.g. a parameter carrying
a high amount of uncertainty because it represents multiple (unresolved)
real-world processes at once.
We provide both a Matlab implementation, which we used for tests with
benchmark problems and a C++ implementation that is suitable to be used to
calibrate the parameters of computaional demanding models (BGC ocean models)
on HPC platforms. The latter is most suitable, if the model running time is
approximately the same for different parameter vectors.

# Matlab Implementation (Folder MATLAB)

The Matlab package consists of the files
- *cmaes.m*: CMA-ES code by N. Hansen (2012), version 3.61.beta from
  http://www.cmap.polytechnique.fr/~nikolaus.hansen/cmaes_inmatlab.html
- *rcmaes.m*: R-CMA-ES code (cf. Algorithm 2 in [1])
- *hybridRanking.m*: Used to select and sort the better half of samples
  (for the applied mirrored sampling, this yields the same as
  simply selecting the better sample from each mirrored pair)
- *testf.m*: Implements a couple of bechnmark functions
- *expectedFitness.m*: Approximates function expectation w.r.t. random component
- *rcmaesVScmaes.m*: Compares R-CMA-ES and CMA-ES w.r.t. a single test instance
- *testRcmaes.m*: Testbed comparing R-CMA-ES and CMA-ES for 10 test instances
  using 5 test functions

Currently, in order to invoke "cmaes.m", the five test functions from "testf.m"
that are used in the testbed are implemented solely as
- *linear1.m*
- *sphere.m*
- *rosenbrock.m*
- *griewank.m*
- *rastrigin.m*

# C++ Implementation (Folder CPP)

## Approach

Optimization is done by chain jobs alternating RCMAES optimization steps with
multiple parallel model runs. The number lambda of model runs corresponds to the
number of sampled parameter vectors in each CMAES iteration (population size).
In each iteration, the serial RCMAES job is started using "serial.job." It
samples lambda parameter vectors from a probability distribution and writes
them to files called "parameters_i_j.txt" where i is the current iteration
number and j stands for the session number in {1,...,lambda}.
All parallel model evaluations are started using "parallel.job".
NOTE: Each model run is supposed to read a parameter file and to write the
results to corresponding files called "fitness_i_j.txt".
After termination of the model runs, a new RCMAES optimization step is started
reading fitness values from the files, updating probability distributions and
sampling new parameter vectors.
After the i-th iteration, RCMAES stores necessary data for the next iteration
in a file "algVars_i.txt" which is read in iteration i+1.
Operational settings and an iteration history is recorded in the file
"nIter.txt", which is accessed by both job files.

## USAGE

### Generating Executables

The package consists of the files
- *eigen-3.4.0.tar.gz*: A copy of the "Eigen" algebra package
- *rcmaes.cpp*: The main source of the CMAES algorithm
- *auxiliaries.cpp/hpp*: Some auxiliary functions required by CMAES
- *eigenreq.hpp*: Path to required EIGEN algebra package and some defs
- *testfunctions.cpp*: Collection of testfunctions for optimization
- *Makefile*: Generates executables
- *nIter0.dat*: Sample interface file (with operational settings)
- *serial.job*: Serial job file (starting RCMAES optimizer)
- *parallel.job*: Parallel job file (starting parallel model runs)
- *unchained.sh*: Shell script doing both, optimizer iterations
  and objective function evaluations

In order to generate executables, the "Eigen" algebra package must be available
on the system. It can be downloaded from
  http://eigen.tuxfamily.org
as ".tar.gz" file and simply extracted to some place on your system without any
installation process (cf. http://eigen.tuxfamily.org/dox/GettingStarted.html).
The Makefile must be adapted to contain the right path, i.e., you have to
change the Makefile line spelling
  EIGEN = -I eigen-3.4.0/Eigen/
to contain your own path to the Eigen folder.
Version 3.4.0 within the folder. You may type "tar -xvf eigen-3.4.0.tar.gz".
The RCMAES executables can be generated by typing "make".

### Operational settings

Operational parameters must be set in the file "nIter.txt" which will be
accessed by the jobfiles/shell script and also serves as iteration history.
The file "nIter0.dat" is a template that can be copied to "nIter.txt".

Currently, also the file "parallel.job" which runs multiple model evaluations
in parallel has to be modified. Namely, the walltime, the number of nodes and
the number of processors per node have to be adapted to fit your requirements
(according to the requirements for a single model simulation, multiplied with
the number of sessions set in "nIter.txt").
The settings in the template "nIter0.dat" are for the "Himmelblau" testfunction
of dimension 2 using 50 RCMAES iterations with 10 samples per iteration.

### Running Optimization

On, e.g.,  HLRN one optimization is started typing "sbatch serial.job".
It is possible to continue a terminated optimization by increasing the
iteration number in nIter.txt and doing "sbatch serial.job" again.
Also in the case of an iterruption caused by technical problems, optimization
can be continued at the point it has been interrupted. In that case, the last
iteration in "nIter.txt" must be deleted if the fitness values of that iteration
have not all been written to the corresponding files!

# Additional Material

The tarball "MOPS.tar" contains the Transport Matrix Method (TMM) code of the
ocean biogeochemical model MOPS as applied in [1] and originally applied in [2].
Recent TMM code can be found under "https://github.com/samarkhatiwala/tmm".
You also need to download additional files from the GIT repositories branch
that is archived under "https://zenodo.org/record/1246300",
and generate the required transport matrices, and forcing and geometry data
according to steps 2), 3), and 5) of their "README.txt".
You also need to install PETSc according to step 1) of their "README.txt"
in order to compile the MOPS code.
After unpacking the folder "MOPS" can be compiled by typing "make rmops".
Our files "nIter.txt", "serial.job", and "parallel.job" in our folder
"CPP/optimizingMOPS" serve as template for the calibration of MOPS
using RCMAES and the "slurm" workload manager, as done by [1] and [2].

# LITERATURE

[1] V. Sauerland, C. von Hallern, I. Kriest, and J. Getzlaff (2023).
    A CMA-ES Algorithm Allowing for Random Parameters in Model Calibration.
    Submitted to AGU Journal of Advances in Modelling Earth Systems (JAMES).
[2] I. Kriest, V. Sauerland, S. Khatiwala, A. Srivastav, and A. Oschlies (2017).
    Calibrating a global three-dimensional biogeochemical ocean model (MOPS-1.0).
    Geoscientific Model Development 10(1):127-154. DOI 10.5194/gmd-10-127-2017. 
[3] N. Hansen (2016). The CMA Evolution Strategy: A Tutorial. arXiv.
[4] CMAES website: https://www.lri.fr/~hansen/cmaesintro.html
