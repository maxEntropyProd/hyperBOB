# hyperBOB
Optimization using MPI parallel Latin hypercube sampling and BOBYQA

hyperBOB is a Fortran subroutine that uses MPI to run embarrassingly parallel instances of the derivative-free box-bounded optimizer BOBYQA. Each instance of BOBYQA will be called using a different initial condition, where the initial conditions are obtained from Latin hypersquares. All available processes in MPI_COMM_WORLD will be used, so the more MPI processes available, the greater the coverage.  There is no guarantee it will locate a global optimum, which is highly dependent on the number of processes used, the dimension of the function being sampled, and the number of local optima in the hypersurface.  It will return the best solution found, but all solutions are available in xMat_hyperBOB and fVal_hyperBOB.

hyperBOB_initialize must be called first to allocate xMat_hyperBOB and fVal_hyperBOB.

Once xMat_hyperBOB and fVal_hyperBOB are no longer needed, they can be deallocated by calling hyperBOB_cleanUp.

A simple test program is included in EggholderFcnTest.f90

<a href="https://zenodo.org/badge/latestdoi/260558296"><img src="https://zenodo.org/badge/260558296.svg" alt="DOI"></a>

