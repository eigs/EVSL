### Updates
- Subspace iteration with rational filter
- Interfaces for some subroutines related to matvec and ASigmaBSol are modified to allow block method

### To-do list
- The new method ratsi is only tested on a linux machine with Intel MKL and Pardiso, tests need to be done in other environments and other direct solvers.
- For tests, RSI/LapRSI, RSI/MMRSI, GEN/MMRSI are all running correctly. Other tests may have segmentation fault due to interface changes, the test driver routines need to be modified.
- Subspace iteration fails to find all eigenvalues in the interval when the number of eigenvalues is greatly underestimated by the DOS estimation.
