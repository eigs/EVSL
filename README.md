

## EVSL:  EigenValues Slicing Library (Version 1.1.1)

 ![Alt](evsl_logo3.png "EVSL-logo")

~~~
                          ChebLanNr, ChebLanTr, ChebSI, RatLanNr and RatLanTr
               Polynomial and Rational Filtered Lanczos and subspace iteration algorithms
                                 For Symmetric Eigenvalue problems
~~~

Welcome to EVSL. EVSL is a  C library for computing the eigenvalues of
a  symmetric  matrix that  are  located  in  a given  interval.   EVSL
implements  a polynomial  filtered Lanczos  method (thick  restart, no
restart)  a  rational  filtered  Lanczos  method  (thick  restart,  no
restart),  and a  polynomial filtered  subspace iteration  for solving
standard  eigenvalue  problems  A  u  =  &lambda;  u  and  generalized
eigenvalue problems A u = &lambda;  B u.  EVSL also provides tools for
spectrum slicing, i.e., the technique  of subdividing a given interval
into  p smaller  subintervals and  computing the  eigenvalues in  each
subinterval independently. This is achieved by estimating the spectral
density of the  matrix (or matrix pair for  generalized problems) with
either the Kernel  Polynomial method (KPM) or a  Lanczos based method.
This release does not yet offer full parallel implementations (trivial
openMP test programs are available  among the test drivers).  Versions
1.1.x  have much  added functionality  relative to  version 1.0.  Most
notably:

+ EVSL now handles both standard and generalized problems
+ Matrix-free mechanism  [whereby user provides own matvec function]
+ Interface for FORTRAN users
+ Spectrum slicing by KPM as well as Lanczos. Also effective
   spectrum slicing for Ax=&lambda; Bx
+ Interfaces to various direct solvers that are needed by the code
+ Recently added: Generalized eigenvalue problems without factorizations
  [using polynomials in B]

For questions/feedback send e-mail to Yousef Saad [saad@umn.edu]

-----------------------------------------------------------------------

###  HANDBOOK AND DOCUMENTATION

-----------------------------------------------------------------------

An article describing EVSL is available from
[SIAM Journal on Scientific Computing](https://epubs.siam.org/doi/abs/10.1137/18M1170935), which
can be cited as follows [bibtex entry]
```
@article{doi:10.1137/18M1170935,
author = {Li, R. and Xi, Y. and Erlandson, L. and Saad, Y.},
title = {The Eigenvalues Slicing Library (EVSL): Algorithms, Implementation, and Software},
journal = {SIAM Journal on Scientific Computing},
volume = {41},
number = {4},
pages = {C393-C415},
year = {2019},
doi = {10.1137/18M1170935},
URL = {https://doi.org/10.1137/18M1170935},
eprint = {https://doi.org/10.1137/18M1170935}}
```

A detailed documentation of the code can be found in the `Documentation'
folder of the package in EVSL_1.1.x

-----------------------------------------------------------------------

###  INSTALLATION

-----------------------------------------------------------------------

**Library**: To build the EVSL library on Unix-like operating system, such as Linux and macOS,
 the simplest way is to run, in the sub-directory `EVSL_1.1.1`,
```
./configure
```
which, upon successful completion, creates file `makefile.in`, EVSL configuration header file `./INC/EVSL_config.h`, file
`config.status` containing information of the current configuration and file `config.log` containing  messages produced by compilers while running configure,
then followed by
```
make
```
which compiles EVSL with default options and installs the library in `./EVSL/include` and `./EVSL/lib`.
There are many options to `configure` and `make` to customize installation directories, compiler selections, compile and load flags, etc. For example,
```
./configure --prefix=evsl-install-dir
make install
```
copies the built EVSL libraries and include files into `evsl-install-dir`, and
```
./configure --with-openmp
```
enables EVSL to use OpenMP threading. For a complete list of options and the descriptions, type
```
./configure --help
```
**The user may directly modify `makefile.in` and `./INC/EVSL_config.h` for customizations.**

**Test programs**:
   In the directories under TESTS/, you will find a number of directories that contain sample drivers
   of different functions in EVSL for spectral slicing, computing eigenvalues, and computing DOS
   (density of states). To compile, type
```
make test
```
in sub-directory `EVSL_1.1.1`.

-----------------------------------------------------------------------

###  BLAS and LAPACK

-----------------------------------------------------------------------
**Note that EVSL can be built without external BLAS and LAPACK being provided** and
can also optionally use external BLAS and LAPACK libraries provided by the user, by doing, e.g.,
```
./configure --with-blas-lib=/usr/local/lib/libblas.a --with-lapack-lib=/usr/local/lib/liblapack.a
```
or
```
./configure --with-blas-lib-name=blas --with-blas-lib-dir=/usr/local/lib --with-lapack-lib-name=lapack
            --with-lapack-lib-dir=/usr/lib
```
On macOS, configure EVSL to use Apple Accelerate framework by
```
./configure --with-framework-accelerate
```


-----------------------------------------------------------------------

###  DIRECT SOLVERS

-----------------------------------------------------------------------

**CXSparse**
   CXSparse is the default direct linear solver used in EVSL for the linear systems
   with complex symmetric matrix (A-SIGMA I) or (A-SIGMA B) that arise in rational filtering methods and
   linear systems with SPD matrix B in generalized eigenvalue problems.
   A (modified) copy of CXSparse is included in `EVSL_1.1.1/DIRECT/CXSparse`. Using alternative high performance direct solvers (see below) can often yield better overall performance of EVSL.

>  NOTE: CXSparse, which is  distributed with EVSL, is Copyrighted by
   Timothy Davis. As noted above much better performance can be
   achieved by other existing direct solvers.
   Refer to CXSparse package for its License. [http://faculty.cse.tamu.edu/davis/suitesparse.html]


**SuiteSparse**:
   EVSL can use UMFPACK in SuiteSparse to solve linear systems with (A-SIGMA I) or (A-SIGMA B),
   and  CHOLMOD to solve linear systems with B. To configure EVSL with SuiteSparse, run
   ```
   ./configure --with-suitesparse-dir=DIR --with-blas-lib-name=LIB --with-lapack-lib-name=LIB
               --with-metis-dir=DIR
   ```

>  NOTE:  SuiteSparse is NOT distributed with EVSL, and it is Copyrighted by Timothy Davis.
>  Refer to the SuiteSparse package for its license. [http://faculty.cse.tamu.edu/davis/suitesparse.html]

**Pardiso**:
   Pardiso in Intel MKL [https://software.intel.com/en-us/mkl-developer-reference-fortran-intel-mkl-pardiso-parallel-direct-sparse-solver-interface] is also supported. To configure EVSL with MKL Pardiso, run
   ```
   ./configure --with-mkl-pardiso --with-intel-mkl
   ```
   ensuring that `MKLROOT` is set correctly.

**Interface**:
   Interface functions to call CXSparse, SuiteSparse and Pardiso from EVSL are
   provided in `EVSL_1.1.1/DIRECT/Interface`.
   Other solvers can also be used by writing an interface of the same type.


###  RATIONAL FILTERING

-----------------------------------------------------------------------
  Rational  filtering  requires  solving  linear  systems  (where  the
  coefficient matrix  is the  original matrix  shifted by  a **complex**
  shift).  A linear solver routine  must be provided.

  After  having  computed  the  rational  filter  by calling

```
int find_ratf(double *intv, ratparams *rat)
```
  the users  can call
```
int SetASigmaBSol(ratparams *rat, int i, SolFuncC func, void *data)
```
  to set a linear solver for the ith pole of the rational filter, where
  `func` is a function pointer  for solving  linear systems  with  pole `i`, of  which the
  coefficient matrix  is `A - s_i I` or `A - s_i B` for generalized eigenvalue problems with `s_i = rat->zk[i]`
  being the *complex* shift, and  `data`  points to the data that are needed  by `func`.

  Note that `func` must be of the following prototype
```
void SolFuncC(int n, double *br, double *bz, double *xr, double *xz, void *data);
```
  where `n`  is the size  of the system,  `br`, `bz` are  the right-hand
  side (real and  imaginary parts of complex vector),  `xr`, `xz` will
  hold the  solution (complex vector) on return,  and `data` contains  all the
  data  needed  for  the  solver.

  Once `SetASigmaBSol` is executed, rational filtering Lanczos methods
  should be ready to use.

-----------------------------------------------------------------------

###  MATRIX-FREE EIGEN-SOLVERS

-----------------------------------------------------------------------
  All the eigensolvers  in EVSL can be used in  `matrix-free' mode. In
  this  mode,  users need  only  to  provide their  own  matrix-vector
  product function of the following prototype:

      void MVFunc(double *x, double *y, void *data);

  where y  = A *  x and data  is the pointer  to the associated  data to
  perform the  matvec. The  `(void *)`  argument is  to provide  a uniform
  interface  to all  user-specific  functions. For  a particular  Matvec
  function, one can pack all data  needed by this function into a struct
  and pass the  pointer of this struct  to EVSL (after casting  it to `(void
  *)`). This  function needs to  be passed to EVSL  as well, so  EVSL can
  call  it  to  perform  all matvecs.    The user can also pass needed
  matrices in the standard CSR format to EVSL, in which case EVSL will use
  its internal   MATVEC routine. This can be set by

      SetAMatrix(csrMat *A)
      SetBMatrix(csrMat *B)

  for the matrices A and B


  In  TESTS/LAP,  an example  of  matvec  functions for  2D/3D  Laplacian
  matrices is  provided, where the  matrix is not explicitly  formed but
  5pt/7pt stencil is used instead to perform the matvecs.
  In  this example, a struct for matvec is first defined:

      typedef struct _lapmv_t {
        int  nx,  ny,  nz;
        double  *stencil;
      } lapmv_t;

  and the matvec function is implemented as:

      void Lap2D3DMatvec(double  *x, double  *y,  void  *data) {
        lapmv_t *lapmv = (lapmv_t *) data;
        int nx = lapmv->nx;
        int ny = lapmv->ny;
        int nz = lapmv->nz;
        double *stencil = lapmv->stencil;
        ...
      }

  in  which   the  pointer  is  first   cast  and  all  the   data  is
  unpacked. Once these are ready, they can be passed to EVSL by calling

      SetAMatvec(n, &Lap2D3DMatvec, (void*) &lapmv) and
      SetBMatvec(n, &Lap2D3DMatvec, (void*) &lapmv)

  to set the matvec routines for A and B respectively,
  where the first input is the size of the "matrix", the second input is
  the  function pointer  and the  third one  is the  data pointer.  Once
  `SetMatvecFunc`  is  called,  EVSL  will  use  the  registered  matvec
  function to perform all matvecs with A.

  Users should first create a function  wrapper of the above type for an
  external matvec routine. Then, following  the steps in the example, it
  will be straightforward to use it in EVSL.

-----------------------------------------------------------------------

###  GENERALIZED EIGENVALUE PROBLEMS

-----------------------------------------------------------------------
  For solving the generalized eigenproblem A x = &lambda; B  x, the users
  must also provide a solver  for the B matrix by calling

      SetBSol(SolFuncR func, void *data).

  To indicate to EVSL that a generalized eigenvalue problem is being solve,
  one must call

      SetGenEig()

  since by default, EVSL assumes that a  standard eigenvalue problem is being
  solved even when a B matrix is provided. Call the function

      SetStdEig()

  for solving standard eigenvalue problem

  To implement spectrum slicing, the current version of EVSL requires system
  solutions with the matrix L', where B=L*L' is the Cholesky factorization of B.
  Call the function

    SetLTSol(SolFuncR func, void *data)

  to set the solver function for L'.

 Version  1.1.1   includes  drivers  for  solving   certain  types  of
 generalized eigenvalue  problems by  polynomial filtering  by methods
 that avoid matrix factorizations.

-----------------------------------------------------------------------

###  CUDA GPU Support

-----------------------------------------------------------------------
 All the eigensolvers in EVSL can be accelerated by Nvidia CUDA GPUs. To enable the GPU support,
 configure EVSL with
 ```
 ./configure --with-cuda
 ```
 with `CUDA_HOME` and `EVSL_CUDA_SM` (the default is 60) set properly. 
 All the test drivers in `TESTS` should work on GPUs with the exception of `Fortran/LapPLanN_MatFree.f90` that requires
 writing CUDA fortran kernels (with PGI CUDA Fortran). It should be clear from these drivers that how the user needs 
 to handle the GPU memory for the matrices and vectors that are given to EVSL solvers.

-----------------------------------------------------------------------

###  Initialization and Completion

-----------------------------------------------------------------------
* Use `EVSLStart()` and `EVSLFinish()` before and after any call to the EVSL functions



-----------------------------------------------------------------------

### DETAILED DESCRIPTION OF CONTENTS

-----------------------------------------------------------------------

 * INC
   - EVSL_config.h     : configuration header file (created from CONFIG/EVSL_config.h.in)
   - evsl.h            : user-level function prototypes and constant definitions
   - internal_header.h : internal header file for SRC/
   - struct.h          : miscellaneous structs used in evsl

 * SRC
   - blas/         :  required C blas routines
   - lapack/       :  required C lapack routines
   - cheblanNr.c   :  Polynomial Filtered no-restart Lanczos
   - cheblanTr.c   :  Polynomial Filtered thick restart Lanczos
   - chebpoly.c    :  Computing and applying polynomial filters
   - chebsi.c      :  Polynomial Filtered Subspace iteration
   - dos_utils.c   :  Miscellaneous functions used for DOS based functions.
   - dumps.c       :  Miscellaneous functions for I/O and for debugging
   - evsl.c        :  Set EVSL solver options and data
   - evsl_f90.c    :  Fortran interface
   - evsl_memory.c :  EVSL memory management
   - exDOS.c       :  Exact DOS given eigenvalues
   - lanTrbounds.c :  A more robust alg. to compute bounds of spectrum based on TR Lanczos
   - lanbounds.c   :  Lanczos alg. to compute bounds of spectrum
   - landos.c      :  Lanczos based DOS algorithm for the standard problem
   - landosG.c     :  Lanczos based DOS algorithm for generalized and standard problems
   - misc_la.c     :  Miscellaneous linear algebra functions
   - ratfilter.c   :  Computing and applying rational filters
   - ratlanNr.c    :  Rational Filtered no-restart Lanczos
   - ratlanTr.c    :  Rational Filtered thick restart Lanczos
   - simpson.c     :  Simpson integrator
   - spmat.c       :  Sparse matrix routines
   - spslice.c     :  Spectrum slicing functions for Kernel Polynomial Method
   - spslice2.c    :  Spectrum slicing functions for Lanczos
   - stats.c       :  Various statistics
   - timing.c      :  Timer
   - vect.c        :  Vector operations

 * TESTS/          : Test drivers

   * Fortran/      : Fortran test drivers

   * PLanR           : Polynomial Filtered Thick Restart Lanczos test drivers
     * LapPLanR.c    : Polynomial filtering T-R Lanczos (Lapliacian)
     * MMPLanR.c     : Polynomial filtering T-R Lanczos (Matrix Market)
     * MMPLanR_omp.c : Polynomial filtering T-R Lanczos (parallelized with OMP for slices) (Matrix Market)

   * RLanR           : Rational Filtered Thick Restart Lanczos test drivers
     * LapRLanR.c    : Rational filtering T-R Lanczos (Laplacian)
     * MMRLanR.c     : Rational filtering T-R Lanczos (Matrix Market)

   * PLanN                : Polynomial Filtered No-Restart Lanczos test drivers
     * MMPLanN.c          : Polynomial filtering non-restart Lanczos (Matrix Market)
     * LapPLanN.c         : Polynomial filtering non-restart Lanczos (Laplacian)
     * LapPLanN_MatFree.c : "matrix-free" version: not forming matrix but passing mat-vec function (Laplacian)

   * RLanN           : Rational Filtered No-Restart Lanczos test drivers
     * LapRLanN.c    : Rational filtering non-restart Lanczos (Laplacian)
     * MMRLanN.c     : Rational filtering non-restart Lanczos (Matrix Market)

   * GEN             : Test drivers for the generalized eigenvalue problem
     * MMsimple.c    : A stripped down (no slicing) drivier based on MMPLanN.c
     * MMPLanN.c     : Polynomial filtering non-restart Lanczos (Matrix Market)
     * MMPLanR.c     : Polynomial filtering T-R Lanczos (Matrix Market)
     * MMRLanN.c     : Rational filtering non-restart Lanczos (Matrix Market)
     * MRLanR.c      : Rational filtering T-R Lanczos (Matrix Market)

   * PSI             : Test drivers for polynomial filter subspace iteration
     * LapPSI.c      : Polynomial filtering subspace iterations (Laplacian)
     * MMPSI.c       : Polynomial filtering subspace iterations (Matrix Market)

   * COMMON          : Routines common to the test drivers
      * io.c         : parse command-line input parameters
      * lapl.c       : Build Laplacian matrices and compute the exact eigenvalues of Laplacians
      * mmio.c       : IO routines for the matrix market format

   * COMMON_GEN      : Routines common to the test drivers for the generalized eigenvalue problem
      * io.c         : parse command-line input parameters
      * lapl.c       : Build Laplacian matrices and compute the exact eigenvalues of Laplacians
      * mmio.c       : IO routines for the matrix market format

   * TESTS/Landos    : test drivers for the lanDOS related functions.
      * LanDos.c     : DOS for standard eigenvalue problem  using Lanczos
      * LanDosG.c    : DOS for generalized eigenvalue problem  using Lanczos

 * DIRECT
   - CXSparse/                    : a copy of modified CXSparse
   - Interface/evsl_suitesparse.c : suitesparse UMFPACK and CHOLMOD interface
   - Interface/evsl_cxsparse.c    : cxsparse interface
   - Interface/evsl_pardiso.c     : pardiso interface
   - Interface/evsl_direct.h      : header for all direct solver interface



-----------------------------------------------------------------------
