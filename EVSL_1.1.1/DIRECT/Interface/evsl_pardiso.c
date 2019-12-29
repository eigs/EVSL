#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
//#include "mkl.h"
#include "mkl_pardiso.h"
//#include "mkl_types.h"
#include "evsl_direct.h"

/**
 * @file evsl_pardiso.c
 * @brief Definitions used for pardiso interface
 */


typedef struct _BSolDataDirect {
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  void *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  /* upper triangular matrix */
  csrMat U;
#ifdef EVSL_USING_CUDA_GPU
  /* workspace for solve, of size 2*n, only for GPU runs */
  double *w;
#endif
} BSolDataDirect;


typedef struct _ASBSolDataDirect {
  /* Internal solver memory pointer pt, */
  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
  /* or void *pt[64] should be OK on both architectures */
  void *pt[64];
  /* Pardiso control parameters. */
  MKL_INT iparm[64];
  /* upper triangular matrix */
  MKL_INT n, *ia, *ja, maxfct, mnum;
  MKL_Complex16 *a;
  /* workspace for solve, of size n */
  MKL_Complex16 *b, *x;
#ifdef EVSL_USING_CUDA_GPU
  /* workspace for solve, of size 2*n, only for GPUs */
  double *w;
#endif
} ASBSolDataDirect;

/** @brief Setup the B-sol by computing the Cholesky factorization of B
 *
 * @param B         matrix B
 * @param data Struct which will be initialized
 * */
int SetupBSolDirect(csrMat *B, void **data) {
  double tms = evsl_timer();
  BSolDataDirect *Bsol_data;
  Bsol_data = evsl_Malloc(1, BSolDataDirect);

  /* PAPT = LLT for symmetric positive-definite matrices */
  MKL_INT mtype = 2;       /* Real SPD matrix */
  /* Maximum number of factors with identical sparsity structure that
   * must be kept in memory at the same time. In most applications this
   * value is equal to 1. It is possible to store several different
   * factorizations with the same nonzero structure at the same time in
   * the internal data structure management of the solver.
   * pardiso can process several matrices with an identical matrix sparsity
   * pattern and it can store the factors of these matrices at the same
   * time. Matrices with a different sparsity structure can be kept in
   * memory with different memory address pointers pt.*/
  MKL_INT maxfct = 1;
  /* Indicates the actual matrix for the solution phase. With this scalar
   * you can define which matrix to factorize. The value must be:
   * 1 <= mnum <= maxfct */
  MKL_INT mnum = 1;
  /* Controls the execution of the solver. Usually it is a two- or
   * three-digit integer. The first digit indicates the starting phase
   * of execution and the second digit indicates the ending phase. */
  MKL_INT phase;
  /* Message level information. If msglvl = 0 then pardiso generates no
   * output, if msglvl = 1 the solver prints statistical information
   * to the screen */
  MKL_INT msglvl = 0;
  /* Initialize error flag */
  MKL_INT error = 0;
  /* Double dummy */
  double ddum;
  /* Integer dummy. */
  MKL_INT idum;
  /* Number of right hand sides. */
  MKL_INT nrhs = 1;

  /* extract upper diag part of B */
  csrMat *U = &Bsol_data->U;
  triuCsr(B, U);

  MKL_INT  n  = U->nrows;
  MKL_INT *ia = U->ia;
  MKL_INT *ja = U->ja;
  double  *a  = U->a;

  /* This function initializes Intel MKL PARDISO internal address pointer
   * pt with zero values (as needed for the very first call of pardiso)
   * and sets default iparm values in accordance with the matrix type. */
#if 1

  pardisoinit(Bsol_data->pt, &mtype, Bsol_data->iparm);

#else
  int i;
  for ( i = 0; i < 64; i++ )
  {
    Bsol_data->iparm[i] = 0;
  }
  for ( i = 0; i < 64; i++ )
  {
    Bsol_data->pt[i] = 0;
  }

  Bsol_data->iparm[0] = 1;         /* No solver default */
  Bsol_data->iparm[1] = 2;         /* Fill-in reordering from METIS */
  Bsol_data->iparm[3] = 0;         /* No iterative-direct algorithm */
  Bsol_data->iparm[4] = 0;         /* No user fill-in reducing permutation */
  Bsol_data->iparm[5] = 0;         /* Write solution into x */
  Bsol_data->iparm[6] = 0;         /* Not in use */
  Bsol_data->iparm[7] = 2;         /* Max numbers of iterative refinement steps */
  Bsol_data->iparm[8] = 0;         /* Not in use */
  Bsol_data->iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
  Bsol_data->iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
  Bsol_data->iparm[11] = 0;        /* Not in use */
  Bsol_data->iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric).
                                      Try iparm[12] = 1 in case of inappropriate accuracy */
  Bsol_data->iparm[13] = 0;        /* Output: Number of perturbed pivots */
  Bsol_data->iparm[14] = 0;        /* Not in use */
  Bsol_data->iparm[15] = 0;        /* Not in use */
  Bsol_data->iparm[16] = 0;        /* Not in use */
  Bsol_data->iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
  Bsol_data->iparm[18] = -1;       /* Output: Mflops for LU factorization */
  Bsol_data->iparm[19] = 0;        /* Output: Numbers of CG Iterations */
#endif

  /* Zero-based indexing of columns and rows.*/
  Bsol_data->iparm[34] = 1;
  /* Max numbers of iterative refinement steps.
   * 0: The solver automatically performs two steps of iterative refinement
   * when perturbed pivots are obtained during the numerical factorization.*/
  Bsol_data->iparm[7] = 0;
  /* Matrix checker */
  Bsol_data->iparm[26] = 1;

  /* --------------------------------------------------------------- */
  /* Reordering and Symbolic Factorization. This step also allocates */
  /* all memory that is necessary for the factorization.             */
  /* --------------------------------------------------------------- */
  phase = 11;
  pardiso(Bsol_data->pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
          &idum, &nrhs, Bsol_data->iparm, &msglvl, &ddum, &ddum, &error);

  if ( error != 0 )
  {
    printf ("\nERROR during symbolic factorization: %d", error);
    exit (1);
  }

  /* printf ("\nReordering completed ... "); */
  printf ("\nNumber of nonzeros in B %d", B->ia[B->nrows]);
  printf ("\nPeak      Memory on symbolic  fact.           %d KB", Bsol_data->iparm[14]);
  printf ("\nPermanent Memory on symbolic  fact.           %d KB", Bsol_data->iparm[15]);
  printf ("\n          Memory on numerical fact. and solve %d KB", Bsol_data->iparm[16]);
  printf ("\nNumber of nonzeros in factors = %d", Bsol_data->iparm[17]);
  printf ("\nNumber of factorization MFLOPS = %d\n", Bsol_data->iparm[18]);

  /* --------------------------------------------------------------- */
  /* Numerical factorization.                                        */
  /* --------------------------------------------------------------- */
  phase = 22;
  pardiso(Bsol_data->pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
          &idum, &nrhs, Bsol_data->iparm, &msglvl, &ddum, &ddum, &error);

  if ( error != 0 )
  {
    printf ("\nERROR during numerical factorization: %d", error);
    exit (2);
  }
  /*
  printf ("\nFactorization completed ... ");
  */

#ifdef EVSL_USING_CUDA_GPU
  Bsol_data->w = evsl_Malloc(2*n, double);
#endif

  *data = (void *) Bsol_data;

  double tme = evsl_timer();
  evslstat.t_setBsv += tme - tms;

  return 0;
}

/** @brief Solver function of B
 *
 * */
void BSolDirect(double *b, double *x, void *data) {
  BSolDataDirect *Bsol_data = (BSolDataDirect *) data;

  MKL_INT maxfct = 1;
  MKL_INT mnum = 1;
  MKL_INT mtype = 2;       /* Real SPD matrix */
  MKL_INT msglvl = 0;
  MKL_INT error = 0;
  /* Integer dummy. */
  MKL_INT idum;
  /* Number of right hand sides. */
  MKL_INT nrhs = 1;

  csrMat *U = &Bsol_data->U;
  MKL_INT  n  = U->nrows;
  MKL_INT *ia = U->ia;
  MKL_INT *ja = U->ja;
  double  *a  = U->a;

#ifdef EVSL_USING_CUDA_GPU
  double *w = Bsol_data->w;
  double *x_device = x;
  double *b_host = w;
  double *x_host = w + n;
  evsl_memcpy_device_to_host(b_host, b, n*sizeof(double));
  b = b_host;
  x = x_host;
#endif

  /* ------------------------------------------------------------- */
  /*  Back substitution and iterative refinement.                  */
  /* ------------------------------------------------------------- */
  /* 33: Solve, iterative refinement */
  MKL_INT phase = 33;
  pardiso(Bsol_data->pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
          &idum, &nrhs, Bsol_data->iparm, &msglvl, b, x, &error);
  if ( error != 0 )
  {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }

  /*
  double *r;
  r = evsl_Malloc(n, double);
  char uplo = 'U';
  mkl_cspblas_dcsrsymv(&uplo, &n, a, ia, ja, x, r);
  double res = 0.0, res0 = 0.0;
  int i;
  for (i = 0; i < n; i++) {
    res += (r[i] - b[i]) * (r[i] - b[i]);
    res0 += b[i] * b[i];
  }
  res = sqrt (res) / sqrt (res0);
  if (res > 1e-12) {
    printf ("\nRelative residual = %e", res);
  }
  */
  /*printf ("\nSolve completed ... ");*/

#ifdef EVSL_USING_CUDA_GPU
  evsl_memcpy_host_to_device(x_device, x, n*sizeof(double));
#endif
}

/** @brief Solver function of L^{T}
 *  x = L^{-T}*b
 * */
void LTSolDirect(double *b, double *x, void *data) {
  BSolDataDirect *Bsol_data = (BSolDataDirect *) data;

  MKL_INT maxfct = 1;
  MKL_INT mnum = 1;
  MKL_INT mtype = 2;       /* Real SPD matrix */
  MKL_INT msglvl = 0;
  MKL_INT error = 0;
  /* Integer dummy. */
  MKL_INT idum;
  /* Number of right hand sides. */
  MKL_INT nrhs = 1;

  csrMat *U = &Bsol_data->U;
  MKL_INT  n  = U->nrows;
  MKL_INT *ia = U->ia;
  MKL_INT *ja = U->ja;
  double  *a  = U->a;

#ifdef EVSL_USING_CUDA_GPU
  double *w = Bsol_data->w;
  double *x_device = x;
  double *b_host = w;
  double *x_host = w + n;
  evsl_memcpy_device_to_host(b_host, b, n*sizeof(double));
  b = b_host;
  x = x_host;
#endif

  /* 333: like phase=33, but only backward substitution */
  MKL_INT phase = 333;
  pardiso(Bsol_data->pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja,
          &idum, &nrhs, Bsol_data->iparm, &msglvl, b, x, &error);
  if ( error != 0 )
  {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }
  /*printf ("\nSolve completed ... ");*/

#ifdef EVSL_USING_CUDA_GPU
  evsl_memcpy_host_to_device(x_device, x, n*sizeof(double));
#endif
}

/** @brief Free solver data
 * */
void FreeBSolDirectData(void *data) {
  BSolDataDirect *Bsol_data = (BSolDataDirect *) data;

  MKL_INT maxfct = 1;
  MKL_INT mnum = 1;
  MKL_INT mtype = 2;       /* Real SPD matrix */
  MKL_INT msglvl = 0;
  MKL_INT error = 0;
  /* Double dummy */
  double ddum;
  /* Integer dummy. */
  MKL_INT idum;
  /* Number of right hand sides. */
  MKL_INT nrhs = 1;

  csrMat *U = &Bsol_data->U;
  MKL_INT  n  = U->nrows;
  MKL_INT *ia = U->ia;
  MKL_INT *ja = U->ja;

  /* -------------------------------------------------------------------- */
  /* .. Termination and release of memory. */
  /* -------------------------------------------------------------------- */
  MKL_INT phase = -1;           /* Release internal memory. */
  pardiso(Bsol_data->pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja,
          &idum, &nrhs, Bsol_data->iparm, &msglvl, &ddum, &ddum, &error);

  if ( error != 0 )
  {
    printf ("\nERROR during termination: %d", error);
    exit (3);
  }
  /*printf ("\nTermination completed ... ");*/

  free_csr(U);
#ifdef EVSL_USING_CUDA_GPU
  evsl_Free(Bsol_data->w);
#endif
  evsl_Free(data);
}

/** @brief setup Pardiso solver for A - SIGMA B
 *
 * The setup invovles shifting the matrix A - SIGMA B and factorizing
 * the shifted matrix
 * The solver function and the data will be saved data
 * Generally speaking, each pole can have a different solver
 *
 * @param A        matrix A
 * @param BB       matrix B, if NULL, it means B is identity
 * @param num      the number of SIGMA's
 * @param zk       array of SIGMA's of length num
 * @param data    all data that are needed for solving the system
 * */
int SetupASIGMABSolDirect(csrMat *A, csrMat *BB, int num,
                          EVSL_Complex *zk, void **data) {
  double tms = evsl_timer();
  int i, j, nrow, /* ncol, */ nnzUB, nnzUC, *map;
  csrMat *B, UC, eye, UA, UB;
  /* the shifted matrix
   * C = A - s * B */
  MKL_INT *UCp, *UCi, n;
  double *UCx;
  /* complex values of C */
  MKL_Complex16 *UCz;
  ASBSolDataDirect *ASBdata, *ASBdata0=NULL;

  MKL_INT mtype = 6;       /* Complex and symmetric matrix */
  MKL_INT maxfct = num;
  MKL_INT mnum;
  MKL_INT phase;
  MKL_INT msglvl = 0;
  MKL_INT error = 0;
  double ddum;
  MKL_INT idum;
  MKL_INT nrhs = 1;

  nrow = A->nrows;
  /* ncol = A->ncols; */

  if (BB) {
    B = BB;
  } else {
    /* if B==NULL, B=I, standard e.v. prob */
    speye(nrow, &eye);
    B = &eye;
  }

  /* upper triangular parts of A and B */
  triuCsr(A, &UA);
  triuCsr(B, &UB);

  /* nnz in the upper tirangular part of B */
  nnzUB = UB.ia[nrow];
  /* map from nnz in UB to nnz in UC, useful for multi-poles */
  map = evsl_Malloc(nnzUB, int);

  /* NOTE: if the matrix entries need to be sorted.
   * The matadd routine will  guarantee this
   * C = A + 0.0 * B
   * This actually computes the pattern of A + B
   * and also can guarantee C has sorted rows
   * map is the mapping from nnzB to nnzC */
  matadd(1.0, 0.0, &UA, &UB, &UC, NULL, map);

  /* arrays for C */
  n   = nrow;
  UCp = UC.ia;
  UCi = UC.ja;
  UCx = UC.a;
  nnzUC = UCp[nrow];

  /* pole loop
   * for each pole we shift with B and factorize */
  for (i = 0; i < num; i++) {
    /* data for pole i */
    ASBdata = evsl_Malloc(1, ASBSolDataDirect);

    /* the complex shift for pole i */
    double zkr = creal(zk[i]);
    double zkc = cimag(zk[i]);

    UCz = evsl_Malloc(nnzUC, MKL_Complex16);
    /* make the values of UC complex */
    for (j = 0; j < nnzUC; j++) {
      UCz[j].real = UCx[j];
      UCz[j].imag = 0.0;
    }

    /* shift with UB */
    for (j = 0; j < nnzUB; j++) {
      int p = map[j];
      double v = UB.a[j];
      if (UCi[p] != UB.ja[j]) {
         fprintf(stderr, "error %s %d\n", __FILE__, __LINE__);
         return 1;
      }
      if (!(p >= 0 && p < nnzUC)) {
         fprintf(stderr,"error %s %d\n", __FILE__, __LINE__);
         return 2;
      }
      UCz[p].real -= zkr * v;
      UCz[p].imag -= zkc * v;
    }

    mnum = i + 1;
    /* only do init and symbolic factorization once */
    if (i == 0) {

      pardisoinit(ASBdata->pt, &mtype, ASBdata->iparm);

      /* Zero-based indexing of columns and rows.*/
      ASBdata->iparm[34] = 1;
      /* Max numbers of iterative refinement steps.
       * 0: The solver automatically performs two steps of iterative refinement
       * when perturbed pivots are obtained during the numerical factorization.*/
      ASBdata->iparm[7] = 0;
      /* Matrix checker */
      ASBdata->iparm[26] = 1;

      /* --------------------------------------------------------------- */
      /* Reordering and Symbolic Factorization. This step also allocates */
      /* all memory that is necessary for the factorization.             */
      /* --------------------------------------------------------------- */
      phase = 11;
      pardiso(ASBdata->pt, &maxfct, &mnum, &mtype, &phase, &n, UCz, UCp, UCi,
              &idum, &nrhs, ASBdata->iparm, &msglvl, &ddum, &ddum, &error);

      if ( error != 0 )
      {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
      }

      ASBdata0 = ASBdata;
    } else {
      /* use the pointer and parm for pole-0 */
      memcpy(ASBdata->pt,    ASBdata0->pt,    64*sizeof(void*));
      memcpy(ASBdata->iparm, ASBdata0->iparm, 64*sizeof(MKL_INT));
    }

    /* --------------------------------------------------------------- */
    /* Numerical factorization.                                        */
    /* --------------------------------------------------------------- */
    phase = 22;
    pardiso(ASBdata->pt, &maxfct, &mnum, &mtype, &phase, &n, UCz, UCp, UCi,
            &idum, &nrhs, ASBdata->iparm, &msglvl, &ddum, &ddum, &error);

    if ( error != 0 )
    {
      printf ("\nERROR during numerical factorization: %d", error);
      exit (2);
    }

    printf ("\nNumber of nonzeros in triu(A - sig*B)         %d",    UCp[n]);
    printf ("\nPeak      Memory on symbolic  fact.           %d KB", ASBdata->iparm[14]);
    printf ("\nPermanent Memory on symbolic  fact.           %d KB", ASBdata->iparm[15]);
    printf ("\n          Memory on numerical fact. and solve %d KB", ASBdata->iparm[16]);
    printf ("\nNumber of nonzeros in factors                 %d",    ASBdata->iparm[17]);
    printf ("\nNumber of factorization MFLOPS                %d\n",  ASBdata->iparm[18]);

    /* save the data */
    ASBdata->maxfct = maxfct;
    ASBdata->mnum = mnum;
    ASBdata->n = n;
    ASBdata->ia = UCp;
    ASBdata->ja = UCi;
    ASBdata->a = UCz;
    ASBdata->b = evsl_Malloc(nrow, MKL_Complex16);
    ASBdata->x = evsl_Malloc(nrow, MKL_Complex16);
#ifdef EVSL_USING_CUDA_GPU
    ASBdata->w = evsl_Malloc(2*nrow, double);
#endif
    data[i] = ASBdata;
  } /* for (i=...)*/

  evsl_Free(map);
  evsl_Free(UC.a);
  free_csr(&UA);
  free_csr(&UB);
  if (!BB) {
    free_csr(&eye);
  }

  double tme = evsl_timer();
  evslstat.t_setASigBsv += tme - tms;

  return 0;
}

/**
 * @brief complex linear solver routine passed to evsl
 *
 * @param n       size of the system
 * @param br,bi   vectors of length n, complex right-hand side (real and imaginary)
 * @param data    all data that are needed for solving the system
 *
 * @param[out] xr,xz     vectors of length n, complex solution (real and imaginary)
 *
 * @warning: This function MUST be of this prototype
 *
 *------------------------------------------------------------------*/
void ASIGMABSolDirect(int n, double *br, double *bi, double *xr,
                      double *xz, void *data) {

  ASBSolDataDirect *sol_data = (ASBSolDataDirect *) data;
  int i;
  MKL_INT mtype = 6;       /* Complex and symmetric matrix */
  MKL_INT msglvl = 0;
  MKL_INT error = 0;
  /* Integer dummy. */
  MKL_INT idum;
  /* Number of right hand sides. */
  MKL_INT nrhs = 1;

  if (n != sol_data->n) {
     fprintf(stderr, "error %s %d\n", __FILE__, __LINE__);
     exit(-1);
  }
  MKL_Complex16 *b = sol_data->b;
  MKL_Complex16 *x = sol_data->x;

#ifdef EVSL_USING_CUDA_GPU
  double *r_host = sol_data->w;
  double *i_host = sol_data->w + n;
  evsl_memcpy_device_to_host(r_host, br, n*sizeof(double));
  evsl_memcpy_device_to_host(i_host, bi, n*sizeof(double));
  br = r_host;
  bi = i_host;
#endif

  /* copy rhs */
  for (i = 0; i < n; i++) {
    b[i].real = br[i];
    b[i].imag = bi[i];
  }

  /* ------------------------------------------------------------- */
  /*  Back substitution and iterative refinement.                  */
  /* ------------------------------------------------------------- */
  /* 33: Solve, iterative refinement */
  MKL_INT phase = 33;
  pardiso(sol_data->pt, &sol_data->maxfct, &sol_data->mnum, &mtype, &phase,
          &sol_data->n, sol_data->a, sol_data->ia, sol_data->ja,
          &idum, &nrhs, sol_data->iparm, &msglvl, b, x, &error);

  if ( error != 0 )
  {
    printf ("\nERROR during solution: %d", error);
    exit (3);
  }

  /*
  MKL_Complex16 *r;
  r = evsl_Malloc(n, MKL_Complex16);
  char uplo = 'U';
  mkl_cspblas_zcsrsymv(&uplo, &sol_data->n, sol_data->a, sol_data->ia, sol_data->ja, x, r);
  double res = 0.0, res0 = 0.0;
  for (i = 0; i < n; i++) {
    res += (r[i].real - b[i].real) * (r[i].real - b[i].real);
    res += (r[i].imag - b[i].imag) * (r[i].imag - b[i].imag);
    res0 += b[i].real * b[i].real;
    res0 += b[i].imag * b[i].imag;
  }
  res = sqrt (res) / sqrt (res0);
  //if (res > 1e-12) {
    printf ("\nRelative residual = %e\n", res);
  //}

  exit(0);
  */

#ifdef EVSL_USING_CUDA_GPU
  double *xr_device = xr;
  double *xz_device = xz;
  xr = r_host;
  xz = i_host;
#endif

  /* copy sol */
  for (i = 0; i < n; i++) {
    xr[i] = x[i].real;
    xz[i] = x[i].imag;
  }

#ifdef EVSL_USING_CUDA_GPU
  evsl_memcpy_host_to_device(xr_device, xr, n*sizeof(double));
  evsl_memcpy_host_to_device(xz_device, xz, n*sizeof(double));
#endif
}

/**
 * @brief free the data needed by Pardiso
 */
void FreeASIGMABSolDirect(int num, void **data) {
  int i;
  MKL_INT mtype = 6;       /* Complex and symmetric matrix */
  MKL_INT msglvl = 0;
  MKL_INT phase = -1;           /* Release internal memory. */
  MKL_INT error = 0;
  /* Integer dummy. */
  MKL_INT idum;
  /* Double dummy */
  double ddum;
  MKL_INT nrhs = 1;

  for (i=0; i<num; i++) {

    ASBSolDataDirect *sol_data = (ASBSolDataDirect *) data[i];

    if (i == 0) {
      /* ----------------------------------------------------------------- */
      /* Termination and release of memory.                                */
      /* ----------------------------------------------------------------- */
      pardiso(sol_data->pt, &sol_data->maxfct, &sol_data->mnum, &mtype, &phase,
              &sol_data->n, sol_data->a, sol_data->ia, sol_data->ja,
              &idum, &nrhs, sol_data->iparm, &msglvl, &ddum, &ddum, &error);

      if ( error != 0 )
      {
        printf ("\nERROR during termination: %d", error);
        exit (3);
      }

      evsl_Free(sol_data->ia);
      evsl_Free(sol_data->ja);
    }

    evsl_Free(sol_data->a);
    evsl_Free(sol_data->b);
    evsl_Free(sol_data->x);
#ifdef EVSL_USING_CUDA_GPU
    evsl_Free(sol_data->w);
#endif
    evsl_Free(sol_data);
  }
}

