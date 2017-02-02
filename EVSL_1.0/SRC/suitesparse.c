#include <stdio.h>
#include <stdlib.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"
#include "cholmod.h"
#include "umfpack.h"

void umfpack_solvefunc(int n, double *br, double *bz, double *xr, double *xz,
                       void *data) {
  /*-------------------------------------------------------------------------
   * complex linear solver routine passed to evsl
   * NOTE: This function MUST be of this prototype
   * INPUT:
   *   n: size of the system
   *   br, bz: vectors of length n, complex right-hand side (real and imaginary)
   *   data: all data that are needed for solving the system
   * OUTPUT:
   *   xr, xz: vectors of length n, complex solution (real and imaginary)
   *-------------------------------------------------------------------------*/
  void* Numeric = data;
  double Control[UMFPACK_CONTROL]; 
  umfpack_zl_defaults(Control);
  Control[UMFPACK_IRSTEP] = 0; // no iterative refinement for umfpack 
  umfpack_zl_solve(UMFPACK_A, NULL, NULL, NULL, NULL, xr, xz, br, bz, 
                   Numeric, Control, NULL); 
}

/* set default solver */
int set_ratf_solfunc_default(csrMat *A, ratparams *rat) {
  int i, j, n, nnz, nnz2=0, status, *diag;
  SuiteSparse_long *Ap, *Ai;
  double *Ax, *Az;
  void *Symbolic=NULL, *Numeric=NULL;

  n = A->nrows;
  nnz = A->ia[n];

  /* save the position of diag entry of each row */
  Malloc(diag, n, int);
  /* UMFPACK matrix */
  Malloc(Ap, n+1, SuiteSparse_long);
  /* allocate nnz+n spaces, in the worst case: A has no nonzero diag entry */
  Malloc(Ai, nnz+n, SuiteSparse_long);
  Malloc(Ax, nnz+n, double);
  Malloc(Az, nnz+n, double);

  /* copy A to Ap, Ai, Ax, Az
   * we consider general cases where A may not have full diagonal
   * this code will handle these cases */
  Ap[0] = 0;
  for (i=0; i<n; i++) {
    int rowi_has_diag = 0;
    for (j=A->ia[i]; j<A->ia[i+1]; j++) {
      int col = A->ja[j];
      double val = A->a[j];
      /* copy value and col idx */
      Ai[nnz2] = col;
      Ax[nnz2] = val;
      Az[nnz2] = 0.0;
      /* mark diagonal entry */
      if (col == i) {
        diag[i] = nnz2;
        rowi_has_diag = 1;
      }
      nnz2++;
    }
    /* if diag not found, we need to create it */
    if (!rowi_has_diag) {
      Ai[nnz2] = i;
      Ax[nnz2] = 0.0;
      Az[nnz2] = 0.0;
      diag[i] = nnz2;
      nnz2++;
    }
    /* done with row i */
    Ap[i+1] = nnz2;
  }

  /* for each pole we shift the diagonal and factorize */
  for (i=0; i<rat->num; i++) {
    /* the complex shift for pole i */
    double zkr = creal(rat->zk[i]);
    double zkc = cimag(rat->zk[i]);
    /* shift the diagonal */
    for (j=0; j<n; j++) {
      int dj = diag[j];
      Ax[dj] -= zkr;
      Az[dj] -= zkc;
    }
    
    /* only do symbolic factorization once */
    if (i==0) {
      /* Symbolic Factorization */
      status = umfpack_zl_symbolic (n, n, Ap, Ai, Ax, Az, &Symbolic, NULL, NULL);
      if (status < 0) {
        printf("umfpack_zl_symbolic failed, %d\n", status);
        return 1;
      }
    }
    /* Numerical Factorization */
    status = umfpack_zl_numeric(Ap, Ai, Ax, Az, Symbolic, &Numeric, NULL, NULL);
    if (status < 0) {
      printf("umfpack_zl_numeric failed and exit, %d\n", status);
      return 1;
    }

    /* set solver pointer and data */
    rat->solshift[i] = umfpack_solvefunc;
    rat->solshiftdata[i] = Numeric;
    
    /* shift the diagonal back */
    for (j=0; j<n; j++) {
      int dj = diag[j];
      Ax[dj] += zkr;
      Az[dj] += zkc;
    }
  } // for (i=0; i<rat->num,...)

  /* free the symbolic fact */
  if (Symbolic) {
    umfpack_zl_free_symbolic(&Symbolic);
  }

  free(diag);
  free(Ap);
  free(Ai);
  free(Ax);
  free(Az);

  return 0;
}

void free_rat_default_sol(ratparams *rat) {
  int i;
  if (rat->use_default_solver) {
    for (i=0; i<rat->num; i++) {
      umfpack_zl_free_numeric(&rat->solshiftdata[i]);
    }
  }
}

