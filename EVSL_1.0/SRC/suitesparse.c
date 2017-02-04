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

/* @brief Create cholmod_sparse matrix just as a wrapper of a csrMat 
 * @warning cholmod_sparse is a CSC format. But since A is symmetric, it is the same */
cholmod_sparse* csrMat_to_cholmod_sparse(csrMat *A, int stype) {
  cholmod_sparse *B = NULL;
  Malloc(B, 1, cholmod_sparse);
  B->nrow = A->nrows;
  B->ncol = A->ncols;
  B->nzmax = A->ia[A->nrows];
  /* column pointers */
  B->p = A->ia;
  /* row indices */
  B->i = A->ja;
  B->nz = NULL;
  B->x = A->a;
  B->z = NULL;
  B->stype = stype;
  B->itype = CHOLMOD_INT;
  B->xtype = CHOLMOD_REAL;
  B->dtype = CHOLMOD_DOUBLE;
  B->sorted = 0;
  B->packed = 1;

  return B;
}

typedef struct _default_LBdata {
  cholmod_common cc;
  cholmod_factor *LB;
  cholmod_dense *Y, *E, *W;
} default_LBdata;

void vector_to_cholmod_dense(int nrow, int ncol, double *v, int ldv,
                             cholmod_dense *x) {
  x->nrow = nrow;
  x->ncol = ncol;
  x->nzmax = nrow * ncol;
  x->d = ldv;
  x->x = v;
  x->z = NULL;
  x->xtype = CHOLMOD_REAL;
  x->dtype = CHOLMOD_DOUBLE;
}

/*
 * soltype = 1 : x = P' * L' \ B
 *         = 2 : x = L \ P * B
 */ 
void cholmod_sol_combine(int soltype, double *b, int nb, int ldb,
                         double *x, int nx, int ldx, void *data) {
  int err1=1, err2=1, n;
  default_LBdata *LBdata = (default_LBdata *) data;
  cholmod_factor *LB = LBdata->LB;
  cholmod_common *cc = &LBdata->cc;
  n = LB->n;
  /* cholmod solve work space */
  cholmod_dense **Y, **E;
  Y = &LBdata->Y;
  E = &LBdata->E;
  cholmod_dense *oldY, *oldE;
  oldY = LBdata->Y;
  oldE = LBdata->E;
  /* convert to cholmod dense */
  cholmod_dense B, X, *Xptr;
  vector_to_cholmod_dense(n, nb, b, ldb, &B);
  vector_to_cholmod_dense(n, nb, x, ldx, &X);
  Xptr = &X;
  cholmod_dense *Wptr = LBdata->W;
  /* solve with LB */
  if (1 == soltype) {
    /* W = L' \ B */
    err1 = cholmod_solve2(5, LB, &B, NULL, &Wptr, NULL, Y, E, cc);
    /* X = P' * W */
    err2 = cholmod_solve2(8, LB, Wptr, NULL, &Xptr, NULL, Y, E, cc);
  } else if (2 == soltype) {
    /* W = P * B */
    err1 = cholmod_solve2(7, LB, &B, NULL, &Wptr, NULL, Y, E, cc);
    /* X = L \ W */
    err2 = cholmod_solve2(4, LB, Wptr, NULL, &Xptr, NULL, Y, E, cc);
  }
  /* W, X should have NOT been reallocated */
  CHKERR(Xptr != &X);
  CHKERR(Wptr != LBdata->W);
  /* Y and E should be only alloced once */
  if (oldY) { CHKERR(oldY != *Y); }
  if (oldE) { CHKERR(oldE != *E); }

  if (err1 != 1 || err2 != 1) {
    printf("cholmod solve2 error err1 %d, err2 %d\n", err1, err2);
  }
}

void default_LSol(double *x, double *y, void *data) {
  cholmod_sol_combine(2, x, y, data);
}

void default_LTSol(double *x, double *y, void *data) {
  cholmod_sol_combine(1, x, y, data);
}

int factor_Bmatrix_default(csrMat *B) {
  int n = B->nrows;
  cholmod_sparse *Bcholmod;
  default_LBdata *LBdata;

  /* unset B just in case it was not freed */
  //if (evsldata.hasB && evsldata.isDefaultLB) {
  //  free_Bfactor_default();
  //}

  Malloc(LBdata, 1, default_LBdata);
  cholmod_common *cc = &LBdata->cc;
  /* start CHOLMOD */
  cholmod_start(cc);
  /* force to have LL factor */
  cc->final_asis = 0;
  cc->final_ll = 1;
  /* convert matrix. 
   * stype=1 means the upper triangular part of B will be accessed */
  Bcholmod = csrMat_to_cholmod_sparse(B, 1);
  /* check common and the matrix */
  cholmod_check_common(cc);
  cholmod_check_sparse(Bcholmod, cc);
  /* symbolic and numeric fact */
  LBdata->LB = cholmod_analyze(Bcholmod, cc);
  cholmod_factorize(Bcholmod, LBdata->LB, cc);
  /* check the factor */
  CHKERR(LBdata->LB->is_ll == 0);
  cholmod_check_factor(LBdata->LB, cc);
  /* set NULL to workspace Y and E,
   * will be allocated at the first call of solve */
  LBdata->Y = NULL;
  LBdata->E = NULL;
  /* allocate workspace W */
  LBdata->W = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, cc);
  /* save the struct to global variable */
  evsldata.LBdata = (void *) LBdata;
  /* free the matrix wrapper */
  free(Bcholmod);

  evsldata.LBsolv = default_LSol;
  evsldata.LBTsolv = default_LTSol;
  
  return 0;
}

void free_Bfactor_default() {
  default_LBdata *LBdata = (default_LBdata *) evsldata.LBdata;
  cholmod_factor *LB = LBdata->LB;
  cholmod_common *cc = &LBdata->cc;
  cholmod_dense *Y = LBdata->Y;
  cholmod_dense *E = LBdata->E;
  cholmod_dense *W = LBdata->W;
  cholmod_free_factor(&LB, cc);
  cholmod_free_dense(&Y, cc);
  cholmod_free_dense(&E, cc);
  cholmod_free_dense(&W, cc);
  cholmod_finish(cc);
}

