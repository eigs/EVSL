#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  csrMat R;
  int *perm;
  double *work;
} default_LBdata;

/*
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
*/

/*
 * soltype = 1 : x = P' * L' \ b
 *         = 2 : x = L \ P * b
 */ 
void default_Lsol_combine(int soltype, double *b, double *x, void *data) {
  int n;
  default_LBdata *LBdata = (default_LBdata *) data;
  csrMat *R = &LBdata->R;
  int *p = LBdata->perm;
  double *w = LBdata->work;
  n = R->nrows;

  /* solve with L */
  if (1 == soltype) {
    /* w = L' \ b */
    tri_sol_upper('N', R, b, w);
    /* x = P' * w */
    vec_iperm(n, p, w, x);
  } else if (2 == soltype) {
    /* w = P * b */
    vec_perm(n, p, b, w);
    /* x = L \ w */
    tri_sol_upper('T', R, w, x);
  }
}

void default_LSol(double *x, double *y, void *data) {
  default_Lsol_combine(2, x, y, data);
}

void default_LTSol(double *x, double *y, void *data) {
  default_Lsol_combine(1, x, y, data);
}

int set_default_LBdata(csrMat *B) {
  int i, n = B->nrows, nnzL;
  cholmod_sparse *Bcholmod, *LBmat;
  cholmod_factor *LB;
  default_LBdata *LBdata;

  /* unset B just in case it was not freed */
  //if (evsldata.hasB && evsldata.isDefaultLB) {
  //  free_Bfactor_default();
  //}

  Malloc(LBdata, 1, default_LBdata);
  cholmod_common cm, *cc;
  cc = &cm;
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
  LB = cholmod_analyze(Bcholmod, cc);
  cholmod_factorize(Bcholmod, LB, cc);
  /* check the factor */
  CHKERR(LB->is_ll == 0);
  cholmod_check_factor(LB, cc);
  /* convert factor to sparse matrix [col format as in CHOLMOD]*/
  LBmat = cholmod_factor_to_sparse(LB, cc);
  /* check some fields of LBmat */
  if (!LBmat->packed) {
    printf("error: cholmod_sparse matrix L is not packed\n");
    CHKERR(1);
  }
  CHKERR(n != LBmat->ncol || n != LBmat->nrow);
  /* convert it to csrMat, which is *upper* triangular */
  nnzL = ((int *) (LBmat->p))[n];
  csr_resize(n, n, nnzL, &LBdata->R);
  LBdata->R.nrows = n;
  LBdata->R.ncols = n;
  memcpy(LBdata->R.ia, LBmat->p, (n+1)*sizeof(int));
  memcpy(LBdata->R.ja, LBmat->i, nnzL*sizeof(int));
  memcpy(LBdata->R.a, LBmat->x, nnzL*sizeof(double));
  /* R should be sorted */
  if (!LBmat->sorted) {
    printf("cholmod_sparse L was not sorted, so we do the sorting\n");
    /* make rows of R have increasing col ids */
    sortrow(&LBdata->R);
  }
  /* check diag of */
  if (check_full_diag('U', &LBdata->R)) {
    printf("error: R has zero diag entry!\n");
    return 1;
  }
  /* copy the perm array */
  int *cholmod_perm = (int*) LB->Perm;
  if (cholmod_perm) {
    Malloc(LBdata->perm, n, int);
    for (i=0; i<n; i++) {
      LBdata->perm[i] = cholmod_perm[i];
    }
  } else {
    LBdata->perm = NULL;
  }
  /* allocate workspace */
  Malloc(LBdata->work, n, double);
  /* save the struct to global variable */
  evsldata.LB_func_data = (void *) LBdata;
  evsldata.LB_solv = default_LSol;
  evsldata.LBT_solv = default_LTSol;
  /* free the matrix wrapper */
  free(Bcholmod);
  /* free the factor */
  cholmod_free_factor(&LB, cc);
  /* free sparse matrix */
  cholmod_free_sparse(&LBmat, cc);
  /* finish cholmod */
  cholmod_finish(cc);
  
  return 0;
}

void free_default_LBdata() {
  default_LBdata *LBdata = (default_LBdata *) evsldata.LB_func_data;
  free_csr(&LBdata->R);
  free(LBdata->perm);
  free(LBdata->work);
}

