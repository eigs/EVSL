#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"
#include "evsl_suitesparse.h"

/** @file suitesparse.c
 *  @brief Default solver function for solving shifted system and factoring B
 *
 *  This file contains the default solver function for solving shifted system 
 *  with A - SIGMA*I or A - SIGMA*B, and for factoring B matrix.
 *  The default solvers are UMFPACK and CHOLMOD from SuiteSparse
 *
 */

/**
 * @brief default complex linear solver routine passed to evsl
 * 
 * @param n       size of the system
 * @param br,bz   vectors of length n, complex right-hand side (real and imaginary)
 * @data          all data that are needed for solving the system
 * 
 * @param[out] xr,xz     vectors of length n, complex solution (real and imaginary)
 *
 * @warning: This function MUST be of this prototype
 *
 *------------------------------------------------------------------*/
void ASIGMABSolSuiteSparse(int n, double *br, double *bz, double *xr, 
                           double *xz, void *data) {
  void* Numeric = data;
  double Control[UMFPACK_CONTROL]; 
  umfpack_zl_defaults(Control);
  Control[UMFPACK_IRSTEP] = 0; // no iterative refinement for umfpack
  umfpack_zl_solve(UMFPACK_A, NULL, NULL, NULL, NULL, xr, xz, br, bz, 
                   Numeric, Control, NULL);
}

/** @brief setup the default solver for A - SIGMA B 
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
 * */
int SetupASIGMABSolSuiteSparse(csrMat *A, csrMat *BB, int num,
                               complex double *zk, void **data) {
  int i, j, nrow, ncol, nnzB, nnzC, *map, status;
  csrMat *B, C, eye;
  /* UMFPACK matrix for the shifted matrix 
   * C = A - s * B */
  SuiteSparse_long *Cp, *Ci;
  double *Cx, *Cz, zkr1;
  void *Symbolic=NULL, *Numeric=NULL;
  
  nrow = A->nrows;
  ncol = A->ncols;

  if (BB) {
    B = BB;
  } else {
    /* if B==NULL, B=I, standard e.v. prob */
    speye(nrow, &eye);
    B = &eye;
  }

  nnzB = B->ia[nrow];
  /* map from nnz in B to nnz in C, useful for multi-poles */
  Malloc(map, nnzB, int);
  /* NOTE: SuiteSparse requires that the matrix must be sorted.
   * The matadd routine can guarantee this
   * C = A + 0.0 * B 
   * This actually computes the pattern of A + B
   * and also can guarantee C has sorted rows 
   * map is the mapping from nnzB to nnzC */
  matadd(1.0, 0.0, A, B, &C, NULL, map);

  nnzC = C.ia[nrow];
  /* malloc and copy to SuiteSparse matrix */
  Malloc(Cp, nrow+1, SuiteSparse_long);
  Malloc(Ci, nnzC, SuiteSparse_long);
  Calloc(Cz, nnzC, double);
  for (i=0; i<nrow+1; i++) {
    Cp[i] = C.ia[i];
  }
  for (i=0; i<nnzC; i++) {
    Ci[i] = C.ja[i];
  }
  Cx = C.a;
  /* pole loop
   * for each pole we shift with B and factorize */
  zkr1 = 0.0;
  for (i=0; i<num; i++) {
    /* the complex shift for pole i */
    double zkr = creal(zk[i]);
    double zkc = cimag(zk[i]);

    // shift B
    for (j=0; j<nnzB; j++) {
      int p = map[j];
      double v = B->a[j];
      CHKERR(Ci[p] != B->ja[j]);
      Cx[p] -= (zkr - zkr1) * v;
      Cz[p] = -zkc * v;
    }

    /* only do symbolic factorization once */
    if (i==0) {
      /* Symbolic Factorization */
      status = umfpack_zl_symbolic(nrow, ncol, Cp, Ci, Cx, Cz, &Symbolic, 
                                   NULL, NULL);
      if (status < 0) {
        printf("umfpack_zl_symbolic failed, %d\n", status);
        return 1;
      }
    }
    /* Numerical Factorization */
    status = umfpack_zl_numeric(Cp, Ci, Cx, Cz, Symbolic, &Numeric, NULL, NULL);
    if (status < 0) {
      printf("umfpack_zl_numeric failed and exit, %d\n", status);
      return 1;
    }
    /* save the data */
    data[i] = Numeric;
    /* for the next shift */
    zkr1 = zkr;
  } /* for (i=...)*/
  
  /* free the symbolic fact */
  if (Symbolic) {
    umfpack_zl_free_symbolic(&Symbolic);
  }

  free(map);
  free(Cp);
  free(Ci);
  free(Cz);
  free_csr(&C);
  if (!BB) {
    free_csr(&eye);
  }

  return 0;
}

/**
 * @brief free the data needed by UMFpack
 */ 
void FreeASIGMABSolSuiteSparse(int num, void **data) {
  int i;
  for (i=0; i<num; i++) {
    umfpack_zl_free_numeric(&data[i]);
  }
}




/**
 * @brief Create cholmod_sparse matrix just as a wrapper of a csrMat
 * This will be useful for the B matrix, since it will be factored
 * @warning cholmod_sparse is a CSC format. But since B is symmetric, it is the same */
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

/** @brief convert an array to cholmod dense (no copying)
 *
 * */
cholmod_dense* arr_to_cholmod_dense(int m, int n, double *A, int lda) {
  cholmod_dense *B = NULL;
  Malloc(B, 1, cholmod_dense);
  B->nrow = m;
  B->ncol = n;
  B->nzmax = m*n;
  B->d = lda;
  B->x = A;
  B->z = NULL;
  B->xtype = CHOLMOD_REAL;
  B->dtype = CHOLMOD_DOUBLE;

  return B;
}

/** @brief Solver function of B (XXX) need to optimize
 *
 * */
void BSolSuiteSparse(double *b, double *x, void *data) {
  int n;

  BSolDataSuiteSparse *Bsol_data = (BSolDataSuiteSparse *) data;
  cholmod_dense *cholmod_x, *cholmod_b;
  cholmod_factor *LB;
  cholmod_common *cc;

  LB = Bsol_data->LB;
  cc = &Bsol_data->cm;
  n = LB->n;
  cholmod_b = arr_to_cholmod_dense(n, 1, b, n);
  cholmod_x = cholmod_solve(CHOLMOD_A, LB, cholmod_b, cc);
  memcpy(x, cholmod_x->x, n*sizeof(double));

  free(cholmod_b);
  cholmod_free_dense(&cholmod_x, cc);
}

/** @brief Setup the B-sol by computing the Cholesky factorization of B
 *
 * @param B       matrix B
 * */
int SetupBSolSuiteSparse(csrMat *B, BSolDataSuiteSparse *Bsol_data) {
  cholmod_sparse *Bcholmod;
  cholmod_common *cc = &Bsol_data->cm;
  cholmod_factor *LB;

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
  cholmod_check_factor(LB, cc);
  /* save the struct to global variable */
  Bsol_data->LB = LB;
  /* check the factor */
  CHKERR(LB->is_ll == 0);
  /* free the matrix wrapper */
  free(Bcholmod);

  return 0;
}

int FreeBSolSuiteSparseData(BSolDataSuiteSparse *data) {
  cholmod_common *cc = &data->cm;
  /* free the factor */
  cholmod_free_factor(&data->LB, cc);
  /* finish cholmod */
  cholmod_finish(cc);
  return 0;
}

/** @brief Solver function of L^{T} 
 *  w = L^{-T}v
 * */
void LTSolSuiteSparse(double *v, double *w, void *data) {
  int n;

  BSolDataSuiteSparse *Bsol_data = (BSolDataSuiteSparse *) data;
  cholmod_dense *cholmod_w, *cholmod_w2, *cholmod_b;
  cholmod_factor *LB;
  cholmod_common *cc;

  LB = Bsol_data->LB;
  cc = &Bsol_data->cm;
  n = LB->n;
  cholmod_b = arr_to_cholmod_dense(n, 1, v, n);
  cholmod_w = cholmod_solve(5, LB, cholmod_b, cc);
  cholmod_w2 = cholmod_solve(8, LB, cholmod_w, cc); 
  memcpy(w, cholmod_w2->x, n*sizeof(double));
  free(cholmod_b);
  cholmod_free_dense(&cholmod_w, cc);
  cholmod_free_dense(&cholmod_w2, cc);
}



