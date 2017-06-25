#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "def.h"
#include "struct.h"
#include "CXSparse/Include/cs.h"
#include "internal_proto.h"
#include "evsl_direct.h"

/**
 * @file evsl_cxsparse.h
 * @brief Definitions used for cxsparse interface
 */

typedef struct _BSolDataDirect {
  /* problem size */
  int n;
  /* symbolic factor */
  cs_dis *S;
  /* numeric factor */
  cs_din *N;
  /* workspace for solve, of size n */
  double *w;
} BSolDataDirect;

typedef struct _ASBSolDataDirect {
  /* problem size */
  int n;
  /* symbolic factor */
  cs_cis *S;
  /* numeric factor */
  cs_cin *N;
  /* workspace for solve, of size n */
  cs_complex_t *b, *x;
} ASBSolDataDirect;


/* true for off-diagonal entries */
static int dropdiag_di (int i, int j, double aij, void *other) { return (i != j) ;}

/* C = A + triu(A,1)' */
static cs_di *make_sym_di (cs_di *A)
{
    cs_di *AT, *C ;
    AT = cs_di_transpose (A, 1) ;          /* AT = A' */
    cs_di_fkeep (AT, &dropdiag_di, NULL) ;    /* drop diagonal entries from AT */
    C = cs_di_add (A, AT, 1, 1) ;          /* C = A+AT */
    cs_di_spfree (AT) ;
    return (C) ;
}

/** @brief Setup the B-sol by computing the Cholesky factorization of B
 *
 * @param B         matrix B
 * @param Bsol_data Struct which will be initialized 
 * */
int SetupBSolDirect(csrMat *B, void **data) {
  BSolDataDirect *Bsol_data;
  Malloc(Bsol_data, 1, BSolDataDirect);
  int i,j,k, n, nnz, *csp, *csi;
  double *csx;
  /* csparse csc matrix and the symmetrized one */
  cs_di Bcs, *Bsym;
  /* csparse factors: symbolic and numeric */
  cs_dis *S;
  cs_din *N;
  
  n = B->nrows;
  nnz = B->ia[n];
  /* allocate memory for Bcs */
  Malloc(csp, n+1, int);
  Malloc(csi, nnz, int);
  Malloc(csx, nnz, double);
  csp[0] = 0;
  /* copy the lower part of B into cs */
  /* nnz in cs */
  k = 0;
  for (i=0; i<n; i++) {
    for (j=B->ia[i]; j<B->ia[i+1]; j++) {
      int c = B->ja[j];
      if (i >= c) {
        csi[k] = c;
        csx[k] = B->a[j];
        k++;
      }
    }
    csp[i+1] = k;
  }
  /* can view this as csc of upper part of B */
  Bcs.nzmax = k;
  Bcs.m = n;
  Bcs.n = n;
  Bcs.p = csp;
  Bcs.i = csi;
  Bcs.x = csx;
  Bcs.nz = -1;
  /* symmetrize B */
  Bsym = make_sym_di(&Bcs);
  if (!Bsym) {
    return -1;
  }
  /* factorization of B */
  S = cs_di_schol(1, Bsym) ;            /* ordering and symbolic analysis */
  N = cs_di_chol(Bsym, S) ;             /* numeric Cholesky factorization */

  if (!(S && N)) {
    return -2;
  }

  //printf("Lnz %d %d\n", N->L->nzmax, N->L->p[N->L->n]);

  /* free */
  cs_di_spfree(Bsym);
  free(csp);
  free(csi);
  free(csx);

  Bsol_data->n = n;
  Bsol_data->S = S;
  Bsol_data->N = N;
  Malloc(Bsol_data->w, n, double);

  *data = (void *) Bsol_data;
  
  return 0;
}

/** @brief Solver function of B
 *
 * */
void BSolDirect(double *b, double *x, void *data) {
  BSolDataDirect *Bsol_data = (BSolDataDirect *) data;
  cs_dis *S = Bsol_data->S;
  cs_din *N = Bsol_data->N;
  int n = Bsol_data->n;
  double *w = Bsol_data->w;

  cs_di_ipvec (S->pinv, b, w, n) ;   /* w = P*b */
  cs_di_lsolve (N->L, w) ;           /* w = L\w */
  cs_di_ltsolve (N->L, w) ;          /* w = L'\w */
  cs_di_pvec (S->pinv, w, x, n) ;    /* x = P'*w */
}

/** @brief Solver function of L^{T} 
 *  x = L^{-T}*b
 * */
void LTSolDirect(double *b, double *x, void *data) {
  BSolDataDirect *Bsol_data = (BSolDataDirect *) data;
  cs_dis *S = Bsol_data->S;
  cs_din *N = Bsol_data->N;
  int n = Bsol_data->n;
  double *w = Bsol_data->w;

  memcpy(w, b, n*sizeof(double));    /* w = b */
  cs_di_ltsolve (N->L, w) ;          /* w = L'\w */
  cs_di_pvec (S->pinv, w, x, n) ;    /* x = P'*w */
}

/** @brief Free solver data
 * */
void FreeBSolDirectData(void *data) {
  BSolDataDirect *Bsol_data = (BSolDataDirect *) data;
  cs_dis *S = Bsol_data->S;
  cs_din *N = Bsol_data->N;
  double *w = Bsol_data->w;
  
  cs_di_sfree(S);
  cs_di_nfree(N);
  free(w);
  free(data);
}

/** @brief setup CXsparse solver for A - SIGMA B 
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
                          complex double *zk, void **data) {
  int i, j, nrow, ncol, nnzB, nnzC, *map;
  csrMat *B, C, eye;
  /* the shifted matrix 
   * C = A - s * B */
  int *Cp, *Ci;
  double *Cx, zkr1;
  /* complex values of C */
  cs_complex_t *Cz;
  /* cs matrix complex integer */
  cs_ci Ccs;
  /* csparse factors: symbolic and numeric */
  cs_cis *S = NULL;
  cs_cin *N = NULL;
  ASBSolDataDirect *ASBdata;
  
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
  /* NOTE: if the matrix entries need to be sorted.
   * The matadd routine will  guarantee this
   * C = A + 0.0 * B 
   * This actually computes the pattern of A + B
   * and also can guarantee C has sorted rows 
   * map is the mapping from nnzB to nnzC */
  matadd(1.0, 0.0, A, B, &C, NULL, map);
  /* arrays for C */
  nnzC = C.ia[nrow];
  Cp = C.ia;
  Ci = C.ja;
  Cx = C.a;
  /* complex values of C */
  Malloc(Cz, nnzC, cs_complex_t);
  for (i=0; i<nnzC; i++) {
    Cz[i] = Cx[i] + 0.0 * I;
  }
  /* since C is complex symmetric, so its CSR and CSC are the same 
   * put them in a cs matrix */
  Ccs.nzmax = nnzC;
  Ccs.m = nrow;
  Ccs.n = ncol;
  Ccs.p = Cp;
  Ccs.i = Ci;
  Ccs.x = Cz;
  Ccs.nz = -1;
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
      double r = creal(Cz[p]) - (zkr - zkr1) * v;
      Cz[p] = r - zkc * v * I;
    }

    /* only do symbolic factorization once */
    if (i == 0) {
      /* Symbolic Factorization */
      /* order == 1: sym, 0: not QR */
      S = cs_ci_sqr(1, &Ccs, 0);     /* ordering and symbolic analysis */
      if (!S) {
        return -1;
      }
    }
    /* Numerical Factorization */
    /*  1.0 results in conventional partial pivoting. */
    /*  0.0 results in no pivoting. */
    N = cs_ci_lu(&Ccs, S, 0.0); 
    if (!N) {
      return -2;
    }

    //printf("Lnz %d %d\n", N->L->nzmax, N->L->p[N->L->n]);
    //printf("Unz %d %d\n", N->U->nzmax, N->U->p[N->U->n]);
    //printf("%e\n", (N->L->nzmax + N->U->nzmax + 0.0) / nnzC);
    
    /* save the data */
    Malloc(ASBdata, 1, ASBSolDataDirect);
    ASBdata->n = nrow;
    ASBdata->S = S;
    ASBdata->N = N;
    Malloc(ASBdata->b, nrow, cs_complex_t);
    Malloc(ASBdata->x, nrow, cs_complex_t);
    data[i] = ASBdata;

    /* for the next shift */
    zkr1 = zkr;
  } /* for (i=...)*/
  
  free(map);
  free(Cz);
  free_csr(&C);
  if (!BB) {
    free_csr(&eye);
  }

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
  CHKERR(n != sol_data->n);
  cs_cis *S = sol_data->S;
  cs_cin *N = sol_data->N;
  cs_complex_t *b = sol_data->b;
  cs_complex_t *x = sol_data->x;
  /* copy rhs */
  for (i=0; i<n; i++) {
    b[i] = br[i] + bi[i] * I;
  }
  cs_ci_ipvec (N->pinv, b, x, n) ;       /* x = b(p) */
  cs_ci_lsolve (N->L, x) ;               /* x = L\x */
  cs_ci_usolve (N->U, x) ;               /* x = U\x */
  cs_ci_ipvec (S->q, x, b, n) ;          /* b(q) = x */
  /* copy sol */
  for (i=0; i<n; i++) {
    xr[i] = creal(b[i]);
    xz[i] = cimag(b[i]);
  }
}

/**
 * @brief free the data needed by CXSparse
 */ 
void FreeASIGMABSolDirect(int num, void **data) {
  int i;
  for (i=0; i<num; i++) {
    ASBSolDataDirect *soldata = (ASBSolDataDirect *) data[i];
    if (i == 0) {
      cs_ci_sfree(soldata->S);
    }
    cs_ci_nfree(soldata->N);
    free(soldata->b);
    free(soldata->x);
    free(soldata);
  }
}

