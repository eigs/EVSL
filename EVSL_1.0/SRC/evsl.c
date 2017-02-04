#include <stdlib.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/* global variable evslData, which is guaranteed to be initialized */
evslData evsldata;

void SetMatvecFunc(int n, MVFunc func, void *data) {
  evsldata.Amatvec.n = n;
  evsldata.Amatvec.func = func;
  evsldata.Amatvec.data = data;
}

void UnsetMatvecFunc() {
  evsldata.Amatvec.n = -1;
  evsldata.Amatvec.func = NULL;
  evsldata.Amatvec.data = NULL;
}

int SetRhsMatrix(csrMat *B) {
  int err;
#ifdef EVSL_WITH_SUITESPARSE
  err = factor_Bmatrix_default(B);
  evsldata.hasB = 1;
  evsldata.isDefaultLB = 1;
  /* malloc some workspace */
  Malloc(evsldata.LBwork, 2*B->nrows, double);
#else
  printf("error: EVSL was not compiled with SuiteSparse, ");
  printf("so the current version cannot solve generalized e.v. problem\n");
  err = -1;
#endif
  return err;
}

void UnsetRhsMatrix() {
#ifdef EVSL_WITH_SUITESPARSE
  if (evsldata.hasB && evsldata.isDefaultLB) {
    free_Bfactor_default();
    free(evsldata.LBdata);
  }
  free(evsldata.LBwork);
#endif
  evsldata.hasB = 0;
  evsldata.isDefaultLB = 0;
  evsldata.LBsolv = NULL;
  evsldata.LBTsolv = NULL;
  evsldata.LBwork = NULL;
}

int matvec_genev(csrMat *A, double *x, double *y) {
  /* if B is not set, so just y = A * x */
  if (!evsldata.hasB) {
    matvec_A(A, x, y);
    return 0;
  }
  /* for gen e.v, y = L \ A / L' *x */
  double *w = evsldata.LBwork;
  evsldata.LBTsolv(x, y, evsldata.LBdata);
  matvec_A(A, y, w);
  evsldata.LBsolv(w, y, evsldata.LBdata);
  return 0;
}

