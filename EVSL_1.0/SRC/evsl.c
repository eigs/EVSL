#include <stdlib.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/* global variable evslData, which is guaranteed to be initialized */
evslData evsldata;

void EVSLStart() {
  evsldata.Amatvec.n = -1;
  evsldata.Amatvec.func = NULL;
  evsldata.Amatvec.data = NULL;
  evsldata.hasB = 0;
  evsldata.isDefaultLB = 0;
  evsldata.LB_mult = NULL;
  evsldata.LBT_mult = NULL;
  evsldata.LB_solv = NULL;
  evsldata.LBT_solv = NULL;
  evsldata.LB_func_data = NULL;
  evsldata.matvec_gen_work = NULL;
}

void EVSLFinish() {
  if (evsldata.hasB && evsldata.isDefaultLB) {
    free_default_LBdata();
    free(evsldata.LB_func_data);
  }
  if (evsldata.matvec_gen_work) {
    free(evsldata.matvec_gen_work);
  }
}

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
  err = set_default_LBdata(B);
  evsldata.hasB = 1;
  evsldata.isDefaultLB = 1;
  /* alloc some workspace */
  Malloc(evsldata.matvec_gen_work, 2*B->nrows, double);
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
    free_default_LBdata();
    free(evsldata.LB_func_data);
  }
#endif
  evsldata.hasB = 0;
  evsldata.isDefaultLB = 0;
  evsldata.LB_mult = NULL;
  evsldata.LBT_mult = NULL;
  evsldata.LB_solv = NULL;
  evsldata.LBT_solv = NULL;
  evsldata.LB_func_data = NULL;
  if (evsldata.matvec_gen_work) {
    free(evsldata.matvec_gen_work);
    evsldata.matvec_gen_work = NULL;
  }
}

/* matvec routine */
int matvec_genev(csrMat *A, double *x, double *y) {
  /* if B is not set, so just y = A * x */
  if (!evsldata.hasB) {
    matvec_A(A, x, y);
    return 0;
  }
  /* for gen e.v, y = L \ A / L' *x */
  double *w = evsldata.matvec_gen_work;
  evsldata.LBT_solv(x, y, evsldata.LB_func_data);
  matvec_A(A, y, w);
  evsldata.LB_solv(w, y, evsldata.LB_func_data);
  return 0;
}

