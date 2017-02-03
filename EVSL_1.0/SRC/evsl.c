#include <stdlib.h>
#include "struct.h"
#include "internal_proto.h"

/* global variable evslData, which is guaranteed to be initialized */
evslData evsldata;

void SetMatvecFunc(int n, matvecFunc func, void *data) {
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
#else
  printf("error: EVSL was not compiled with the default solver \n");
  err = -1;
#endif
  return err;
}

void UnsetRhsMatrix() {
#ifdef EVSL_WITH_SUITESPARSE
  if (evsldata.hasB && evsldata.isDefaultLB) {
    free_Bfactor_default();
  }
#endif
  evsldata.hasB = 0;
  evsldata.isDefaultLB = 0;
}
