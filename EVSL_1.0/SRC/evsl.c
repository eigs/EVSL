#include <stdlib.h>
#include "struct.h"

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

