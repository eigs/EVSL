#ifndef EVSL_SSP_H
#define EVSL_SSP_H

#include "cholmod.h"
#include "umfpack.h"

typedef struct _BSolDataSuiteSparse {
  cholmod_factor *LB;
  cholmod_common cm;
} BSolDataSuiteSparse;

/* solver routine with SuiteSparse */
void BSolSuiteSparse(double *b, double *x, void *data);

int SetupBSolSuiteSparse(csrMat *B, BSolDataSuiteSparse *data);

int FreeBSolSuiteSparseData(BSolDataSuiteSparse *data);

void ASIGMABSolSuiteSparse(int n, double *br, double *bz, double *xr, 
                           double *xz, void *data);

int SetupASIGMABSolSuiteSparse(csrMat *A, csrMat *BB, int num,
                               complex double *zk, void **data);

void FreeASIGMABSolSuiteSparse(int num, void **data);

void LTSolSuiteSparse(double *v, double *w, void *data);

#endif

