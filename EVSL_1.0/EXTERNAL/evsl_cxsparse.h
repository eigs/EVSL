#ifndef EVSL_CXSPARSE_H
#define EVSL_CSSPARSE_H

#include "CXSparse/Include/cs.h"

/**
 * @file evsl_cxsparse.h
 * @brief Definitions used for cxsparse interface
 */

typedef struct _BSolDataCXSparse {
  /* problem size */
  int n;
  /* symbolic factor */
  cs_dis *S;
  /* numeric factor */
  cs_din *N;
  /* workspace for solve, of size n */
  double *w;
} BSolDataCXSparse;

int SetupBSolCXSparse(csrMat *B, BSolDataCXSparse *data);
void BSolCXSparse(double *b, double *x, void *data);
void LTSolCXSparse(double *b, double *x, void *data);
void FreeBSolCXSparseData(BSolDataCXSparse *data);


typedef struct _ASBSolDataCXSparse {
  /* problem size */
  int n;
  /* symbolic factor */
  cs_cis *S;
  /* numeric factor */
  cs_cin *N;
  /* workspace for solve, of size n */
  cs_complex_t *b, *x;
} ASBSolDataCXSparse;

int SetupASIGMABSolCXSparse(csrMat *A, csrMat *BB, int num,
                            complex double *zk, void **data);
void ASIGMABSolCXSparse(int n, double *br, double *bi, double *xr, 
                        double *xz, void *data);
void FreeASIGMABSolCXSparse(int num, void **data);


#endif
