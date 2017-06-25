#ifndef EVSL_DIRECT_H
#define EVSL_DIRECT_H

/**
 * @file evsl_direct.h
 * @brief Definitions used for direct solver interface
 * 
 * Note that this file is meant to be the header file
 * of both evsl_suitesparse.c and evsl_cxsparse.c.
 * If more direct solver options will be added later,
 * this should serve them as well.
 */

/* functions for B solve */
int  SetupBSolDirect(csrMat *B, void **data);
void BSolDirect(double *b, double *x, void *data);
void LTSolDirect(double *b, double *x, void *data);
void FreeBSolDirectData(void *data);

/* functions for A-SIGMA*B solve */
int  SetupASIGMABSolDirect(csrMat *A, csrMat *BB, int num,
                          complex double *zk, void **data);
void ASIGMABSolDirect(int n, double *br, double *bi, double *xr, 
                      double *xz, void *data);
void FreeASIGMABSolDirect(int num, void **data);

#endif

