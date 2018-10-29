#ifndef EVSL_DIRECT_H
#define EVSL_DIRECT_H

#include "evsl.h"

/**
 * @file evsl_direct.h
 * @brief Definitions used for direct solver interface
 *
 * Note that this file is meant to be the header file
 * of evsl_suitesparse.c, evsl_pardiso.c, and evsl_cxsparse.c.
 * If more direct solver options will be added later,
 * this should serve them as well.
 */

/* functions for B solve */
int  SetupBSolDirect(csrMat *B, void **data);
void BSolDirect(double *b, double *x, void *data);
void LTSolDirect(double *b, double *x, void *data);
void FreeBSolDirectData(void *data);

/* functions for A-SIGMA*B solve */
int  SetupASIGMABSolDirect(csrMat *A, csrMat *BB, int num, EVSL_Complex *zk, void **data);
void ASIGMABSolDirect(int n, double *br, double *bi, double *xr, double *xz, void *data);
void FreeASIGMABSolDirect(int num, void **data);

/* evsl_direct_f90.c */
#ifdef __cplusplus
extern "C" {
#endif
void EVSLFORT(setup_bsol_direct,SETUP_BSOL_DIRECT)(uintptr_t *Bf90, uintptr_t *Bsoldataf90);
void EVSLFORT(free_bsol_direct,FREE_BSOL_DIRECT)(uintptr_t *Bsolf90);
void EVSLFORT(setup_asigmabsol_direct,SETUP_ASIGMABSOL_DIRECT)(uintptr_t *Af90, int *flagB, uintptr_t *Bf90, uintptr_t *ratf90, uintptr_t *solshiftf90);
void EVSLFORT(free_asigmabsol_direct,FREE_ASIGMABSOL_DIRECT)(uintptr_t *ratf90, uintptr_t *solshiftf90);
#ifdef __cplusplus
}
#endif

#endif

