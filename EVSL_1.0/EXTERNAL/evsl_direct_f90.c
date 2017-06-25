#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"
#include "evsl_direct.h"
/** 
 * @file evsl_cxsparse_f90.c
 * @brief Fortran interface definitions
 ***/

/** @brief Fortran interface for SetupBSolDirect and
 * also SetBsol and SetLTSol
 * @param[in] Bf90: CSR matrix of B
 * @param[out] Bsoldataf90: data pointer for Bsol and LTsol
 */
void EVSLFORT(setup_bsol_direct)(uintptr_t *Bf90, 
                                 uintptr_t *Bsoldataf90) {
  /* cast csr pointer of B */
  csrMat *B = (csrMat *) (*Bf90);
  void *Bsol;
  /* setup B sol and LT sol*/
  SetupBSolDirect(B, &Bsol);
  SetBSol(BSolDirect, Bsol);
  SetLTSol(LTSolDirect, Bsol); 
  /* cast pointer for output */
  *Bsoldataf90 = (uintptr_t) Bsol;
}

/** @brief Fortran interface for FreeBSolDirectData */
void EVSLFORT(free_bsol_direct)(uintptr_t *Bsolf90) {
  /* cast pointer */
  void *Bsol = (void *) (*Bsolf90);
  FreeBSolDirectData(Bsol);
}

/** @brief Fortran interface for SetupASIGMABSolDirect 
 * @param[in] Af90: CSR matrix of A
 * @param[in] flagB: if B is present (gen. eig. problem)
 * @param[in] Bf90: CSR matrix of B
 * @param[in] ratf90: rational filter
 * @param[out] solshiftf90: pointer of solshift array
 */
void EVSLFORT(setup_asigmabsol_direct)(uintptr_t *Af90,
                                       int *flagB,
                                       uintptr_t *Bf90,
                                       uintptr_t *ratf90, 
                                       uintptr_t *solshiftf90) {
  /* cast csr pointer of A */
  csrMat *A = (csrMat *) (*Af90);
  /* cast csr pointer of B */
  csrMat *B = (*flagB) ? (csrMat *) (*Bf90) : NULL;
  /* cast pointer */
  ratparams *rat = (ratparams *) (*ratf90);
  /* allocate and setup solshiftdata */
  void **solshiftdata = (void **) malloc(rat->num*sizeof(void *));
  SetupASIGMABSolDirect(A, B, rat->num, rat->zk, solshiftdata);
  /* cast pointer for output */
  *solshiftf90 = (uintptr_t) solshiftdata;
}

/** @brief Fortran interface for SetASigmaBSol with Direct solve
 * @param[in,out] ratf90: pointer of rational filter
 * @param[out] solshiftf90: pointer of solshift array
 */
void EVSLFORT(set_asigmabsol_direct)(uintptr_t *ratf90, 
                                     uintptr_t *solshiftf90) {
  /* cast pointers */
  ratparams *rat = (ratparams *) (*ratf90);
  void **solshiftdata = (void **) (*solshiftf90);

  SetASigmaBSol(rat, NULL, ASIGMABSolDirect, solshiftdata);
}

/** @brief Fortran interface for FreeASIGMABSolDirect
 * @param ratf90: pointer of rational filter [in/out] 
 * @param solshiftf90: pointer of solshift array [in/out] 
 */
void EVSLFORT(free_asigmabsol_direct)(uintptr_t *ratf90, 
                                      uintptr_t *solshiftf90) {
  /* cast pointers */
  ratparams *rat = (ratparams *) (*ratf90);
  void **solshiftdata = (void **) (*solshiftf90);
  FreeASIGMABSolDirect(rat->num, solshiftdata);
  free(solshiftdata);
}

