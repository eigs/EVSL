#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"
#include "evsl_suitesparse.h"

/** @brief Fortran interface for SetBSol and SetLTSol 
 * @warning Will use evsldata.B, so must be used after setting B */
void EVSLFORT(setup_bsol_suitesparse)(uintptr_t *Bsolf90) {
  BSolDataSuiteSparse *Bsol;
  Malloc(Bsol, 1, BSolDataSuiteSparse);
  SetupBSolSuiteSparse(evsldata.B, Bsol);
  
  SetBSol(BSolSuiteSparse, (void *) Bsol);
  SetLTSol(LTSolSuiteSparse, (void *) Bsol);
  
  *Bsolf90 = (uintptr_t) Bsol;
}

/** @brief Fortran interface for FreeBSolSuiteSparseData */
void EVSLFORT(free_bsol_suitesparse)(uintptr_t *Bsolf90) {
  /* cast pointer */
  BSolDataSuiteSparse *Bsol = (BSolDataSuiteSparse *) (*Bsolf90);
  FreeBSolSuiteSparseData(Bsol);
  free(Bsol);
}

/** @brief Fortran interface for SetupASIGMABSolSuiteSparse 
 * @warning Will use evsldata.A and evsldata.B, 
 * so must be used after setting A and B 
 * @param[in] ratf90: pointer of rational filter
 * @param[out] solshiftf90: pointer of solshift array
 */
void EVSLFORT(setup_asigmabsol_suitesparse)(uintptr_t *ratf90, 
                                            uintptr_t *solshiftf90) {
  /* cast pointer */
  ratparams *rat = (ratparams *) (*ratf90);
  /* allocate and setup solshiftdata */
  void **solshiftdata = (void **) malloc(rat->num*sizeof(void *));
  SetupASIGMABSolSuiteSparse(evsldata.A, evsldata.B, rat->num, rat->zk, 
                             solshiftdata);
  /* cast pointer in output */
  *solshiftf90 = (uintptr_t) solshiftdata;
}

/** @brief Fortran interface for SetASigmaBSol
 * @param[in/out] ratf90: pointer of rational filter
 * @param[out] solshiftf90: pointer of solshift array
 */
void EVSLFORT(evsl_setasigmabsol)(uintptr_t *ratf90, 
                                  uintptr_t *solshiftf90) {
  /* cast pointers */
  ratparams *rat = (ratparams *) (*ratf90);
  void **solshiftdata = (void **) (*solshiftf90);

  SetASigmaBSol(rat, NULL, ASIGMABSolSuiteSparse, solshiftdata);
}

/** @brief Fortran interface for FreeASIGMABSolSuiteSparse
 * @param[in/out] ratf90: pointer of rational filter
 * @param[in/out] solshiftf90: pointer of solshift array
 */
void EVSLFORT(free_asigmabsol_suitesparse)(uintptr_t *ratf90, 
                                           uintptr_t *solshiftf90) {
  /* cast pointers */
  ratparams *rat = (ratparams *) (*ratf90);
  void **solshiftdata = (void **) (*solshiftf90);
  FreeASIGMABSolSuiteSparse(rat->num, solshiftdata);
  free(solshiftdata);
}

