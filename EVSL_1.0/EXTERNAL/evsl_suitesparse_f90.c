#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"
#include "evsl_suitesparse.h"

void setup_bsol_suitesparse_f90_(size_t *Bsolf90) {
  BSolDataSuiteSparse *Bsol;
  Malloc(Bsol, 1, BSolDataSuiteSparse);
  SetupBSolSuiteSparse(evsldata.B, Bsol);
  
  SetBSol(BSolSuiteSparse, (void *) Bsol);
  SetLTSol(LTSolSuiteSparse); 
  
  *Bsolf90 = (size_t) Bsol;
}

void free_bsol_suitesparse_f90_(size_t *Bsolf90) {
  BSolDataSuiteSparse *Bsol = (BSolDataSuiteSparse *) (*Bsolf90);
  FreeBSolSuiteSparseData(Bsol);
  free(Bsol);
}

void setup_asigmabsol_suitesparse_f90_(size_t *ratf90, size_t *solshiftf90) {
  ratparams *rat = (ratparams *) (*ratf90);
  
  void **solshiftdata = (void **) malloc(rat->num*sizeof(void *));
  SetupASIGMABSolSuiteSparse(evsldata.A, evsldata.B, rat->num, rat->zk, 
                             solshiftdata);
  *solshiftf90 = (size_t) solshiftdata;
}

void evsl_setasigmabsol_f90_(size_t *ratf90, size_t *solshiftf90) {
  ratparams *rat = (ratparams *) (*ratf90);
  void **solshiftdata = (void **) (*solshiftf90);
  SetASigmaBSol(rat, NULL, ASIGMABSolSuiteSparse, solshiftdata);
}

void free_asigmabsol_suitesparse_f90_(size_t *ratf90, size_t *solshiftf90) {
  ratparams *rat = (ratparams *) (*ratf90);
  void **solshiftdata = (void **) (*solshiftf90);
  FreeASIGMABSolSuiteSparse(rat->num, solshiftdata);
  free(solshiftdata);
}

