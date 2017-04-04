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
