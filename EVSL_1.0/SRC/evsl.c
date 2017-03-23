#include <stdlib.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/** \brief global variable of EVSL
 *
 * global variable is guaranteed to be initialized
 * */
evslData evsldata;

/**
 * @brief Initialize evslData 
 *
 * */
int EVSLStart() {
  evsldata.ifGenEv = 0;
  evsldata.A = NULL;
  evsldata.B = NULL;
  evsldata.Amv = NULL;
  evsldata.Bmv = NULL;
  evsldata.Bsol = NULL;
  evsldata.LTsol = NULL;
  return 0;
}

/**
 * @brief Finish EVSL.
 *
 * */
int EVSLFinish() {
  if (evsldata.Amv) {
    free(evsldata.Amv);
  }
  if (evsldata.Bmv) {
    free(evsldata.Bmv);
  }
  if (evsldata.Bsol) {
    free(evsldata.Bsol);
  }
  if (evsldata.LTsol) {
    free(evsldata.LTsol);
  }
  return 0;
}

/** 
 * @brief Set the matrix A
 * 
 * */
int SetAMatrix(csrMat *A) {
  evsldata.A = A;
  
  return 0;
}

/**
 * @brief Set the B matrix.
 * 
 * */
int SetBMatrix(csrMat *B) {
  evsldata.B = B;

  return 0;
}

/**
 * @brief Set the user-input matvec routine and the associated data for A.
 * Save them in evsldata
 * @warning Once this matvec func is set, matrix A will be ignored even it
 * is provided
 * */
int SetAMatvec(int n, MVFunc func, void *data) {
  if (!evsldata.Amv) {
    Calloc(evsldata.Amv, 1, externalMatvec);
  }
  evsldata.Amv->n = n;
  evsldata.Amv->func = func;
  evsldata.Amv->data = data;
  
  return 0;
}

/**
 * @brief Set the user-input matvec routine and the associated data for B.
 * Save them in evsldata
 * @warning Once this matvec func is set, matrix B will be ignored even it
 * is provided
 * */
int SetBMatvec(int n, MVFunc func, void *data) {
  if (!evsldata.Bmv) {
    Calloc(evsldata.Bmv, 1, externalMatvec);
  }
  evsldata.Bmv->n = n;
  evsldata.Bmv->func = func;
  evsldata.Bmv->data = data;
  
  return 0;
}


/**
 * @brief Unset the matvec routine and the associated data for A
 * */
int UnsetAMatvec() {
  if (evsldata.Amv) {
    free(evsldata.Amv);
  }
  evsldata.Amv = NULL;
  
  return 0;
}

/**
 * @brief Unset the matvec routine and the associated data for B
 * */
int UnsetBMatvec() {
  if (evsldata.Bmv) {
    free(evsldata.Bmv);
  }
  evsldata.Bmv = NULL;
  
  return 0;
}

/**
 * @brief Set the solve routine and the associated data for B
 * */
int SetBSol(SolFuncR func, void *data) {
  if (!evsldata.Bsol) {
    Calloc(evsldata.Bsol, 1, BSolType);
  }

  evsldata.Bsol->func = func;
  evsldata.Bsol->data = data;

  return 0;
}


/**
 * @brief Set the problem to standard eigenvalue problem
 * 
 * */
int SetStdEig() {
  evsldata.ifGenEv = 0;

  return 0;
}

/**
 * @brief Set the problem to generalized eigenvalue problem
 * 
 * */
int SetGenEig() {
  evsldata.ifGenEv = 1;

  return 0;
}

/**
 * @brief Set the solve routine and the associated data for A-SIGMA*B
 * if func == NULL, set all functions to be allf
 * */
int SetASigmaBSol(ratparams *rat, SolFuncC *func, SolFuncC allf, void **data) {
  int i,num;
  num = rat->num;
  Calloc(rat->ASIGBsol, num, ASIGMABSolType);

  for (i=0; i<num; i++) {
    rat->ASIGBsol[i].func = func ? func[i] : allf;
    rat->ASIGBsol[i].data = data[i];
  }

  return 0;
}


/**
 * @brief Set the solve routine for LT
 * */
int SetLTSol(SolFuncR func) {
  if (!evsldata.LTsol) {
    Calloc(evsldata.LTsol, 1, LTSolType);
  }
  evsldata.LTsol->func = func;
  return 0;
}

