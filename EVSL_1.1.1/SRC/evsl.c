#include "internal_header.h"

/**
 * @file evsl.c
 * @brief Set EVSL solver options and data
 */

/** \brief global variable of EVSL
 *
 * global variable is guaranteed to be initialized
 * */
evslData evsldata;

/** \brief global statistics of EVSL
 *
 * global variable is guaranteed to be initialized
 * */
evslStat evslstat;

/**
 * @brief Initialize evslData
 *
 * */
int EVSLStart() {
  evsldata.n = 0;
  evsldata.ifGenEv = 0;
  evsldata.Amv = NULL;
  evsldata.Bmv = NULL;
  evsldata.Bsol = NULL;
  evsldata.LTsol = NULL;
  evsldata.ds = NULL;

  StatsReset();

  /* do a dummy timer call to improve accuracy of timing */
  evsl_timer();

#ifdef EVSL_USING_CUDA_GPU
  cusparseStatus_t cusparseStat = cusparseCreate(&evsldata.cusparseH);
  CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

  cublasStatus_t cublasStat = cublasCreate(&evsldata.cublasH);
  CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);

  curandStatus_t curandStatus = curandCreateGenerator(&evsldata.curandGen,
                                                      CURAND_RNG_PSEUDO_DEFAULT);
  CHKERR(CURAND_STATUS_SUCCESS != curandStatus);

  //curandSetPseudoRandomGeneratorSeed(evsldata.randGen ,1234ULL);
#endif

  return 0;
}

/**
 * @brief Finish EVSL.
 *
 * Frees parts of the evsldata struct
 * */
int EVSLFinish() {
#ifdef EVSL_USING_CUDA_GPU
  cusparseStatus_t cusparseStat = cusparseDestroy(evsldata.cusparseH);
  CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

  cublasStatus_t cublasStat = cublasDestroy(evsldata.cublasH);
  CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);

  curandStatus_t curandStatus = curandDestroyGenerator(evsldata.curandGen);
  CHKERR(CURAND_STATUS_SUCCESS != curandStatus);
#endif
  if (evsldata.Amv) {
    evsl_Free(evsldata.Amv);
  }
  if (evsldata.Bmv) {
    evsl_Free(evsldata.Bmv);
  }
  if (evsldata.Bsol) {
    evsl_Free(evsldata.Bsol);
  }
  if (evsldata.LTsol) {
    evsl_Free(evsldata.LTsol);
  }
#ifdef EVSL_USING_CUDA_GPU
  evsl_last_device_err();
#endif
  return 0;
}

/**
 * @brief Set the CSR matrix A
 *
 *
 * @param[in] A The matrix to set
 * */
int SetAMatrix(csrMat *A) {
  evsldata.n = A->ncols;
  if (!evsldata.Amv) {
     evsldata.Amv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Amv->func = matvec_csr;
  evsldata.Amv->data = (void *) A;

  return 0;
}

/**
 * @brief Set the B matrix.
 *
 * @param[in] B the matrix to set
 * */
int SetBMatrix(csrMat *B) {
  evsldata.n = B->ncols;
  if (!evsldata.Bmv) {
     evsldata.Bmv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Bmv->func = matvec_csr;
  evsldata.Bmv->data = (void *) B;

  return 0;
}

/**
 * @brief Set the user-input matvec routine and the associated data for A.
 * Save them in evsldata
 * @warning Once this matvec func is set, matrix A will be ignored even it
 * is provided
 *
 * @param[in] n Size of problem
 * @param[in] func Function to use for matvec
 * @param[in] data Data required
 * */
int SetAMatvec(int n, MVFunc func, void *data) {
  evsldata.n = n;
  if (!evsldata.Amv) {
    evsldata.Amv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Amv->func = func;
  evsldata.Amv->data = data;

  return 0;
}

/**
 * @brief Set the user-input matvec routine and the associated data for B.
 * Save them in evsldata
 * @warning Once this matvec func is set, matrix B will be ignored even it
 * is provided
 * @param[in] n Size of problem
 * @param[in] func Function to use for matvec
 * @param[in] data Data required for matvec
 * */
int SetBMatvec(int n, MVFunc func, void *data) {
  evsldata.n = n;
  if (!evsldata.Bmv) {
    evsldata.Bmv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Bmv->func = func;
  evsldata.Bmv->data = data;

  return 0;
}


/**
 * @brief Set the solve routine and the associated data for B
 * @param[in] func Function to use for solve
 * @param[in] data Data for solnve
 * */
int SetBSol(SolFuncR func, void *data) {
  if (!evsldata.Bsol) {
    evsldata.Bsol = evsl_Calloc(1, EVSLBSol);
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
 *
 * @param[in] rat Rational parameters
 * @param[in] func function array
 * @param[in] allf Function to be used for all solves
 * @param[in] data data array
 *
 * */
int SetASigmaBSol(ratparams *rat, int i, SolFuncC func, void *data) {
  int num = rat->num;

  if (rat->ASIGBsol == NULL)
  {
    rat->ASIGBsol = evsl_Calloc(num, EVSLASIGMABSol);
  }

  CHKERR(i < 0 || i >= num);

  rat->ASIGBsol[i].func = func;
  rat->ASIGBsol[i].data = data;

  return 0;
}


/**
 * @brief Set the solve routine for LT
 *
 * @param[in] func Function to use
 * @param[in] data Data to use
 * */
int SetLTSol(SolFuncR func, void *data) {
  if (!evsldata.LTsol) {
    evsldata.LTsol = evsl_Calloc(1, EVSLLTSol);
  }
  evsldata.LTsol->func = func;
  evsldata.LTsol->data = data;
  return 0;
}

/**
 * @brief Set diagonal scaling matrix D
 * @param[in] ds Sets diagonal scaling
 * */
void SetDiagScal(double *ds) {
  evsldata.ds = ds;
}

#ifdef EVSL_USING_CUDA_GPU
/**
 * @brief Set the GPU CSR matrix A
 * @param[in] A The matrix to set
 * */
int SetAMatrix_device_csr(csrMat *A) {
  evsldata.n = A->ncols;
  if (!evsldata.Amv) {
     evsldata.Amv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Amv->func = matvec_cusparse_csr;
  evsldata.Amv->data = (void *) A;

  return 0;
}

/**
 * @brief Set the GPU CSR matrix B
 * @param[in] B The matrix to set
 * */
int SetBMatrix_device_csr(csrMat *B) {
  evsldata.n = B->ncols;
  if (!evsldata.Bmv) {
     evsldata.Bmv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Bmv->func = matvec_cusparse_csr;
  evsldata.Bmv->data = (void *) B;

  return 0;
}

#ifdef EVSL_USING_CUSPARSE_HYB
/**
 * @brief Set the GPU HYB matrix A
 * @param[in] A The matrix to set
 * */
int SetAMatrix_device_hyb(hybMat *A) {
  evsldata.n = A->ncols;
  if (!evsldata.Amv) {
     evsldata.Amv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Amv->func = matvec_cusparse_hyb;
  evsldata.Amv->data = (void *) A;

  return 0;
}

/**
 * @brief Set the GPU HYB matrix B
 * @param[in] B The matrix to set
 * */
int SetBMatrix_device_hyb(hybMat *B) {
  evsldata.n = B->ncols;
  if (!evsldata.Bmv) {
     evsldata.Bmv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Bmv->func = matvec_cusparse_hyb;
  evsldata.Bmv->data = (void *) B;

  return 0;
}
#endif

/**
 * @brief Get the last CUDA device error
 * */
void evsl_last_device_err() {
  cudaError_t cudaerr = cudaGetLastError();
  if (cudaerr != cudaSuccess) {
    printf(" EVSL CUDA error: %s\n",cudaGetErrorString(cudaerr));
    CHKERR(1);
  }
}
#endif

