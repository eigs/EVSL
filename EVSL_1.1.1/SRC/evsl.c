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
 * @brief Set the GPU HYB matrix A
 * @param[in] A The matrix to set
 * */
int SetAMatrix_device(hybMat *A) {
  evsldata.n = A->ncols;
  if (!evsldata.Amv) {
     evsldata.Amv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Amv->func = matvec_cusparse;
  evsldata.Amv->data = (void *) A;

  return 0;
}

/**
 * @brief Set the GPU HYB matrix B
 * @param[in] B The matrix to set
 * */
int SetBMatrix_device(hybMat *B) {
  evsldata.n = B->ncols;
  if (!evsldata.Bmv) {
     evsldata.Bmv = evsl_Calloc(1, EVSLMatvec);
  }
  evsldata.Bmv->func = matvec_cusparse;
  evsldata.Bmv->data = (void *) B;

  return 0;
}

/**
 * @brief create HYB matrix A on GPU
 *
 * @param[in] A The CSR matrix to set [on host].
 * A cusparse HYB matrix will be generated
 * */
int evsl_CreateHybMat(csrMat *A, hybMat *Ahyb) {

  Ahyb->nrows = A->nrows;
  Ahyb->ncols = A->ncols;

  cusparseStatus_t cusparseStat = cusparseCreateMatDescr(&Ahyb->descr);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  cusparseStat = cusparseCreateHybMat(&Ahyb->hyb);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  cusparseSetMatIndexBase(Ahyb->descr, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatType(Ahyb->descr, CUSPARSE_MATRIX_TYPE_GENERAL);

  /* CSR on GPU */
  csrMat Agpu;
  evsl_copy_csr_to_gpu(A, &Agpu);

  /* convert to hyb */
  cusparseStat = cusparseDcsr2hyb(evsldata.cusparseH, A->nrows, A->ncols,
                                  Ahyb->descr, Agpu.a, Agpu.ia, Agpu.ja,
                                  Ahyb->hyb, -1, CUSPARSE_HYB_PARTITION_AUTO);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  evsl_free_csr_gpu(&Agpu);

  return 0;
}

/**
 * @brief Query CUDA device and set device
 *
 * @param[in] set_dev: the device number to set
 * */
void evsl_device_query(int set_dev) {
  int deviceCount, dev;
  cudaGetDeviceCount(&deviceCount);
  printf("=========================================\n");
  if (deviceCount == 0) {
    printf("There is no device supporting CUDA\n");
  }

  for (dev = 0; dev < deviceCount; ++dev) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if (dev == 0) {
      if (deviceProp.major == 9999 && deviceProp.minor == 9999)
        printf("There is no device supporting CUDA.\n");
      else if (deviceCount == 1)
        printf("There is 1 device supporting CUDA\n");
      else
        printf("There are %d devices supporting CUDA\n", deviceCount);
    }
    printf("\nDevice %d: \"%s\"\n", dev, deviceProp.name);
    printf("  Major revision number:          %d\n",
           deviceProp.major);
    printf("  Minor revision number:          %d\n",
           deviceProp.minor);
    printf("  Total amount of global memory:  %.2f GB\n",
           deviceProp.totalGlobalMem/1e9);
  }

  dev = set_dev;
  CHKERR(dev < 0 || dev >= deviceCount);
  cudaSetDevice(dev);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  printf("\nRunning on Device %d: \"%s\"\n", dev, deviceProp.name);
  printf("=========================================\n");
}

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

