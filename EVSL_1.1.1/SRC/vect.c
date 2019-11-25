#include <stdlib.h>
#include <string.h>
#include "internal_header.h"
/**
 * @file vect.c
 * @brief Vector operations
 * @param[in] n Size of vector
 * @param[out] n vector
 */

/**
 * Generates a uniform distributed random vector of length n, [-1, 1]
 */
void rand_double(int n, double *v) {
  int i;
  double t = ((double) RAND_MAX)/2.0;
  for (i=0; i<n; i++) {
    v[i] = (rand() -t)/t;
  }
}

/**
 * Generates a normally distributed random vector of length n
 *
 * Uses the Box-Muller transformation
 * @param[out] v Vector to be populated
 * @param[in] n Number of elements to generate (should be the length of v)
 * */
void randn_double(int n, double *v) {
  const double two_pi = 2.0 * 3.1415926535;
  int i;
  for(i = 0; i < n; i++) {
    static double Z0;
    static double Z1;
    static int regen = 0;//A boolean
    regen = !regen;
    if(!regen)
    {
      v[i] = Z1;
    }

    double U1 = rand() * (1.0 / RAND_MAX);
    double U2 = rand() * (1.0 / RAND_MAX);

    Z0 = sqrt(-2.0 * log(U1)) * cos(two_pi  * U2);
    Z1 = sqrt(-2.0 * log(U1)) * sin(two_pi  * U2);
    v[i] = Z0;
  }
}

/* https://docs.nvidia.com/cuda/curand/group__HOST.html#group__HOST*/

#define CURAND_METHOD 3 /* 1: generate 32-bit quasirandom numbers, and then convert to double
                           2: generate uniformly distributed doubles directly with cuRand
                           3: generate on host and copy to device */

#ifdef EVSL_USING_CUDA_GPU
#if CURAND_METHOD == 1
__global__ void evsl_rand_double_cudakernel(int n, unsigned int *ui, double *d) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      const double t = ((double) 0xFFFFFFFF) / 2.0;
      d[tid] = (ui[tid] - t) / t;
   }
}
#else
__global__ void evsl_rand_double_cudakernel(int n, double *d) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      d[tid] = d[tid] * 2.0 - 1.0;
   }
}
#endif
#endif

/**
 * Generates a uniform distributed random vector of length n, [-1, 1]
 * on device
 */
void rand_double_device(int n, double *v) {
#ifdef EVSL_USING_CUDA_GPU
#if CURAND_METHOD == 1 || CURAND_METHOD == 2
  int bDim = 512;
  int gDim = (n + bDim - 1) / bDim;
#endif
#if CURAND_METHOD == 1
  unsigned int *d_uint = evsl_Malloc_device(n, unsigned int);
  curandStatus_t curandStatus = curandGenerate(evsldata.curandGen, d_uint, n);
  CHKERR(CURAND_STATUS_SUCCESS != curandStatus);
  evsl_rand_double_cudakernel<<<gDim, bDim>>>(n, d_uint, v);
  evsl_Free_device(d_uint);
#elif CURAND_METHOD == 2
  curandStatus_t curandStatus = curandGenerateUniformDouble(evsldata.curandGen, v, n);
  CHKERR(CURAND_STATUS_SUCCESS != curandStatus);
  evsl_rand_double_cudakernel<<<gDim, bDim>>>(n, v);
#else
  double *h_v = evsl_Malloc(n, double);
  rand_double(n, h_v);
  evsl_memcpy_host_to_device(v, h_v, n*sizeof(double));
  evsl_Free(h_v);
#endif
#else
  rand_double(n, v);
#endif
}

/*
  if (CURAND_STATUS_SUCCESS != curandStatus) {
     if (curandStatus == CURAND_STATUS_NOT_INITIALIZED) printf("1 \n");
     if (curandStatus == CURAND_STATUS_PREEXISTING_FAILURE) printf("2 \n");
     if (curandStatus == CURAND_STATUS_LAUNCH_FAILURE) printf("3 \n");
     if (curandStatus == CURAND_STATUS_LENGTH_NOT_MULTIPLE) printf("4 \n");
     if (curandStatus == CURAND_STATUS_DOUBLE_PRECISION_REQUIRED) printf("5 \n");
  }
*/

/**
 * Generates a normally distributed random vector of length n on device
 * on device
 */
void randn_double_device(int n, double *v) {
#ifdef EVSL_USING_CUDA_GPU
#if CURAND_METHOD == 2
  curandStatus_t curandStatus = curandGenerateNormalDouble(evsldata.curandGen, v, n, 0.0, 1.0);
  CHKERR(CURAND_STATUS_SUCCESS != curandStatus);
#else
  double *h_v = evsl_Malloc(n, double);
  randn_double(n, h_v);
  evsl_memcpy_host_to_device(v, h_v, n*sizeof(double));
  evsl_Free(h_v);
#endif
#else
  randn_double(n, v);
#endif
}


/**
 * Sets all elements of v to t
 * @param[in] n Number of elements
 * @param[in] t Value which elements should be set to
 * @param[out] v Vector to set
 * */
void vecset(int n, double t, double *v) {
  int i;
  for (i=0; i<n; i++)
    v[i] = t;
}

/**
 * Creates a vector whose elements are linearly spaced
 * @param[in] a Lower bound
 * @param[in] b Upper bound
 * @param[in] num Number of values
 * @param[out] arr Output vector
 */
void linspace(double a, double b, int num, double *arr){
  double h;
  h = (num==1? 0: (b-a)/(num-1));
  int i;
 //-------------------- careful at the boundaries!
  arr[0] = a;
  arr[num-1] = b;
  for (i=1; i<num-1; i++)
    arr[i] = a+i*h;

}

/**
 * @brief Compares a,b as doubles
 * @param[in] a First value
 * @param[in] b Second value
 * @return -1 if b>a, 0 if a==b, 1 otherwise
 * */
int compare1(const void *a, const void *b) {
  double *aa = (double*) a;
  double *bb = (double*) b;
  if (*aa < *bb) {
    return -1;
  } else if (*aa == *bb) {
    return 0;
  } else {
    return 1;
  }
}
typedef struct _doubleint {
  int i;
  double d;
} doubleint;

/**
 * @brief Compares the doubles of a,b as double/int pairs
 * @param[in] a First value
 * @param[in] b Second value
 * @return -1 if b>a, 0 if a==b, 1 otherwise
 * */
int compare2(const void *a, const void *b) {
  const doubleint *aa = (doubleint*) a;
  const doubleint *bb = (doubleint*) b;
  if (aa->d < bb->d) {
    return -1;
  } else if (aa->d == bb->d) {
    return 0;
  } else {
    return 1;
  }
}
/**
 * @brief Sorts a vector, and potentially indices
 * @param[in] n Number of elements
 * @param[in, out] v Vector to sort
 * @param[in, out] ind Indices to sort
 *
 * */
void sort_double(int n, double *v, int *ind) {
  /* if sorting indices are not wanted */
  if (ind == NULL) {
    qsort(v, n, sizeof(double), compare1);
    return;
  }
  doubleint *vv;
  vv = evsl_Malloc(n, doubleint);
  int i;
  for (i=0; i<n; i++) {
    vv[i].d = v[i];
    vv[i].i = i;
  }
  qsort(vv, n, sizeof(doubleint), compare2);
  for (i=0; i<n; i++) {
    v[i] = vv[i].d;
    ind[i] = vv[i].i;
  }
  evsl_Free(vv);
}

/** @brief y = x(p)
 * @param[in] n Number of points in vector
 * @param[in] p Permutation vector
 * @param[in] x Input vector
 * @param[out] y Output vector
 * */
void vec_perm(int n, int *p, double *x, double *y) {
  if (!p) {
    memcpy(y, x, n*sizeof(double));
  } else {
    int i;
    for (i=0; i<n; i++) {
      y[i] = x[p[i]];
    }
  }
}


/* @brief y(p) = x
 * @param[in] n Number of elements in vector
 * @param[in] p Permutation vector
 * @param[in] x Input vector
 * @param[in] y Output vector */
void vec_iperm(int n, int *p, double *x, double *y) {
  if (!p) {
    memcpy(y, x, n*sizeof(double));
  } else {
    int i;
    for (i=0; i<n; i++) {
      y[p[i]] = x[i];
    }
  }
}

/**
 * BLAS routines that are performed on device (CUBLAS)
 * When EVSL is not configured with CUDA, these are just wrappers to BLAS calls
 */
/*
   if (CUBLAS_STATUS_SUCCESS != cublasStat) {
   if (cublasStat == CUBLAS_STATUS_NOT_INITIALIZED) printf("err 1\n");
   if (cublasStat == CUBLAS_STATUS_ALLOC_FAILED) printf("err 2\n");
   if (cublasStat == CUBLAS_STATUS_INVALID_VALUE) printf("err 3\n");
   if (cublasStat == CUBLAS_STATUS_ARCH_MISMATCH) printf("err 4\n");
   if (cublasStat == CUBLAS_STATUS_MAPPING_ERROR) printf("err 5\n");
   if (cublasStat == CUBLAS_STATUS_EXECUTION_FAILED) printf("err 6\n");
   if (cublasStat == CUBLAS_STATUS_INTERNAL_ERROR) printf("err 7\n");
   if (cublasStat == CUBLAS_STATUS_NOT_SUPPORTED) printf("err 8\n");
   if (cublasStat == CUBLAS_STATUS_LICENSE_ERROR) printf("err 9\n");
   }
*/
double evsl_dnrm2_device(int *n, double *x, int *incr) {
   double t;
#ifdef EVSL_USING_CUDA_GPU
   cublasStatus_t cublasStat = cublasDnrm2(evsldata.cublasH, *n, x, *incr, &t);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   t = evsl_dnrm2(n, x, incr);
#endif
   return t;
}

void evsl_dscal_device(int *n, double *a, double *x, int *incr) {
#ifdef EVSL_USING_CUDA_GPU
   cublasStatus_t cublasStat = cublasDscal(evsldata.cublasH, *n, a, x, *incr);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   evsl_dscal(n, a, x, incr);
#endif
}

double evsl_ddot_device(int *n, double *x, int *incx, double *y, int *incy) {
   double t;
#ifdef EVSL_USING_CUDA_GPU
   cublasStatus_t cublasStat = cublasDdot(evsldata.cublasH, *n, x, *incx, y, *incy, &t);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   t = evsl_ddot(n, x, incx, y, incy);
#endif
   return t;
}

void evsl_daxpy_device(int *n, double *a, double *x, int *incx, double *y, int *incy) {
#ifdef EVSL_USING_CUDA_GPU
   cublasStatus_t cublasStat = cublasDaxpy(evsldata.cublasH, *n, a, x, *incx, y, *incy);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   evsl_daxpy(n, a, x, incx, y, incy);
#endif
}

void evsl_dcopy_device(int *n, double *x, int *incx, double *y, int *incy) {
#ifdef EVSL_USING_CUDA_GPU
   cublasStatus_t cublasStat = cublasDcopy(evsldata.cublasH, *n, x, *incx, y, *incy);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   evsl_dcopy(n, x, incx, y, incy);
#endif
}

void evsl_dgemv_device(const char *trans, int *m, int *n, double *alpha, double *a, int *lda,
                       double *x, int *incx, double *beta, double *y, int *incy) {
#ifdef EVSL_USING_CUDA_GPU
   cublasOperation_t cublas_trans = (*trans == 'T' || *trans == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;
   cublasStatus_t cublasStat = cublasDgemv(evsldata.cublasH, cublas_trans, *m, *n,
                                           alpha, a, *lda, x, *incx, beta, y, *incy);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   evsl_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif
}

void evsl_dgemm_device(const char *transa, const char *transb, int *m, int * n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc)
{
#ifdef EVSL_USING_CUDA_GPU
   cublasOperation_t cublas_transa = (*transa == 'T' || *transa == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;
   cublasOperation_t cublas_transb = (*transb == 'T' || *transb == 't') ? CUBLAS_OP_T : CUBLAS_OP_N;
   cublasStatus_t cublasStat = cublasDgemm(evsldata.cublasH, cublas_transa, cublas_transb, *m, *n, *k,
                                           alpha, a, *lda, b, *ldb, beta, c, *ldc);
   CHKERR(CUBLAS_STATUS_SUCCESS != cublasStat);
#else
   evsl_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
#endif
}

#ifdef EVSL_USING_CUDA_GPU
__global__ void evsl_element_mult_kernel(int n, double *a, double *b) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      b[tid] *= a[tid];
   }
}
#endif

/* b := a .* b */
void evsl_element_mult_device(int n, double *a, double *b) {
#ifdef EVSL_USING_CUDA_GPU
  int bDim = 512;
  int gDim = (n + bDim - 1) / bDim;
  evsl_element_mult_kernel<<<gDim, bDim>>>(n, a, b);
#else
  int i;
  for (i = 0; i < n; i++) {
     b[i] *= a[i];
  }
#endif
}

#ifdef EVSL_USING_CUDA_GPU
__global__ void evsl_element_divide_kernel(int n, double *a, double *b) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      b[tid] /= a[tid];
   }
}
#endif

/* b := a ./ b */
void evsl_element_divide_device(int n, double *a, double *b) {
#ifdef EVSL_USING_CUDA_GPU
  int bDim = 512;
  int gDim = (n + bDim - 1) / bDim;
  evsl_element_divide_kernel<<<gDim, bDim>>>(n, a, b);
#else
  int i;
  for (i = 0; i < n; i++) {
     b[i] /= a[i];
  }
#endif
}

