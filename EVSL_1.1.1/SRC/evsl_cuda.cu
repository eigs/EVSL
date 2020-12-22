#include "internal_header.h"

/**
 * @file cuda_util.cu
 * @brief cuda functions
 */

#ifdef EVSL_USING_CUDA_GPU
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

/* vkp1[i] = t*(vkp1[i]-cc*vk[i]) - vkm1[i]; y[i] += s*vkp1[i]; */
__global__ void evsl_pnav_kernel(int n, int k, double t, double s, double cc, double *vkp1, double *v_cur,
                                 double *v_old, double *y) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      double zi = k > 1 ? v_old[tid] : 0.0;
      double wi = t*(vkp1[tid] - cc*v_cur[tid]) - zi;
      vkp1[tid] = wi;
      y[tid] += s * wi;
   }
}

void evsl_pnav_device(int n, int k, double t, double s, double cc, double *vkp1, double *v_cur,
                      double *v_old, double *y) {
  const int bDim = 512;
  int gDim = (n + bDim - 1) / bDim;
  evsl_pnav_kernel<<<gDim, bDim>>>(n, k, t, s, cc, vkp1, v_cur, v_old, y);
}

/* https://docs.nvidia.com/cuda/curand/group__HOST.html#group__HOST*/

#define CURAND_METHOD 3 /* 1: generate 32-bit quasirandom numbers, and then convert to double
                           2: generate uniformly distributed doubles directly with cuRand
                           3: generate on host and copy to device */

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

/**
 * Generates a uniform distributed random vector of length n, [-1, 1]
 * on device
 */
void rand_double_device(int n, double *v) {
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
#if CURAND_METHOD == 2
  curandStatus_t curandStatus = curandGenerateNormalDouble(evsldata.curandGen, v, n, 0.0, 1.0);
  CHKERR(CURAND_STATUS_SUCCESS != curandStatus);
#else
  double *h_v = evsl_Malloc(n, double);
  randn_double(n, h_v);
  evsl_memcpy_host_to_device(v, h_v, n*sizeof(double));
  evsl_Free(h_v);
#endif
}

__global__ void evsl_element_mult_kernel(int n, double *a, double *b) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      b[tid] *= a[tid];
   }
}

/* b := a .* b */
void evsl_element_mult_device(int n, double *a, double *b) {
  int bDim = 512;
  int gDim = (n + bDim - 1) / bDim;
  evsl_element_mult_kernel<<<gDim, bDim>>>(n, a, b);
}

__global__ void evsl_element_divide_kernel(int n, double *a, double *b) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      b[tid] /= a[tid];
   }
}

/* b := a ./ b */
void evsl_element_divide_device(int n, double *a, double *b) {
  int bDim = 512;
  int gDim = (n + bDim - 1) / bDim;
  evsl_element_divide_kernel<<<gDim, bDim>>>(n, a, b);
}

/* vkp1[i] = t*(vkp1[i]-cc*vk[i]) - vkm1[i]; y[i] += s*vkp1[i]; */
__global__ void evsl_chebAv_kernel(int n, int k, double t, double s, double cc, double *vkp1, double *vk,
                                   double *vkm1, double *y) {
   int tid = threadIdx.x + blockIdx.x * blockDim.x;
   if (tid < n) {
      double zi = k > 1 ? vkm1[tid] : 0.0;
      double wi = t*(vkp1[tid] - cc*vk[tid]) - zi;
      vkp1[tid] = wi;
      y[tid] += s * wi;
   }
}

void evsl_chebAv_device(int n, int k, double t, double s, double cc, double *vkp1, double *vk,
                        double *vkm1, double *y) {
    const int bDim = 512;
    int gDim = (n + bDim - 1) / bDim;
    evsl_chebAv_kernel<<<gDim, bDim>>>(n, k, t, s, cc, vkp1, vk, vkm1, y);
}
#endif
