#include <string.h>
#include "internal_header.h"

/**
 * @file evsl_memory.c
 * @brief EVSL memory management routines
 */

void *_evsl_Malloc(size_t nbytes)
{
   void *ptr = malloc(nbytes);
   if (ptr == NULL && nbytes > 0)
   {
      fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes);
      CHKERR(1);
   }
   return ptr;
}

void *_evsl_Calloc(size_t count, size_t nbytes)
{
   void *ptr = calloc(count, nbytes);
   if (ptr == NULL)
   {
      fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", count * nbytes);
      CHKERR(1);
   }
   return ptr;
}

void *_evsl_Realloc(void *ptr, size_t nbytes)
{
   ptr = realloc(ptr, nbytes);
   if (ptr == NULL && nbytes > 0)
   {
      fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes);
      CHKERR(1);
   }
   return ptr;
}

void _evsl_Free(void *ptr)
{
   if (!ptr)
   {
      return;
   }

   free(ptr);
}

/* GPU memory management, which call the CPU counterparts if not configured with CUDA 
 * i.e., (device == host) */
void *_evsl_Malloc_device(size_t nbytes)
{
#ifdef EVSL_USING_CUDA_GPU
   void *ptr;
   cudaError_t stat = cudaMalloc(&ptr, nbytes);

   /* fprintf(stdout, "evsl alloc device memory %.2f GB\n", nbytes/1024.0/1024.0/1024.0); */

   if (stat != cudaSuccess || (ptr == NULL && nbytes > 0))
   {
      fprintf(stdout, "EVSL Error: out of device memory [%zu bytes asked]\n", nbytes);
      CHKERR(1);
   }
   return ptr;
#else
   return _evsl_Malloc(nbytes);
#endif
}

void *_evsl_Calloc_device(size_t count, size_t nbytes)
{
#ifdef EVSL_USING_CUDA_GPU
   void *ptr = _evsl_Malloc_device(count * nbytes);
   cudaError_t stat = cudaMemset(ptr, 0, count * nbytes);
   CHKERR(stat != cudaSuccess);
   return ptr;
#else
   return _evsl_Calloc(count, nbytes);
#endif
}

void _evsl_Free_device(void *ptr)
{
#ifdef EVSL_USING_CUDA_GPU
   if (!ptr)
   {
      return;
   }
   cudaError_t stat = cudaFree(ptr);
   CHKERR(stat != cudaSuccess);
#else
   _evsl_Free(ptr);
#endif
}

/**
 * @brief Device memory realloc
 * @warning Note that the `old_size' must also be provied
 */
void *_evsl_Realloc_device(void *old_ptr, size_t old_nbytes, size_t new_nbytes)
{
#ifdef EVSL_USING_CUDA_GPU
   void *new_ptr = _evsl_Malloc_device(new_nbytes);
   cudaError_t stat = cudaMemcpy(new_ptr, old_ptr, evsl_min(old_nbytes, new_nbytes), cudaMemcpyDeviceToDevice);
   CHKERR(stat != cudaSuccess);
   _evsl_Free_device(old_ptr);

   return new_ptr;
#else
   return _evsl_Realloc(old_ptr, new_nbytes);
#endif
}

void evsl_memcpy_device_to_host(void *dst, void *src, size_t nbytes) {
#ifdef EVSL_USING_CUDA_GPU
   cudaError_t stat = cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToHost);
   CHKERR(stat != cudaSuccess);
#else
   memcpy(dst, src, nbytes);
#endif
}

void evsl_memcpy_device_to_device(void *dst, void *src, size_t nbytes) {
#ifdef EVSL_USING_CUDA_GPU
   cudaError_t stat = cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToDevice);
   CHKERR(stat != cudaSuccess);
#else
   memcpy(dst, src, nbytes);
#endif
}

void evsl_memcpy_host_to_device(void *dst, void *src, size_t nbytes) {
#ifdef EVSL_USING_CUDA_GPU
   cudaError_t stat = cudaMemcpy(dst, src, nbytes, cudaMemcpyHostToDevice);
   CHKERR(stat != cudaSuccess);
#else
   memcpy(dst, src, nbytes);
#endif
}

void evsl_memset_device(void * ptr, int value, size_t nbytes) {
#ifdef EVSL_USING_CUDA_GPU
   cudaError_t stat = cudaMemset(ptr, value, nbytes);
   CHKERR(stat != cudaSuccess);
#else
   memset(ptr, value, nbytes);
#endif
}
