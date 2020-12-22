#include "lapl.h"

#ifdef EVSL_USING_CUDA_GPU
/* y = A * x
 * MatVec routine provided by users.
 * Example: 7-pt stencil
 * data: pointer to a struct that contains all needed data
 */

__global__ void
Lap2D3DMatvecCUDA(int nx, int ny, int nz, double *x, double *y,
                  double s0, double s1, double s2, double s3, double s4,
                  double s5, double s6) {
   int n, i, j, k, p;
   n = nx * ny * nz;
   p = blockIdx.x * blockDim.x + threadIdx.x;
   if (p < n) {
      i = p;
      k = i / (nx*ny); i -= k*nx*ny;
      j = i / nx; i -= j * nx;

      y[p] = s0 * x[p];
      // x-1, x+1
      if (i>0)    { y[p] += s1 * x[p-1]; }
      if (i<nx-1) { y[p] += s2 * x[p+1]; }
      // y-1, y+1
      if (j>0)    { y[p] += s3 * x[p-nx]; }
      if (j<ny-1) { y[p] += s4 * x[p+nx]; }
      // z-1, z+1
      if (k>0)    { y[p] += s5 * x[p-nx*ny]; }
      if (k<nz-1) { y[p] += s6 * x[p+nx*ny]; }
   }
}

void Lap2D3DMatvec(double *x, double *y, void *data) {
  lapmv_t *lapmv = (lapmv_t *) data;
  int nx = lapmv->nx;
  int ny = lapmv->ny;
  int nz = lapmv->nz;
  int n = nx * ny * nz;
  double *stencil = lapmv->stencil;
  int bdim = 1024;
  int gdim = (n + bdim - 1) / bdim;
  Lap2D3DMatvecCUDA<<<gdim, bdim>>>(nx, ny, nz, x, y, stencil[0], stencil[1],
        stencil[2], stencil[3], stencil[4], stencil[5], stencil[6]);
}
#endif

