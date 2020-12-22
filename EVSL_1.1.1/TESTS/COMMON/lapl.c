#include "lapl.h"

/**-----------------------------------------------------------------------
 *
 * @brief Laplacean Matrix generator
 *
 * @param[in] nx  Number of points in x-direction
 * @param[in] ny  Number of points in y-direction
 * @param[in] nz  Number of points in z-direction
 * @param[out] *Acoo matrix in coordinate format.
 *
 -----------------------------------------------------------------------**/
int lapgen(int nx, int ny, int nz, cooMat *Acoo) {
  int n = nx * ny * nz;
  Acoo->nrows = n;
  Acoo->ncols = n;

  int nzmax = nz > 1 ? 7*n : 5*n;
  Acoo->ir = evsl_Malloc(nzmax, int);
  Acoo->jc = evsl_Malloc(nzmax, int);
  Acoo->vv = evsl_Malloc(nzmax, double);

  int ii, nnz=0;
  for (ii=0; ii<n; ii++) {
    double v = -1.0;
    int i,j,k,jj;
    k = ii / (nx*ny);
    i = (ii - k*nx*ny) / nx;
    j = ii - k*nx*ny - i*nx;

    if (k > 0) {
      jj = ii - nx * ny;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }
    if (k < nz-1) {
      jj = ii + nx * ny;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }

    if (i > 0) {
      jj = ii - nx;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }
    if (i < ny-1) {
      jj = ii + nx;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }

    if (j > 0) {
      jj = ii - 1;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }
    if (j < nx-1) {
      jj = ii + 1;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }

    v = nz > 1 ? 6.0 : 4.0;
    Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = ii;  Acoo->vv[nnz] = v;  nnz++;
  }

  Acoo->nnz = nnz;

  return 0;
}

/**-----------------------------------------------------------------------
 *
 * @brief Exact eigenvalues of Laplacean in interval [a b]
 * @param[in] nx  Number of points in x-direction
 * @param[in] ny  Number of points in y-direction
 * @param[in] nz  Number of points in z-direction
 * @param[in] a  Left bound
 * @param[in] b  Right bound
 * @param[out] m number of eigenvalues found
 * @param[out] **vo pointer to array of eigenvalues found
 *
 *-----------------------------------------------------------------------**/
int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo) {
  double thetax = 0.5 * EVSL_PI / (nx + 1.0);
  double thetay = 0.5 * EVSL_PI / (ny + 1.0);
  double thetaz = 0.5 * EVSL_PI / (nz + 1.0);
  int i, j, k, l=0, n=nx*ny*nz;
  double *v;
  double tx, ty, tz, ev;
  v = evsl_Malloc(n, double);
  for (i=1; i<=nx; i++) {
    tx = sin(i*thetax);
    for (j=1; j<=ny; j++) {
      ty = sin(j*thetay);
      for (k=1; k<=nz; k++) {
        tz = sin(k*thetaz);
        if (1 == nz) {
          tz = 0.0;
        }
        ev = 4*(tx*tx+ty*ty+tz*tz);
        if (ev >= a - EVSL_DBL_EPS_MULT * DBL_EPSILON && ev <= b + EVSL_DBL_EPS_MULT * DBL_EPSILON) {
          v[l++] = ev;
        }
      }
    }
  }
  v = evsl_Realloc(v, l, double);
  sort_double(l, v, NULL);

  *m = l;
  *vo = v;
  return 0;
}

#ifndef EVSL_USING_CUDA_GPU
/**-----------------------------------------------------------------------
 *
 * @brief Matvec with Laplacian in the form of constant stencils
 */
void Lap2D3DMatvec(double *x, double *y, void *data) {
  lapmv_t *lapmv = (lapmv_t *) data;
  int nx = lapmv->nx;
  int ny = lapmv->ny;
  int nz = lapmv->nz;
  double *stencil = lapmv->stencil;
  int i,j,k,p;

  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
        p = k*nx*ny + j*nx + i;
        y[p] = stencil[0] * x[p];
        // x-1, x+1
        if (i>0)    { y[p] += stencil[1] * x[p-1]; }
        if (i<nx-1) { y[p] += stencil[2] * x[p+1]; }
        // y-1, y+1
        if (j>0)    { y[p] += stencil[3] * x[p-nx]; }
        if (j<ny-1) { y[p] += stencil[4] * x[p+nx]; }
        // z-1, z+1
        if (k>0)    { y[p] += stencil[5] * x[p-nx*ny]; }
        if (k<nz-1) { y[p] += stencil[6] * x[p+nx*ny]; }
      }
    }
  }
}
#endif


