#include <math.h>
#include <float.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/**-----------------------------------------------------------------------
 * @brief Laplacean Matrix generator 
 * @param nx  Number of points in x-direction
 * @param ny  Number of points in y-direction
 * @param nz  Number of points in z-direction
 * @param[out] *Acoo matrix in coordinate format. 
 -----------------------------------------------------------------------**/
int lapgen(int nx, int ny, int nz, cooMat *Acoo) {
  int n = nx * ny * nz;
  Acoo->nrows = n;
  Acoo->ncols = n;

  int nzmax = nz > 1 ? 7*n : 5*n;
  Malloc(Acoo->ir, nzmax, int);
  Malloc(Acoo->jc, nzmax, int);
  Malloc(Acoo->vv, nzmax, double);

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
 * @brief Exact eigenvalues of Laplacean in interval [a b]
 * @param nx  Number of points in x-direction
 * @param ny  Number of points in y-direction
 * @param nz  Number of points in z-direction
 * @param[out] m number of eigenvalues found 
 * @param[out] **vo pointer to array of eigenvalues found 
 *-----------------------------------------------------------------------**/
int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo) {
  double thetax = 0.5 * PI / (nx + 1.0);
  double thetay = 0.5 * PI / (ny + 1.0);
  double thetaz = 0.5 * PI / (nz + 1.0);
  int i, j, k, l=0, n=nx*ny*nz;
  double *v;
  double tx, ty, tz, ev;
  Malloc(v, n, double);
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
        if (ev >= a - DBL_EPSILON && ev <= b + DBL_EPSILON) {
          v[l++] = ev;
        }
      }
    }
  }
  Realloc(v, l, double);
  sort_double(l, v, NULL);

  *m = l;
  *vo = v;
  return 0;
}

