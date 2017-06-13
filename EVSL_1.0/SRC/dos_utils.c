#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "blaslapack.h"
#include "def.h"
#include "evsl.h"
#include "internal_proto.h"
#include "struct.h"

double rec(const double a) { return 1.0 / a; }  // Reciprocal

double isqrt(const double a) { return 1.0 / sqrt(a); }  // Inverse square root

/*
 * Extract the square root of diagonal entries of the mass matrix B
 *
 * @param[in] B Matrix to extract the square root of the diagonals
 * @param[out] sqrtdiag preallocated vector of lengeth B.ncols to put the sqrt
 * of
 *  diagonals into
 */
void extractDiag(cooMat *B, double *sqrtdiag) {
  const int nnz = B->nnz;
  int i, row, col;
  for (i = 0; i < nnz; i++) { /* Iterate over all elements looking for diags */
    row = B->ir[i];
    col = B->jc[i];
    if (row == col) {
      sqrtdiag[col] = sqrt(B->vv[i]);
    }
  }
}

/*
 * Initialize the member of BSolDataPol struct for solving B
 */
void SetupBSolPol(csrMat *B, BSolDataPol *data) {
  const int n = B->nrows, mdeg = 200;
  set_pol_def(&data->pol_sol);
  double *mu_sol = (double *)calloc(mdeg, sizeof(double));
  data->pol_sol.mu = mu_sol;
  data->wk = (double *)malloc(3 * n * sizeof(double));
  lsPol(&data->intv[0], mdeg, rec, 1e-5, &data->pol_sol);
}

/*
 * Initialize the member of BSolDataPol struct for solving B^{1/2}
 */
void SetupBsqrtSolPol(csrMat *B, BSolDataPol *data) {
  const int n = B->nrows, mdeg = 200;
  set_pol_def(&data->pol_sol);
  double *mu_sol = (double *)calloc(mdeg, sizeof(double));
  data->pol_sol.mu = mu_sol;
  data->wk = (double *)malloc(3 * n * sizeof(double));
  lsPol(&data->intv[0], mdeg, isqrt, 1e-5, &data->pol_sol);
}

/*
 * Free the BSolDataPol struct
 */
void FreeBSolPolData(BSolDataPol *data) {
  free(data->wk);
  free_pol(&data->pol_sol);
}

/*
 * Setup the function pointer for evsl struct to call B_sol function
 */
void BSolPol(double *b, double *x, void *data) {
  BSolDataPol *Bsol_data = (BSolDataPol *)data;
  double *wk = Bsol_data->wk;
  polparams pol = Bsol_data->pol_sol;
  pnav(pol.mu, pol.deg, pol.cc, pol.dd, b, x, wk);
}

/*
 * Diagonal scaling for A and B such that A(i,j) =
 * A(i,j)/(sqrtdiag(i)*sqrtdiag(j)) and B(i,j) =
 * B(i,j)/(sqrtdiag(i)*sqrtdiag(j))
 *
 * @param[in,out] A Matrix to scale using sqrtdiag
 * @param[in,out] B Matrix to scale using sqrtdiag
 * @param[in] sqrtdiag The vector contating the square root of the diagonal
 *  elements of B obtained via extractDiag
 */
void diagScaling(cooMat *A, cooMat *B, double *sqrtdiag) {
  int i, row, col, nnz;
  double tmp;
  // diagonal scaling for A
  nnz = A->nnz;
  for (i = 0; i < nnz; i++) {
    row = A->ir[i];
    col = A->jc[i];
    tmp = 1.0 / (sqrtdiag[row] * sqrtdiag[col]);
    A->vv[i] = A->vv[i] * tmp;
  }
  // diagonal scaling for B
  nnz = B->nnz;
  for (i = 0; i < nnz; i++) {
    row = B->ir[i];
    col = B->jc[i];
    tmp = 1.0 / (sqrtdiag[row] * sqrtdiag[col]);
    B->vv[i] = B->vv[i] * tmp;
  }
}

/**----------------------------------------------------------------------
 *
 *    Evalutes ffun at the xi's.
 *    Assumes a transformation of original inetrval [a b] into [-1, 1] so:
 *    the xi's are between -1 and 1
 *
 *    @param[in] c Value to increase xi's by
 *    @param[in] h Value to scale xi's by
 *    @param[in] *xi Points for which to evaluate ffun at
 *    @param[in] npts Number of points in xi to evaluate
 *    @param[in] ffun Function to evaluate
 *
 *    @param[out] yi ffun evaluated at xi's
 *
 *
 *----------------------------------------------------------------------*/

int apfun(const double c, const double h, const double *const xi,
          double (*ffun)(double), const int npts, double *yi) {
  int i = 0;
  for (i = 0; i < npts; i++) {
    yi[i] = ffun(c + h * xi[i]);
  }
  return 0;
}

/**
 * Computes y=P(A) v, where pn is a Cheb. polynomial expansion
 *
 * This explicitly calls matvec, so it can be useful for implementing
 * user-specific matrix-vector multiplication.
 *
 * @param pol Struct containing the paramenters and expansion coefficient of
 * the polynomail.
 * @param v input vector
 *
 * @param[out] y p(A)v
 *
 * @b Workspace
 * @param w Work vector of length 3*n [allocate before call]
 * @param v is untouched
 **/
int pnav(double *mu, const int m, const double cc, const double dd, double *v,
         double *y, double *w) {  // Really just ChebAv
  const int n = evsldata.n;
  /*-------------------- pointers to v_[k-1],v_[k], v_[k+1] from w */
  double *vk = w;
  double *vkp1 = w + n;
  double *vkm1 = vkp1 + n;
  /*-------------------- */
  int k, i;
  double t, s, *tmp;
  const double t1 = 1.0 / dd;
  const double t2 = 2.0 / dd;
  /*-------------------- vk <- v; vkm1 <- zeros(n,1) */
  memcpy(vk, v, n * sizeof(double));
  memset(vkm1, 0, n * sizeof(double));
  /*-------------------- special case: k == 0 */
  s = mu[0];
  for (i = 0; i < n; i++) {
    y[i] = s * vk[i];
  }
  /*-------------------- degree loop. k IS the degree */
  for (k = 1; k <= m; k++) {
    /*-------------------- y = mu[k]*Vk + y */
    t = k == 1 ? t1 : t2;
    /*-------------------- */
    s = mu[k];
    matvec_B(vk, vkp1);

    for (i = 0; i < n; i++) {
      vkp1[i] = t * (vkp1[i] - cc * vk[i]) - vkm1[i];
      /*-------------------- for degree 2 and up: */
      y[i] += s * vkp1[i];
    }
    /*-------------------- next: rotate vectors via pointer exchange */
    tmp = vkm1;
    vkm1 = vk;
    vk = vkp1;
    vkp1 = tmp;
  }

  return 0;
}

/**----------------------------------------------------------------------
 *
 *    Finds the least-square polynomial approximation to function ffun
 *    in interval given by intv
 *
 *    @param[in] *intv Stores the interval of interest \\
 *        intv[0:1] Interval of interest
 *    @param[in] maxDeg Maximum degree of the polynomial appoximation
 *    @param[in] ffun Function to generate an approximation for
 *    @param[in] tol Tolerance for approximation
 *
 *    @param[out] pol polparams struct \\
 *      Contains: cc = (a + b) / 2 \\
 *                dd = (b - a) / 2 \\
 *                mu = coefficients \\
 *                deg = number of coefficients
 *
 *
 *----------------------------------------------------------------------*/

int lsPol(const double *const intv, const int maxDeg, double (*ffun)(double),
          const double tol, polparams *pol) {
  const double a = intv[0];
  const double b = intv[1];
  pol->cc = (a + b) / 2;
  pol->dd = (b - a) / 2;
  //------------------------- Number of points for Gauss-Chebyshev
  // integration
  const int npts = maxDeg * 4;

  double *theti;
  Malloc(theti, npts, double);
  int i = 0;
  for (i = 0; i < npts * 2; i += 2) {
    theti[i / 2] = (i + 1) * (PI / (2 * npts));
  }

  double *xi;
  Malloc(xi, npts, double);

  for (i = 0; i < npts; i++) {
    xi[i] = cos(theti[i]);
  }

  double *gi;
  Malloc(gi, npts, double);
  apfun(pol->cc, pol->dd, xi, ffun, npts, gi);
  free(xi);

  // ----------------------- Degree loop
  // ----------------------- Each coefficient generated by gauss-Chebyshev
  // quadature.
  const int ptsNrm = 300;  // Number of points for norm infinity
  Malloc(xi, ptsNrm, double);
  linspace(-1, 1, ptsNrm, xi);

  // Compute f(x) once

  double *yx;
  Malloc(yx, ptsNrm, double);
  apfun(pol->cc, pol->dd, xi, ffun, ptsNrm, yx);

  double *yi;
  Malloc(yi, npts, double);

  double *ya;
  Malloc(ya, ptsNrm, double);

  double na;

  int k = 0;
  const double t1 = 1.0 / npts;
  const double t2 = 2.0 / npts;
  for (k = 0; k < maxDeg; k++) {
    for (i = 0; i < npts; i++) {
      yi[i] = cos(k * theti[i]);
    }
    const double t = k == 0 ? t1 : t2;
    double sum = 0;
    for (i = 0; i < npts; i++) {
      sum += yi[i] * gi[i];
    }
    pol->mu[k] = t * sum;
    chebxPltd(k, pol->mu, ptsNrm, xi, ya);
    double nrm = 0;
    for (i = 0; i < ptsNrm; i++) {
      na = (ya[i] - yx[i]) / yx[i];
      nrm += fabs(na);
    }
    if (nrm < tol) {
      pol->deg = k + 1;
      break;
    }
  }
  free(xi);
  free(yi);
  free(ya);
  free(yx);
  free(theti);
  free(gi);
  return 0;
}
