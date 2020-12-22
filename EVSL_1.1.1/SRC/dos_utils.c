#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "internal_header.h"

/**
 *  @file dos_utils.c
 *  @brief A number of utility functions related to DOS functionality
 */
/**
 *  @brief Reciprocal
 *  @param[in] a Number to take reciprocal of
 *  @return reciprocal of a
 */
double rec(const double a) { return 1.0 / a; }  // Reciprocal

/**
 * @brief Inverse  square root
 * @param[in] a number to take inverse square root of
 * @return inverse square root of a
 * */
double isqrt(const double a) { return 1.0 / sqrt(a); }  // Inverse square root

/**
 * Initalize BSolDataPol: use L-S polynomial to approximate a function 'ffun'
 * @param[in] n Number
 * @param[in] max_deg Max degree of poliynomial
 * @param[in] tol Tolerance to be used
 * @param[in] lmin Minimum eigenvalue
 * @param[in] lmax Maximum eigenvalue
 * @param[in] ffun Function approximate
 * @param[in, out] data Structure to be initialized
 *
 */
void SetupBPol(int n, int max_deg, double tol, double lmin, double lmax,
               double (*ffun)(double), BSolDataPol *data) {
  data->max_deg = max_deg;
  data->intv[0] = lmin;
  data->intv[1] = lmax;
  data->tol = tol;
  data->mu = evsl_Calloc(max_deg, double);
  data->wk = evsl_Malloc_device(3*n, double);
  /* determine degree and compute poly coeff */
  lsPol(ffun, data);
}

/**
 * Initialize the member of BSolDataPol struct for solving B^{-1}
 *
 *
 * @param[in] n Number
 * @param[in] max_deg Max degree of poliynomial
 * @param[in] tol Tolerance to be used
 * @param[in] lmin Minimum eigenvalue
 * @param[in] lmax Maximum eigenvalue
 * @param[in, out] data Structure to be initialized
 */

void SetupPolRec(int n, int max_deg, double tol, double lmin, double lmax,
                 BSolDataPol *data) {
  SetupBPol(n, max_deg, tol, lmin, lmax, rec, data);
}

/**
 * Initialize the member of BSolDataPol struct for solving B^{1/2}
 *
 * @param[in] n Number
 * @param[in] max_deg Max degree of poliynomial
 * @param[in] tol Tolerance to be used
 * @param[in] lmin Minimum eigenvalue
 * @param[in] lmax Maximum eigenvalue
 * @param[in, out] data Structure to be initialized
 */
void SetupPolSqrt(int n, int max_deg, double tol, double lmin, double lmax,
                  BSolDataPol *data) {
  SetupBPol(n, max_deg, tol, lmin, lmax, isqrt, data);
}

/**
 * Free the BSolDataPol struct
 * @param[in, out] data struct to free data
 */
void FreeBSolPolData(BSolDataPol *data) {
  evsl_Free_device(data->wk);
  evsl_Free(data->mu);
}

/*
 * Setup the function pointer for evsl struct to call B_sol function
 *
 * @param[in] b Input vector
 * @param[out] x p(A)b
 * @param[in, out] data Data to be cast to BSolDataPol
 */
void BSolPol(double *b, double *x, void *data) {
  BSolDataPol *pol = (BSolDataPol *)data;
  pnav(pol->mu, pol->deg, pol->cc, pol->dd, b, x, pol->wk);
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
 * @param[in] mu Coefficents of the cheb. polynomial (size m+1)
 * @param[in] cc cc member of pol struct
 * @param[in] dd dd member of pol struct
 * @param[in] m m member of pol struct
 * @param[in] v input vector
 *
 * @param[out] y p(A)v
 *
 * @b Workspace
 * @param[in,out] w Work vector of length 3*n [allocate before call
 **/
int pnav(double *mu, const int m, const double cc, const double dd, double *v,
         double *y, double *w) {  // Really just ChebAv
  int n = evsldata.n;
  /*-------------------- pointers to v_[k-1],v_[k], v_[k+1] from w */
  double *vk = w;
  double *vkp1 = w + n;
  double *vkm1 = vkp1 + n;
  /*-------------------- */
  int k, one=1;
  double t1 = 1.0 / dd;
  double t2 = 2.0 / dd;
  /*-------------------- vk <- v; vkm1 <- zeros(n,1) */
#if 0
  /* LEAVE HERE IT FOR REFERENCE */
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
#else
  /* OPTIMIZATION */
  evsl_memcpy_device_to_device(y, v, n*sizeof(double));
  evsl_dscal_device(&n, mu, y, &one);
  /*-------------------- degree loop. k IS the degree */
  for (k = 1; k <= m; k++) {
    double *v_cur = k == 1 ? v : vk;
    double *v_old = k == 2 ? v : vkm1;
    double t = k == 1 ? t1 : t2;

    matvec_B(v_cur, vkp1);

#ifdef EVSL_USING_CUDA_GPU
    /* fuse 3 blas-1 into 1 kernel */
    evsl_pnav_device(n, k, t, mu[k], cc, vkp1, v_cur, v_old, y);
#else
    double ncc = -cc;
    double dmone = -1.0;
    evsl_daxpy_device(&n, &ncc, v_cur, &one, vkp1, &one);
    evsl_dscal_device(&n, &t, vkp1, &one);
    /*-------------------- for degree 2 and up: */
    if (k > 1) {
      evsl_daxpy_device(&n, &dmone, v_old, &one, vkp1, &one);
    }
    /*-------------------- y = mu[k]*Vk + y */
    evsl_daxpy_device(&n, &mu[k], vkp1, &one, y, &one);
#endif
#endif
    /*-------------------- next: rotate vectors via pointer exchange */
    double *tmp = vkm1;
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
 *    @param[in] ffun Function to generate an approximation for
 *
 *    @param[in,out] pol polparams struct \\
 *      Contains: cc = (a + b) / 2 \\
 *                dd = (b - a) / 2 \\
 *                mu = coefficients \\
 *                max_deg = number of coefficients
 *
 *
 *----------------------------------------------------------------------*/

int lsPol(double (*ffun)(double), BSolDataPol *pol) {
  const double a = pol->intv[0];
  const double b = pol->intv[1];
  pol->cc = (a + b) / 2;
  pol->dd = (b - a) / 2;
  /*------------------------- Number of points for Gauss-Chebyshev integration */
  int maxDeg = pol->max_deg;
  const int npts = maxDeg * 4;
  double *theti;
  theti = evsl_Malloc(npts, double);

  int i = 0;
  for (i = 0; i < npts * 2; i += 2) {
    theti[i / 2] = (i + 1) * (EVSL_PI / (2 * npts));
  }

  double *xi;
  xi = evsl_Malloc(npts, double);

  for (i = 0; i < npts; i++) {
    xi[i] = cos(theti[i]);
  }

  double *gi;
  gi = evsl_Malloc(npts, double);
  apfun(pol->cc, pol->dd, xi, ffun, npts, gi);
  evsl_Free(xi);

  /* ----------------------- Degree loop
   * ----------------------- Each coefficient generated by gauss-Chebyshev quadature */
  const int ptsNrm = 50 * maxDeg;  // Number of points for norm infinity
  xi = evsl_Malloc(ptsNrm, double);
  linspace(-1, 1, ptsNrm, xi);

  /* Compute f(x) once */
  double *yx;
  yx = evsl_Malloc(ptsNrm, double);
  apfun(pol->cc, pol->dd, xi, ffun, ptsNrm, yx);

  double *yi;
  yi = evsl_Malloc(npts, double);

  double *ya;
  ya = evsl_Malloc(ptsNrm, double);

  double na;

  int k = 0;
  const double t1 = 1.0 / npts;
  const double t2 = 2.0 / npts;

  /* degree loop */
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
    if (nrm < pol->tol) {
      k++;
      break;
    }
  }
  /* mu is of size k, so deg = k - 1 */
  pol->deg = k - 1;

  evsl_Free(xi);
  evsl_Free(yi);
  evsl_Free(ya);
  evsl_Free(yx);
  evsl_Free(theti);
  evsl_Free(gi);

  return 0;
}

