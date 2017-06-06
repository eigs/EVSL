#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blaslapack.h"
#include "def.h"
#include "evsl.h"
#include "internal_proto.h"
#include "struct.h"
/** ---------------------------------------------------------------------
 * Contains various functions and utilities used primarily for DOS functions.
 * */

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

int apfun1(const double c, const double h, const double *const xi,
           double (*ffun)(double), const int npts, double *yi) {
  int i = 0;
  for (i = 0; i < npts; i++) {
    yi[i] = ffun(c + h * xi[i]);
  }
  return 0;
}

/**-------------------------------------------------------------------------
 *
 * Evaluates p_n(A) v where p_n is expressed as a combination of Chebyshev
 * polynomials.
 *
 * @param[in] A Matrix
 * @param[in] c
 * @param[in] h
 * @param[in] mu
 * @param[in] v
 *
 * @param[out] v
 *
 * */
#if 0
int pnav(const double c, const double h, const double *const mu,
         const int degMu, double v, double *y) {
  int n = evsldata.n;
  int one = 1;
  double *vm1;
  Malloc(vm1, n, double);

  double *ytmp;
  Malloc(ytmp, n, double);
  double *ytmp2;
  Malloc(ytmp2, n, double);
  // Where to malloc y?
  int i;
  for (i = 0; i < n; i++) {
    y[i] = mu[0] * v;
  }
  double scal;
  int j;
  for (i = 0; i < degMu; i++) {
    scal = 2 / h;
    if (i == 0) {
      scal = 1 / h;
    }
    matvec_A(v, ytmp);
    DCOPY(&n, ytmp2, &one, v, &one);
    DSCAL(&n, &c, ytmp2, &one);
    for (j = 0; j < n; j++) {
      ytmp = ytmp[j] - ytmp2[j];
    }
    DSCAL(&n, &scal, ytmp, &one);
  }
}
#endif
