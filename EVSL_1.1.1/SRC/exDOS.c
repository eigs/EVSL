#include <stdio.h>
#include <string.h>  //for memset
#include <stdlib.h>
#include "internal_header.h"

/**
 * Calculate the exact DOS given eigenvalues
 *
 * @param[in] vals eigenvalues
 * @param[in] n number of eigenvalues
 * @param[in] npts number of points for dos curve
 * @param[in] intv intervals of interests
 * @param[out] x coordinates for dos plot
 * @param[out] y y coordinates for dos plot
 * @note Both x and y are expected to be preallocated.
 */
int exDOS(double *vals, int n, int npts, double *x, double *y, double *intv) {
/* sigma = coefficient for gaussian
   we select sigma so that at end of subinterval
   we have exp[-0.5 [H/sigma]^2] = 1/K. */
  int i,j=0, one=1;
  double h, xi, t, scaling;
  const double lm = intv[2];
  const double lM = intv[3];
  const double kappa = 1.25;
  const int M = evsl_min(n, 30);
  const double H = (lM - lm) / (M - 1);
  const double sigma = H / sqrt(8 * log(kappa));
  const double a = intv[0];
  const double b = intv[1];
  double sigma2 = 2 * sigma * sigma;
/*-------------------- sigma related..
  -------------------- if gaussian smaller than tol ignore point. */
  const double tol = 1e-08;
  double width = sigma * sqrt(-2.0 * log(tol));
  /* ------------------- define width of sub-intervals for DOS */
  linspace(a, b, npts, x);
  h = x[1]-x[0];
  memset(y, 0, npts*sizeof(double));
/*-------------------- scan all relevant eigenvalues and add to its
                       interval [if any] */
  for (i=0; i<n; i++) {
    t = vals[i];
    if (t < a || t > b) {
      continue;
    }
    /*-------------------- first point to consider  */
    j = evsl_max((int) ((t-a-width)/h), 0);
    /*xi = a+j*h;
    ! check!  */
    for  (xi = a+j*h; xi <= evsl_min(t+width, b); xi+=h)
      y[j++] += exp(-(xi-t)*(xi-t)/sigma2);
  }

  if (j == 0) {
    return 1;
  }

  /*-------------------- scale dos - [n eigenvalues in all]
    scaling 2 -- for plots due to the missing 1/sqrt[2*pi*sigma^2] factor.
    However, note that this does not guarantee that
    sum[y] * h = the correct number of eigenvalues
    y = y / sqrt[pi*sigma2] ; */
  scaling = 1.0 / (n*sqrt(sigma2*EVSL_PI));
  evsl_dscal(&npts, &scaling, y, &one);

  return 0;
}

