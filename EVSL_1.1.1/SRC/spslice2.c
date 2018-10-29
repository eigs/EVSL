#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "internal_header.h"

/**
 * @file spslice2.c
 * @brief Spectrum slicing
 */
/**----------------------------------------------------------------------
 *
 *    @brief Interval partitioner based for Lanczos DOS output
 *
 *    @param[in] xi coordinates of interval [a b]
 *    @param[in, out] yi yi[k] = integral of the does from a to xi[k]
 *    @param[in] n_int Number of desired sub-intervals
 *    @param[in] npts number of integration points (length of xi)
 *
 *    @param[out] sli Array of length n_int containing the boundaries of
 *       the intervals. [sli[i], [sli[i+1]] is the i-th interval with sli[0]
 *       = xi[0] and sli[n_int] = xi[npts-1]
 *
 *----------------------------------------------------------------------*/

void spslicer2(double* xi, double* yi, int n_int, int npts, double* sli) {
  /*-------------------- makes a call here to  integration by Simpson */
  double want;
  int k = 0;
  double t;
  int ls = 0;

  //-------------------- in-place integration ydos<--- int ydos..
  simpson(xi, yi, npts);
  //
  t = yi[0];
  want = (yi[npts - 1] - yi[0]) / (double)n_int;
  sli[ls] = xi[k];
  //-------------------- First point - t should be zero actually
  for (k = 1; k < npts; k++) {
    if (yi[k] - t >= want) {
      //-------------------- New interval defined
      ls = ls + 1;
      sli[ls] = xi[k];
      t = yi[k];
    }
  }
  //-------------------- bound for last interval is last point.
  sli[n_int] = xi[npts - 1];
}
