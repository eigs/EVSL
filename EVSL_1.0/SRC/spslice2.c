#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blaslapack.h"
#include "def.h"
#include "internal_proto.h"
#include "string.h"  //for memset
#include "struct.h"

/**----------------------------------------------------------------------
 *
 *    Interval partitioner
 *
 *    @param[in] xi coordinates of interval [a b]
 *    @param[in] yi yi[k] = integral of the does from a to xi[k]
 *    @param[in] n_int Number of desired sub-intervals
 *    @param[in] npts number of integration points (length of xi)
 *
 *    @param[out] sli Array of length n_int containing the boundaries of
 *       the intervals. [sli[i], [sli[i+1]] is the i-th interval with sli[0]
 *       = xi[0] and sli[n_int] = xi[npts-1]
 *
 *
 *----------------------------------------------------------------------*/

void spslicer2(double* xi, double* yi, int n_int, int npts, double* sli) {
  // Assumes yi has already been integrated
  double want = (yi[npts - 1] - yi[0]) / n_int;
  int k = 0;
  double t = yi[0];
  int ls = 0;

  sli[ls] = xi[k];
  // First point - t should be zero actually
  for (k = 1; k < npts; k++) {
    if (yi[k] - t >= want) {
      ls = ls + 1;  // New interval found
      sli[ls] = xi[k];
      t = yi[k];
    }
  }
  sli[n_int - 1] = xi[npts - 1];
}
