#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blaslapack.h"
#include "def.h"
#include "internal_proto.h"
#include "string.h"  //for memset
#include "struct.h"
//#include "vector.h"

/**----------------------------------------------------------------------
 *
 *    Interval partitioner
 *
 *    @param[in] xi coordinates of interval [a b]
 *    @param[in] yi[k] = integral of the does from a to xi[k]
 *    @param[in] n_int Number of desired sub-intervals
 *    @param[in] npts number of integration points (length of xi)
 *
 *    @param[out] sli Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *
 *----------------------------------------------------------------------*/

void spslicer2(double* xi, double* yi, int n_int, int npts, double* sli) {
  // Assume yi has already been integrated
  double want = (yi[npts - 1] - yi[0]) / n_int;
  int k = 0;
  double t = yi[0];
  int ls = 0;
  sli[ls] = xi[k];
  for (k = 1; k < npts; k++) {
    if (yi[k] - t >= want) {
      ls = ls + 1;
      sli[ls] = xi[k];
      t = yi[k];
    }
  }
  sli[n_int] = xi[npts - 1];
}
