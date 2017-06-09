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

int apfun(const double c, const double h, const double *const xi,
           double (*ffun)(double), const int npts, double *yi) {
  int i = 0;
  for (i = 0; i < npts; i++) {
    yi[i] = ffun(c + h * xi[i]);
  }
  return 0;
}

