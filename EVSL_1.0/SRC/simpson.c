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
 *    This function computes the integrals from xi[0] to xi[j] for j=0:npts-1
 *
 *    @param[in] xi npts equally space points
 *    @param[in] yi values of a function f at the xi
 *    @param[in] npts number of sample points
 *
 *    @param[out] si Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Uses simposen's rule completed by trapezoidal
 *    rule for middle points
 *
 *
 *----------------------------------------------------------------------*/

void simpson(double* xi, double* yi, int npts, double* si) {
  double tc = 0;
  double ti = 0;
  si[0] = 0;

  for (int i = 1; i < npts - 1; i += 2) {  // simpson rule loop
    ti = (xi[i + 1] - xi[i - 1]) * (yi[i - 1] + 4 * yi[i] + yi[i + 1]) / 6.0;
    tc = tc + ti;
    si[i + 1] = tc;
  }

  for (int i = 1; i < npts; i += 2) {  // trapezoidal rule for other points
    ti = (xi[i] - xi[i - 1]) * (yi[i] + yi[i - 1]) / 2.0;
    si[i] = si[i - 1] + ti;
  }
}
