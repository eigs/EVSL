#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "internal_header.h"
//#include "vector.h"
/**
 * @file simpson.c
 * @brief Simpson integrater
 */

/**----------------------------------------------------------------------
 *
 *    This function computes the integrals from xi[0] to xi[j] for j=0:npts-1
 *
 *    @param[in] xi npts equally space points
 *    @param[in,out] yi values of a function f at the xi
 *    @param[in] npts number of sample points
 *
 *
 *
 *
 *    @note In-place version.
 *----------------------------------------------------------------------*/

void simpson(double* xi, double* yi, int npts) {
  double tc = 0.0, ti, tm, ysav;
  int i = 1;
  //-------------------- save yi[even_i]
  ysav = yi[0];
  yi[0] = 0;
  while (i < npts - 1) {
    // -------------------- ysav = yi[i-1] saved
    // trapeze portion
    ti = (xi[i] - xi[i - 1]) * (ysav + yi[i]) / 2.0;
    tm = tc + ti;
    // simpsion portion
    ti = (xi[i + 1] - xi[i - 1]) * (ysav + 4 * yi[i] + yi[i + 1]) / 6.0;
    tc += ti;
    // save and replace
    yi[i] = tm;
    ysav = yi[i + 1];
    yi[i + 1] = tc;
    i += 2;
  }
  //-------------------- case when npts is even
  if (i == npts - 1) {
    ti = (xi[i] - xi[i - 1]) * (yi[i] + ysav) / 2.0;
    yi[i] = tc + ti;
  }
}
