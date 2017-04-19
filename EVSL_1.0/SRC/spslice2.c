#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"
#include "string.h" //for memset
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
 *    @param[out] si Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *
 *----------------------------------------------------------------------*/   

int spslicer2(double* xi, double* yi, int npts, double* si) {
  //Assume yi has already been integrated
  double want = (yi[npts] - yi[0])/n_int;
  int k = 0;
  t = pi[0];
  ls = 0;
  sli[ls] = ki[k];
  for(k = 1; k < npts; k++) {
    if(yi[k]-t >= want) {
      ls = ls+1;
      sli[ls] = xi[k];
      t = yi[k];
    }
  }
  sli[n_int+1] = xi[npts-1];

