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
 *    Integrate the dos function at all poiints
 *
 *    @param[in] xi npts equally space points
 *    @param[in] yi values of a function f at the xi
 *    @param[in] npts number of sample points
 *
 *    @param[out] si Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *
 *----------------------------------------------------------------------*/   

int simspon(double* xi, double* yi, int npts, double* si) {
  int tmp = ((npts-1)/2);
  int m = tmp>0?floor(tmp):ceil(tmp); // fix
  int tc = 0;
  si[0] = 0;

  for(int ii = 0; ii < m; ii++) {
    int i = 2*ii;
    double ti = (xi[i+1] - xi[i-1]) * (yi[i-1]+4*yi[i]+yi[i+1])/6.0;
    tc = tc + ti;
    si[i+1] = tc;
  }
  for(int i=1; i < npts; i+=2) {
    ti = (xi[i]-xi[i-1]) * (yi[i] + yi[i-1])/2.0;
    si[i] = si[i-1]+ti;
  }
  return 0;
}

