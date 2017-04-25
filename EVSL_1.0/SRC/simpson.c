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

void simpson(double* xi, double* yi, int npts, double* si) {
  int tmp = ((npts - 1) / 2);
  int m = tmp > 0 ? floor(tmp) : ceil(tmp);  // fix
  double tc = 0;
  double ti = 0;
  si[0] = 0;

  printf("\n npts: %i \n", npts);
  printf("\nx:\n");
  for (int i = 0; i < npts; i++) {
    printf("%f\t", xi[i]);
  }
  printf("\ny:\n");
  for (int i = 0; i < npts; i++) {
    printf("%f, ", yi[i]);
  }
  printf("\nsi:\n");
  for (int i = 0; i < m; i++) {
    printf("%f, ", si[i]);
  }
  printf("\n\n\n");
  printf("\n");
  for (int ii = 1; ii < m; ii++) {
    int i = 2 * ii;
    ti = (xi[i + 1] - xi[i - 1]) * (yi[i - 1] + 4 * yi[i] + yi[i + 1]) / 6.0;
    tc = tc + ti;
    si[i + 1] = tc;
    printf("i: %i,", i);
  }
  printf("si:\n");
  for (int i = 0; i < m; i++) {
    printf("%f, ", si[i]);
  }
  printf("\n\n\n");
  for (int i = 1; i < npts; i += 2) {
    ti = (xi[i] - xi[i - 1]) * (yi[i] + yi[i - 1]) / 2.0;
    si[i] = si[i - 1] + ti;
  }
  printf("si:\n");
  for (int i = 0; i < m; i++) {
    printf("%f, ", si[i]);
  }
  printf("\n\n\n");
}
