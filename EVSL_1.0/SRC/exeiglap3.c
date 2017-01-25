#include <math.h>
#include <float.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo) {
  double thetax = 0.5 * PI / (nx + 1.0);
  double thetay = 0.5 * PI / (ny + 1.0);
  double thetaz = 0.5 * PI / (nz + 1.0);
  int i, j, k, l=0, n=nx*ny*nz;
  double *v;
  double tx, ty, tz, ev;
  Malloc(v, n, double);
  for (i=1; i<=nx; i++) {
    tx = sin(i*thetax);
    for (j=1; j<=ny; j++) {
      ty = sin(j*thetay);
      for (k=1; k<=nz; k++) {
        tz = sin(k*thetaz);
        if (1 == nz) {
          tz = 0.0;
        }
        ev = 4*(tx*tx+ty*ty+tz*tz);
        if (ev >= a - DBL_EPSILON && ev <= b + DBL_EPSILON) {
          v[l++] = ev;
        }
      }
    }
  }
  Realloc(v, l, double);
  sort_double(l, v, NULL);

  *m = l;
  *vo = v;
  return 0;
}

