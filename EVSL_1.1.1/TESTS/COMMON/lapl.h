#include <math.h>
#include <float.h>
#include "evsl.h"

/* datatype for performing matvec for Laplacians */
typedef struct _lapmv_t
{
   int nx, ny, nz;
   double *stencil;
} lapmv_t;

#ifdef __cplusplus
extern "C" {
#endif

int lapgen(int nx, int ny, int nz, cooMat *Acoo);
int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo);
void Lap2D3DMatvec(double *x, double *y, void *data);

#ifdef __cplusplus
}
#endif

