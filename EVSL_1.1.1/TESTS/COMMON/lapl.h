#include <math.h>
#include <float.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

#ifdef __cplusplus
extern "C" {
#endif

int lapgen(int nx, int ny, int nz, cooMat *Acoo);
int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo);
#ifdef __cplusplus
}
#endif
