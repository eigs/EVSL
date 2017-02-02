#ifndef DEF_H
#define DEF_H

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define PI M_PI 
//3.14159265358979323846
#define orthTol 1e-14

#define CHKERR(ierr) assert(!(ierr))
//#define CHKREQ(ierr) { if (ierr) { return (ierr); } }

#define Malloc(base, nmem, type) {\
  (base) = (type*) malloc((nmem)*sizeof(type)); \
  CHKERR((base) == NULL); \
}
#define Calloc(base, nmem, type) {\
  (base) = (type*) calloc((nmem), sizeof(type)); \
  CHKERR((base) == NULL); \
}
#define Realloc(base, nmem, type) {\
  (base) = (type*) realloc((base), (nmem)*sizeof(type)); \
  CHKERR((base) == NULL && nmem > 0); \
}

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/* max number of Gram–Schmidt process in orthogonalization */
#define NGS_MAX 2

#endif
