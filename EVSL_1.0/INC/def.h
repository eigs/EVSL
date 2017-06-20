#ifndef DEF_H
#define DEF_H

#include <stdlib.h>
#include <assert.h>
#include <math.h>

/*! \file def.h
    \brief defs in EVSL
*/

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
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

/*!
  \def max(x,y)
  Computes the maximum of \a x and \a y.
*/
#define max(a, b) ((a) > (b) ? (a) : (b))

/*!
  \def min(x,y)
  Computes the minimum of \a x and \a y.
*/
#define min(a, b) ((a) < (b) ? (a) : (b))

/*! Fortran interface naming convention
 */
#define EVSLFORT(name) name ## _f90_

/*! max number of Gram–Schmidt process in orthogonalization
 */
#define NGS_MAX 2

#endif
