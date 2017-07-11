#ifndef BLAS_LAPACK_H
#define BLAS_LAPACK_H

#include <complex.h>

#ifdef USE_MKL

#define MKL_Complex16 complex
#include "mkl.h"

#else

#define DCOPY    dcopy_
#define DDOT     ddot_
#define DNRM2    dnrm2_
#define DSCAL    dscal_
#define DASUM    dasum_
#define DGEMV    dgemv_
#define DGEMM    dgemm_
#define DAXPY    daxpy_
#define DSTEV    dstev_
#define DSYEV    dsyev_
#define DSTEMR   dstemr_
#define DHSEQR   dhseqr_
#define ZGESV    zgesv_

/**
 * @file blaslapack.h
 * @brief Defs for blaslapack routines
 */

// Fortran logical type
typedef int logical;

void DCOPY(int *n, double *dx, int *incx, double *dy, int *incy);
void DAXPY(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
void DSCAL(int *n,double *a,double *x,int *incx);
double DASUM(int *n,double *x,int *incx);
double DDOT(int *n,double *x,int *incx,double *y,int *incy);
double DNRM2(int *n,double *x,int *incx);
void DGEMM(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
void DGEMV(char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);
void DSTEV(char *jobz, int *n, double *diagonal, double *subdiagonal, double *V, int *ldz, double *work, int *info);
void DSYEV(char* jobz,char* uplo,int* n,double* fa,int* lda,double* w, double* work,int* lwork,int* info);
void DSTEMR(char *jobz, char *range, int *n, double *D, double *E, double *VL, double *VU, int *IL, int *IU, 
	    int *M, double *W, double *Z, int *LDZ, int *NZC, int *ISUPPZ, logical *TRYRAC, double *WORK, 
	    int *LWORK, int *IWORK, int *LIWORK, int *INFO);
void DHSEQR(char* jobz,char* compz,int* n,int* ilo,int* ihi,double* h,int* ldh,double* wr,double* wi,
	    double* z,int* ldz,double* work, int* lwork,int* info);
void ZGESV(int *n, int *nrow, complex double * A, int* m, int* ipiv, complex double *rhs, int* k, int* INFO);

#endif

/* USE SQRT(DDOT()) INSTEAD OF DNRM2, FOR BETTER PERFORMANCE */
#undef DNRM2
#define DNRM2(n, x, incx) sqrt(DDOT(n, x, incx, x, incx))

#endif
