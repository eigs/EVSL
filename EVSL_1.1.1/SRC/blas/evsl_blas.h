#ifndef EVSL_BLAS_H
#define EVSL_BLAS_H

/* use this file outside blas directory to provide headers */

#if !defined(EVSL_USING_EVSL_BLAS)
/* if not using EVSL_BLAS, cast names back to standard BLAS names */
#define evsl_daxpy  daxpy_
#define evsl_dcabs1 dcabs1_
#define evsl_dcopy  dcopy_
#define evsl_ddot   ddot_
#define evsl_dgemm  dgemm_
#define evsl_dgemv  dgemv_
#define evsl_dger   dger_
#define evsl_dnrm2  dnrm2_
#define evsl_dscal  dscal_
#define evsl_dswap  dswap_
#define evsl_dsymv  dsymv_
#define evsl_dsyr2  dsyr2_
#define evsl_dsyr2k dsyr2k_
#define evsl_dtrmm  dtrmm_
#define evsl_dtrmv  dtrmv_
#define evsl_izamax izamax_
#define evsl_lsame  lsame_
#define evsl_xerbla xerbla_
#define evsl_zgemm  zgemm_
#define evsl_zgeru  zgeru_
#define evsl_zscal  zscal_
#define evsl_zswap  zswap_
#define evsl_ztrsm  ztrsm_
#endif

#if defined(__cplusplus) && !defined(EVSL_USING_EVSL_BLAS)
extern "C" {
#endif
int evsl_daxpy(const EVSL_Int *n, const EVSL_Real *da, const EVSL_Real *dx, const EVSL_Int *incx, EVSL_Real *dy, const EVSL_Int *incy);
int evsl_dcopy(const EVSL_Int *n, const EVSL_Real *dx, const EVSL_Int *incx, EVSL_Real *dy, const EVSL_Int *incy);
EVSL_Real evsl_ddot(const EVSL_Int *n, const EVSL_Real *dx, const EVSL_Int *incx, const EVSL_Real *dy, const EVSL_Int *incy);
int evsl_dgemm(const char *transa, const char *transb, const EVSL_Int *m, const EVSL_Int * n, const EVSL_Int *k, const EVSL_Real *alpha, const EVSL_Real *a, const EVSL_Int *lda, const EVSL_Real *b, const EVSL_Int *ldb, const EVSL_Real *beta, EVSL_Real *c__, const EVSL_Int *ldc);
int evsl_dgemv(const char *trans, const EVSL_Int *m, const EVSL_Int *n, const EVSL_Real *alpha, const EVSL_Real *a, const EVSL_Int *lda, const EVSL_Real *x, const EVSL_Int *incx, const EVSL_Real *beta, EVSL_Real *y, const EVSL_Int *incy);
EVSL_Real evsl_dnrm2(const EVSL_Int *n, const EVSL_Real *x, const EVSL_Int *incx);
int evsl_dscal(const EVSL_Int *n, const EVSL_Real *da, EVSL_Real *dx, const EVSL_Int *incx);
int evsl_dcopy(const EVSL_Int *n, const EVSL_Real *dx, const EVSL_Int *incx, EVSL_Real *dy, const EVSL_Int *incy);
#if defined(LAPACK_F2C_INCLUDE)
/* blas needed from evsl lapack only */
int evsl_dtrmm(const char *side, const char *uplo, const char *transa, const char *diag, const EVSL_Int *m, const EVSL_Int *n, const EVSL_Real *alpha, const EVSL_Real *a, const EVSL_Int *lda, EVSL_Real *b, const EVSL_Int *ldb);
int evsl_dsyr2k(const char *uplo, const char *trans, const EVSL_Int *n, const EVSL_Int *k, const EVSL_Real *alpha, const EVSL_Real *a, const EVSL_Int *lda, const EVSL_Real *b, const EVSL_Int *ldb, const EVSL_Real *beta, EVSL_Real *c__, const EVSL_Int *ldc);
int evsl_dsymv(const char *uplo, EVSL_Int *n, const EVSL_Real *alpha, const EVSL_Real *a, const EVSL_Int *lda, const EVSL_Real *x, const EVSL_Int *incx, const EVSL_Real *beta, EVSL_Real *y, const EVSL_Int *incy);
int evsl_dger(const EVSL_Int *m, const EVSL_Int *n, const EVSL_Real *alpha, const EVSL_Real *x, const EVSL_Int *incx, const EVSL_Real *y, const EVSL_Int *incy, EVSL_Real *a, const EVSL_Int *lda);
EVSL_Int lsame_(const char *ca, const char *cb);
int evsl_xerbla(const char *srname, EVSL_Int *info);
int evsl_zgemm(const char *transa, const char *transb, const EVSL_Int *m, const EVSL_Int *n, const EVSL_Int *k, const doublecomplex *alpha, const doublecomplex *a, const EVSL_Int *lda, const doublecomplex *b, const EVSL_Int *ldb, const doublecomplex *beta, doublecomplex *c__, const EVSL_Int *ldc);
int evsl_zgeru(const EVSL_Int *m, const EVSL_Int *n, const doublecomplex *alpha, const doublecomplex *x, const EVSL_Int *incx, const doublecomplex *y, const EVSL_Int *incy, doublecomplex *a, const EVSL_Int *lda);
int evsl_ztrsm(const char *side, const char *uplo, const char *transa, const char *diag, const EVSL_Int *m, const EVSL_Int *n, const doublecomplex *alpha, const doublecomplex *a, const EVSL_Int *lda, doublecomplex *b, const EVSL_Int *ldb);
EVSL_Int evsl_izamax(const EVSL_Int *n, const doublecomplex *zx, const EVSL_Int *incx);
int evsl_zswap(const EVSL_Int *n, doublecomplex *zx, const EVSL_Int *incx, doublecomplex *zy, const EVSL_Int *incy);
int evsl_zscal(const EVSL_Int *n, const doublecomplex *za, doublecomplex *zx, const EVSL_Int *incx);
int evsl_dtrmv(const char *uplo, const char *trans, const char *diag, const EVSL_Int *n, const EVSL_Real *a, const EVSL_Int *lda, EVSL_Real *x, const EVSL_Int *incx);
// dswap is not specified with correct const due to conflicting decl in dsteqr.c
int evsl_dswap(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
int evsl_dsyr2(const char *uplo, const EVSL_Int *n, const EVSL_Real *alpha, const EVSL_Real *x, const EVSL_Int *incx, const EVSL_Real *y, const EVSL_Int *incy, EVSL_Real *a, const EVSL_Int *lda);
#endif
#if defined(__cplusplus) && !defined(EVSL_USING_EVSL_BLAS)
}
#endif

#endif
