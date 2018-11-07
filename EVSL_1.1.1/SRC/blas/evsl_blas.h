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
int evsl_daxpy(EVSL_Int *n, EVSL_Real *da, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
int evsl_dcopy(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
EVSL_Real evsl_ddot(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
int evsl_dgemm(const char *transa, const char *transb, EVSL_Int *m, EVSL_Int * n, EVSL_Int *k, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *b, EVSL_Int *ldb, EVSL_Real *beta, EVSL_Real *c__, EVSL_Int *ldc);
int evsl_dgemv(const char *trans, EVSL_Int *m, EVSL_Int *n, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *x, EVSL_Int *incx, EVSL_Real *beta, EVSL_Real *y, EVSL_Int *incy);
EVSL_Real evsl_dnrm2(EVSL_Int *n, EVSL_Real *x, EVSL_Int *incx);
int evsl_dscal(EVSL_Int *n, EVSL_Real *da, EVSL_Real *dx, EVSL_Int *incx);
int evsl_dcopy(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
#if defined(LAPACK_F2C_INCLUDE)
/* blas needed from evsl lapack only */
int evsl_dtrmm(const char *side, const char *uplo, const char *transa, const char *diag, EVSL_Int *m, EVSL_Int *n, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *b, EVSL_Int *ldb);
int evsl_dsyr2k(const char *uplo, const char *trans, EVSL_Int *n, EVSL_Int *k, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *b, EVSL_Int *ldb, EVSL_Real *beta, EVSL_Real *c__, EVSL_Int *ldc);
int evsl_dsymv(const char *uplo, EVSL_Int *n, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *x, EVSL_Int *incx, EVSL_Real *beta, EVSL_Real *y, EVSL_Int *incy);
int evsl_dger(EVSL_Int *m, EVSL_Int *n, EVSL_Real *alpha, EVSL_Real *x, EVSL_Int *incx, EVSL_Real *y, EVSL_Int *incy, EVSL_Real *a, EVSL_Int *lda);
EVSL_Int lsame_(const char *ca, const char *cb);
int evsl_xerbla(const char *srname, EVSL_Int *info);
int evsl_zgemm(const char *transa, const char *transb, EVSL_Int *m, EVSL_Int *n, EVSL_Int *k, doublecomplex *alpha, doublecomplex *a, EVSL_Int *lda, doublecomplex *b, EVSL_Int *ldb, doublecomplex *beta, doublecomplex *c__, EVSL_Int *ldc);
int evsl_zgeru(EVSL_Int *m, EVSL_Int *n, doublecomplex *alpha, doublecomplex *x, EVSL_Int *incx, doublecomplex *y, EVSL_Int *incy, doublecomplex *a, EVSL_Int *lda);
int evsl_ztrsm(const char *side, const char *uplo, const char *transa, const char *diag, EVSL_Int *m, EVSL_Int *n, doublecomplex *alpha, doublecomplex *a, EVSL_Int *lda, doublecomplex *b, EVSL_Int *ldb);
EVSL_Int evsl_izamax(EVSL_Int *n, doublecomplex *zx, EVSL_Int *incx);
int evsl_zswap(EVSL_Int *n, doublecomplex *zx, EVSL_Int *incx, doublecomplex *zy, EVSL_Int *incy);
int evsl_zscal(EVSL_Int *n, doublecomplex *za, doublecomplex *zx, EVSL_Int *incx);
int evsl_dtrmv(const char *uplo, const char *trans, const char *diag, EVSL_Int *n, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *x, EVSL_Int *incx);
int evsl_dswap(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
int evsl_dsyr2(const char *uplo, EVSL_Int *n, EVSL_Real *alpha, EVSL_Real *x, EVSL_Int *incx, EVSL_Real *y, EVSL_Int *incy, EVSL_Real *a, EVSL_Int *lda);
#endif
#if defined(__cplusplus) && !defined(EVSL_USING_EVSL_BLAS)
}
#endif

#endif

