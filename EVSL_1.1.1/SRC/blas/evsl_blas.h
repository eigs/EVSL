#ifndef EVSL_BLAS_H
#define EVSL_BLAS_H

/* use this file outside blas directory to provide headers */

#ifndef EVSL_USING_EVSL_BLAS
/* if not using EVSL_BLAS, cast them back standard BLAS names */
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
#ifdef __cplusplus
extern "C" {
#endif
#endif

int evsl_daxpy(EVSL_Int *n, EVSL_Real *da, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
int evsl_dcopy(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
EVSL_Real evsl_ddot(EVSL_Int *n, EVSL_Real *dx, EVSL_Int *incx, EVSL_Real *dy, EVSL_Int *incy);
int evsl_dgemm(const char *transa, const char *transb, EVSL_Int *m, EVSL_Int * n, EVSL_Int *k, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *b, EVSL_Int *ldb, EVSL_Real *beta, EVSL_Real *c__, EVSL_Int *ldc);
int evsl_dgemv(const char *trans, EVSL_Int *m, EVSL_Int *n, EVSL_Real *alpha, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *x, EVSL_Int *incx, EVSL_Real *beta, EVSL_Real *y, EVSL_Int *incy);
EVSL_Real evsl_dnrm2(EVSL_Int *n, EVSL_Real *x, EVSL_Int *incx);
int evsl_dscal(EVSL_Int *n, EVSL_Real *da, EVSL_Real *dx, EVSL_Int *incx);

#ifndef EVSL_USING_EVSL_BLAS
#ifdef __cplusplus
}
#endif
#endif

#endif
