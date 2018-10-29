#ifndef EVSL_LAPACK_H
#define EVSL_LAPACK_H

/* use this file outside lapack directory to provide headers */

#ifndef EVSL_USING_EVSL_LAPACK
#define evsl_dstev  dstev_
#define evsl_dsyev  dsyev_
#define evsl_dstemr dstemr_
#define evsl_zgesv  zgesv_
#ifdef __cplusplus
extern "C" {
#endif
#endif

int evsl_dstev(const char *jobz, EVSL_Int *n, EVSL_Real *d__, EVSL_Real *e, EVSL_Real *z__, EVSL_Int *ldz, EVSL_Real *work, EVSL_Int *info);
int evsl_dsyev(const char *jobz, const char *uplo, EVSL_Int *n, EVSL_Real *a, EVSL_Int *lda, EVSL_Real *w, EVSL_Real *work, EVSL_Int *lwork, EVSL_Int *info);
int evsl_dstemr(const char *jobz, const char *range, EVSL_Int *n, EVSL_Real *d__, EVSL_Real *e, EVSL_Real *vl, EVSL_Real *vu, EVSL_Int *il, EVSL_Int *iu, EVSL_Int *m, EVSL_Real *w, EVSL_Real *z__, EVSL_Int *ldz, EVSL_Int *nzc, EVSL_Int *isuppz, EVSL_Int *tryrac, EVSL_Real *work, EVSL_Int *lwork, EVSL_Int *iwork, EVSL_Int *liwork, EVSL_Int *info);
int evsl_zgesv(EVSL_Int *n, EVSL_Int *nrhs, EVSL_Complex *a, EVSL_Int *lda, EVSL_Int *ipiv, EVSL_Complex *b, EVSL_Int *ldb, EVSL_Int *info);

#ifndef EVSL_USING_EVSL_LAPACK
#ifdef __cplusplus
}
#endif
#endif

#endif

