#ifndef LAPACK_F2C_INCLUDE
#define LAPACK_F2C_INCLUDE

#include "../blas/f2c.h"

#ifndef EVSL_USING_EVSL_BLAS
/* have these so that evsl's lapack can link external blas */
#define evsl_daxpy daxpy_
#define evsl_dcabs1 dcabs1_
#define evsl_dcopy dcopy_
#define evsl_ddot ddot_
#define evsl_dgemm dgemm_
#define evsl_dgemv dgemv_
#define evsl_dger dger_
#define evsl_dnrm2 dnrm2_
#define evsl_dscal dscal_
#define evsl_dswap dswap_
#define evsl_dsymv dsymv_
#define evsl_dsyr2 dsyr2_
#define evsl_dsyr2k dsyr2k_
#define evsl_dtrmm dtrmm_
#define evsl_dtrmv dtrmv_
#define evsl_izamax izamax_
#define evsl_lsame lsame_
#define evsl_xerbla xerbla_
#define evsl_zgemm zgemm_
#define evsl_zgeru zgeru_
#define evsl_zscal zscal_
#define evsl_zswap zswap_
#define evsl_ztrsm ztrsm_
#endif

/* add prefix evsl_ to all lapack subroutines,
 * so can live together with external lapack */
#define disnan_   evsl_disnan
#define dlae2_    evsl_dlae2
#define dlaebz_   evsl_dlaebz
#define dlaev2_   evsl_dlaev2
#define dlaisnan_ evsl_dlaisnan
#define dlamch_   evsl_dlamch
#define dlaneg_   evsl_dlaneg
#define dlanst_   evsl_dlanst
#define dlansy_   evsl_dlansy
#define dlapy2_   evsl_dlapy2
#define dlar1v_   evsl_dlar1v
#define dlarfb_   evsl_dlarfb
#define dlarf_    evsl_dlarf
#define dlarfg_   evsl_dlarfg
#define dlarft_   evsl_dlarft
#define dlarnv_   evsl_dlarnv
#define dlarra_   evsl_dlarra
#define dlarrb_   evsl_dlarrb
#define dlarrc_   evsl_dlarrc
#define dlarrd_   evsl_dlarrd
#define dlarre_   evsl_dlarre
#define dlarrf_   evsl_dlarrf
#define dlarrj_   evsl_dlarrj
#define dlarrk_   evsl_dlarrk
#define dlarrr_   evsl_dlarrr
#define dlarrv_   evsl_dlarrv
#define dlartg_   evsl_dlartg
#define dlaruv_   evsl_dlaruv
#define dlascl_   evsl_dlascl
#define dlaset_   evsl_dlaset
#define dlasq2_   evsl_dlasq2
#define dlasq5_   evsl_dlasq5
#define dlasq6_   evsl_dlasq6
#define dlasr_    evsl_dlasr
#define dlasrt_   evsl_dlasrt
#define dlassq_   evsl_dlassq
#define dlatrd_   evsl_dlatrd
#define dlazq3_   evsl_dlazq3
#define dlazq4_   evsl_dlazq4
#define dorg2l_   evsl_dorg2l
#define dorg2r_   evsl_dorg2r
#define dorgql_   evsl_dorgql
#define dorgqr_   evsl_dorgqr
#define dorgtr_   evsl_dorgtr
#define dstemr_   evsl_dstemr
#define dsteqr_   evsl_dsteqr
#define dsterf_   evsl_dsterf
#define dstev_    evsl_dstev
#define dsyev_    evsl_dsyev
#define dsytd2_   evsl_dsytd2
#define dsytrd_   evsl_dsytrd
#define ieeeck_   evsl_ieeeck
#define ilaenv_   evsl_ilaenv
#define iparmq_   evsl_iparmq
#define zgesv_    evsl_zgesv
#define zgetf2_   evsl_zgetf2
#define zgetrf_   evsl_zgetrf
#define zgetrs_   evsl_zgetrs
#define zlaswp_   evsl_zlaswp

#endif
