#ifndef BLAS_F2C_INCLUDE
#define BLAS_F2C_INCLUDE

/* do not use this file outside blas or lapack directory */

#include "EVSL_config.h"
#include <stdio.h>
#include <complex.h>
#include <math.h>

#define integer       EVSL_Int
#define doublereal    EVSL_Real
#define real          float
/* #define doublecomplex EVSL_Complex */
typedef struct { doublereal r, i; } doublecomplex; /* this must be compatible with EVSL_Complex */
#define logical       EVSL_Int
#define flag          EVSL_LongInt
#define ftnlen        EVSL_LongInt
#define ftnint        EVSL_LongInt

#define TRUE_  (1)
#define FALSE_ (0)

#define VOID void

/*external read, write*/
typedef struct
{
  flag cierr;
  ftnint ciunit;
  flag ciend;
  char *cifmt;
  ftnint cirec;
} cilist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#if !defined(LAPACK_F2C_INCLUDE) || defined(EVSL_USING_EVSL_BLAS)
/* 1. for evsl blas:
 *      add prefix evsl_ to all blas subroutines,
 *      so can live together with external blas
 * 2. for evsl lapack:
 *      if using evsl blas, cast name to evsl_name
 *      so that evsl lapack can link external blas
 */
#define daxpy_  evsl_daxpy
#define dcabs1_ evsl_dcabs1
#define dcopy_  evsl_dcopy
#define ddot_   evsl_ddot
#define dgemm_  evsl_dgemm
#define dgemv_  evsl_dgemv
#define dger_   evsl_dger
#define dnrm2_  evsl_dnrm2
#define dscal_  evsl_dscal
#define dswap_  evsl_dswap
#define dsymv_  evsl_dsymv
#define dsyr2_  evsl_dsyr2
#define dsyr2k_ evsl_dsyr2k
#define dtrmm_  evsl_dtrmm
#define dtrmv_  evsl_dtrmv
#define izamax_ evsl_izamax
#define lsame_  evsl_lsame
#define xerbla_ evsl_xerbla
#define zgemm_  evsl_zgemm
#define zgeru_  evsl_zgeru
#define zscal_  evsl_zscal
#define zswap_  evsl_zswap
#define ztrsm_  evsl_ztrsm
#endif

/* f2c */
#define f__cabs evsl_f__cabs
#define d_cnjg  evsl_d_cnjg
#define d_imag  evsl_d_imag
#define d_sign  evsl_d_sign
#define i_nint  evsl_i_nint
#define pow_di  evsl_pow_di
#define s_cmp   evsl_s_cmp
#define s_copy  evsl_s_copy
#define sig_die evsl_sig_die
#define z_abs   evsl_z_abs
#define z_div   evsl_z_div

#endif
