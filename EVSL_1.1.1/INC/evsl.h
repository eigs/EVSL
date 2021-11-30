/**
 * @file evsl.h
 * @brief This file contains function prototypes, constant and macro definitions for EVSL users
*/

#ifndef EVSL_H
#define EVSL_H

#include <stdio.h>
#include <stdint.h>
#include "EVSL_config.h"

#ifdef EVSL_USING_CUDA_GPU
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse.h>
#include <curand.h>
#endif

#include "struct.h"

#define EVSL_PI 3.14159265358979323846

/** for comparing floating-point numbers (double)
 *     for "a >= b", compare "a + EVSL_DBL_EPS_MULT * DBL_EPSILON >= b"
 *     for "a <= b", compare "a - EVSL_DBL_EPS_MULT * DBL_EPSILON >= b"
 */
#define EVSL_DBL_EPS_MULT 10

/*!
  \def max(x,y)
  Computes the maximum of \a x and \a y.
*/
#define evsl_max(a, b) ((a) > (b) ? (a) : (b))

/*!
  \def min(x,y)
  Computes the minimum of \a x and \a y.
*/
#define evsl_min(a, b) ((a) < (b) ? (a) : (b))

/* memory management */
/**
 * @brief A malloc wrapper which provides basic error checking
 * @param count Number of elements to be allocated
 * @param type  Type of elements to be allocated
 */
#define evsl_Malloc(count, type) \
(\
 /* printf("[%s:%d] MALLOC %ld bytes\n", __FILE__,__LINE__, (size_t)(sizeof(type) * (count))) ,*/ \
 (type *) _evsl_Malloc( (size_t)(sizeof(type) * (count)) ) \
)

/**
 * @brief A calloc wrapper which provides basic error checking
 * @param count Number of elements to be allocated
 * @param type Type of elements to be allocated
 */
#define evsl_Calloc(count, type) \
(\
 /* printf("[%s:%d] CALLOC %ld bytes\n", __FILE__,__LINE__, (size_t)(sizeof(type) * (count))) ,*/ \
 (type *) _evsl_Calloc( (size_t)(count), (size_t)sizeof(type) ) \
)

/**
 * @brief A realloc wrapper which provides basic error checking
 * @param count Number of elements to be allocated
 * @param type Type of elements to be allocated
 */
#define evsl_Realloc(ptr, count, type) \
(\
 /* printf("[%s:%d] ReALLOC %ld bytes\n", __FILE__,__LINE__, (size_t)(sizeof(type) * (count))) ,*/ \
 (type *) _evsl_Realloc( (void *)ptr, (size_t)(sizeof(type) * (count)) ) \
)

/**
 * @brief A free wrapper which sets NULL
 * @param ptr Pointer to a memory block previously allocated with evsl_Malloc, evsl_Calloc or evsl_Realloc
 */
#define evsl_Free(ptr) \
(\
 _evsl_Free((void *)ptr), \
 ptr = NULL \
)

/**
 * @brief GPU memory management
 */

/**
 * @brief A cuda-malloc wrapper which provides basic error checking
 * @param count Number of elements to be allocated
 * @param type Type of elements to be allocated
 */
#define evsl_Malloc_device(count, type) \
(\
 /* printf("[%s:%d] DEVICE_MALLOC %ld bytes\n", __FILE__,__LINE__, (size_t)(sizeof(type) * (count))) ,*/ \
 (type *) _evsl_Malloc_device( (size_t)(sizeof(type) * (count)) ) \
)

/**
 * @brief A cuda-calloc wrapper which provides basic error checking
 * @param count Number of elements to be allocated
 * @param type Type of elements to be allocated
 */
#define evsl_Calloc_device(count, type) \
(\
 /* printf("[%s:%d] DEVICE_CALLOC %ld bytes\n", __FILE__,__LINE__, (size_t)(sizeof(type) * (count))) ,*/ \
 (type *) _evsl_Calloc_device( (size_t)(count), (size_t)sizeof(type) ) \
)

/**
 * @brief A cuda-realloc wrapper which provides basic error checking
 * @param old_count Number of elements previously allocated
 * @param new_count Number of elements to be allocated
 * @param old_type Type of elements previous allocated
 * @param new_type Type of elements to be allocated
 */
#define evsl_Realloc_device(ptr, old_count, old_type, new_count, new_type) \
(\
 /* printf("[%s:%d] DEVICE_ReALLOC %ld bytes\n", __FILE__,__LINE__, (size_t)(sizeof(type) * (count))) ,*/ \
 (new_type *) _evsl_Realloc_device( (void *)ptr, (size_t)(sizeof(old_type) * (old_count)),  \
                                                 (size_t)(sizeof(new_type) * (new_count)) ) \
)

/**
 * @brief A cuda-free wrapper which sets NULL
 * @param ptr Pointer to a memory block previously allocated with evsl_Malloc_device,
 *        evsl_Calloc_device or evsl_Realloc_device
 */
#define evsl_Free_device(ptr) \
(\
 _evsl_Free_device((void *)ptr), \
 ptr = NULL \
)

#ifdef __cplusplus
extern "C" {
#endif

/*- - - - - - - - - cheblanNr.c */
int ChebLanNr(double *intv, int maxit, double tol, double *vinit, polparams *pol, int *nevOut, double **lamo, double **Wo, double **reso, FILE *fstats);

/*- - - - - - - - - cheblanTr.c */
int ChebLanTr(int lanm, int nev, double *intv, int maxit, double tol, double *vinit, polparams *pol, int *nev2, double **vals, double **W, double **resW, FILE *fstats);

/*- - - - - - - - - chebpoly.c */
void set_pol_def(polparams *pol);
int find_pol(double *intv, polparams *pol);
void free_pol(polparams *pol);
/* int ChebAv(polparams *pol, double *v, double *y, double *w); */

/*- - - - - - - - - chebsi.c */
int ChebSI(int nev, double *intv, int maxit, double tol, double *vinit, polparams *pol, int *nevo, double **lamo, double **Yo, double **reso, FILE *fstats);

/*- - - - - - - - - - dos_utils.c */
void SetupBPol(int n, int max_deg, double tol, double lmin, double lmax, double (*ffun)(double), BSolDataPol *data);
void SetupPolRec(int n, int max_deg, double tol, double lmin, double lmax, BSolDataPol *data);
void SetupPolSqrt(int n, int max_deg, double tol, double lmin, double lmax, BSolDataPol *data);
void FreeBSolPolData(BSolDataPol *data);
void BSolPol(double *b, double *x, void *data);
void extrDiagCsr(csrMat *B, double *d);
void diagScalCsr(csrMat *A, double *d);

/*- - - - - - - - - evsl_memory.c */
void *_evsl_Malloc(size_t nbytes);
void *_evsl_Calloc(size_t count, size_t nbytes);
void *_evsl_Realloc(void *ptr, size_t nbytes);
void _evsl_Free(void *ptr);
void *_evsl_Malloc_device(size_t nbytes);
void *_evsl_Calloc_device(size_t count, size_t nbytes);
void *_evsl_Realloc_device(void *old_ptr, size_t old_nbytes, size_t new_nbytes);
void _evsl_Free_device(void *ptr);
void evsl_memcpy_device_to_host(void *dst, void *src, size_t nbytes);
void evsl_memcpy_device_to_device(void *dst, void *src, size_t nbytes);
void evsl_memcpy_host_to_device(void *dst, void *src, size_t nbytes);
void evsl_memset_device(void * ptr, int value, size_t nbytes);

/*- - - - - - - - - exDos.c */
int exDOS(double *vals, int n, int npts, double *x, double *y, double *intv);

/*- - - - - - - - - lanbounds.c */
int LanBounds(int msteps, double *v, double *lmin, double *lmax);

/*- - - - - - - - - - landos.c */
// Computes the density of states (DOS, or spectral density)
int LanDos(const int nvec, int msteps, const int npts, double *xdos, double *ydos, double *neig, const double *const intv);

//*- - - - - - - - -  landosG.c - Generalized lanDOS */
int LanDosG(const int nvec, int msteps, const int npts, double *xdos, double *ydos, double *neig, const double *const intv);

/*- - - - - - - - - lanTrbounds.c */
int LanTrbounds(int lanm, int maxit, double tol, double *vinit, int bndtype, double *lammin, double *lammax, FILE *fstats);

/*- - - -- - - - - - misc_la.c */
int scalEigVec(int n, int nev, double *Y, double* sqrtdiag);

/*- - - - - - - - - ratfilter.c */
void set_ratf_def(ratparams *rat);
int find_ratf(double *intv, ratparams *rat);
int set_ratf_solfunc(ratparams *rat, csrMat *A, csrMat *B, SolFuncC *funcs, void **data);
void free_rat(ratparams *rat);

/*- - - - - - - - - ratlanNr.c */
int RatLanNr(double *intv, int maxit, double tol, double *vinit, ratparams *rat, int *nevOut, double **lamo, double **Wo, double **reso, FILE *fstats);

/*- - - - - - - - - ratlanTr.c */
int RatLanTr(int lanm, int nev, double *intv, int maxit, double tol, double *vinit, ratparams *rat, int *nev2, double **vals, double **W, double **resW, FILE *fstats);

/*- - - - - - - - - spmat.c */
void csr_copy(csrMat *A, csrMat *B, int allocB);
// convert a COO matrix to a CSR matrix
int cooMat_to_csrMat(int cooidx, cooMat *coo, csrMat *csr);
// free a CSR
void free_csr(csrMat *csr);
// free a COO
void free_coo(cooMat *coo);
// matvec y = A*x
void matvec_csr(double *x, double *y, void *data);
/* add two csr matrices */
int matadd(double alp, double bet, csrMat *A, csrMat *B, csrMat *C, int *mapA, int *mapB);
/* sparse identity */
int speye(int n, csrMat *A);
/* extract upper triangular part of A */
void triuCsr(csrMat *A, csrMat *U);
#ifdef EVSL_USING_CUDA_GPU
void matvec_cusparse_csr(double *x, double *y, void *data);
void evsl_create_csr_gpu(csrMat *Acpu, csrMat *Agpu);
void evsl_free_csr_gpu(csrMat *csr);
#ifdef EVSL_USING_CUSPARSE_HYB
void matvec_cusparse_hyb(double *x, double *y, void *data);
int evsl_create_hybMat(csrMat *A, hybMat *hyb);
void evsl_free_hybMat(hybMat *hyb);
#endif
#endif

int arrays_copyto_csrMat(int nrow, int ncol, int *ia, int *ja, double *a, csrMat *A);
/*- - - - - - - - - evsl.c */
/* set an external matvec function */
int SetAMatvec(int n, MVFunc func, void *data);
int SetBMatvec(int n, MVFunc func, void *data);
/* unset an external matvec function */
int UnsetAMatvec();
int UnsetBMatvec();
int SetAMatrix(csrMat *A);
int SetBMatrix(csrMat *B);
int SetBSol(SolFuncR func, void *data);
int SetLTSol(SolFuncR func, void *data);
int SetASigmaBSol(ratparams *rat, int i, SolFuncC func, void *data);
int SetStdEig();
int SetGenEig();
/* start EVSL */
int EVSLStart();
/* finalize EVSL */
int EVSLFinish();
void SetDiagScal(double *ds);
#ifdef EVSL_USING_CUDA_GPU
#ifdef EVSL_USING_CUSPARSE_HYB
int SetAMatrix_device_hyb(hybMat *A);
int SetBMatrix_device_hyb(hybMat *B);
#endif
int SetAMatrix_device_csr(csrMat *A);
int SetBMatrix_device_csr(csrMat *B);
void evsl_device_query(int dev);
void evsl_last_device_err();
void evsl_pnav_device(int n, int k, double t, double s, double cc, double *vkp1, double *v_cur, double *v_old, double *y);
void evsl_chebAv_device(int n, int k, double t, double s, double cc, double *vkp1, double *vk, double *vkm1, double *y);
#endif

/*- - - - - - - - - spslicer.c */
int spslicer(double *sli, double *mu, int Mdeg, double *intv, int n_int, int npts);
int kpmdos(int Mdeg, int damping, int nvec, double *ab, double *mu, double *ecnt);

/*- - - - - - - - - - spslicer2.c */
void spslicer2(double *xi, double *yi, int n_int, int npts, double *sli);

/*- - - - - - - - - timing.c */
double evsl_timer();

/*- - - - - - - - - vect.c */
// generate a random vector
void rand_double(int n, double *v);
void rand_double_device(int n, double *v);
// generate a normally distributed random vector
void randn_double(int n, double *v);
void randn_double_device(int n, double *v);
// sort a vector
void sort_double(int n, double *v, int *ind);
void linspace(double a, double b, int num, double *arr);

double evsl_dnrm2_device(int *n, double *x, int *incr);
void evsl_dscal_device(int *n, double *a, double *x, int *incr);
double evsl_ddot_device(int *n, double *x, int *incx, double *y, int *incy);
void evsl_daxpy_device(int *n, double *a, double *x, int *incx, double *y, int *incy);
void evsl_dcopy_device(int *n, double *x, int *incx, double *y, int *incy);
void evsl_dgemv_device(const char *trans, int *m, int *n, double *alpha, double *a, int *lda,
                       double *x, int *incx, double *beta, double *y, int *incy);
void evsl_dgemm_device(const char *transa, const char *transb, int *m, int * n, int *k, 
                       double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void evsl_element_mult_device(int n, double *a, double *b);
void evsl_element_divide_device(int n, double *a, double *b);

/*- - - - - - - - - stats.c */
void StatsPrint(FILE *fstats);
void StatsReset();

/* evsl_f90.c: */
void EVSLFORT(evsl_start,EVSL_START)();
void EVSLFORT(evsl_finish,EVSL_FINISH)();
void EVSLFORT(evsl_coo2csr,EVSL_COO2CSR)(int *n, int *nnz, int *ir, int *jc, double *vv, uintptr_t *csrf90);
void EVSLFORT(evsl_arr2csr,EVSL_ARR2CSR)(int *n, int *ia, int *ja, double *a, uintptr_t *csrf90);
void EVSLFORT(evsl_arr2devicecsr,EVSL_ARR2DEVICECSR)(int *n, int *ia, int *ja, double *a, uintptr_t *csrf90);
void EVSLFORT(evsl_free_csr,EVSL_FREE_CSR)(uintptr_t *csrf90);
void EVSLFORT(evsl_free_devicecsr,EVSL_FREE_DEVICECSR)(uintptr_t *csrf90);
void EVSLFORT(evsl_seta_csr,EVSL_SETA_CSR)(uintptr_t *Af90);
void EVSLFORT(evsl_seta_devicecsr,EVSL_SETA_DEVICECSR)(uintptr_t *Af90);
void EVSLFORT(evsl_setb_csr,EVSL_SETB_CSR)(uintptr_t *Bf90);
void EVSLFORT(evsl_setb_devicecsr,EVSL_SETB_DEVICECSR)(uintptr_t *Bf90);
void EVSLFORT(evsl_setamv,EVSL_SETAMV)(int *n, void *func, void *data);
void EVSLFORT(evsl_setbmv,EVSL_SETBMV)(int *n, void *func, void *data);
void EVSLFORT(evsl_setbsol,EVSL_SETBSOL)(void *func, void *data);
void EVSLFORT(evsl_setltsol,EVSL_SETLTSOL)(void *func, void *data);
void EVSLFORT(set_asigmabsol,SET_ASIGMABSOL)(uintptr_t *ratf90, int *i, void *func, void *data);
void EVSLFORT(evsl_set_geneig,EVSL_SET_GENEIG)();
void EVSLFORT(evsl_lanbounds,EVSL_LANBOUNDS)(int *nsteps, double *lmin, double *lmax);
void EVSLFORT(evsl_kpm_spslicer,EVSL_KPM_SPSLICER)(int *Mdeg, int *nvec, double *xintv, int *nslices, double *sli, int *evint);
void EVSLFORT(evsl_find_pol,EVSL_FIND_POL)(double *xintv, double *thresh_int, double *thresh_ext, uintptr_t *polf90);
void EVSLFORT(evsl_free_pol,EVSL_FREE_POL)(uintptr_t *polf90);
void EVSLFORT(evsl_find_rat,EVSL_FIND_RAT)(int *defrat, int *num, int *pw, int *method, double *beta, double *bar, double *intv, uintptr_t *ratf90);
void EVSLFORT(evsl_free_rat,EVSL_FREE_RAT)(uintptr_t *ratf90);
void EVSLFORT(evsl_cheblantr,EVSL_CHEBLANTR)(int *mlan, int *nev, double *xintv, int *max_its, double *tol, uintptr_t *polf90);
void EVSLFORT(evsl_cheblannr,EVSL_CHEBLANNR)(double *xintv, int *max_its, double *tol, uintptr_t *polf90);
void EVSLFORT(evsl_ratlannr,EVSL_RATLANNR)(double *xintv, int *max_its, double *tol, uintptr_t *ratf90);
void EVSLFORT(evsl_ratlantr,EVSL_RATLANTR)(int *lanm, int *nev, double *xintv, int *max_its, double *tol, uintptr_t *ratf90);
void EVSLFORT(evsl_get_nev,EVSL_GET_NEV)(int *nev);
void EVSLFORT(evsl_copy_result,EVSL_COPY_RESULT)(double *val, double *vec);

#ifdef __cplusplus
}
#endif

#endif
