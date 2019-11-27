#include <string.h>
#include "internal_header.h"

/**
 * @file evsl_f90.c
 * @brief Definitions used by fortran interface
 */
/** global variables that hold results from EVSL
 * evsl_copy_result_f90 will copy results from these vars and reset them
 * */
int evsl_nev_computed=0, evsl_n=0;
double *evsl_eigval_computed=NULL, *evsl_eigvec_computed=NULL;

/** @brief Fortran interface for EVSLStart */
void EVSLFORT(evsl_start,EVSL_START)() {
  EVSLStart();
}

/** @brief Fortran interface for EVSLFinish */
void EVSLFORT(evsl_finish,EVSL_FINISH)() {
  EVSLFinish();
}

/** @brief Fortran interface to convert a COO matrix to CSR
 * the pointer of CSR will be returned
 * @param[in] *n   : size of A
 * @param[in] *nnz : nnz of A
 * @param[in] *ir, *jc, *vv : COO triplets
 * @param[out] csrf90 : CSR pointer
 */
void EVSLFORT(evsl_coo2csr,EVSL_COO2CSR)(int *n, int *nnz, int *ir, int *jc,
                                         double *vv, uintptr_t *csrf90) {
  cooMat coo;
  coo.nrows = *n;
  coo.ncols = *n;
  coo.nnz = *nnz;
  coo.ir = ir;
  coo.jc = jc;
  coo.vv = vv;
  csrMat *csr;
  csr = evsl_Malloc(1, csrMat);
  cooMat_to_csrMat(1, &coo, csr);
  *csrf90 = (uintptr_t) csr;
}

/** @brief Fortran interface to return a CSR struct from (ia,ja,a)
 * the pointer of CSR will be returned
 * @warning : no memory allocation inside the CSR
 * @param[in] *n   : size of A
 * @param[in] *ia, *ja, *a : CSR array
 * @param[out] csrf90 : CSR pointer
 */
void EVSLFORT(evsl_arr2csr,EVSL_ARR2CSR)(int *n, int *ia, int *ja,
                                         double *a, uintptr_t *csrf90) {
  csrMat *csr;
  csr = evsl_Malloc(1, csrMat);
  /* does not own the data */
  csr->owndata = 0;
  csr->nrows = *n;
  csr->ncols = *n;
  csr->ia = ia;
  csr->ja = ja;
  csr->a = a;
  *csrf90 = (uintptr_t) csr;
}

void EVSLFORT(evsl_arr2devicecsr,EVSL_ARR2DEVICECSR)(int *n, int *ia, int *ja,
                                                     double *a, uintptr_t *csrf90) {
  csrMat *csr;
  csr = evsl_Malloc(1, csrMat);
  /* does not own the data */
  csr->owndata = 0;
  csr->nrows = *n;
  csr->ncols = *n;
  csr->ia = ia;
  csr->ja = ja;
  csr->a = a;
#ifdef EVSL_USING_CUDA_GPU
  csrMat *csr_gpu;
  csr_gpu = evsl_Malloc(1, csrMat);
  evsl_create_csr_gpu(csr, csr_gpu);
  free_csr(csr);
  evsl_Free(csr);
  *csrf90 = (uintptr_t) csr_gpu;
#else
  *csrf90 = (uintptr_t) csr;
#endif
}

/** @brief Fortran interface to free a CSR matrix
 * @param[in] csrf90 : CSR pointer
 */
void EVSLFORT(evsl_free_csr,EVSL_FREE_CSR)(uintptr_t *csrf90) {
  csrMat *csr = (csrMat *) (*csrf90);
  free_csr(csr);
  evsl_Free(csr);
}

void EVSLFORT(evsl_free_devicecsr,EVSL_FREE_DEVICECSR)(uintptr_t *csrf90) {
  csrMat *csr = (csrMat *) (*csrf90);
#ifdef EVSL_USING_CUDA_GPU
  evsl_free_csr_gpu(csr);
#else
  free_csr(csr);
#endif
  evsl_Free(csr);
}

/** @brief Fortran interface to set matrix A from a CSR matrix
 * @param[in] Af90 : CSR pointer of A
 */
void EVSLFORT(evsl_seta_csr,EVSL_SETA_CSR)(uintptr_t *Af90) {
  csrMat *A = (csrMat *) (*Af90);
  SetAMatrix(A);
}

void EVSLFORT(evsl_seta_devicecsr,EVSL_SETA_DEVICECSR)(uintptr_t *Af90) {
  csrMat *A = (csrMat *) (*Af90);
#ifdef EVSL_USING_CUDA_GPU
  SetAMatrix_device_csr(A);
#else
  SetAMatrix(A);
#endif
}

/** @brief Fortran interface to set matrix B from a CSR matrix
 * @param[in] Bf90 : CSR pointer of B
 */
void EVSLFORT(evsl_setb_csr,EVSL_SETB_CSR)(uintptr_t *Bf90) {
  csrMat *B = (csrMat *) (*Bf90);
  SetBMatrix(B);
}

void EVSLFORT(evsl_setb_devicecsr,EVSL_SETB_DEVICECSR)(uintptr_t *Bf90) {
  csrMat *B = (csrMat *) (*Bf90);
#ifdef EVSL_USING_CUDA_GPU
  SetBMatrix_device_csr(B);
#else
  SetBMatrix(B);
#endif
}

/** @brief Fortran interface for SetAMatvec
 * @param[in] n : size of A
 * @param[in] func : function pointer
 * @param[in] data : associated data
 */
void EVSLFORT(evsl_setamv,EVSL_SETAMV)(int *n, void *func, void *data) {
  SetAMatvec(*n, (MVFunc) func, data);
}

/** @brief Fortran interface for SetBMatvec
 * @param[in] n : size of B
 * @param[in] func : function pointer
 * @param[in] data : associated data
 */
void EVSLFORT(evsl_setbmv,EVSL_SETBMV)(int *n, void *func, void *data) {
  SetBMatvec(*n, (MVFunc) func, data);
}

/** @brief Fortran interface for SetBSol
 * @param[in] func: func pointer of Bsol
 * @param[in] data: data pointer of Bsol
 */
void EVSLFORT(evsl_setbsol,EVSL_SETBSOL)(void *func, void *data) {
  SetBSol((SolFuncR) func, data);
}

/** @brief Fortran interface for SetLTSol
 * @param[in] func: func pointer of LTsol
 * @param[in] data: data pointer of LTsol
 */
void EVSLFORT(evsl_setltsol,EVSL_SETLTSOL)(void *func, void *data) {
  SetLTSol((SolFuncR) func, data);
}

/** @brief Fortran interface for SetASigmaBSol
 * @param[in,out] ratf90: pointer of rational filter
 * @param[in]           i: pole index
 * @param[in]        func: func pointer of A-siB
 * @param[in]        data: data pointer of A-siB
 */
void EVSLFORT(set_asigmabsol,SET_ASIGMABSOL)(uintptr_t *ratf90,
                                             int       *i,
                                             void      *func,
                                             void      *data) {
  /* cast pointers */
  ratparams *rat = (ratparams *) (*ratf90);
  SetASigmaBSol(rat, *i, (SolFuncC) func, data);
}

/** @brief Fortran interface for SetGenEig */
void EVSLFORT(evsl_set_geneig,EVSL_SET_GENEIG)() {
  SetGenEig();
}

/** @brief Fortran interface for evsl_lanbounds
 * @param[in] nsteps: number of steps
 * @param[out] lmin: lower bound
 * @param[out] lmax: upper bound
 * */
void EVSLFORT(evsl_lanbounds,EVSL_LANBOUNDS)(int *nsteps, double *lmin, double *lmax) {
  int n = evsldata.n;
  double *vinit;
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);
  LanTrbounds(50, *nsteps, 1e-10, vinit, 1, lmin, lmax, NULL);
  //LanBounds(*nsteps, vinit, lmin, lmax);
  evsl_Free_device(vinit);
}

/** @brief Fortran interface for kpmdos and spslicer
 * @param[in] Mdeg: degree
 * @param[in] nvec: number of vectors to use
 * @param[in] xintv:   an array of length 4  \n
 *                 [intv[0] intv[1]] is the interval of desired eigenvalues
 *                 that must be cut (sliced) into n_int  sub-intervals \n
 *                 [intv[2],intv[3]] is the global interval of eigenvalues
 *                 it must contain all eigenvalues of A \n
 * @param[in] nslices: number of slices
 * @param[out] sli: slices [of size nslices+1]
 * @param[out] evint: estimated ev count per slice
 */
void EVSLFORT(evsl_kpm_spslicer,EVSL_KPM_SPSLICER)(int *Mdeg, int *nvec, double *xintv,
                                                   int *nslices, double *sli, int *evint) {
  int npts;
  double *mu, ecount;
  mu = evsl_Malloc(*Mdeg+1, double);
  kpmdos(*Mdeg, 1, *nvec, xintv, mu, &ecount);
  fprintf(stdout, " estimated eig count in interval: %.15e \n",ecount);

  npts = 10 *ecount;
  int ierr = spslicer(sli, mu, *Mdeg, xintv, *nslices, npts);
  if (ierr) {
    printf("spslicer error %d\n", ierr);
    exit(1);
  }
  *evint = (int) (1 + ecount / (*nslices));
  evsl_Free(mu);
}

/** @brief Fortran interface for find_pol
 * @param[in] xintv Intervals of interest
 * @param[in]  thresh_int Threshold for accepting interior values
 * @param[in]  thresh_ext Threshold for accepting exterior values
 * @param[out] polf90 : pointer of pol
 * @warning: The pointer will be cast to uintptr_t
 *
 * uintptr_t: Integer type capable of holding a value converted from
 * a void pointer and then be converted back to that type with a value
 * that compares equal to the original pointer
 */
void EVSLFORT(evsl_find_pol,EVSL_FIND_POL)(double *xintv, double *thresh_int,
                                           double *thresh_ext, uintptr_t *polf90) {
  polparams *pol;
  pol = evsl_Malloc(1, polparams);
  set_pol_def(pol);
  pol->damping = 2;
  pol->thresh_int = *thresh_int;
  pol->thresh_ext = *thresh_ext;
  pol->max_deg  = 300;
  find_pol(xintv, pol);

  fprintf(stdout, " polynomial deg %d, bar %e gam %e\n",
          pol->deg,pol->bar, pol->gam);

  *polf90 = (uintptr_t) pol;
}

/** @brief Fortran interface for free_pol
 * @param[out] polf90 : pointer of pol
 *
 * */
void EVSLFORT(evsl_free_pol,EVSL_FREE_POL)(uintptr_t *polf90) {
  /* cast pointer */
  polparams *pol = (polparams *) (*polf90);
  free_pol(pol);
  evsl_Free(pol);
}

/** @brief Fortran interface for find_rat
 * @param[in]  defrat: if true (1) use default rational filter, input parameters are ignored
 * @param[in]  num, pw, method, beta, bar: rational filter parameters
 * @param[in]  intv: interval of interest
 * @param[out] ratf90 : pointer of rat
 * @warning: The pointer will be cast to uintptr_t
 * uintptr_t: Integer type capable of holding a value converted from
 * a void pointer and then be converted back to that type with a value
 * that compares equal to the original pointer
 */
void EVSLFORT(evsl_find_rat,EVSL_FIND_RAT)(int *defrat, int *num, int *pw, int *method, double *beta,
                                           double *bar, double *intv, uintptr_t *ratf90) {
  ratparams *rat;
  rat = evsl_Malloc(1, ratparams);
  set_ratf_def(rat);
  if (*defrat == 0) {
     set_ratf_def(rat);
     rat->num    = *num;
     rat->pw     = *pw;
     rat->method = *method;
     rat->beta   = *beta;
     rat->bar    = *bar;
  }
  find_ratf(intv, rat);

  *ratf90 = (uintptr_t) rat;
}

/** @brief Fortran interface for free_rat */
void EVSLFORT(evsl_free_rat,EVSL_FREE_RAT)(uintptr_t *ratf90) {
  ratparams *rat = (ratparams *) (*ratf90);
  free_rat(rat);
  evsl_Free(rat);
}

/** @brief Fortran interface for ChebLanTr
 *  the results will be saved in the internal variables
 */
void EVSLFORT(evsl_cheblantr,EVSL_CHEBLANTR)(int *mlan, int *nev, double *xintv, int *max_its,
                                             double *tol, uintptr_t *polf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;

  n = evsldata.n;
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);

  /* cast pointer */
  polparams *pol = (polparams *) (*polf90);

  ierr = ChebLanTr(*mlan, *nev, xintv, *max_its, *tol, vinit,
                   pol, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("ChebLanTr error %d\n", ierr);
  }

  evsl_Free_device(vinit);
  if (res) {
    evsl_Free(res);
  }
  /* save pointers to the global variables */
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

/** @brief Fortran interface for ChebLanNr
 *  the results will be saved in the internal variables
 */
void EVSLFORT(evsl_cheblannr,EVSL_CHEBLANNR)(double *xintv, int *max_its, double *tol,
                                             uintptr_t *polf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;

  n = evsldata.n;
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);

  /* cast pointer */
  polparams *pol = (polparams *) (*polf90);

  ierr = ChebLanNr(xintv, *max_its, *tol, vinit,
                   pol, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("ChebLanNr error %d\n", ierr);
  }

  evsl_Free_device(vinit);
  if (res) {
    evsl_Free(res);
  }
  /* save pointers to the global variables */
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

/** @brief Fortran interface for RatLanNr
 *  the results will be saved in the internal variables
 */
void EVSLFORT(evsl_ratlannr,EVSL_RATLANNR)(double *xintv, int *max_its, double *tol,
                                           uintptr_t *ratf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;

  n = evsldata.n;
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);

  /* cast pointer */
  ratparams *rat = (ratparams *) (*ratf90);

  ierr = RatLanNr(xintv, *max_its, *tol, vinit,
                  rat, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("RatLanNr error %d\n", ierr);
  }

  evsl_Free_device(vinit);
  if (res) {
    evsl_Free(res);
  }
  /* save pointers to the global variables */
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

/** @brief Fortran interface for RatLanNr
 *  the results will be saved in the internal variables
 */
void EVSLFORT(evsl_ratlantr,EVSL_RATLANTR)(int *lanm, int *nev, double *xintv,
                                           int *max_its, double *tol, uintptr_t *ratf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;

  n = evsldata.n;
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);

  /* cast pointer */
  ratparams *rat = (ratparams *) (*ratf90);

  ierr = RatLanTr(*lanm, *nev, xintv, *max_its, *tol, vinit,
                  rat, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("RatLanNr error %d\n", ierr);
  }

  evsl_Free_device(vinit);
  if (res) {
    evsl_Free(res);
  }
  /* save pointers to the global variables */
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

/** @brief Get the number of last computed eigenvalues
 */
void EVSLFORT(evsl_get_nev,EVSL_GET_NEV)(int *nev) {
  *nev = evsl_nev_computed;
}

/** @brief copy the computed eigenvalues and vectors
 * @warning: after this call the internal saved results will be freed
 */
void EVSLFORT(evsl_copy_result,EVSL_COPY_RESULT)(double *val, double *vec) {
  memcpy(val, evsl_eigval_computed, evsl_nev_computed*sizeof(double));
  evsl_memcpy_device_to_host(vec, evsl_eigvec_computed, evsl_nev_computed*evsl_n*sizeof(double));
  /* reset global variables */
  evsl_nev_computed = 0;
  evsl_n = 0;
  evsl_Free(evsl_eigval_computed);
  evsl_Free_device(evsl_eigvec_computed);
  evsl_eigval_computed = NULL;
  evsl_eigvec_computed = NULL;
}

