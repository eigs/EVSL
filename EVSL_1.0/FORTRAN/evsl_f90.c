#include <string.h>
#include <stdint.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/** global variables that hold results from EVSL
 * evsl_copy_result_f90 will copy results from these vars and reset them
 * */
int evsl_nev_computed=0, evsl_n=0;
double *evsl_eigval_computed=NULL, *evsl_eigvec_computed=NULL;

/** @brief Fortran interface for EVSLStart */
void EVSLFORT(evsl_start)() {
  EVSLStart();
}

/** @brief Fortran interface for EVSLFinish */
void EVSLFORT(evsl_finish)() {
  EVSLFinish();
}

/** @brief Fortran interface to set matrix A as a COO matrix
 * The matrix will be converted to CSR and given to EVSL 
 * EVSL will take care of deallocation of the CSR matrix
 * @param[in] *n   : size of A
 * @param[in] *nnz : nnz of A
 * @param[in] *ir, *jc, *vv : COO triplets
 */ 
void EVSLFORT(evsl_seta_coo)(int *n, int *nnz, int *ir, int *jc, double *vv) {
  cooMat coo;
  coo.nrows = *n;
  coo.ncols = *n;
  coo.nnz = *nnz;
  coo.ir = ir;
  coo.jc = jc;
  coo.vv = vv;
  csrMat *csr;
  Malloc(csr, 1, csrMat);
  cooMat_to_csrMat(1, &coo, csr);
  evsldata.A = csr;
  /* evsl owns this csr and will free it after done */
  evsldata.ifOwnA = 1;
}

/** @brief Fortran interface to set matrix B as a COO matrix
 * The matrix will be converted to CSR and given to EVSL 
 * EVSL will take care of deallocation of the CSR matrix
 * @param[in] n   : size of B
 * @param[in] nnz : nnz of B
 * @param[in] ir, jc, vv : COO triplets
 */ 
void EVSLFORT(evsl_setb_coo)(int *n, int *nnz, int *ir, int *jc, double *vv) {
  cooMat coo;
  coo.nrows = *n;
  coo.ncols = *n;
  coo.nnz = *nnz;
  coo.ir = ir;
  coo.jc = jc;
  coo.vv = vv;
  csrMat *csr;
  Malloc(csr, 1, csrMat);
  cooMat_to_csrMat(1, &coo, csr);
  evsldata.B = csr;
  /* evsl owns this csr and will free it after done */
  evsldata.ifOwnB = 1;
}

/** @brief Fortran interface for SetAMatvec 
 * @param[in] n : size of A
 * @param[in] func : function pointer 
 * @param[in] data : associated data
 */
void EVSLFORT(evsl_setamv)(int *n, void *func, void *data) {
  SetAMatvec(*n, (MVFunc) func, data);
}

/** @brief Fortran interface for SetBMatvec 
 * @param[in] n : size of B
 * @param[in] func : function pointer 
 * @param[in] data : associated data
 */
void EVSLFORT(evsl_setbmv)(int *n, void *func, void *data) {
  SetBMatvec(*n, (MVFunc) func, data);
}

/** @brief Fortran interface for SetBSol
 * @param[in] func: func pointer of Bsol
 * @param[in] data: data pointer of Bsol
 */
void EVSLFORT(evsl_setbsol)(void *func, void *data) {
  SetBSol((SolFuncR) func, data);
}

/** @brief Fortran interface for SetLTSol
 * @param[in] func: func pointer of LTsol
 * @param[in] data: data pointer of LTsol
 */
void EVSLFORT(evsl_setltsol)(void *func, void *data) {
  SetLTSol((SolFuncR) func, data);
}

/** @brief Fortran interface for SetGenEig */
void EVSLFORT(evsl_set_geneig)() {
  SetGenEig();
}

/** @brief Fortran interface for evsl_lanbounds 
 * @param[in] nstpes: number of steps
 * @param[out] lmin: lower bound
 * @param[out] lmax: upper bound
 * */
void EVSLFORT(evsl_lanbounds)(int *nsteps, double *lmin, double *lmax) {
  int n;
  double *vinit;
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  Malloc(vinit, n, double);
  rand_double(n, vinit);
  LanTrbounds(50, *nsteps, 1e-10, vinit, 1, lmin, lmax, NULL);
  //LanBounds(*nsteps, vinit, lmin, lmax);
  free(vinit);
}

/** @brief Fortran interface for kpmdos and spslicer
 * @param[in] Mdeg: degree
 * @param[in] nvec: number of vectors to use
 * @param[in] intv:   an array of length 4  \n
 *                 [intv[0] intv[1]] is the interval of desired eigenvalues 
 *                 that must be cut (sliced) into n_int  sub-intervals \n
 *                 [intv[2],intv[3]] is the global interval of eigenvalues 
 *                 it must contain all eigenvalues of A \n
 * @param[in] nslices: number of slices
 * @param[out] sli: slices [of size nslices+1]
 * @param[out] evint: estimated ev count per slice
 */ 
void EVSLFORT(evsl_kpm_spslicer)(int *Mdeg, int *nvec, double *xintv,
                                int *nslices, double *sli, int *evint) {
  int npts;
  double *mu, ecount;
  Malloc(mu, *Mdeg+1, double);
  kpmdos(*Mdeg, 1, *nvec, xintv, mu, &ecount);
  fprintf(stdout, " estimated eig count in interval: %.15e \n",ecount);
  
  npts = 10 *ecount;
  int ierr = spslicer(sli, mu, *Mdeg, xintv, *nslices, npts);
  if (ierr) {
    printf("spslicer error %d\n", ierr);
    exit(1);
  }
  *evint = (int) (1 + ecount / (*nslices));
  free(mu);
}

/** @brief Fortran interface for find_pol 
 * @param[out] polf90 : pointer of pol
 * @warning: The pointer will be cast to uintptr_t
 * 
 * uintptr_t: Integer type capable of holding a value converted from 
 * a void pointer and then be converted back to that type with a value 
 * that compares equal to the original pointer
 */
void EVSLFORT(evsl_find_pol)(double *xintv, double *thresh_int, 
                             double *thresh_ext, uintptr_t *polf90) {
  polparams *pol;
  Malloc(pol, 1, polparams);
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

/** @brief Fortran interface for free_pol */
void EVSLFORT(evsl_free_pol)(uintptr_t *polf90) {
  /* cast pointer */
  polparams *pol = (polparams *) (*polf90);
  free_pol(pol);
  free(pol);
}

/** @brief Fortran interface for find_rat
 * @param[out] ratf90 : pointer of rat
 * @warning: The pointer will be cast to uintptr_t
 * 
 * uintptr_t: Integer type capable of holding a value converted from 
 * a void pointer and then be converted back to that type with a value 
 * that compares equal to the original pointer
 */
void EVSLFORT(evsl_find_rat)(double *intv, uintptr_t *ratf90) {
  int pow = 2;
  int num = 1;
  double beta = 0.01;
  ratparams *rat;
  Malloc(rat, 1, ratparams);
  set_ratf_def(rat);
  rat->pw = pow;
  rat->num = num;
  rat->beta = beta;
  find_ratf(intv, rat);

  *ratf90 = (uintptr_t) rat;
}

/** @brief Fortran interface for free_rat */
void EVSLFORT(evsl_free_rat)(uintptr_t *ratf90) {
  ratparams *rat = (ratparams *) (*ratf90);
  free_rat(rat);
  free(rat);
}

/** @brief Fortran interface for ChebLanTr
 *  the results will be saved in the internal variables
 */
void EVSLFORT(evsl_cheblantr)(int *mlan, int *nev, double *xintv, int *max_its,
                              double *tol, uintptr_t *polf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;
 
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  Malloc(vinit, n, double);
  rand_double(n, vinit);

  /* cast pointer */
  polparams *pol = (polparams *) (*polf90);

  ierr = ChebLanTr(*mlan, *nev, xintv, *max_its, *tol, vinit,
                   pol, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("ChebLanTr error %d\n", ierr);
  }

  free(vinit);
  if (res) {
    free(res);
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
void EVSLFORT(evsl_cheblannr)(double *xintv, int *max_its, double *tol, 
                              uintptr_t *polf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;
 
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  Malloc(vinit, n, double);
  rand_double(n, vinit);

  /* cast pointer */
  polparams *pol = (polparams *) (*polf90);

  ierr = ChebLanNr(xintv, *max_its, *tol, vinit,
                   pol, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("ChebLanNr error %d\n", ierr);
  }

  free(vinit);
  if (res) {
    free(res);
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
void EVSLFORT(evsl_ratlannr)(double *xintv, int *max_its, double *tol, 
                             uintptr_t *ratf90) {
  int n, nev2, ierr;
  double *lam, *Y, *res;
  FILE *fstats = stdout;
  double *vinit;
 
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  Malloc(vinit, n, double);
  rand_double(n, vinit);

  /* cast pointer */
  ratparams *rat = (ratparams *) (*ratf90);

  ierr = RatLanNr(xintv, *max_its, *tol, vinit,
                  rat, &nev2, &lam, &Y, &res, fstats);

  if (ierr) {
    printf("RatLanNr error %d\n", ierr);
  }

  free(vinit);
  if (res) {
    free(res);
  }
  /* save pointers to the global variables */
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}


/** @brief Get the number of last computed eigenvalues 
 */
void EVSLFORT(evsl_get_nev)(int *nev) {
  *nev = evsl_nev_computed;
}

/** @brief copy the computed eigenvalues and vectors
 * @warning: after this call the internal saved results will be freed
 */
void EVSLFORT(evsl_copy_result)(double *val, double *vec) {
  memcpy(val, evsl_eigval_computed, evsl_nev_computed*sizeof(double));
  memcpy(vec, evsl_eigvec_computed, evsl_nev_computed*evsl_n*sizeof(double));
  /* reset global variables */
  evsl_nev_computed = 0;
  evsl_n = 0;
  free(evsl_eigval_computed);
  free(evsl_eigvec_computed);
  evsl_eigval_computed = NULL;
  evsl_eigvec_computed = NULL;
}




