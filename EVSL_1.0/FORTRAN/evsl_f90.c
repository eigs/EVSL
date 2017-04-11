#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/* global variable to hold results from EVSL */
int evsl_nev_computed=0, evsl_n=0;
double *evsl_eigval_computed=NULL, *evsl_eigvec_computed=NULL;

void evsl_start_f90_() {
  EVSLStart();
}

void evsl_finish_f90_() {
  EVSLFinish();
}

void evsl_seta_coo_f90_(int *n, int *nnz, int *ir, int *jc, double *vv) {
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
  evsldata.ifOwnA = 1;
}

void evsl_setb_coo_f90_(int *n, int *nnz, int *ir, int *jc, double *vv) {
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
  evsldata.ifOwnB = 1;
  //savemat(csr, "B.mtx");
  //exit(0);
}

void evsl_set_geneig_f90_() {
  SetGenEig();
}

void evsl_lanbounds_f90_(int *nsteps, double *lmin, double *lmax) {
  int n;
  double *vinit;
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  Malloc(vinit, n, double);
  rand_double(n, vinit);
  LanBounds(*nsteps, vinit, lmin, lmax);
  free(vinit);
}

void evsl_kpm_spslicer_f90_(int *Mdeg, int *nvec, double *xintv,
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

void evsl_find_pol_f90_(double *xintv, double *thresh_int, double *thresh_ext, 
                        size_t *polf90) {
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

  *polf90 = (size_t) pol;
}

void evsl_cheblantr_f90_(int *mlan, int *nev, double *xintv, int *max_its,
                         double *tol, size_t *polf90) {
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
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

void evsl_cheblannr_f90_(double *xintv, int *max_its, double *tol, size_t *polf90) {
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
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

void evsl_free_pol_f90_(size_t *polf90) {
  polparams *pol = (polparams *) (*polf90);
  free_pol(pol);
  free(pol);
}

void evsl_get_nev_f90_(int *nev) {
  *nev = evsl_nev_computed;
}

void evsl_copy_result_f90_(double *val, double *vec) {
  memcpy(val, evsl_eigval_computed, evsl_nev_computed*sizeof(double));
  memcpy(vec, evsl_eigvec_computed, evsl_nev_computed*evsl_n*sizeof(double));

  evsl_nev_computed = 0;
  evsl_n = 0;
  free(evsl_eigval_computed);
  free(evsl_eigvec_computed);
  evsl_eigval_computed = NULL;
  evsl_eigvec_computed = NULL;
}

void evsl_setamv_f90_(int *n, void *func, void *data) {
  SetAMatvec(*n, (MVFunc) func, data);
}

void evsl_find_rat_f90_(double *intv, size_t *ratf90) {
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

  *ratf90 = (size_t) rat;
}

void evsl_ratlannr_f90_(double *xintv, int *max_its, double *tol, size_t *ratf90) {
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
  evsl_nev_computed = nev2;
  evsl_n = n;
  evsl_eigval_computed = lam;
  evsl_eigvec_computed = Y;
}

void evsl_free_rat_f90_(size_t *ratf90) {
  ratparams *rat = (ratparams *) (*ratf90);
  free_rat(rat);
  free(rat);
}

