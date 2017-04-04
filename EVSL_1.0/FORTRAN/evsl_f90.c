#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

/* global variable to hold results from EVSL */
int evsl_nev_computed, evsl_n;
double *evsl_eigval_computed, *evsl_eigvec_computed;

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
  npts = 10 *ecount;
  spslicer(sli, mu, *Mdeg, xintv, *nslices, npts);
  *evint = (int) (1 + ecount / (*nslices));
  /*
  printf("%.15e, %.15e, %.15e, %.15e\n", xintv[0], xintv[1], xintv[2], xintv[3]);
  int j;
  for (j=0; j<*nslices;j++) {
    printf(" %2d: [% .15e , % .15e]\n", j+1, sli[j],sli[j+1]);
  }
  */
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
}

