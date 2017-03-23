/*
  This file contains function prototypes and constant definitions internally used in EVSL
*/
#ifndef INTERNAL_PROTO_H
#define INTERNAL_PROTO_H

#include <stdio.h>
#include <complex.h>
#include "evsl.h"

/*- - - - - - - - - chebpoly.c */
//
int chebxCoefd(int m, double gam, int damping, double *mu);
//
int dampcf(int m, int damping, double *jac);
//
int chebxPltd(int m, double *mu, int n, double *xi, double *yi);
//
int ChebAv(polparams *pol, double *v, double *y, double *w);

int ChebAv0(polparams *pol, double *v, double *y, double *w);
//
void chext(polparams *pol, double aIn, double bIn);


/*- - - - - - - - - dump.c */
//
void savemat(csrMat *A, const char *fn);
//
void savedensemat(double *A, int lda, int m, int n, const char *fn);
//
void save_vec(int n, double *x, const char fn[]);

/*- - - - - - - - - evsl.c */
//

/*- - - - - - - - - misc_la.c */
//
int SymmTridEig(double *eigVal, double *eigVec, int n, const double *diag, const double *sdiag);
//
int SymmTridEigS(double *eigVal, double *eigVec, int n, double vl, double vu, int *nevO, const double *diag, const double *sdiag);
//
void SymEigenSolver(int n, double *A, int lda, double *Q, int ldq, double *lam);
//
void CGS_DGKS(int n, int k, int i_max, double *Q, double *v, double *nrmv, double *w);
//
void CGS_DGKS2(int n, int k, int i_max, double *Z, double *Q, double *v, double *w);
//
void orth(double *V, int n, int k, double *Vo, double *work);


/*- - - - - - - - - ratfilter.c */
//
void contQuad(int method, int n, complex double* zk);
//
void ratf2p2(int n, int *mulp, complex double *zk, complex double* alp, int m, double *z, double *x);
//
void pfe2(complex double s1, complex double s2, int k1, int k2, complex double* alp, complex double* bet);
//
complex double integg2(complex double s1, complex double s2, complex double* alp, int k1, complex double* bet, int k2, double a, double b);
//
void weights(int n, complex double* zk, int* pow, double lambda, complex double* omega);
//
int scaleweigthts(int n, double a, double b, complex double *zk, int* pow, complex double* omegaM);

/*- - - - - - - - - ratlanNr.c */
void RatFiltApply(int n, ratparams *rat, double *b, double *x, double *w3);

/*- - - - - - - - - spmat.c */
// matvec: y = A * x
int matvec_A(double *x, double *y);
// matvec: y = B * x
int matvec_B(double *x, double *y);

// memory allocation/reallocation for a CSR matrix
void csr_resize(int nrow, int ncol, int nnz, csrMat *csr);
//
void sortrow(csrMat *A);
//
int check_tri_full_diag(char type, csrMat *A);
//
int tri_sol_upper(char trans, csrMat *R, double *b, double *x);
//
int matadd(double alp, double bet, csrMat *A, csrMat *B, csrMat *C,
           int *mapA, int *mapB);
/* sparse identity */
int speye(int n, csrMat *A);

/*- - - - - - - - - suitesparse.c */
#ifdef EVSL_WITH_SUITESPARSE
int set_ratf_solfunc_default(csrMat *A, csrMat *B, ratparams *rat);
void free_rat_default_sol(ratparams *rat);
int set_default_LBdata(csrMat *B);
void free_default_LBdata();
void default_Bsol(double *b, double *x, void *data);
int set_Bsol_default(csrMat *B);
#endif

/*- - - - - - - - - timing.c */
int time_seeder();


/*- - - - - - - - - vect.c */
void vecset(int n, double t, double *v); 
void vec_perm(int n, int *p, double *x, double *y);
void vec_iperm(int n, int *p, double *x, double *y);

/*- - - - - - - - - - check if an interval is valid */
static inline int check_intv(double *intv, FILE *fstats) {
  /* intv[4]: ( intv[0], intv[1] ) is the inteval of interest
   *          ( intv[2], intv[3] ) is the spectrum bounds
   * return   0: ok
   *        < 0: interval is invalid
   */
  double a=intv[0], b=intv[1], lmin=intv[2], lmax=intv[3];
  if (a >= b) {
    fprintf(fstats, " error: invalid interval (%e, %e)\n", a, b);
    return -1;
  }
  
  if (a >= lmax || b <= lmin) {
    fprintf(fstats, " error: interval (%e, %e) is outside (%e %e) \n", 
            a, b, lmin, lmax);
    return -2;
  } 
  
  return 0;
}

#endif
