/*
  This file contains function prototypes and constant definitions for EVSL
*/

#ifndef EVSL_H
#define EVSL_H

#include <complex.h>
#include "struct.h"

/*- - - - - - - - - cheblanNr.c */
int ChebLanNr(csrMat *A, double *intv, int maxit, double tol, double *vinit,
              polparams *pol, int *nevOut, double **lamo, double **Wo,
              double **reso, FILE *fstats);

/*- - - - - - - - - cheblanTr.c */
int ChebLanTr(csrMat *A, int lanm, int nev, double *intv, int maxit, double tol,
              double *vinit, polparams *pol, int *nev2, double **vals,
              double **W, double **resW, FILE *fstats);

/*- - - - - - - - - chebpoly.c */
//
void set_pol_def(polparams *pol);
//
int find_pol(double *intv, polparams *pol);
//
void free_pol(polparams *pol);

/*- - - - - - - - - chebsi.c */
int ChebSI(csrMat *A, int nev, double *intv, int maxit, double tol,
           double *vinit, polparams *pol, int *nevo, double **lamo, double **Yo,
           double **reso, FILE *fstats);

/*- - - - - - - - - lanbounds.c */
int LanBounds(csrMat *A, int msteps, double *v, double *lmin, double *lmax);

/*- - - - - - - - - ratfilter.c */
//
void set_ratf_def(ratparams *rat);
//
int find_ratf(double *intv, ratparams *rat);
//
int set_ratf_solfunc(ratparams *rat, csrMat *A, linSolFunc *funcs, void **data);
//
void free_rat(ratparams *rat);

/*- - - - - - - - - ratlanNr.c */
//
int RatLanNr(csrMat *A, double *intv, ratparams *, int maxit, double tol,
             double *vinit, int *nevOut, double **lamo, double **Wo,
             double **reso, FILE *fstats);

/*- - - - - - - - - ratlanTr.c */
//
int RatLanTr(csrMat *A, int lanm, int nev, double *intv, ratparams *, int maxit,
             double tol, double *vinit, int *nev2, double **lamo, double **Yo,
             double **reso, FILE *fstats);

/*- - - - - - - - - spmat.c */
// convert a COO matrix to a CSR matrix
int cooMat_to_csrMat(int cooidx, cooMat *coo, csrMat *csr);
// free a CSR
void free_csr(csrMat *csr);
// free a COO
void free_coo(cooMat *coo);
// matvec y = A*x
int matvec_csr(csrMat *A, double *x, double *y);

/*- - - - - - - - - evsl.c */
/* set an external matvec function */
void SetMatvecFunc(int n, MVFunc func, void *data);
/* unset an external matvec function */
void UnsetMatvecFunc();
/* set matrix B */
int SetRhsMatrix(csrMat *B);
/* unset matrix B */
void UnsetRhsMatrix();
/* start EVSL */
void EVSLStart();
/* finalize EVSL */
void EVSLFinish();

/*- - - - - - - - - spslicer.c */
//
int spslicer(double *sli, double *mu, int Mdeg, double *intv, int n_int,
             int npts);
//
int kpmdos(csrMat *A, int Mdeg, int damping, int nvec, double *ab, double *mu,
           double *ecnt);

/*- - - - - - - - - timing.c */
//
double cheblan_timer();

/*- - - - - - - - - vect.c */
// generate a random vector
void rand_double(int n, double *v);
// generate a normally distributed random vector
void randn_double(int n, double *v);
// sort a vector
void sort_double(int n, double *v, int *ind);
//
void linspace(double a, double b, int num, double *arr);

/*- - - - - - - - - - landos.c */
// Computes the density of states (DOS, or spectral density)
int LanDos(csrMat *A, const int nvec, int msteps, const int npts, double *xdos,
           double *ydos, double *neig, const double *const intv);

/*- - - - - - - - - - simpson.c */
void simpson(double *xi, double *yi, int npts, double *si);

/*- - - - - - - - - - spslicer2.c */
void spslicer2(double *xi, double *yi, int n_int, int npts, double *sli);

#endif
