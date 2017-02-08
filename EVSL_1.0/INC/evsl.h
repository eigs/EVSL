/*
  This file contains function prototypes and constant definitions for EVSL
*/

#ifndef EVSL_H
#define EVSL_H

#include <complex.h>
#include "struct.h"

/*- - - - - - - - - cheblanNr.c */
int ChebLanNr(csrMat *A, double *intv, int maxit, double tol, double *vinit, polparams *pol, int *nevOut,  double **lamo, double **Wo, double **reso, FILE *fstats);


/*- - - - - - - - - cheblanTr.c */
int ChebLanTr(csrMat *A, int lanm, int nev, double *intv, int maxit, double tol, double *vinit, polparams *pol, int *nev2, double **vals, double **W, double **resW, FILE *fstats);


/*- - - - - - - - - chebpoly.c */
//
void set_pol_def(polparams *pol);
//
int find_pol(double *intv, polparams *pol);
//
void free_pol(polparams *pol);

/*- - - - - - - - - chebsi.c */
int ChebSI(csrMat *A, int nev, double *intv, int maxit, double tol, double *vinit, polparams *pol, int *nevo, double **lamo, double **Yo, double **reso, FILE *fstats);


/*- - - - - - - - - exeiglap3.c */
int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo);


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
int RatLanNr(csrMat *A, double *intv, ratparams *, int maxit, double tol, double *vinit, int *nevOut, double **lamo, double **Wo, double **reso, FILE *fstats);


/*- - - - - - - - - ratlanTr.c */
//
int RatLanTr(csrMat *A, int lanm, int nev, double *intv, ratparams*, int maxit, double tol, double *vinit, int *nev2, double **lamo, double **Yo, double **reso, FILE *fstats);

/*- - - - - - - - - spmat.c */
// convert a COO matrix to a CSR matrix
int cooMat_to_csrMat(int cooidx, cooMat *coo, csrMat *csr);
// free a CSR
void free_csr(csrMat *csr);
// free a COO
void free_coo(cooMat *coo);
// generate a 2D/3D Laplacian matrix in COO format
int lapgen(int nx, int ny, int nz, cooMat *Acoo);
// set an external matvec function
void SetMatvecFunc(int n, matvecFunc func, void *data);
// unset an external matvec function
void UnsetMatvecFunc();

/*- - - - - - - - - spslicer.c */
//
int spslicer(double *sli, double *mu, int Mdeg, double *intv, int n_int,  int npts);
//
int kpmdos(csrMat *A, int Mdeg, int damping, int nvec, double *ab, double *mu, double *ecnt);


/*- - - - - - - - - timing.c */
//
double cheblan_timer();


/*- - - - - - - - - vect.c */
// generate a random vector
void rand_double(int n, double *v);
// sort a vector
void sort_double(int n, double *v, int *ind);
//
void linspace(double a, double b, int num, double *arr);

#endif
