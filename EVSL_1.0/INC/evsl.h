/*
  This file contains function prototypes and constant definitions for EVSL
*/

#ifndef EVSL_H
#define EVSL_H

#include <complex.h>
#include "struct.h"

/*- - - - - - - - - cheblanNr.c */
int ChebLanNr(double *intv, int maxit, double tol, double *vinit, 
              polparams *pol, int *nevOut,  double **lamo, double **Wo, 
              double **reso, FILE *fstats);

/*- - - - - - - - - cheblanTr.c */
int ChebLanTr(int lanm, int nev, double *intv, int maxit, 
              double tol, double *vinit, polparams *pol, int *nev2, 
              double **vals, double **W, double **resW, FILE *fstats);

/*- - - - - - - - - chebpoly.c */
//
void set_pol_def(polparams *pol);
//
int find_pol(double *intv, polparams *pol);
//
void free_pol(polparams *pol);

/*- - - - - - - - - chebsi.c */
int ChebSI(int nev, double *intv, int maxit, double tol, 
           double *vinit, polparams *pol, int *nevo, double **lamo, 
           double **Yo, double **reso, FILE *fstats);


/*- - - - - - - - - lanbounds.c */
int LanBounds(int msteps, double *v, double *lmin, double *lmax);


/*- - - - - - - - - ratfilter.c */
//
void set_ratf_def(ratparams *rat);
//
int find_ratf(double *intv, ratparams *rat);
//
int set_ratf_solfunc(ratparams *rat, csrMat *A, csrMat *B, SolFuncC *funcs, 
                     void **data);
//
void free_rat(ratparams *rat);

/*- - - - - - - - - ratlanNr.c */
//
int RatLanNr(double *intv, int maxit, double tol, double *vinit, 
             ratparams *rat, int *nevOut, double **lamo, double **Wo, 
             double **reso, FILE *fstats);

/*- - - - - - - - - ratlanTr.c */
//
int RatLanTr(int lanm, int nev, double *intv, int maxit, 
             double tol, double *vinit, ratparams *rat, int *nev2, 
             double **vals, double **W, double **resW, FILE *fstats);

/*- - - - - - - - - spmat.c */
// convert a COO matrix to a CSR matrix
int cooMat_to_csrMat(int cooidx, cooMat *coo, csrMat *csr);
// free a CSR
void free_csr(csrMat *csr);
// free a COO
void free_coo(cooMat *coo);
// matvec y = A*x
int matvec(char trans, csrMat *A, double *x, double *y);

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
int SetLTSol(SolFuncR func);
int SetASigmaBSol(ratparams *rat, SolFuncC *func, SolFuncC allf, void **data);
int SetStdEig();
int SetGenEig();
/* start EVSL */
int EVSLStart();
/* finalize EVSL */
int EVSLFinish();

/*- - - - - - - - - spslicer.c */
//
int spslicer(double *sli, double *mu, int Mdeg, double *intv, int n_int,  int npts);
//
int kpmdos(int Mdeg, int damping, int nvec, double *ab, double *mu, 
           double *ecnt);

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
