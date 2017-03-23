//============================================================================
// Routines for computing eigenvalues of a symmetric tridiagonal matrix.
// They are wrappers of the LAPACK routine DSTEV() or sstev_()
//============================================================================

#include <stdio.h>
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

/**----------------------------------------------------------------------- 
 *  @brief compute all eigenvalues and eigenvectors of a symmetric tridiagonal
 *  matrix  
 *  @param n                The  dimension of the symmetric tridiagonal  matrix
 *  @param diag[],sdiag[]   Define the symmetric tridiagonal  matrix:  the
 *          diagonal elements are diag[0,...,n-1]  in order and the subdiagonal
 *          elements are sdiag[0,...,n-2] in order  
 *  @param[out] eigVal The output vector of length n containing all eigenvalues
 *          in ascending order 
 *  @param[out] eigVec The output n-by-n matrix with columns as eigenvectors,
 *          in the order as elements in eigVal 
 *  @return The flag returned by the
 *  LAPACK routine DSTEV() (if double  is double) or stev_() (if double
 *  is float)
 * --------------------------------------------------------------------- */

int SymmTridEig(double *eigVal, double *eigVec, int n, 
                const double *diag, const double *sdiag) {
  char jobz = 'V';  // compute eigenvalues and eigenvectors
  int nn = n;
  int ldz = n;
  int info;  // output flag
  // copy diagonal and subdiagonal elements to alp and bet
  double *alp = eigVal;
  double *bet;
  Malloc(bet, n-1, double);
  memcpy(alp, diag, n*sizeof(double));
  memcpy(bet, sdiag, (n-1)*sizeof(double));
  // allocate storage for computation
  double *sv = eigVec;
  double *work;
  Malloc(work, 2*n-2, double);
  DSTEV(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
  // free memory
  free(bet);
  free(work);
  // return info
  return info;
}

/**----------------------------------------------------------------------- 
 *  @brief compute  eigenvalues and  eigenvectors of  a symmetric  tridiagonal
 *  matrix in a slice
 *  @param n The  dimension of  the  symmetric tridiagonal  matrix
 *  @param diag[],sdiag[]  define  the   symmetric  tridiagonal  matrix.  
 *  @param[out] eigVal Total number of eigenvalues found.
 *  @param[out] eigVec The first M elements contain teh selected eigenvalues in
 *  ascending oredr
 *  @param[in] vl If range='V', The lower bound of the interval to be searched
 *  for eigen values.
 *  @param[in] vu If  range='V', the upper bound of the interval to be searched
 *  for eigenvalues.
 *  @param[in] nevO If range='I', the index of the smallest eigen value to be
 *  returned.
 *
 *  This
 *  routine  computes selected  eigenvalues/vectors as  specified by  a
 *  range of values. This is a wrapper to the LAPACK routine DSTEMR().
 * ----------------------------------------------------------------------- */
int SymmTridEigS(double *eigVal, double *eigVec, int n, double vl, double vu,
                 int *nevO, const double *diag, const double *sdiag) {
  char jobz = 'V';  // compute eigenvalues and eigenvectors
  char range = 'V'; // compute eigenvalues in an interval

  // this does not use mwlapack for mex files
  int info;  
  //int idum = 0;
  //-------------------- isuppz not needed  
  int *isuppz;
  Malloc(isuppz, 2*n, int);
  //-------------------- real work array
  double *work;
  int lwork = 18*n;
  Malloc(work, lwork, double);
  //-------------------- int work array
  int *iwork;
  int liwork = 10*n;
  Calloc(iwork, liwork, int);
  //-------------------- copy diagonal + subdiagonal elements
  //                     to alp and bet
  double *alp; 
  double *bet;
  Malloc(bet, n, double);
  Malloc(alp, n, double);
  //
  memcpy(alp, diag, n*sizeof(double));
  if (n > 1) {
    memcpy(bet, sdiag, (n-1)*sizeof(double));
  }

  //-------------------- allocate storage for computation
  logical tryrac = 1;
  double t0 = vl;
  double t1 = vu;

  DSTEMR(&jobz, &range, &n, alp, bet, &t0, &t1, NULL, NULL, nevO, 
         eigVal, eigVec, &n, &n, isuppz, &tryrac, work, &lwork, 
         iwork, &liwork, &info);

  if (info) {
    printf("dstemr_ error %d\n", info);
  }

  //-------------------- free memory
  free(bet);
  free(alp);
  free(work);
  free(iwork);
  free(isuppz);
  //
  return info;
}


/**- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     @brief interface to   LAPACK SYMMETRIC EIGEN-SOLVER 
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
void SymEigenSolver(int n, double *A, int lda, double *Q, int ldq, double *lam) {
  /* compute eigenvalues/vectors of A that n x n, symmetric
   * eigenvalues saved in lam
   * eigenvectors saved in Q */
  char jobz='V';/* want eigenvectors */
  char uplo='U';/* use upper triangular part of the matrix */
  /*   copy A to Q */
  int i,j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Q[j+i*ldq] = A[j+i*lda];
    }
  }
  /*   first call: workspace query */
  double work_size;
  int lwork = -1, info;
  DSYEV(&jobz, &uplo, &n, Q, &ldq, lam, &work_size, &lwork, &info);
  if (info != 0) {
    fprintf(stdout, "DSYEV error [query call]: %d\n", info);
    exit(0);
  }
  /*   second call: do the computation */
  lwork = (int) work_size;
  double *work;
  Malloc(work, lwork, double);
  DSYEV(&jobz, &uplo, &n, Q, &ldq, lam, work, &lwork, &info);
  if (info != 0) {
    fprintf(stdout, "DSYEV error [comp call]: %d\n", info);
    exit(0);
  }
  free(work);
}

/**
 * @brief Classical GS reortho with Daniel, Gragg, Kaufman, Stewart test
 **/
void CGS_DGKS(int n, int k, int i_max, double *Q, double *v, double *nrmv, double *w) {
  double eta = 1.0 / sqrt(2.0);
  int i, one=1;
  char cT = 'T', cN = 'N';
  double done=1.0, dmone=-1.0, dzero=0.0;
  double old_nrm = DNRM2(&n, v, &one);
  double new_nrm = 0.0;
  for (i=0; i<i_max; i++) {
    DGEMV(&cT, &n, &k, &done,  Q, &n, v, &one, &dzero, w, &one);
    DGEMV(&cN, &n, &k, &dmone, Q, &n, w, &one, &done,  v, &one);
    new_nrm = DNRM2(&n, v, &one);
    if (new_nrm > eta * old_nrm) {
      break;
    }
    old_nrm = new_nrm;
  }

  if (nrmv) {
    *nrmv = new_nrm;
  }
}

/**
 * @brief Classical GS reortho. No test. just do i_max times 
 **/
void CGS_DGKS2(int n, int k, int i_max, double *Z, double *Q, 
               double *v, double *w) {
  int i, one=1;
  char cT = 'T', cN = 'N';
  double done=1.0, dmone=-1.0, dzero=0.0;
  for (i=0; i<i_max; i++) {
    DGEMV(&cT, &n, &k, &done,  Q, &n, v, &one, &dzero, w, &one);
    DGEMV(&cN, &n, &k, &dmone, Z, &n, w, &one, &done,  v, &one);
  }
}

//  max number of reorthogonalizations 
#define NGS_MAX 2
/**
 * @brief Orthogonalize columns of n-by-k matrix V 
 * @param n number of rows in V
 * @param V Matrix which columns are to be orthogonalized
 * @param k number of columns in V
 * @param[out] Vo Output matrix
 * @param work work 
 */
void orth(double *V, int n, int k, double *Vo, double *work) {
  int i;
  int one=1;
  int nk = n*k;
  DCOPY(&nk, V, &one, Vo, &one);
  double tt = DDOT(&n, Vo, &one, Vo, &one);
  double nrmv = sqrt(tt);
  double t = 1.0 / nrmv;
  DSCAL(&n, &t, Vo, &one); 
  for (i = 1; i < k; i++) {
    int istart = i*n;   
    CGS_DGKS(n, i, NGS_MAX, Vo, Vo+istart, &nrmv, work);
    t = 1.0 / nrmv;
    DSCAL(&n, &t, Vo+istart, &one); 
  }
}

