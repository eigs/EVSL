//============================================================================
// Routines for computing eigenvalues of a symmetric tridiagonal matrix.
// They are wrappers of the LAPACK routine dstev() or sstev_()
//============================================================================

#include <stdio.h>
#include <string.h>  // for memcpy, strcmp, strncmp, etc.
#include "internal_header.h"

/* in Gram-Schmidt, use BLAS-2 GEMV or BLAS-1 DOT/AXPY */
#define USE_DGEMV 1

/**
 * @file misc_la.c
 * @brief Miscellaneous linear algebra functions
 */
/**-----------------------------------------------------------------------
 *  @brief compute all eigenvalues and eigenvectors of a symmetric tridiagonal
 *  matrix
 *  @param[in] n                The  dimension of the symmetric tridiagonal  matrix
 *  @param[in] diag   Define the symmetric tridiagonal  matrix:  the
 *          diagonal elements are diag[0,...,n-1]
 *  @param[in] sdiag Subdiagonal elements
 *  @param[out] eigVal The output vector of length n containing all eigenvalues
 *          in ascending order
 *  @param[out] eigVec The output n-by-n matrix with columns as eigenvectors,
 *          in the order as elements in eigVal. If NULL, then no eigenvector
 *          will be computed
 *  @return The flag returned by the
 *  LAPACK routine dstev() (if double) or sstev_() (if float)
 * --------------------------------------------------------------------- */

int SymmTridEig(double *eigVal, double *eigVec, int n,
                const double *diag, const double *sdiag) {
  double tms = evsl_timer();
  // compute eigenvalues and eigenvectors or eigvalues only
  const char jobz = eigVec ? 'V' : 'N';
  const int nn = n;
  const int ldz = n;
  int info;  // output flag
  // copy diagonal and subdiagonal elements to alp and bet
  double *alp = eigVal;
  double *bet;
  bet = evsl_Malloc(n-1, double);
  memcpy(alp, diag, n*sizeof(double));
  memcpy(bet, sdiag, (n-1)*sizeof(double));
  // allocate storage for computation
  double *sv = eigVec;
  double *work = NULL;
  if (jobz == 'V') {
    work = evsl_Malloc(2*n-2, double);
  }
  evsl_dstev(&jobz, &nn, alp, bet, sv, &ldz, work, &info);
  // free memory
  evsl_Free(bet);
  if (work) {
    evsl_Free(work);
  }
  if (info) {
    printf("evsl_dstev ERROR: INFO %d\n", info);
    save_vec(n, diag, "alp");
    save_vec(n-1, sdiag, "bet");
    exit(0);
  }
  double tme = evsl_timer();
  evslstat.t_eig += tme - tms;
  // return info
  return info;
}

/**-----------------------------------------------------------------------
 *  @brief compute  eigenvalues and  eigenvectors of  a symmetric  tridiagonal
 *  matrix in a slice
 *  @param[in] n The  dimension of  the  symmetric tridiagonal  matrix
 *  @param[in] diag Diagonal elements
 *  @param[in] sdiag Sub-diagonal elements
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
 *  range of values. This is a wrapper to the LAPACK routine dstemr().
 * ----------------------------------------------------------------------- */
int SymmTridEigS(double *eigVal, double *eigVec, int n, double vl, double vu,
                 int *nevO, const double *diag, const double *sdiag) {
  const double tms = evsl_timer();
  const char jobz = 'V';  // compute eigenvalues and eigenvectors
  const char range = 'V'; // compute eigenvalues in an interval

  // this does not use mwlapack for mex files
  int info;
  //int idum = 0;
  //-------------------- isuppz not needed
  int *isuppz;
  isuppz = evsl_Malloc(2*n, int);
  //-------------------- real work array
  double *work;
  const int lwork = 18*n;
  work = evsl_Malloc(lwork, double);
  //-------------------- int work array
  int *iwork;
  const int liwork = 10*n;
  iwork = evsl_Calloc(liwork, int);
  //-------------------- copy diagonal + subdiagonal elements
  //                     to alp and bet
  double *alp;
  double *bet;
  bet = evsl_Malloc(n, double);
  alp = evsl_Malloc(n, double);
  //
  memcpy(alp, diag, n*sizeof(double));
  if (n > 1) {
    memcpy(bet, sdiag, (n-1)*sizeof(double));
  }

  //-------------------- allocate storage for computation
  //logical tryrac = 1;
  int tryrac = 1;
  double t0 = vl;
  double t1 = vu;

  evsl_dstemr(&jobz, &range, &n, alp, bet, &t0, &t1, NULL, NULL, nevO,
              eigVal, eigVec, &n, &n, isuppz, &tryrac, work, &lwork,
              iwork, &liwork, &info);

  if (info) {
    printf("dstemr_ error %d\n", info);
  }

  //-------------------- free memory
  evsl_Free(bet);
  evsl_Free(alp);
  evsl_Free(work);
  evsl_Free(iwork);
  evsl_Free(isuppz);

  double tme = evsl_timer();
  evslstat.t_eig += tme - tms;
  //
  return info;
}


/**- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *     @brief interface to   LAPACK SYMMETRIC EIGEN-SOLVER
 *     @param[in] n Size of problem
 *     @param[in] A Matrix
 *     @param[in] lda Leading dimension
 *     @param[out] Q Eigenvectors
 *     @param[in] ldq Leading dimension q
 *     @param[out] lam Eigenvalues
 *- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
void SymEigenSolver(int n, const double *A, int lda, double *Q, int ldq, double *lam) {
  double tms = evsl_timer();
  /* compute eigenvalues/vectors of A that n x n, symmetric
   * eigenvalues saved in lam: the eigenvalues in ascending order
   * eigenvectors saved in Q */
  const char jobz='V';/* want eigenvectors */
  const char uplo='U';/* use upper triangular part of the matrix */
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
  evsl_dsyev(&jobz, &uplo, &n, Q, &ldq, lam, &work_size, &lwork, &info);
  if (info != 0) {
    fprintf(stdout, "DSYEV error [query call]: %d\n", info);
    exit(0);
  }
  /*   second call: do the computation */
  lwork = (int) work_size;
  double *work;
  work = evsl_Malloc(lwork, double);
  evsl_dsyev(&jobz, &uplo, &n, Q, &ldq, lam, work, &lwork, &info);
  if (info != 0) {
    fprintf(stdout, "DSYEV error [comp call]: %d\n", info);
    exit(0);
  }
  evsl_Free(work);
  double tme = evsl_timer();
  evslstat.t_eig += tme - tms;
}

/**
 * @brief Classical GS reortho with Daniel, Gragg, Kaufman, Stewart test
 * @param[in] n Number of rows in Q
 * @param[in] k Number of cols in Q
 * @param[in] i_max Number iterations
 * @param[in] Q matrix
 * @param[out] nrmv
 * @param[out] v Output
 * @param[out] w Output
 **/
void CGS_DGKS(int n, int k, int i_max, const double *Q, double *v, double *nrmv, double *w) {
  const double tms = evsl_timer();
  const double eta = 1.0 / sqrt(2.0);
  int i, one=1;
#if USE_DGEMV
  char cT = 'T', cN = 'N';
  double done=1.0, dmone=-1.0, dzero=0.0;
#endif
  double old_nrm = evsl_dnrm2(&n, v, &one);
  double new_nrm = 0.0;

  for (i=0; i<i_max; i++) {
#if USE_DGEMV
    evsl_dgemv(&cT, &n, &k, &done,  Q, &n, v, &one, &dzero, w, &one);
    evsl_dgemv(&cN, &n, &k, &dmone, Q, &n, w, &one, &done,  v, &one);
#else
    int j;
    for (j=0; j<k; j++) {
       const double t = -evsl_ddot(&n, &Q[j*n], &one, v, &one);
       evsl_daxpy(&n, &t, &Q[j*n], &one, v, &one);
    }
#endif
    new_nrm = evsl_dnrm2(&n, v, &one);
    if (new_nrm > eta * old_nrm) {
      break;
    }
    old_nrm = new_nrm;
  }

  if (nrmv) {
    *nrmv = new_nrm;
  }
  const double tme = evsl_timer();
  evslstat.t_reorth += tme - tms;
}

/**
 * @brief Classical GS reortho. No test. just do i_max times
 * used in generalized ev problems
 *
 * @param[in] n Number of rows in Z,Q
 * @param[in] k Number of cols in Z,Q
 * @param[in] i_max Number iterations
 * @param[in] Z matrix
 * @param[in] Q matrix
 * @param[out] v Output
 * @param[out] w Output
 **/
void CGS_DGKS2(int n, int k, int i_max, const double *Z, const double *Q,
               double *v, double *w) {
  const double tms = evsl_timer();
  int i;
  const int one=1;
#if USE_DGEMV
  const char cT = 'T';
  const char cN = 'N';
  const double done=1.0;
  const double dmone=-1.0;
  const double dzero=0.0;
#endif
  for (i=0; i<i_max; i++) {
#if USE_DGEMV
    evsl_dgemv(&cT, &n, &k, &done,  Q, &n, v, &one, &dzero, w, &one);
    evsl_dgemv(&cN, &n, &k, &dmone, Z, &n, w, &one, &done,  v, &one);
#else
    int j;
    for (j=0; j<k; j++) {
       const double t = -evsl_ddot(&n, &Q[j*n], &one, v, &one);
       evsl_daxpy(&n, &t, &Z[j*n], &one, v, &one);
    }
#endif
  }
  const double tme = evsl_timer();
  evslstat.t_reorth += tme - tms;
}

//  max number of reorthogonalizations
#define NGS_MAX 2
/**
 * @brief Orthogonalize columns of n-by-k matrix V
 * @param[in] n number of rows in V
 * @param[in] V Matrix which columns are to be orthogonalized
 * @param[in] k number of columns in V
 * @param[out] Vo Output matrix
 * @param[in, out] work work
 *
 * @warning Aliasing happens in call to CGS_DGKS
 */
void orth(const double *V, int n, int k, double *Vo, double *work) {
  int i;
  const int one=1;
  const int nk = n*k;
  evsl_dcopy(&nk, V, &one, Vo, &one);
  const double tt = evsl_ddot(&n, Vo, &one, Vo, &one);
  double nrmv = sqrt(tt);
  double t = 1.0 / nrmv;
  evsl_dscal(&n, &t, Vo, &one);
  for (i = 1; i < k; i++) {
    const int istart = i*n;
    CGS_DGKS(n, i, NGS_MAX, Vo, Vo+istart, &nrmv, work);
    t = 1.0 / nrmv;
    evsl_dscal(&n, &t, Vo+istart, &one);
  }
}

