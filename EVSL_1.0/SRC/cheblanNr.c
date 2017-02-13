#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"
/**-----------------------------------------------------------------------
 *  @brief Chebyshev polynomial filtering Lanczos process [NON-restarted version]
 *
 *  @param A        Matrix of size n x n
 * 
 *  @param intv     An array of length 4 \n
 *          [intv[0], intv[1]] is the interval of desired eigenvalues \n
 *          [intv[2], intv[3]] is the global interval of all eigenvalues \n
 *          Must contain all eigenvalues of A
 *  
 *  @param maxit    Max number of outer Lanczos steps  allowed --[max dim of Krylov 
 *          subspace]
 *  
 *  @param tol       
 *          Tolerance for convergence. The code uses a stopping criterion based
 *          on the convergence of the restricted trace. i.e., the sum of the
 *          eigenvalues of T_k that  are in the desired interval. This test  is
 *          rather simple since these eigenvalues are between `bar' and  1.0.
 *          We want the relative error on  this restricted  trace to be  less
 *          than  tol.  Note that the test  performed on filtered matrix only
 *          - *but* the actual residual norm associated with the original
 *          matrix A is returned
 *  
 *  @param vinit     Initial  vector for Lanczos -- [optional]
 * 
 *  @param pol       A struct containing the parameters of the polynomial..
 *  This is set up by a call to find_deg prior to calling chenlanNr 
 * 
 *  @b Modifies:
 *  @param[out] nevOut    Number of eigenvalues/vectors computed
 *  @param[out] Wo        A set of eigenvectors  [n x nevOut matrix]
 *  @param[out] reso      Associated residual norms [nev x 1 vector]
 *  @param[out] lamo      Lambda computed
 *  @param[out] fstats    File stream which stats are printed to
 *
 *  @return Returns 0 on success (or if check_intv() is non-positive),  -1
 *  if |gamB| < 1
 *
 *
 * @warning memory allocation for Wo/lamo/reso within this function 
 **/
int ChebLanNr(csrMat *A, double *intv, int maxit, double tol, double *vinit, 
    polparams *pol, int *nevOut,  double **lamo, double **Wo, 
    double **reso, FILE *fstats) {
  /*-------------------- for stats */
  double tm,  tmv=0.0, tr0, tr1, tall;
  double *y, flami; 
  //-------------------- to report timings/
  tall = cheblan_timer();
  int i, k, kdim;
  // handle case where fstats is NULL. Then no output. Needed for openMP.
  int do_print = 1;   
  if (fstats == NULL){
    do_print = 0;
  }  
  /*--------------------   frequently used constants  */
  char cN = 'N';   
  int one = 1;
  double done=1.0,dzero=0.0;
  /*--------------------   Ntest = when to start testing convergence */
  int Ntest = 30; 
  /*--------------------   how often to test */
  int cycle = 20; 
  /* size of the matrix */
  int n;
  /* if users provided their own matvec function, input matrix A will be ignored */
  if (evsldata.Amatvec.func) {
    n = evsldata.Amatvec.n;
  } else {
    n = A->nrows;
  }
  maxit = min(n, maxit);
  /*--------------------polynomial filter  approximates the delta
    function centered at 'gamB'. 
    bar: a bar value to threshold Ritz values of p(A)         */
  double bar = pol->bar;
  double gamB = pol->gam;
  /*-------------------- interval [aa, bb], used for testing only at end */
  if (check_intv(intv, fstats) < 0) {
    *nevOut = 0;
    *lamo = NULL; *Wo = NULL; *reso = NULL;
    return 0;
  }
  double aa = intv[0];
  double bb = intv[1];
  int deg = pol->deg;
  if (do_print) 
    fprintf(fstats, " ** Cheb Poly of deg = %d, gam = %.15e, bar: %.15e\n", 
        deg, gamB, bar);
  /*-------------------- gamB must be within [-1, 1] */
  if (gamB > 1.0 || gamB < -1.0) {
    fprintf(stdout, "gamB error %.15e\n", gamB);
    return -1;
  }
  /*-----------------------------------------------------------------------* 
   * *Non-restarted* Lanczos iteration 
   *-----------------------------------------------------------------------*/
  if (do_print) 
    fprintf(fstats, " ** Cheb-LanNr \n");
  /*-------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *dT, *eT;
  Malloc(V, n*(maxit+1), double);
  /*-------------------- diag. subdiag of Tridiagional matrix */
  Malloc(dT, maxit, double);
  Malloc(eT, maxit, double);
  double *W, *Lam, *res, *EvalT, *EvecT/*, *FY*/;
  /*-------------------- Lam, W: the converged (locked) Ritz values*/
  Malloc(Lam, maxit, double);         // holds computed Ritz values
  Malloc(res, maxit, double);         // residual norms (w.r.t. ro(A))
  Malloc(EvalT, maxit, double);       // eigenvalues of tridia. matrix  T
  Malloc(EvecT, maxit*maxit, double); // Eigen vectors of T  
  //Malloc(FY, maxit*maxit, double);    // coeffs of Ritz vectors in Lan basis
  /*-------------------- nev = current number of converged e-pairs 
    nconv = converged eigenpairs from looking at Tk alone */
  int nev, nconv = 0;
  /*-------------------- nmv counts  matvecs */
  int nmv = 0;
  /*-------------------- copy initial vector to V(:,1)   */
  DCOPY(&n, vinit, &one, V, &one);
  /*--------------------  normalize it     */
  double t, t1, t2, nt, res0;
  t = DDOT(&n, V, &one, V, &one); // add a test here.
  t = 1.0 / sqrt(t);
  DSCAL(&n, &t, V, &one);
  /*-------------------- u  is just a pointer. wk == work space */
  double *u, *wk;
  Malloc(wk, 3*n, double);
  /*-------------------- for ortho test */
  double wn = 0.0;
  int nwn = 0;
  /*-------------------- for stopping test [restricted trace]*/
  tr0 = 0;
  /*-------------------- lanczos vectors updated by rotating pointer*/
  /*-------------------- pointers to Lanczos vectors */
  double *vold, *v, *w;
  /*--------------------  Lanczos recurrence coefficients */
  double alpha, nalpha, beta=0.0, nbeta, resi;
  int count = 0;
  // ---------------- main Lanczos loop 
  for (k=0; k<maxit; k++) {
    /*-------------------- quick reference to V(:,k-1) when k>0*/
    vold = k > 0 ? V+(k-1)*n : NULL; 
    /*-------------------- a quick reference to V(:,k) */
    v = &V[k*n];
    /*--------------------   next Lanczos vector V(:,k+1)*/
    w = v + n;
    /*-------------------- compute   w = p[(A-cc)/dd] * v */
    /*  orthgonlize against the locked ones first */
    tm = cheblan_timer();
    ChebAv(A, pol, v, w, wk);
    tmv += cheblan_timer() - tm;
    nmv += deg;
    /*   w = w - beta*vold */
    if (vold) {
      nbeta = -beta;
      DAXPY(&n, &nbeta, vold, &one, w, &one);
    }
    /*--------------------   alpha = w'*v */
    alpha = DDOT(&n, v, &one, w, &one);
    dT[k] = alpha;       //   T(k,k) = alpha 
    wn += fabs(alpha);
    /*   w = w - alpha*v */
    nalpha = -alpha;
    DAXPY(&n, &nalpha, v, &one, w, &one);
    /*--------------------   FULL reortho to all previous Lan vectors */
    /*   w = w - V(:,1:k)*V(:,1:k)'*w */
    /*   beta = norm(w) */
    CGS_DGKS(n, k+1, NGS_MAX, V, w, &beta, wk);
    wn += 2.0 * beta;
    nwn += 3;
    //vold = v;
    /*-------------------- lucky breakdown test */
    if (beta*nwn < orthTol*wn) {
      if (do_print) 
        fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", k, beta);
      rand_double(n, w);
      beta = DDOT(&n, w, &one, w, &one);
      beta = sqrt(beta);
    }
    eT[k] = beta;
    /*   w = w / beta --------------------*/
    t = 1.0 / beta;
    DSCAL(&n, &t, w, &one);
    /*--------------------  test for Ritz vectors */
    if ( (k < Ntest || (k-Ntest) % cycle != 0) && k != maxit-1 ) {
      continue;
    }
    /*--------------------   diagonalize  T(1:k,1:k)       */
    /*                       vals in EvalT, vecs in EvecT  */
    kdim = k+1;
#if 1
    //-------------------- THIS uses dsetv
    SymmTridEig(EvalT, EvecT, kdim, dT, eT);
    count = kdim;
#else
    //-------------------- THIS uses dstemr:
    double vl = bar - DBL_EPSILON, vu = 2.0;  // needed by SymmTridEigS
    SymmTridEigS(EvalT, EvecT, kdim, vl, vu, &count, dT, eT);
#endif
    tr1 = 0;     // restricted trace -- used for convergence test
    /*-------------------- get residual norms and check acceptance of
      Ritz values for p(A). nconv records number of eigenvalues whose
      residual for p(A) is smaller than tol. */
    nconv = 0;
    for (i=0; i<count; i++) {
      flami = EvalT[i];
      if (flami >= bar) tr1+= flami;
      // the last row of EvecT: EvecT[i*kdim+kdim-1]
      if (beta*fabs(EvecT[(i+1)*kdim-1]) < tol) nconv++;
    }

    if (do_print) {
      fprintf(fstats, "k %4d:   nMV %8d, nconv %4d  tr1 %21.15e\n",
              k, nmv, nconv,tr1);
    }
    //-------------------- simple test because all eigenvalues
    // are between gamB and ~1.
    if (fabs(tr1-tr0)<tol*fabs(tr1)) {
      break;
    }
    tr0 = tr1;
  }

  //-------------------- done == compute Ritz vectors --
  Malloc(W,nconv*n, double);       // holds computed Ritz vectors
  //
  nev = 0;
  for (i=0; i<count;i++) {
    flami = EvalT[i];
    //-------------------- reject eigenvalue if rho(lam)<bar
    if (flami < bar) 
      continue;
    y = &EvecT[i*kdim];
    //-------------------- make sure to normalize
    t = DNRM2(&kdim, y, &one);  t = 1.0 / t; 
    DSCAL(&kdim, &t, y, &one);
    //-------------------- residual norm 
    resi = beta*fabs(y[kdim-1]);
    if (resi > tol)
      continue;
    //-------------------- compute Ritz vectors.
    u = &W[nev*n];  
    DGEMV(&cN, &n, &kdim, &done, V, &n, y, &one, &dzero, u, &one);
    /*--------------------   w = A*u        */
    matvec_genev(A, u, wk);
    nmv ++;
    /*--------------------   Ritzval: t = (y'*w)/(y'*y) */
    t1 = DDOT(&n, u, &one, u, &one);  // should be one
    t2 = DDOT(&n, wk, &one, u, &one);
    t  = t2 / t1;
    /*--------------------  if lambda (==t) is in [a,b] */
    if (t < aa - DBL_EPSILON || t > bb + DBL_EPSILON) 
      continue;
    /*-------------------- compute residual wrt A for this pair */
    nt = -t;
    /*-------------------- w = w - t*y */
    DAXPY(&n, &nt, u, &one, wk, &one);
    /*--------------------   res0 = norm(w) */
    res0 = DNRM2(&n, wk, &one); 
    /*--------------------   accept (t, y) */
    Lam[nev] = t;
    res[nev] = res0;
    nev++;
  }

  /* for generalized eigenvalue problem: L' \ W */
  if (evsldata.hasB) {
    for (i=0; i<nev; i++) {
      evsldata.LBT_solv(W+i*n, wk, evsldata.LB_func_data);
      DCOPY(&n, wk, &one, W+i*n, &one);
    }
  }

  /*-------------------- Done.  output : */
  *nevOut = nev;
  *lamo = Lam;
  *Wo = W;
  *reso = res;
  /*-------------------- free arrays */
  free(V);
  free(dT);
  free(eT);
  free(EvalT);
  free(EvecT);
  /*free(FY);*/
  free(wk);
  /*-------------------- record stats */
  tall = cheblan_timer() - tall;
  /*-------------------- print stat */
  if (do_print){
    fprintf(fstats, "------This slice consumed: \n");
    fprintf(fstats, "Matvecs :        %d\n", nmv);
    fprintf(fstats, "total  time :    %.2f\n", tall);
    fprintf(fstats, "matvec time :    %.2f\n", tmv);
    fprintf(fstats,"======================================================\n");
  }
  return 0;
}


