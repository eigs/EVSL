#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"
/**
 * @brief Chebyshev polynomial filtering Lanczos process [Thick restart version]
 *
 * @param A         Matrix of size n x n
 * @param lanm      Dimension of Krylov subspace [restart dimension]
 * @param nev       Estimate of number of eigenvalues in the interval --
 *         ideally nev == exact number or a little larger.  This is not used
 *         for testing convergence but it helps make decisions as to when to
 *         test convergence ChebLanTr attempts to compute *all* eigenvalues in
 *         the interval and stops only when no more eigenvalyes are left. The
 *         convergenve test is a very simple one based on the residual norm for
 *         the filtered matrix 
 *
 * @param intv   an array of length 4  \n
 *         [intv[0], intv[1]] is the interval of desired eigenvalues \n
 *         [intv[2], intv[3]] is the global interval of all eigenvalues \n
 *         it must contain all eigenvalues of A
 * 
 * @param maxit     max Num of outer Lanczos iterations (restarts) allowed -- 
 *         Each restart may or use the full lanm lanczos steps or fewer.
 * 
 * @param tol       tolerance for convergence. stop when ||res||< tol
 * @param vinit     initial  vector for Lanczos -- [optional]
 * @param pol       a struct containing the parameters of the polynomial. This 
 *         is set up by a call to find_deg prior to calling chenlanNr 
 *
 * @b Modifies:
 * @param[out]  nev2     Number of eigenvalues/vectors computed
 * @param[out] W         A set of eigenvectors  [n x nev2 matrix]
 * @param[out] vals      Associated eigenvalues [nev2 x 1 vector]
 * @param[out] resW      Associated residual norms [nev x 1 vector]
 * @param[out] fstats    File stream which stats are printed to
 *
 * @return Returns 0 on success (or if check_intv() is non-positive), -1 if
 * gamB is outside [-1, 1], and 2 if there are no eigenvalues found.
 *
 *
 * @warning memory allocation for W/vals/resW within this function 
 *
 **/
int ChebLanTr(csrMat *A, int lanm, int nev, double *intv, int maxit, 
    double tol, double *vinit, polparams *pol, int *nev2, 
    double **vals, double **W, double **resW, FILE *fstats) {
  /*-------------------- for stats */
  double tm, tall=0.0, tmv=0.0;
  double tolP = tol;
  tall = cheblan_timer();
  int do_print = 1;
  // handle case where fstats is NULL. Then no output. Needed for openMP.
  if (fstats == NULL){
    do_print = 0;
  }  
  /*--------------------- size of A */
  int n = A->nrows;
  /*--------------------- adjust lanm and maxit */
  lanm = min(lanm, n);
  int lanm1=lanm+1;
  // if use full lanczos, should not do more than n iterations
  if (lanm == n) { maxit = min(maxit, n); }
  /*-------------------- this is needed to increment the space when we
    discover more than nev eigenvalues in interval */
  double nevInc = 0.2;   // add 1  + 20% each time it is needed
  /*-------------------- if we have at least nev/ev_frac good candidate 
    eigenvalues from p(A) == then we restart to lock them in */
  int evFrac = 2;
  /*--------------------   some constants frequently used */
  //char cT='T';
  char cN = 'N';
  int one = 1;
  double done=1.0,dmone=-1.0,dzero=0.0;
  /*--------------------   Ntest = when to start testing convergence */
  int Ntest = min(lanm, nev+50);
  /*--------------------   how often to test */
  int cycle = 30; 
  int i, ll, count, last_count, jl, last_jl;
  /*-----------------------------------------------------------------------
   -----------------------------------------------------------------------*/
  if (check_intv(intv, fstats) < 0) {
    *nev2 = 0;
    *vals = NULL; *W = NULL; *resW = NULL;
    return 0;
  }
  double aa = intv[0];
  double bb = intv[1];
  /*-------------------- Assumption: polynomial pol computed before calling cheblanTr
  pol.  approximates the delta function centered at 'gamB'
  bar: a bar value to threshold Ritz values of p(A)
  */
  int deg = pol->deg;
  double gamB=pol->gam, bar=pol->bar;
  /*-------------------- gamB must be within [-1, 1] */
  if (gamB > 1.0 || gamB < -1.0) {
    fprintf(fstats, "gamB error %.15e\n", gamB);
    return -1;
  }
  /*  save mu in file */ 
  //save_vec(deg+1, mu, "OUT/mu.mtx");
  /*-----------------------------------------------------------------------* 
   * *thick restarted* Lanczos step 
   *-----------------------------------------------------------------------*/
  if (do_print) 
    fprintf(fstats, " Cheb-LanTR, dim %d  cycle %d  \n",lanm, cycle);
  /*   the min number of steps to be performed for each innter loop, 
       must be >= 1 - if not, it means that the Krylov dim is too small */
  int min_inner_step = 5;
  /*-------------------- it = number of Lanczos steps */
  int it = 0;
  /*-------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *T;
  Malloc(V, n*lanm1, double);
  /*-------------------- T must be zeroed out initially*/
  Calloc(T, lanm1*lanm1, double);
  /*-------------------- Lam, Y: the converged (locked) Ritz values/vectors 
res: related residual norms */
  double *Y, *Lam, *res;
  Malloc(Y, n*nev, double);
  Malloc(Lam, nev, double);
  Malloc(res, nev, double);
  /*-------------------- lock =  number of locked vectors */
  int lock = 0;
  /*-------------------- trlen = dim. of thick restart set */
  int trlen = 0, prtrlen=-1;
  /*-------------------- nmv counts  matvecs */
  int nmv = 0;
  /*-------------------- Ritz values and vectors of p(A) */
  double *Rval, *Rvec, *resi;
  Malloc(Rval, lanm, double);
  Malloc(resi, lanm, double);
  Malloc(Rvec, n*lanm, double);
  /*-------------------- Eigen vectors of T */
  double *EvecT;
  Malloc(EvecT, lanm1*lanm1, double);
  /*-------------------- s used by TR (spike of 1st block in Tm)*/
  double *s;
  Malloc(s, lanm, double);
  /*-------------------- copy initial vector to V(:,1)   */
  DCOPY(&n, vinit, &one, V, &one);
  /*--------------------  normalize it     */
  double t = DDOT(&n, V, &one, V, &one);
  t = 1.0 / sqrt(t);
  DSCAL(&n, &t, V, &one);
  /*-------------------- alloc some work space */
  double *work;
  int work_size = 3*n;
  Malloc(work, work_size, double);
  /*-------------------- main (restarted Lan) outer loop */
  while (it < maxit) {
    /*-------------------- for ortho test */
    double wn = 0.0;
    int nwn = 0;
    /*  beta */
    double beta = 0.0, r, res0; 
    /*  start with V(:,k) */
    int k = trlen;
    int k1 = k+1;
    /* ! add a test if dimension exceeds (m+1) 
     * (trlen + 1) + min_inner_step <= lanm + 1 */
    if (k1+min_inner_step > lanm1) {
      fprintf(fstats, "Krylov dim too small for this problem. Try a larger dim\n");
      exit(1);
    }
    /*--------------------  thick restart special step */
    if (trlen > 0) {
      /*    a quick reference to V(:,k) */
      double *v = &V[k*n];
      /*    next Lanczos vector */
      double *w = v + n;
      /*    w = p[(A-cc)/dd] * v     */
      tm = cheblan_timer();
      ChebAv(A, pol, v, w, work);
      tmv += cheblan_timer() - tm;
      nmv += deg;
      if (lock > 0) {
        /* orthgonlize against locked vectors first, w = w - Y*Y'*w */
        CGS_DGKS(n, lock, NGS_MAX, Y, w, NULL, work);
      }
      /*-------------------- restart with 'trlen' Ritz values/vectors
        T = diag(Rval(1:trlen)) */
      for (i=0; i<trlen; i++) {
        T[i*lanm1+i] = Rval[i];
        wn += fabs(Rval[i]);
      }
      /*   s(k) = V(:,k)'*w */
      s[k] = DDOT(&n, v, &one, w, &one);
      /*   w = w - V(:,1:k)*s(1:k) */
      DGEMV(&cN, &n, &k1, &dmone, V, &n, s, &one, &done, w, &one);
      /*-------------------- expand T matrix to k-by-k, arrow-head shape
        T = [T, s(1:k-1)] then T = [T; s(1:k)'] */
      for (i=0; i<k; i++) {
        T[trlen*lanm1+i] = s[i];
        T[i*lanm1+trlen] = s[i];
        wn += 2.0 * fabs(s[i]);
      }
      T[trlen*lanm1+trlen] = s[k];
      wn += fabs(s[k]);
      /*   beta = norm(w) */
      beta = DDOT(&n, w, &one, w, &one);
      beta = sqrt(beta);
      /*   T(k+1,k) = beta;    then   T(k,k+1) = beta; */
      T[k*lanm1+k] = beta;
      T[k1*lanm1+k] = beta;
      wn += 2.0 * beta;
      nwn += 3*k1;
      /*   beta ~ 0 */
      if (beta*nwn < orthTol*wn) {
        rand_double(n, w);
        beta = DDOT(&n, w, &one, w, &one);
        beta = sqrt(beta);
      }
      /*   w = w / beta */
      double ibeta = 1.0 / beta;
      DSCAL(&n, &ibeta, w, &one);
    } 
    /*-------------------- Done with TR step. Rest of Lanczos step */
    //-------------------- reset Ntest at each restart. 
    Ntest = max(20,nev-lock+10);
    last_count = 0;  last_jl = 0;
    /*   regardless of trlen,  *(k+1)* is the current number of Lanczos 
         vectors in V */
    /*  pointer to the previous Lanczos vector */
    double *vold = k > 0 ? V+(k-1)*n : NULL;
    /*------------------------------------------------------*/
    /*------------------ Lanczos inner loop ----------------*/
    /*------------------------------------------------------*/
    while (k < lanm && it < maxit) {
      k++;
      /*   a quick reference to V(:,k) */
      double *v = &V[(k-1)*n];
      /*   next Lanczos vector */
      double *w = v + n;
      /*   w = p[(A-cc)/dd] * v */
      tm = cheblan_timer();
      ChebAv(A, pol, v, w, work);
      tmv += cheblan_timer() - tm;
      nmv += deg;
      it++;
      /*-------------------- orthgonalize vs locked ones first */
      if (lock > 0) {
        /*--------------------   w = w - Y*Y'*w */
        CGS_DGKS(n, lock, NGS_MAX, Y, w, NULL, work);
      }
      /*--------------------  w = w - beta*vold */
      if (vold) {
        double nbeta = -beta;
        DAXPY(&n, &nbeta, vold, &one, w, &one);
      }
      /*--------------------   alpha = w'*v */
      double alpha = DDOT(&n, v, &one, w, &one);
      /*--------------------   T(k,k) = alpha */
      T[(k-1)*lanm1+(k-1)] = alpha;
      wn += fabs(alpha);
      /*--------------------   w = w - alpha*v */
      double nalpha = -alpha;
      DAXPY(&n, &nalpha, v, &one, w, &one);
      /*   FULL reortho to all previous Lan vectors */
      /*   w = w - V(:,1:k)*V(:,1:k)'*w */
      /*   beta = norm(w) */
      CGS_DGKS(n, k, NGS_MAX, V, w, &beta, work);
      /*--------------------  T(k,k+1) = T(k+1,k) = beta */
      T[k*lanm1+(k-1)] = beta;
      T[(k-1)*lanm1+k] = beta;
      wn += 2.0 * beta;
      nwn += 3;
      /*   vold = v */
      vold = v;
      /*-------------------- lucky breakdown  test */
      if (beta*nwn < orthTol*wn) {
        if (do_print) {
          fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", it, beta);
        }
        rand_double(n, w);
        beta = DNRM2(&n, w, &one);
        //beta = DDOT(&n, w, &one, w, &one);  beta = sqrt(beta);
      }
      /*--------------------   w = w / beta */
      double ibeta = 1.0 / beta;
      DSCAL(&n, &ibeta, w, &one);
      /*-------------------- Restarting test */
      k1 = k-trlen-Ntest;
      if ( ((k1>=0) && (k1 % cycle == 0)) || (k == lanm) || it == maxit) {
        /*--------------------   solve eigen-problem for T(1:k,1:k) */
        /*                       vals in Rval, vecs in EvecT */
        //printf("k %d, trlen %d, Ntest %d, its %d\n", k, trlen, Ntest, it);
        SymEigenSolver(k, T, lanm1, EvecT, lanm1, Rval);
        /*-------------------- max dim reached-break from the inner loop */
        if (k == lanm || it == maxit) {
          break;
        }
        count = 0;
        /*-------------------- get residual norms and check acceptance of
          Ritz values for p(A). */
        //--------------------count e.vals in interval + those that convergeed
        jl = 0;
        for (i=0; i<k; i++) {
          if (Rval[i]>= bar) {
            jl++;
            r = fabs(beta*EvecT[i*lanm1+(k-1)]);
            resi[i] = r;
            if (r < tolP) {
              count++;
            }
          }
        }
        //-------------------- testing (partial) convergence for restart
        if (do_print)
          fprintf(fstats,"  --> testing conv k %4d, it %4d, count %3d  jl %3d trlen %3d\n",
              k, it, count, jl, trlen);
        /*-------------------- enough good candidates 1st time -> break */
        //if ((count*evFrac >= nev-lock) && (prtrlen==-1))
        if (count*evFrac >= nev-lock)
          break;
        /*-------------------- count & jl unchanged since last test --> break */
        if ((count<=last_count) && (jl<=last_jl))
          break;
        last_count = count;  
        last_jl = jl;
      }   //                 if -- [Restarting test] block
      /*-------------------- end of inner (Lanczos) loop - Next: restart*/        
    }     //                 while (k<mlan) loop

    /*--------------------   TWO passes to select good candidates */
    /*                       Pass-1: based on if ``p(Ritzvalue) > bar'' */	    
    jl = 0;
    for (i=0; i<k; i++) {
      /*--------------------   if this Ritz value is higher than ``bar'' */
      if (Rval[i] >= bar) {
        // move good eigenvectors/vals to front
        if (i != jl) {
          DCOPY(&k, EvecT+i*lanm1, &one, EvecT+jl*lanm1, &one);
          Rval[jl] = Rval[i];
          resi[jl] = resi[i];
        }
        jl++;
      }
    }
    /*   Compute the Ritz vectors: Rvec(:,1:jl) = V(:,1:k) * EvecT(:,1:jl) */
    DGEMM(&cN,&cN,&n,&jl,&k,&done,V,&n,EvecT,&lanm1,&dzero,Rvec,&n);
    /*--------------------  Pass-2: check if Ritz vals of A are in [a,b] */
    /*                      number of Ritz values in [a,b] */
    ll = 0;
    /*-------------------- trlen = # Ritz vals that will  go to TR set */
    prtrlen = trlen; 
    trlen = 0;
    for (i=0; i<jl; i++) {
      double *y = Rvec + i*n;
      double *w = work;
      /*--------------------   normalize just in case. */
      t = DNRM2(&n, y, &one); 
      // return  code 2 --> zero eigenvector found 
      if (t == 0.0)
        return 2;
      t = 1.0/t; 
      DSCAL(&n, &t, y, &one);
      /*--------------------   w = A*y */
      //-- matvec
      matvec(A, y, w);
      nmv ++;
      /*--------------------   Ritzval: t3 = (y'*w)/(y'*y) */
      //-- Rayleigh quotient 
      double t3 = DDOT(&n, y, &one, w, &one);
      /*--------------------  if lambda (==t3) is in [a,b] */
      if (t3 >= aa - DBL_EPSILON && t3 <= bb + DBL_EPSILON) {
        ll++;
        /*-------------------- compute residual wrt A for this pair */
        double nt3 = -t3;
        /*-------------------- w = w - t3*y */
        DAXPY(&n, &nt3, y, &one, w, &one);
        /*--------------------   res0 = norm(w) */
        res0 = DNRM2(&n, w, &one); 
        r = resi[i] ;
        /*-------------------- test res. of this Ritz pair against tol */
        if (r < tol) {
          //-------------------- check if need to realloc
          if (lock >= nev){
            nev += 1 + (int) (nev*nevInc);
            if (do_print) 
              fprintf(fstats, "-- More eigval found: realloc space for %d evs\n", nev);
            Realloc(Y, nev*n, double);
            Realloc(Lam, nev, double);
            Realloc(res, nev, double);
          }
          /*--------------------   accept (t3, y) */
          DCOPY(&n, y, &one, Y+lock*n, &one);
          Lam[lock] = t3;
          res[lock] = res0;
          lock++;
        } else {
          /*-------------------- restart; move Ritz pair for TR to front */
          Rval[trlen] = Rval[i];
          DCOPY(&n, y, &one, V+trlen*n, &one);
          /* special vector for TR that is the bottom row of 
           * eigenvectors of Tm */
          s[trlen] = beta * EvecT[i*lanm1+(k-1)];
          trlen ++;
        }
      }
    }/* end of 2nd pass */

    /*-------Note:   jl = #evs that passed the 1st test,
     *       ll = #evs that passed the 1st and the 2nd tests.
     *       These ll evs are either locked (accepted) 
     *       or put into trlen (candidates for next restarts)
     *       step when trlen = 0 last restart and ll=0 this time.
     *       this is a sort of confirmation that nothing is left. 
     *       another test may be added later to make it more rigorous.
     */       
    if (do_print) 
      fprintf(fstats,"it %4d:  nMV %7d, k %3d, jl %3d, ll %3d, lock %3d, trlen %3d\n",
          it, nmv, k, jl, ll, lock, trlen);
    /*-------------------- TESTs for stopping */
    if ((prtrlen == 0) && (ll==0)) {
      /*--------------------  It looks like: Nothing left to compute  */
      if (do_print) {
        fprintf(fstats, "--------------------------------------------------\n");
        fprintf(fstats, " --> No ev.s left to be computed\n");
      }
      break;
    }
    if (it == maxit){
      if (do_print) {
        fprintf(fstats, "--------------------------------------------------\n");
        fprintf(fstats, " --> Max its reached without convergence\n");
      }
      break;
    }
    /*-------------------- prepare to restart.  First zero out all T */
    memset(T, 0, lanm1*lanm1*sizeof(double));
    /* move starting vector vector V(:,k+1);  V(:,trlen+1) = V(:,k+1) */
    DCOPY(&n, V+k*n, &one, V+trlen*n, &one);
  }      /* outer loop (it) */

  if (do_print) {
    fprintf(fstats, "     Number of evals found = %d\n", lock);
    fprintf(fstats, "--------------------------------------------------\n");
  }
  /*-------------------- Done.  output : */
  *nev2 = lock;
  *vals = Lam;
  *W = Y;
  *resW = res;
  /*-------------------- free arrays */
  free(V);
  free(T);
  free(Rval);
  free(resi);
  free(EvecT);
  free(Rvec);
  free(s);
  free(work);
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

