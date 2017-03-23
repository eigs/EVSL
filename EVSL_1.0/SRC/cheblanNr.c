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
 *  @param vinit  initial  vector for Lanczos -- [optional]
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
 * ------------------------------------------------------------ */
int ChebLanNr(double *intv, int maxit, double tol, double *vinit, 
              polparams *pol, int *nevOut, double **lamo, double **Wo, 
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
  /* if users provided their own matvec function, matrix A will be ignored */
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  maxit = min(n, maxit);
  /*-------------------- polynomial filter  approximates the delta
                         function centered at 'gamB'. 
                         bar: a bar value to threshold Ritz values of p(A) */
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
  if (do_print) {
    fprintf(fstats, " ** Cheb Poly of deg = %d, gam = %.15e, bar: %.15e\n", 
            deg, gamB, bar);
  }
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
  double *V, *dT, *eT, *Z;
  Malloc(V, n*(maxit+1), double);
  if (evsldata.ifGenEv) {
    /* storage for Z = B * V */
    Malloc(Z, n*(maxit+1), double);
  } else {
    /* Z and V are the same */
    Z = V;
  }
  /*-------------------- diag. subdiag of Tridiagional matrix */
  Malloc(dT, maxit, double);
  Malloc(eT, maxit, double);
  double *Rvec, *Lam, *res, *EvalT, *EvecT;
  /*-------------------- Lam, Rvec: the converged (locked) Ritz values vecs*/
  Malloc(Lam, maxit, double);         // holds computed Ritz values
  Malloc(res, maxit, double);         // residual norms (w.r.t. ro(A))
  Malloc(EvalT, maxit, double);       // eigenvalues of tridia. matrix  T
  Malloc(EvecT, maxit*maxit, double); // Eigen vectors of T
  /*-------------------- nev = current number of converged e-pairs 
                         nconv = converged eigenpairs from looking at Tk alone */
  int nev, nconv = 0;
  /*-------------------- nmv counts  matvecs */
  int nmv = 0;
  /*-------------------- copy initial vector to V(:,1) */
  DCOPY(&n, vinit, &one, V, &one);
  /*--------------------  normalize it */
  double t, nt, res0;
  if (evsldata.ifGenEv) {
    /* B norm */
    matvec_B(V, Z);
    t = 1.0 / sqrt(DDOT(&n, V, &one, Z, &one));
    DSCAL(&n, &t, Z, &one);
  } else {
    /* 2-norm */
    t = 1.0 / DNRM2(&n, V, &one); // add a test here.
  }
  /* unit B-norm or 2-norm */
  DSCAL(&n, &t, V, &one);
  /*-------------------- u  is just a pointer. wk == work space */
  double *u, *wk, *w2;
  int wk_size = evsldata.ifGenEv ? 4*n : 3*n;
  Malloc(wk, wk_size, double);
  w2 = wk + n;
  /*-------------------- for ortho test */
  double wn = 0.0;
  int nwn = 0;
  /*-------------------- for stopping test [restricted trace]*/
  tr0 = 0;
  /*-------------------- lanczos vectors updated by rotating pointer*/
  /*-------------------- pointers to Lanczos vectors */
  double *zold, *z, *znew;
  double *v, *vnew;
  /*--------------------  Lanczos recurrence coefficients */
  double alpha, nalpha, beta=0.0, nbeta, resi;
  int count = 0;
  // ---------------- main Lanczos loop 
  for (k=0; k<maxit; k++) {
    /*-------------------- quick reference to Z(:,k-1) when k>0*/
    zold = k > 0 ? Z+(k-1)*n : NULL;
    /*-------------------- a quick reference to V(:,k) */
    v = &V[k*n];
    /*-------------------- a quick reference to Z(:,k) */
    z = &Z[k*n];
    /*-------------------- next Lanczos vector V(:,k+1)*/
    vnew = v + n;
    /*-------------------- next Lanczos vector Z(:,k+1)*/
    znew = z + n;
    /*-------------------- compute w = p[(A-cc)/dd] * v */
    tm = cheblan_timer();
    /*------------------ NOTE: z is used!!! [TODO: FIX ME] */
    ChebAv(pol, z, znew, wk);
    tmv += cheblan_timer() - tm;
    nmv += deg;
    /*------------------ znew = znew - beta*zold */
    if (zold) {
      nbeta = -beta;
      DAXPY(&n, &nbeta, zold, &one, znew, &one);
    }
    /*-------------------- alpha = znew'*v */
    alpha = DDOT(&n, v, &one, znew, &one);
    /*-------------------- T(k,k) = alpha */
    dT[k] = alpha;
    wn += fabs(alpha);
    /*-------------------- znew = znew - alpha*z */
    nalpha = -alpha;
    DAXPY(&n, &nalpha, z, &one, znew, &one);
    /*-------------------- FULL reortho to all previous Lan vectors */
    if (evsldata.ifGenEv) {
      /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
      CGS_DGKS2(n, k+1, NGS_MAX, Z, V, znew, wk);
      /* vnew = B \ znew */
      evsldata.Bsol->func(znew, vnew, evsldata.Bsol->data);
      /*-------------------- beta = (vnew, znew)^{1/2} */
      beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
    } else {
      /*   vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
      /*   beta = norm(vnew) */
      CGS_DGKS(n, k+1, NGS_MAX, V, vnew, &beta, wk);
    }
    /*-------------------- T(k+1,k) = alpha */
    eT[k] = beta;
    wn += 2.0 * beta;
    nwn += 3;
    /*-------------------- lucky breakdown test */
    if (beta*nwn < orthTol*wn) {
      if (do_print) {
        fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", k, beta);
      }
      /*------------------ generate a new init vector */
      rand_double(n, vnew);
      if (evsldata.ifGenEv) {
      /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
        CGS_DGKS2(n, k+1, NGS_MAX, V, Z, vnew, wk);        
        matvec_B(vnew, znew);
        beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
      } else {
        CGS_DGKS(n, k+1, NGS_MAX, V, vnew, &beta, wk);      
      }
    }
    /*-------------------- vnew = vnew / beta */
    t = 1.0 / beta;
    DSCAL(&n, &t, vnew, &one);
    if (evsldata.ifGenEv) {
      /*-------------------- znew = znew / beta */
      DSCAL(&n, &t, znew, &one);
    }
    /*---------------------- test for Ritz vectors */
    if ( (k < Ntest || (k-Ntest) % cycle != 0) && k != maxit-1 ) {
      continue;
    }
    /*---------------------- diagonalize  T(1:k,1:k)       
                             vals in EvalT, vecs in EvecT  */
    kdim = k+1;
#if 1
    /*-------------------- THIS uses dsetv */
    SymmTridEig(EvalT, EvecT, kdim, dT, eT);
    count = kdim;
#else
    /*-------------------- THIS uses dstemr */
    double vl = bar - DBL_EPSILON, vu = 2.0;  /* needed by SymmTridEigS */
    SymmTridEigS(EvalT, EvecT, kdim, vl, vu, &count, dT, eT);
#endif
    /*-------------------- restricted trace: used for convergence test */
    tr1 = 0;  
    /*-------------------- get residual norms and check acceptance of Ritz 
     *                     values for p(A). nconv records number of 
     *                     eigenvalues whose residual for p(A) is smaller 
     *                     than tol. */
    nconv = 0;
    for (i=0; i<count; i++) {
      flami = EvalT[i];
      if (flami + DBL_EPSILON >= bar) {
        tr1+= flami;
        /*---------------- the last row of EvecT: EvecT[i*kdim+kdim-1] */
        if (beta*fabs(EvecT[(i+1)*kdim-1]) < tol) {
          nconv++;
        }
      }
    }

    if (do_print) {
      fprintf(fstats, "k %4d:   nMV %8d, nconv %4d  tr1 %21.15e\n",
              k, nmv, nconv,tr1);
    }
    /* -------------------- simple test because all eigenvalues
                            are between gamB and ~1. */
    if (fabs(tr1-tr0) < tol*fabs(tr1)) {
      break;
    }
    tr0 = tr1;
  } /* end of the main loop */

  /*-------------------- done == compute Ritz vectors */
  Malloc(Rvec, nconv*n, double);       // holds computed Ritz vectors

  nev = 0;
  for (i=0; i<count; i++) {
    flami = EvalT[i];
    //-------------------- reject eigenvalue if rho(lam)<bar
    if (flami < bar) {
      continue;
    }
    y = &EvecT[i*kdim];
    /*-------------------- make sure to normalize */
    /*
    t = DNRM2(&kdim, y, &one);
    t = 1.0 / t;
    DSCAL(&kdim, &t, y, &one);
    */
    /*-------------------- residual norm for transformed Pb. */
    resi = beta*fabs(y[kdim-1]);
    if (resi > tol) {
      continue;
    }
    /*-------------------- compute Ritz vectors */
    u = &Rvec[nev*n];
    DGEMV(&cN, &n, &kdim, &done, V, &n, y, &one, &dzero, u, &one);
    /*-------------------- normalize u */
    if (evsldata.ifGenEv) {
      /* B-norm, w2 = B*u */
      matvec_B(u, w2);
      t = sqrt(DDOT(&n, u, &one, w2, &one)); /* should be one */
    } else {
      /* 2-norm */
      t = DNRM2(&n, u, &one); /* should be one */
    }
    /*-------------------- return code 2 --> zero eigenvector found */
    if (t == 0.0) {
      return 2;
    }
    /*-------------------- scal u */
    t = 1.0 / t;
    DSCAL(&n, &t, u, &one);
    if (evsldata.ifGenEv) {
      /*------------------ w2 = B*u */
      DSCAL(&n, &t, w2, &one);
    }
    /*-------------------- w = A*u */
    matvec_A(u, wk);
    nmv ++;
    /*-------------------- Ritz val: t = (u'*w)/(u'*u)
                                     t = (u'*w)/(u'*B*u) */
    t = DDOT(&n, wk, &one, u, &one);
    /*-------------------- if lambda (==t) is in [a,b] */
    if (t < aa - DBL_EPSILON || t > bb + DBL_EPSILON) {
      continue;
    }
    /*-------------------- compute residual wrt A for this pair */
    nt = -t;
    if (evsldata.ifGenEv) {
      /*-------------------- w = w - t*B*u */
      DAXPY(&n, &nt, w2, &one, wk, &one);
      /*-------------------- B norm */
      matvec_B(wk, w2);
      res0 = sqrt(DDOT(&n, wk, &one, w2, &one));
    } else {
      /*-------------------- w = w - t*u */
      DAXPY(&n, &nt, u, &one, wk, &one);
      /*-------------------- res0 = norm(w) */
      res0 = DNRM2(&n, wk, &one); 
    }
    /*--------------------   accept (t, y) */
    Lam[nev] = t;
    res[nev] = res0;
    nev++;
  }

  /*-------------------- Done.  output : */
  *nevOut = nev;
  *lamo = Lam;
  *Wo = Rvec;
  *reso = res;
  /*-------------------- free arrays */
  free(V);
  free(dT);
  free(eT);
  free(EvalT);
  free(EvecT);
  free(wk);
  if (evsldata.ifGenEv) {
    free(Z);
  }
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


