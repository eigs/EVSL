#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "blaslapack.h"
#include "def.h"
#include "evsl.h"
#include "internal_proto.h"
#include "string.h"  //for memset
#include "struct.h"

#define USE_RAND_NUM 0
double fsqrt(const double a) { return sqrt(1 / a); }

double rec(const double a) { return 1 / a; }

void ones(int n, double* v) { int i = 0; for(i = 0; i < n; i++) { v[i] = 1; }}


/**
 * @brief @b Computes y=P(A) y, where pn is a Cheb. polynomial expansion 
 * 
 * This explicitly calls matvec_B, so it can be useful for implementing
 * user-specific matrix-vector multiplication.
 *
 * @param mu Vector containing Chebyshev Coefficients
 * @param cc 
 * @param dd 
 * @param v input vector
 * @param v input vector
 * @param v input vector
 *
 * @param[out] y p(A)v
 *
 * @b Workspace
 * @param w Work vector of length 3*n [allocate before call]
 * @param v is untouched
 **/
int pnav(double* mu, const int m, const double cc, const double dd, double *v, double *y, double *w) {//Really just ChebAv
  printf("eln: %i", evsldata.n);
  int n = evsldata.n;
  printf("N: %i \n", n);
  printf("mu[0]: %f", mu[0]);
  printf("v[0]: %f", v[0]);
  printf("dd: %f", dd);
  printf("cc: %f", cc);
  printf("m: %i", m);
  /*-------------------- pointers to v_[k-1],v_[k], v_[k+1] from w */
  double *vk   = w;
  double *vkp1 = w+n;
  double *vkm1 = vkp1+n;
  double *w2 = evsldata.ifGenEv ? vkm1 + n : NULL;
  /*-------------------- */
  int k, i;
  double t, s, *tmp, t1= 1.0 / dd, t2 = 2.0 / dd; 
  /*-------------------- vk <- v; vkm1 <- zeros(n,1) */
  memcpy(vk, v, n*sizeof(double));
  memset(vkm1, 0, n*sizeof(double));
  /*-------------------- special case: k == 0 */
  s = mu[0];
  printf("s: %f \n", s);
  printf("vk0: %f \n", vk[0]);
  for (i=0; i<n; i++) {
    y[i] = s*vk[i];
  }
  printf("y0: %f \n", y[0]);
  printf("y1: %f \n", y[1]);
  printf("y2: %f \n", y[2]);
  /*-------------------- degree loop. k IS the degree */
  for (k=1; k<=m; k++) {
    /*-------------------- y = mu[k]*Vk + y */    
    t = k == 1 ? t1 : t2; 
    /*-------------------- */    
    s = mu[k];
      matvec_B(vk, vkp1);
    printf("vk0: %f \n", vk[0]);
    printf("vkp0: %f \n", vkp1[0]);
    printf("N: %i \n", n);

    for (i=0; i<n; i++) {
      vkp1[i] = t*(vkp1[i]-cc*vk[i]) - vkm1[i];
      /*-------------------- for degree 2 and up: */
      y[i] += s*vkp1[i];
    }
    /*-------------------- next: rotate vectors via pointer exchange */
    tmp = vkm1;
    vkm1 = vk;
    vk = vkp1;
    vkp1 = tmp;
  }
    printf("End y[0]: %f \n",  y[0]);
    printf("End y[%i]: %f \n", n, y[0]);

  return 0;
}





/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density)
 *
 *    @param[in] *A  -- used through calls to matvec_A
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] degB  degree used to approximate B
 *    @param[in] npts number of sample points used for he DOS curve
 *    @param[in] *intv Stores the  intervals of interest \\
 *    intv[3:4] = [a b] = interval where DOS is to be computed \\
 *    intv[5:6] = [lambda_min(B), lambda_max(B)] 
 *    Used only if polynomial
 *    approximation is used instead of Cholesky
 *
 *    @param[in] tau Tolerance used for the  approximation
 *
 *    @param[out] xdos Length-npts long vector, x-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] ydos Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] neig == estimated number of eigenvalues
 *
 *
 *----------------------------------------------------------------------*/


int LanDosG2(const int nvec, int msteps, const int degB, int npts, double *xdos, double *ydos,
           double *neig, const double *const intv, const double tau) {
  //--------------------

  int maxit = msteps;
  const double tol = tau;
  int n = evsldata.n;

  const int mdeg = 200;

  polparams pol_sqr;
  polparams pol_sol;


  double c_sqr, c_sol, h_sqr, h_sol;
  double *mu_sol, *mu_sqr;
  double *mu_sol_tmp, *mu_sqr_tmp;
  int *deg_sqr_tmp, *deg_sol_tmp;

  Malloc(pol_sqr.mu, mdeg, double);
  Malloc(pol_sol.mu, mdeg, double);
  Malloc(deg_sqr_tmp, 1, int);
  Malloc(deg_sol_tmp, 1, int);

  double* vinit;
  Malloc(vinit, n, double);
#if USE_RAND_NUM
  printf("Using rand numbers \n");
  randn_double(n,vinit);
#else
  ones(n,vinit);
  printf("Using predetermined numbers \n");
  printf("vinit[0]: %f \n", vinit[0]);
#endif


  if (degB <= 0) {
    printf("LanDos with degB <= 0 not yet implemented");
    exit(-1);
  }
  else {
    lsPol2(&intv[4], mdeg, fsqrt, tau, &pol_sqr);
    lsPol2(&intv[4], mdeg, rec, tau, &pol_sol);

  }
  const int deg_sol = deg_sol_tmp[0];
  const int deg_sqr = deg_sqr_tmp[0];
  free(deg_sol_tmp);
  free(deg_sqr_tmp);

  double tm,  tmv=0.0, tr0, tr1, tall;
  double *y, flami; 
  //-------------------- to report timings/
  tall = cheblan_timer();
  int i, k, kdim;
  // handle case where fstats is NULL. Then no output. Needed for openMP.
  int do_print = 1;   
#if 0
  if (fstats == NULL){
    do_print = 0;
  }  
#endif
  /*--------------------   frequently used constants  */
  char cN = 'N';   
  int one = 1;
  double done=1.0,dzero=0.0;
  /*--------------------   Ntest = when to start testing convergence */
  int Ntest = 30; 
  /*--------------------   how often to test */
  int cycle = 20; 
  /* size of the matrix */
  maxit = min(n, maxit);
  /*-------------------- Rational filter with pole at ((a+b)/2,(b-a)/2) with 
    multiplicity pow, bar value equals 1/2        */
  /*-------------------- a, b, used for testing only at end */
  double bar = 0.5;
  
#if 0
  if (check_intv(intv, fstats) < 0) {
    *nevOut = 0;
    *lamo = NULL; *Wo = NULL; *reso = NULL;
    return 0;
  }
#endif
  double aa = intv[0];
  double bb = intv[1];
  /*-----------------------------------------------------------------------* 
   * *Non-restarted* Lanczos iteration 
   *-----------------------------------------------------------------------*/
#if 0
  if (do_print) {
    fprintf(fstats, " ** Rat-LanNr \n");
  }
#endif
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
  //Malloc(EvecT, maxit*maxit, double); // Eigen vectors of T
  /*-------------------- nev = current number of converged e-pairs 
    nconv = converged eigenpairs from looking at Tk alone */
  int nev, nconv = 0;
  /*-------------------- nmv counts  matvecs */
  int nmv = 0;
  /*-------------------- nsv counts  solves */
  int nsv = 0;
  /*-------------------- u  is just a pointer. wk == work space */
  double *u, *wk, *w2, *vrand = NULL;
  int wk_size = evsldata.ifGenEv ? 6*n : 4*n;
  Malloc(wk, wk_size, double);
  w2 = wk + n;
  /*-------------------- copy initial vector to Z(:,1) */
#if 1
  /* Filter the initial vector */
  tm = cheblan_timer();
  printf("vinit[0]: %f \n", vinit[0]);
  pnav(pol_sqr.mu, pol_sqr.deg, pol_sqr.cc, pol_sqr.dd, vinit, V, wk); //Ish
  for(i = 0; i < pol_sqr.deg; i++) {
    printf("V[i]: %f \n", V[i]);
  }
  tmv += cheblan_timer() - tm;
#if 0
  nsv += deg;
#endif
  Malloc(vrand, n, double);
  if(evsldata.ifGenEv){
    DCOPY(&n, V, &one, Z, &one);
#if 0
    nmv += (deg-1);
#endif
  }
#else
  DCOPY(&n, vinit, &one, Z, &one);
#endif  
  /*--------------------  normalize it */
  double t, nt, res0;
  if (evsldata.ifGenEv) {
    /* B norm */
    matvec_B(Z, V);
    t = 1.0 / sqrt(DDOT(&n, V, &one, Z, &one));
    DSCAL(&n, &t, Z, &one);
    nmv++;
  } else {
    /* 2-norm */
    t = 1.0 / DNRM2(&n, V, &one); // add a test here.
  }
  /* unit B^{-1}-norm or 2-norm */
  DSCAL(&n, &t, V, &one);
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
  double alpha, nalpha, beta=0.0, nbeta;
  int count = 0;
  // ---------------- main Lanczos loop 
  //
  // ##Should be good through about here
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
    /*-------------------- compute w = rat * v */
    tm = cheblan_timer();
    pnav(pol_sol.mu, pol_sol.deg, pol_sol.cc, pol_sol.dd, v, znew, wk); //Ish
    tmv += cheblan_timer() - tm;
#if 0
    if (evsldata.ifGenEv) {
      nmv += (deg-1);
    }
    nsv += deg;
#endif
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
      /* -------------- NOTE: B-matvec 
       *                vnew = B * znew */
      matvec_B(znew, vnew);
      /*-------------------- beta = (vnew, znew)^{1/2} */
      beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
      nmv++;
    } else {
      /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
      /* beta = norm(vnew) */
      CGS_DGKS(n, k+1, NGS_MAX, V, vnew, &beta, wk);
    }
    wn += 2.0 * beta;
    nwn += 3;
    /*-------------------- lucky breakdown test */
    if (beta*nwn< orthTol*wn) {
#if 0
      if (do_print) {
        fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", k, beta);
      }
#endif
#if 1      
      /*------------------ generate a new init vector in znew */
#if USE_RAND_NUM
      rand_double(n, vrand);
#else
      ones(n, vrand);
#endif
      /* Filter the initial vector */
      tm = cheblan_timer();
      pnav(pol_sol.mu, pol_sol.deg, pol_sol.cc, pol_sol.dd, vrand, znew, wk); //Ish
      tmv += cheblan_timer() - tm;
#if 0
      nmv += (deg-1);
#endif
#else 
#if USE_RAND_NUM
      rand_double(n, znew);
#else
      ones(n,vrand);
#endif
#endif            
      if (evsldata.ifGenEv) {
	/* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
        CGS_DGKS2(n, k+1, NGS_MAX, Z, V, znew, wk);
	/* -------------- NOTE: B-matvec */
        matvec_B(znew, vnew);
        nmv ++;
        beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
	/*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        DSCAL(&n, &t, vnew, &one);
	/*-------------------- znew = znew / beta */
        DSCAL(&n, &t, znew, &one);
        beta = 0.0;
#if 0
        nsv += deg;        
#endif
      } else {
	/* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
	/* beta = norm(vnew) */
        CGS_DGKS(n, k+1, NGS_MAX, V, vnew, &beta, wk);
	/*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        DSCAL(&n, &t, vnew, &one);
        beta = 0.0;
      }
    } else {
      /*-------------------- vnew = vnew / beta */
      t = 1.0 / beta;
      DSCAL(&n, &t, vnew, &one);
      if (evsldata.ifGenEv) {
	/*-------------------- znew = znew / beta */
	DSCAL(&n, &t, znew, &one);
      }    
    }
    /*-------------------- T(k+1,k) = beta */
    eT[k] = beta;
    /*---------------------- test for Ritz vectors */
    if ( (k < Ntest || (k-Ntest) % cycle != 0) && k != maxit-1 ) {
      continue;
    }
    /*---------------------- diagonalize  T(1:k,1:k)       
      vals in EvalT, vecs in EvecT  */
    kdim = k+1;
#if 1
    /*-------------------- THIS uses dsetv, do not need eig vector */    
    SymmTridEig(EvalT, NULL, kdim, dT, eT);
    count = kdim;
#else
    /*-------------------- THIS uses dstemr */
    double vl = bar - DBL_EPSILON, vu = 10000.0;  /* needed by SymmTridEigS */
    SymmTridEigS(EvalT, EvecT, kdim, vl, vu, &count, dT, eT);
#endif
    /*-------------------- restricted trace: used for convergence test */
    tr1 = 0;  
    /*-------------------- get residual norms and check acceptance of Ritz 
     *                     values for r(A). nconv records number of 
     *                     eigenvalues whose residual for r(A) is smaller 
     *                     than tol. */
    nconv = 0;
    for (i=0; i<count; i++) {
      flami = EvalT[i];
      if (flami + DBL_EPSILON >= bar) {
        tr1+= flami;
        nconv++;
      }
    }

#if 0
    if (do_print) {
      fprintf(fstats, "k %4d:   # sols %8d, nconv %4d  tr1 %21.15e\n",
              k, nsv, nconv,tr1);
    }
#endif
    /* -------------------- simple test because all eigenvalues
       are between gamB and ~1. */
    if (fabs(tr1-tr0) < tol*fabs(tr1)) {
      break;
    }
    tr0 = tr1;
  } /* end of the main loop */

  /*-------------------- compute eig vals and vector */    
  Malloc(EvecT, kdim*kdim, double); // Eigen vectors of T
  SymmTridEig(EvalT, EvecT, kdim, dT, eT);
  
  /*-------------------- done == compute Ritz vectors */
  Malloc(Rvec, nconv*n, double);       // holds computed Ritz vectors

  nev = 0;
  for (i=0; i<count; i++) {
    flami = EvalT[i];
    //-------------------- reject eigenvalue if rho(lam)<bar
    if (fabs(flami) < bar) {
      continue;
    }
    y = &EvecT[i*kdim];
    /*-------------------- make sure to normalize */
    /*
      t = DNRM2(&kdim, y, &one);
      t = 1.0 / t;
      DSCAL(&kdim, &t, y, &one);
    */
    /*-------------------- compute Ritz vectors 
     *                     NOTE: use Z for gen e.v */
    u = &Rvec[nev*n];
    DGEMV(&cN, &n, &kdim, &done, Z, &n, y, &one, &dzero, u, &one);
    /*-------------------- normalize u */
    if (evsldata.ifGenEv) {
      /* B-norm, w2 = B*u */
      matvec_B(u, w2);
      nmv ++;
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
    /*-------------------- scal B*u */
    if (evsldata.ifGenEv) {
      /*-------------------- w2 = B*u */
      DSCAL(&n, &t, w2, &one);
    }
    /*-------------------- w = A*u */
    matvec_A(u, wk);
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
      /*-------------------- 2 norm of res */
      res0 = DNRM2(&n, wk, &one);
    } else {
      /*-------------------- w = w - t*u */
      DAXPY(&n, &nt, u, &one, wk, &one);
      /*-------------------- res0 = norm(w) */
      res0 = DNRM2(&n, wk, &one); 
    }
    /*--------------------   accept (t, y) */
    if (res0 < tol) {
       Lam[nev] = t;
       res[nev] = res0;
       nev++;
    }    
  }

  /*-------------------- Done.  output : */
#if 0
  *nevOut = nev;
  *lamo = Lam;
  *Wo = Rvec;
  *reso = res;
#endif
  /*-------------------- free arrays */
  free(V);
  free(dT);
  free(eT);
  free(EvalT);
  free(EvecT);
  free(wk);
  if (vrand) {
    free(vrand);
  }
  if (evsldata.ifGenEv) {
    free(Z);
  }
  /*-------------------- record stats */
  tall = cheblan_timer() - tall;
  /*-------------------- print stat */
#if 0
  if (do_print){
    fprintf(fstats, "------This slice consumed: \n");
    fprintf(fstats, "# of solves :        %d\n", nsv);
    if (evsldata.ifGenEv) {
      fprintf(fstats, "# of Matvec with B:  %d\n", nmv);
    }
    fprintf(fstats, "total time  :        %.2f\n", tall);
    fprintf(fstats, "solve time  :        %.2f\n", tmv);
    fprintf(fstats,"======================================================\n");
  }
#endif
  for(i = 0; i < n; i++) {
    printf("LanDos y[%i] : %f", i, y[i]);
    } 
  return 0;
}