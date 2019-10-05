#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "internal_header.h"

/**
 * @file cheblanNr.c
 * @brief Polynomial Filtered no-restart Lanczos
 */

/**
 * if filter the initial vector
 */
#define FILTER_VINIT 1

/**
 *  @brief Chebyshev polynomial filtering Lanczos process [NON-restarted version]
 *
 *  @param[in] intv     An array of length 4 \n
 *          [intv[0], intv[1]] is the interval of desired eigenvalues \n
 *          [intv[2], intv[3]] is the global interval of all eigenvalues \n
 *          Must contain all eigenvalues of A
 *
 *  @param[in] maxit    Max number of outer Lanczos steps  allowed --[max dim of Krylov
 *          subspace]
 *
 *  @param[in] tol
 *          Tolerance for convergence. The code uses a stopping criterion based
 *          on the convergence of the restricted trace. i.e., the sum of the
 *          eigenvalues of T_k that  are in the desired interval. This test  is
 *          rather simple since these eigenvalues are between `bar' and  1.0.
 *          We want the relative error on  this restricted  trace to be  less
 *          than  tol.  Note that the test  performed on filtered matrix only
 *          - *but* the actual residual norm associated with the original
 *          matrix A is returned
 *
 *  @param[in] vinit  initial  vector for Lanczos -- [optional]
 *
 *  @param[in] pol       A struct containing the parameters of the polynomial..
 *  This is set up by a call to find_deg prior to calling chenlanNr
 *
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
 *
 * ------------------------------------------------------------ */
int ChebLanNr(const double *intv, int maxit, double tol, const double *vinit,
              const polparams *pol, int *nevOut, double **lamo, double **Wo,
              double **reso, FILE *fstats) {
  //-------------------- to report timings/
  double tall, tm1 = 0.0, tt;
  tall = evsl_timer();
  const int ifGenEv = evsldata.ifGenEv;
  double tr0, tr1;
  double *y, flami;
  int i, k, kdim;
  /* handle case where fstats is NULL. Then no output. Needed for openMP */
  int do_print = 1;
  if (fstats == NULL) {
    do_print = 0;
  }
  /*-------------------- frequently used constants  */
  const char cN = 'N';
  const int one = 1;
  const double done = 1.0;
  const double dzero = 0.0;
  /*-------------------- Ntest = when to start testing convergence */
  const int Ntest = 30;
  /*-------------------- how often to test */
  const int cycle = 20;
  /* size of the matrix */
  const int n = evsldata.n;
  /* max num of its */
  maxit = evsl_min(n, maxit);
  /*-------------------- Caveat !!!: To prevent integer overflow, we save
   *                     another size_t version of n
   *                     Note n itself may not overflow but something like
   *                     (maxit * n) is much more likely
   *                     When using size_t n, in all the multiplications
   *                     integer is promoted to size_t */
  size_t n_l = n;
  /*-------------------- polynomial filter  approximates the delta
                         function centered at 'gamB'.
                         bar: a bar value to threshold Ritz values of p(A) */
  const double bar = pol->bar;
  const double gamB = pol->gam;
  /*-------------------- interval [aa, bb], used for testing only at end */
  if (check_intv(intv, fstats) < 0) {
    *nevOut = 0;
    *lamo = NULL; *Wo = NULL; *reso = NULL;
    return 0;
  }
  const double aa = intv[0];
  const double bb = intv[1];
  const int deg = pol->deg;
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
  if (do_print) {
    fprintf(fstats, " ** Cheb-LanNr \n");
  }
  /*-------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *dT, *eT, *Z;
  V = evsl_Malloc(n_l*(maxit+1), double);
  if (ifGenEv) {
    /* storage for Z = B * V */
    Z = evsl_Malloc(n_l*(maxit+1), double);
  } else {
    /* Z and V are the same */
    Z = V;
  }
  /*-------------------- diag. subdiag of Tridiagional matrix */
  dT = evsl_Malloc(maxit, double);
  eT = evsl_Malloc(maxit, double);
  double *Rvec, *Lam, *res, *EvalT, *EvecT;
  /*-------------------- Lam, Rvec: the converged (locked) Ritz values vecs*/
  Lam = evsl_Malloc(maxit, double);       // holds computed Ritz values
  res = evsl_Malloc(maxit, double);       // residual norms (w.r.t. ro(A))
  EvalT = evsl_Malloc(maxit, double);       // eigenvalues of tridiag matrix T
  /*-------------------- nev = current number of converged e-pairs
                         nconv = converged eigenpairs from looking at Tk alone */
  int nev, nconv = 0;
  /*-------------------- u is just a pointer. wk == work space */
  double *u, *wk, *w2, *vrand = NULL;
  size_t wk_size = ifGenEv ? 4*n_l : 3*n_l;
  wk = evsl_Malloc(wk_size, double);
  w2 = wk + n;
  /*-------------------- copy initial vector to Z(:,1) */
#if FILTER_VINIT
  /*-------------------- compute w = p[(A-cc)/dd] * v */
  /*-------------------- Filter the initial vector */
  ChebAv(pol, vinit, V, wk);
  vrand = evsl_Malloc(n, double);
#else
  /*-------------------- copy initial vector to V(:,1) */
  evsl_dcopy(&n, vinit, &one, V, &one);
#endif
  /*-------------------- normalize it */
  double t, nt, res0;
  if (ifGenEv) {
    /* B norm */
    matvec_B(V, Z);
    t = 1.0 / sqrt(evsl_ddot(&n, V, &one, Z, &one));
    evsl_dscal(&n, &t, Z, &one);
  } else {
    /* 2-norm */
    t = 1.0 / evsl_dnrm2(&n, V, &one); // add a test here.
  }
  /* unit B-norm or 2-norm */
  evsl_dscal(&n, &t, V, &one);
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
  for (k=0; k<maxit; k++) {
    /*-------------------- quick reference to Z(:,k-1) when k>0 */
    zold = k > 0 ? Z+(k-1)*n_l : NULL;
    /*-------------------- a quick reference to V(:,k) */
    v = &V[k*n_l];
    /*-------------------- a quick reference to Z(:,k) */
    z = &Z[k*n_l];
    /*-------------------- next Lanczos vector V(:,k+1)*/
    vnew = v + n;
    /*-------------------- next Lanczos vector Z(:,k+1)*/
    znew = z + n;
    /*-------------------- compute w = p[(A-cc)/dd] * v */
    /*-------------------- NOTE: z is used!!! [TODO: FIX ME] */
    ChebAv(pol, z, znew, wk);
    /*------------------ znew = znew - beta*zold */
    if (zold) {
      nbeta = -beta;
      evsl_daxpy(&n, &nbeta, zold, &one, znew, &one);
    }
    /*-------------------- alpha = znew'*v */
    alpha = evsl_ddot(&n, v, &one, znew, &one);
    /*-------------------- T(k,k) = alpha */
    dT[k] = alpha;
    wn += fabs(alpha);
    /*-------------------- znew = znew - alpha*z */
    nalpha = -alpha;
    evsl_daxpy(&n, &nalpha, z, &one, znew, &one);
    /*-------------------- FULL reortho to all previous Lan vectors */
    if (ifGenEv) {
      /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
      CGS_DGKS2(n, k+1, NGS_MAX, Z, V, znew, wk);
      /* vnew = B \ znew */
      solve_B(znew, vnew);
      /*-------------------- beta = (vnew, znew)^{1/2} */
      beta = sqrt(evsl_ddot(&n, vnew, &one, znew, &one));
    } else {
      /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
      /* beta = norm(vnew) */
      CGS_DGKS(n, k+1, NGS_MAX, V, vnew, &beta, wk);
    }
    wn += 2.0 * beta;
    nwn += 3;
    /*-------------------- lucky breakdown test */
    if (beta*nwn < orthTol*wn) {
      if (do_print) {
        fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", k, beta);
      }
#if FILTER_VINIT
      /*------------------ generate a new init vector */
      rand_double(n, vrand);
      /*------------------  Filter the initial vector*/
      ChebAv(pol, vrand, vnew, wk);
#else
      rand_double(n, vnew);
#endif
      if (ifGenEv) {
        /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
        CGS_DGKS2(n, k+1, NGS_MAX, V, Z, vnew, wk);
        matvec_B(vnew, znew);
        beta = sqrt(evsl_ddot(&n, vnew, &one, znew, &one));
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        evsl_dscal(&n, &t, vnew, &one);
	/*-------------------- znew = znew / beta */
        evsl_dscal(&n, &t, znew, &one);
        beta = 0.0;
      } else {
        /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
        /* beta = norm(vnew) */
        CGS_DGKS(n, k+1, NGS_MAX, V, vnew, &beta, wk);
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        evsl_dscal(&n, &t, vnew, &one);
        beta = 0.0;
      }
    } else {
      /*-------------------- vnew = vnew / beta */
      t = 1.0 / beta;
      evsl_dscal(&n, &t, vnew, &one);
      if (ifGenEv) {
        /*-------------------- znew = znew / beta */
        evsl_dscal(&n, &t, znew, &one);
      }
    }
    /*-------------------- T(k+1,k) = beta */
    eT[k] = beta;
#if 0
    /*-------------------- re-allocate memory if maxit is smaller than # of eigs */
    if (k == maxit-1) {
      maxit = 1 + (int) (maxit * 1.5);
      V = evsl_Realloc(V, (maxit+1)*n_l, double);
      if (ifGenEv) {
        Z = evsl_Realloc(Z, (maxit+1)*n_l, double);
      } else {
        /* make sure Z == V since V may be changed in the re-alloc above */
        Z = V;
      }
      dT = evsl_Realloc(dT,    maxit, double);
      eT = evsl_Realloc(eT,    maxit, double);
      Lam = evsl_Realloc(Lam,   maxit, double);
      rest = evsl_Realloc(res,   maxit, double);
      EvalT = evsl_Realloc(EvalT, maxit, double);
    }
#endif
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
    //Used for convergence test
    for (i=0; i<count; i++) {
      flami = EvalT[i];
      if (flami + EVSL_DBL_EPS_MULT * DBL_EPSILON >= bar) {
        tr1+= flami;
        nconv++;
      }
    }

    if (do_print) {
      fprintf(fstats, "k %4d: nconv %4d  tr1 %21.15e\n",
              k, nconv,tr1);
    }
    /* -------------------- simple test because all eigenvalues
                            are between gamB and ~1. */
    if ( (fabs(tr1-tr0) < 2e-12) || (fabs(tr1)+fabs(tr0) < 1e-10) ) {
      break;
    }
    tr0 = tr1;
  } /* end of the main loop */

  if (k >= maxit) {
     fprintf(fstats, "The max number of iterations [%d] has been reached. The eigenvalues computed may not have converged\n", maxit);
  }

  /*-------------------- compute eig vals and vector */
  size_t kdim_l = kdim; /* just in case if kdim is > 65K */
  EvecT = evsl_Malloc(kdim_l*kdim_l, double); // Eigen vectors of T
  SymmTridEig(EvalT, EvecT, kdim, dT, eT);

  tt = evsl_timer();
  /*-------------------- done == compute Ritz vectors */
  Rvec = evsl_Malloc(nconv*n_l, double);       // holds computed Ritz vectors

  nev = 0;
  for (i=0; i<count; i++) {
    flami = EvalT[i];
    //-------------------- reject eigenvalue if rho(lam)<bar
    if (flami < bar) {
      continue;
    }
    y = &EvecT[i*kdim_l];
    /*-------------------- make sure to normalize */
    /*
    t = evsl_dnrm2(&kdim, y, &one);
    t = 1.0 / t;
    evsl_dscal(&kdim, &t, y, &one);
    */
    /*-------------------- compute Ritz vectors */
    u = &Rvec[nev*n_l];
    evsl_dgemv(&cN, &n, &kdim, &done, V, &n, y, &one, &dzero, u, &one);
    /*-------------------- normalize u */
    if (ifGenEv) {
      /* B-norm, w2 = B*u */
      matvec_B(u, w2);
      t = sqrt(evsl_ddot(&n, u, &one, w2, &one)); /* should be one */
    } else {
      /* 2-norm */
      t = evsl_dnrm2(&n, u, &one); /* should be one */
    }
    /*-------------------- return code 2 --> zero eigenvector found */
    if (t == 0.0) {
      return 2;
    }
    /*-------------------- scal u */
    t = 1.0 / t;
    evsl_dscal(&n, &t, u, &one);
    /*-------------------- scal B*u */
    if (ifGenEv) {
      /*------------------ w2 = B*u */
      evsl_dscal(&n, &t, w2, &one);
    }
    /*-------------------- w = A*u */
    matvec_A(u, wk);
    /*-------------------- Ritz val: t = (u'*w)/(u'*u)
                                     t = (u'*w)/(u'*B*u) */
    t = evsl_ddot(&n, wk, &one, u, &one);
    /*-------------------- if lambda (==t) is in [a,b] */
    if (t < aa - EVSL_DBL_EPS_MULT * DBL_EPSILON || t > bb + EVSL_DBL_EPS_MULT * DBL_EPSILON) {
      continue;
    }
    /*-------------------- compute residual wrt A for this pair */
    nt = -t;
    if (ifGenEv) {
      /*-------------------- w = w - t*B*u */
      evsl_daxpy(&n, &nt, w2, &one, wk, &one);
    } else {
      /*-------------------- w = w - t*u */
      evsl_daxpy(&n, &nt, u, &one, wk, &one);
    }
    /*-------------------- res0 = 2-norm(wk) */
    res0 = evsl_dnrm2(&n, wk, &one);
    /*--------------------   accept (t, y) */
    if (res0 < tol) {
      Lam[nev] = t;
      res[nev] = res0;
      nev++;
    }
  }
  tm1 = evsl_timer() - tt;

  /*-------------------- Done.  output : */
  *nevOut = nev;
  *lamo = Lam;
  *Wo = Rvec;
  *reso = res;
  /*-------------------- free arrays */
  evsl_Free(V);
  evsl_Free(dT);
  evsl_Free(eT);
  evsl_Free(EvalT);
  evsl_Free(EvecT);
  evsl_Free(wk);
  if (vrand) {
    evsl_Free(vrand);
  }
  if (ifGenEv) {
    evsl_Free(Z);
  }
  /*-------------------- record stats */
  tall = evsl_timer() - tall;

  evslstat.t_iter = tall;
  evslstat.t_ritz = tm1;

  return 0;
}

