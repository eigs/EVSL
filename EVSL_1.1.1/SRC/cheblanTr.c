#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "internal_header.h"

/**
 * @file cheblanTr.c
 * @brief Polynomial Filtered thick restart Lanczos
 */
/**
 * if filter the initial vector
 */
#define FILTER_VINIT 1

/**
 * @brief Chebyshev polynomial filtering Lanczos process [Thick restart version]
 *
 * @param[in] lanm      Dimension of Krylov subspace [restart dimension]
 * @param[in] nev       Estimate of number of eigenvalues in the interval --
 *         ideally nev == exact number or a little larger.  This is not used
 *         for testing convergence but it helps make decisions as to when to
 *         test convergence ChebLanTr attempts to compute *all* eigenvalues in
 *         the interval and stops only when no more eigenvalyes are left. The
 *         convergenve test is a very simple one based on the residual norm for
 *         the filtered matrix
 *
 * @param[in] intv   an array of length 4  \n
 *         [intv[0], intv[1]] is the interval of desired eigenvalues \n
 *         [intv[2], intv[3]] is the global interval of all eigenvalues \n
 *         it must contain all eigenvalues of A
 *
 * @param[in] maxit  max Num of outer Lanczos iterations (restarts) allowed --
 *         Each restart may or use the full lanm lanczos steps or fewer.
 *
 * @param[in] tol       tolerance for convergence. stop when ||res||< tol
 * @param[in] vinit     initial  vector for Lanczos -- [optional]
 * @param[in] pol       a struct containing the parameters of the polynomial. This
 *                  is set up by a call to find_deg prior to calling chenlanTr
 *
 * @param[out] nev2     Number of eigenvalues/vectors computed
 * @param[out] vals     Associated eigenvalues [nev2 x 1 vector]
 * @param[out] W        A set of eigenvectors  [n x nev2 matrix]
 *                      of unit 2-norm for standard eig prob
 *                      of unit B-norm for generalized eig prob
 * @param[out] resW     Associated residual norms [nev x 1 vector]
 *                      2-norm for standard eig prob
 *                      B-norm for generalized eig prob
 * @param[out] fstats   File stream which stats are printed to
 *
 * @return Returns 0 on success (or if check_intv() is non-positive), -1 if
 * gamB is outside [-1, 1], and 2 if there are no eigenvalues found.
 *
 *
 * @warning memory allocation for W/vals/resW within this function
 *
 **/
int ChebLanTr(int lanm, int nev, double *intv, int maxit,
              double tol, double *vinit, polparams *pol, int *nev2,
              double **vals, double **W, double **resW, FILE *fstats) {
  const int ifGenEv = evsldata.ifGenEv;
  /*-------------------- for stats */
  double tall=0.0;
  //double tolP = tol;
  double tr, last_tr;
  tall = evsl_timer();
  int do_print = 1;
  /* handle case where fstats is NULL. Then no output. Needed for openMP. */
  if (fstats == NULL) {
    do_print = 0;
  }
  /* size of the matrix */
  int n = evsldata.n;
  size_t n_l = n;
  /*--------------------- adjust lanm and maxit */
  lanm = evsl_min(lanm, n);
  int lanm1=lanm+1;
  size_t lanm1_l = lanm1;
  /*  if use full lanczos, should not do more than n iterations */
  if (lanm == n) {
    maxit = evsl_min(maxit, n);
  }
  /*-------------------- this is needed to increment the space when we
                         discover more than nev eigenvalues in interval */
  double nevInc = 0.2;   /* add 1  + 20% each time it is needed */
  /*-------------------- if we have at least nev/ev_frac good candidate
                         eigenvalues from p(A) == then we restart to lock them in */
  //int evFrac = 2;
  /*--------------------   some constants frequently used */
  /* char cT='T'; */
  char cN = 'N';
  int one = 1;
  double done=1.0, dmone=-1.0, dzero=0.0;
  /*-------------------- Ntest = when to start testing convergence */
  int Ntest = evsl_min(lanm, nev+50);
  /*--------------------   how often to test */
  int cycle = 50;
  int i, ll, /* count, last_count,*/ jl, last_jl;
  /*-----------------------------------------------------------------------
    -----------------------------------------------------------------------*/
  /* check if the given interval is valid */
  if (check_intv(intv, fstats) < 0) {
    *nev2 = 0;
    *vals = NULL; *W = NULL; *resW = NULL;
    return 0;
  }
  double aa = intv[0];
  double bb = intv[1];
  /*-------------------- Assumption: polynomial pol computed before calling cheblanTr
                         pol.  approximates the delta function centered at 'gamB'
                         bar: a bar value to threshold Ritz values of p(A) */
  //int deg = pol->deg;
  double gamB=pol->gam, bar=pol->bar;
  /*-------------------- gamB must be within [-1, 1] */
  if (gamB > 1.0 || gamB < -1.0) {
    if (fstats) {
      fprintf(fstats, "gamB error %.15e\n", gamB);
    }
    return -1;
  }
  /*-----------------------------------------------------------------------*
   * *thick restarted* Lanczos step
   *-----------------------------------------------------------------------*/
  if (do_print) {
    fprintf(fstats, " Cheb-LanTR, dim %d  cycle %d  \n",lanm, cycle);
  }
  /*--------------------- the min number of steps to be performed for
   *                      each innter loop, must be >= 1 - if not,
   *                      it means that the Krylov dim is too small */
  int min_inner_step = 5;
  /*-------------------- it = number of Lanczos steps */
  int it = 0;
  /*-------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *T;
  V = evsl_Malloc_device(n_l*lanm1, double);
  /*-------------------- for gen eig prob, storage for Z = B * V */
  double *Z;
  if (ifGenEv) {
    Z = evsl_Malloc_device(n_l*lanm1, double);
  } else {
    Z = V;
  }
  /*-------------------- T must be zeroed out initially */
  T = evsl_Calloc(lanm1_l*lanm1_l, double);
  /*-------------------- Lam, Y: the converged (locked) Ritz values/vectors
                         res: related residual norms */
  double *Y, *Lam, *res;
  Y = evsl_Malloc_device(n_l*nev, double);
  Lam = evsl_Malloc(nev, double);
  res = evsl_Malloc(nev, double);
  double *BY = NULL;
  /*-------------------- for gen eig prob, storage for B*Y */
  if (ifGenEv) {
    BY = evsl_Malloc_device(n_l*nev, double);
  }
  /*-------------------- lock =  number of locked vectors */
  int lock = 0;
  /*-------------------- trlen = dim. of thick restart set */
  int trlen = 0, prtrlen=-1;
  /*-------------------- Ritz values and vectors of p(A) */
  double *Rval, *Rvec, *resi, *BRvec = NULL;
  Rval = evsl_Malloc(lanm, double);
  resi = evsl_Malloc(lanm, double);
  Rvec = evsl_Malloc_device(n_l*lanm, double);
  if (ifGenEv) {
    BRvec = evsl_Malloc_device(n_l*lanm, double);
  }
  /*-------------------- Eigen vectors of T */
  double *EvecT, *EvecT_device;
  EvecT = evsl_Malloc(lanm1_l*lanm1_l, double);
#ifdef EVSL_USING_CUDA_GPU
  EvecT_device = evsl_Malloc_device(lanm1_l*lanm1_l, double);
#else
  EvecT_device = EvecT;
#endif
  /*-------------------- s used by TR (the ``spike'' of 1st block in Tm)*/
  double *s;
  s = evsl_Malloc(lanm, double);
  /*-------------------- alloc some work space */
  double *work, *vrand = NULL;
  size_t work_size = ifGenEv ? 4*n_l : 3*n_l;
  work = evsl_Malloc_device(work_size, double);
#if FILTER_VINIT
  /*------------------  Filter the initial vector*/
  ChebAv(pol, vinit, V, work);
  vrand = evsl_Malloc_device(n, double);
#else
  /*-------------------- copy initial vector to V(:,1)   */
  evsl_dcopy_device(&n, vinit, &one, V, &one);
#endif
  /*-------------------- normalize it */
  double t;
  if (ifGenEv) {
    /* B norm */
    matvec_B(V, Z);
    t = 1.0 / sqrt(evsl_ddot_device(&n, V, &one, Z, &one));
    /* z = B*v */
    evsl_dscal_device(&n, &t, Z, &one);
  } else {
    /* 2-norm */
    t = 1.0 / evsl_dnrm2_device(&n, V, &one);
  }
  /* unit B-norm or 2-norm */
  evsl_dscal_device(&n, &t, V, &one);
  /*-------------------- main (restarted Lan) outer loop */
  while (it < maxit) {
    /*-------------------- for ortho test */
    double wn = 0.0;
    int nwn = 0;
    /*  beta */
    double beta = 0.0, r, res0;
    /*  start with V(:,k) */
    int k = trlen > 0 ? trlen + 1 : 0;
    /* ! add a test if dimension exceeds (m+1)
     * (trlen + 1) + min_inner_step <= lanm + 1 */
    if (k+min_inner_step > lanm1) {
      if (fstats) {
        fprintf(fstats, "Krylov dim too small for this problem. Try a larger dim\n");
      }
      exit(1);
    }
    /*-------------------- thick restart special step */
    if (trlen > 0) {
      int k1 = k-1;
      /*------------------ a quick reference to V(:,k) */
      double *v = &V[k1*n_l];
      double *z = &Z[k1*n_l];
      /*------------------ next Lanczos vector */
      double *vnew = v + n;
      double *znew = z + n;
      /*------------------ znew = p[(A-cc)/dd] * v */
      /*------------------ NOTE: z is used!!! [TODO: FIX ME] */
      ChebAv(pol, z, znew, work);
      /*------------------ deflation */
      if (lock > 0) {
        if (ifGenEv) {
          /* orthogonalize against locked vectors first, w = w - B*Y*Y'*w */
          CGS_DGKS2(n, lock, NGS_MAX, BY, Y, znew, work);
        } else {
          /* orthogonalize against locked vectors first, w = w - Y*Y'*w */
          CGS_DGKS(n, lock, NGS_MAX, Y, vnew, NULL, work);
        }
      }
      /*-------------------- restart with 'trlen' Ritz values/vectors
                             T = diag(Rval(1:trlen)) */
      for (i=0; i<trlen; i++) {
        T[i*lanm1_l+i] = Rval[i];
        wn += fabs(Rval[i]);
      }
      /*--------------------- s(k) = V(:,k)'* znew */
      s[k1] = evsl_ddot_device(&n, v, &one, znew, &one);
      /*--------------------- znew = znew - Z(:,1:k)*s(1:k) */
      evsl_dgemv_device(&cN, &n, &k, &dmone, Z, &n, s, &one, &done, znew, &one);
      /*-------------------- expand T matrix to k-by-k, arrow-head shape
                             T = [T, s(1:k-1)] then T = [T; s(1:k)'] */
      for (i=0; i<k1; i++) {
        T[trlen*lanm1_l+i] = s[i];
        T[i*lanm1_l+trlen] = s[i];
        wn += 2.0 * fabs(s[i]);
      }
      T[trlen*lanm1_l+trlen] = s[k1];
      wn += fabs(s[k1]);
      if (ifGenEv) {
        /*-------------------- vnew = B \ znew */
        solve_B(znew, vnew);
        /*-------------------- beta = (vnew, znew)^{1/2} */
        beta = sqrt(evsl_ddot_device(&n, vnew, &one, znew, &one));
      } else {
        /*-------------------- beta = norm(w) */
        beta = evsl_dnrm2_device(&n, vnew, &one);
      }
      wn += 2.0 * beta;
      nwn += 3*k;
      /*   beta ~ 0 */
      if (beta*nwn < orthTol*wn) {
        if (do_print) {
          fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", k, beta);
        }
#if FILTER_VINIT
        /* filter random vector */
        rand_double_device(n, vrand);
        ChebAv(pol, vrand, vnew, work);
#else
        rand_double_device(n, vnew);
#endif
        if (ifGenEv) {
          /* orthogonalize against locked vectors first, v = v - Y*(B*Y)'*v */
          CGS_DGKS2(n, lock, NGS_MAX, Y, BY, vnew, work);
          /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
          CGS_DGKS2(n, k, NGS_MAX, V, Z, vnew, work);
          matvec_B(vnew, znew);
          beta = sqrt(evsl_ddot_device(&n, vnew, &one, znew, &one));
          double ibeta = 1.0 / beta;
          evsl_dscal_device(&n, &ibeta, vnew, &one);
          evsl_dscal_device(&n, &ibeta, znew, &one);
          beta = 0.0;
        } else {
          /* orthogonalize against locked vectors first, w = w - Y*Y'*w */
          CGS_DGKS(n, lock, NGS_MAX, Y, vnew, NULL, work);
          /*   vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
          /*   beta = norm(w) */
          CGS_DGKS(n, k, NGS_MAX, V, vnew, &beta, work);
          double ibeta = 1.0 / beta;
          evsl_dscal_device(&n, &ibeta, vnew, &one);
          beta = 0.0;
        }
      } else {
        /*------------------- w = w / beta */
        double ibeta = 1.0 / beta;
        evsl_dscal_device(&n, &ibeta, vnew, &one);
        if (ifGenEv) {
          evsl_dscal_device(&n, &ibeta, znew, &one);
        }
      }
      /*------------------- T(k+1,k) = beta; T(k,k+1) = beta; */
      T[k1*lanm1_l+k] = beta;
      T[k*lanm1_l+k1] = beta;
    } /* if (trlen > 0) */
    /*-------------------- Done with TR step. Rest of Lanczos step */
    /*-------------------- reset Ntest at each restart. */
    Ntest = evsl_max(20,nev-lock+10);
    /*last_count = 0;*/  last_jl = 0;  last_tr = 0.0;
    /*-------------------- regardless of trlen, *(k+1)* is the current
     *                     number of Lanczos vectors in V */
    /*-------------------- pointer to the previous Lanczos vector */
    double *zold = k > 0 ? Z+(k-1)*n_l : NULL;
    /*------------------------------------------------------*/
    /*------------------ Lanczos inner loop ----------------*/
    /*------------------------------------------------------*/
    while (k < lanm && it < maxit) {
      k++;
      /*---------------- a quick reference to V(:,k) */
      double *v = &V[(k-1)*n_l];
      double *z = &Z[(k-1)*n_l];
      /*---------------- next Lanczos vector */
      double *vnew = v + n;
      double *znew = z + n;
      /*------------------ znew = p[(A-cc)/dd] * v */
      /*------------------ NOTE: z is used!!! [TODO: FIX ME] */
      ChebAv(pol, z, znew, work);
      it++;
      /*-------------------- deflation: orthogonalize vs locked ones first */
      if (lock > 0) {
        if (ifGenEv) {
          /* orthogonalize against locked vectors first, znew = znew - B*Y*Y'*znew */
          CGS_DGKS2(n, lock, NGS_MAX, BY, Y, znew, work);
        } else {
          /*--------------------   vnew = vnew - Y*Y'*vnew */
          CGS_DGKS(n, lock, NGS_MAX, Y, vnew, NULL, work);
        }
      }
      /*-------------------- znew = znew - beta*zold */
      if (zold) {
        double nbeta = -beta;
        evsl_daxpy_device(&n, &nbeta, zold, &one, znew, &one);
      }
      /*-------------------- alpha = znew'*v */
      double alpha = evsl_ddot_device(&n, v, &one, znew, &one);
      /*-------------------- T(k,k) = alpha */
      T[(k-1)*lanm1_l+(k-1)] = alpha;
      wn += fabs(alpha);
      /*-------------------- znew = znew - alpha*z */
      double nalpha = -alpha;
      evsl_daxpy_device(&n, &nalpha, z, &one, znew, &one);
      /*-------------------- FULL reortho to all previous Lan vectors */
      if (ifGenEv) {
        /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
        CGS_DGKS2(n, k, NGS_MAX, Z, V, znew, work);
        /*-------------------- vnew = B \ znew */
        solve_B(znew, vnew);
        /*-------------------- beta = (vnew, znew)^{1/2} */
        beta = sqrt(evsl_ddot_device(&n, vnew, &one, znew, &one));
      } else {
        /*   vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
        /*   beta = norm(w) */
        CGS_DGKS(n, k, NGS_MAX, V, vnew, &beta, work);
      }
      wn += 2.0 * beta;
      nwn += 3;
      /*-------------------- zold = z */
      zold = z;
      /*-------------------- lucky breakdown test */
      if (beta*nwn < orthTol*wn) {
        if (do_print) {
          fprintf(fstats, "it %4d: Lucky breakdown, beta = %.15e\n", it, beta);
        }
#if FILTER_VINIT
        /* filter random vector */
        rand_double_device(n, vrand);
        ChebAv(pol, vrand, vnew, work);
#else
        /* generate a new init vector*/
        rand_double_device(n, vnew);
#endif
        if (ifGenEv) {
          /* orthogonalize against locked vectors first, w = w - Y*(B*Y)'*w */
          CGS_DGKS2(n, lock, NGS_MAX, Y, BY, vnew, work);
          /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
          CGS_DGKS2(n, k, NGS_MAX, V, Z, vnew, work);
          matvec_B(vnew, znew);
          beta = sqrt(evsl_ddot_device(&n, vnew, &one, znew, &one));
          double ibeta = 1.0 / beta;
          evsl_dscal_device(&n, &ibeta, vnew, &one);
          evsl_dscal_device(&n, &ibeta, znew, &one);
          beta = 0.0;
        } else {
          /* orthogonalize against locked vectors first, w = w - Y*Y'*w */
          CGS_DGKS(n, lock, NGS_MAX, Y, vnew, NULL, work);
          /*   vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
          /*   beta = norm(w) */
          CGS_DGKS(n, k, NGS_MAX, V, vnew, &beta, work);
          double ibeta = 1.0 / beta;
          evsl_dscal_device(&n, &ibeta, vnew, &one);
          beta = 0.0;
        }
      } else {
        /*---------------------- vnew = vnew / beta */
        double ibeta = 1.0 / beta;
        evsl_dscal_device(&n, &ibeta, vnew, &one);
        if (ifGenEv) {
          /*-------------------- znew = znew / beta */
          evsl_dscal_device(&n, &ibeta, znew, &one);
        }
      }
      /*-------------------- T(k,k+1) = T(k+1,k) = beta */
      T[k*lanm1_l+(k-1)] = beta;
      T[(k-1)*lanm1_l+k] = beta;
      /*---------------------- Restarting test */
      int k1 = k-trlen-Ntest;
      if ( ((k1>=0) && (k1 % cycle == 0)) || (k == lanm) || it == maxit) {
        /*-------------------- solve eigen-problem for T(1:k,1:k)
                               vals in Rval, vecs in EvecT */
        /* printf("k %d, trlen %d, Ntest %d, its %d\n", k, trlen, Ntest, it); */
        SymEigenSolver(k, T, lanm1, EvecT, lanm1, Rval);
        /*-------------------- max dim reached-break from the inner loop */
        if (k == lanm || it == maxit) {
          break;
        }
#if 1
        /* we change the convergence test to be simple:
         * we test if the sum and the number of the Ritz values that are >= bar no longer change */
        jl = 0;
        tr = 0.0;
        for (i=0; i<k; i++) {
          if (Rval[i] + EVSL_DBL_EPS_MULT * DBL_EPSILON >= bar) {
            jl++;
            tr += Rval[i];
          }
        }
        if (fabs(tr-last_tr) <= 2e-12 && jl == last_jl) {
          if (fstats) {
            fprintf(fstats,"break: [it %d, k %d]: last_tr %.15e, tr %.15e, last_jl %d  jl %d\n",
                    it, k, last_tr, tr, last_jl, jl);
          }
          break;
        }
        if (fstats) {
          fprintf(fstats,"[it %d, k %d] testing: last_tr %.15e, tr %.15e, last_jl %d  jl %d\n",
                  it, k, last_tr, tr, last_jl, jl);
        }
        last_tr = tr;
        last_jl = jl;
#else
        count = 0;
        /*-------------------- get residual norms and check acceptance of
                               Ritz values for p(A). */
        /*-------------------- count e.vals in interval + those that convergeed */
        jl = 0;
        for (i=0; i<k; i++) {
          if (Rval[i]>= bar) {
            jl++;
            /* for standard e.v prob, this is the 2-norm of A*y-lam*y
             * for gen e.v prob, this is the B-norm of B^{-1}*A*y-lam*y */
            r = fabs(beta*EvecT[i*lanm1+(k-1)]);
            resi[i] = r;
            if (r < tolP) {
              count++;
            }
          }
        }
        /*-------------------- testing (partial) convergence for restart */
        if (do_print) {
          fprintf(fstats,"inner: testing conv k %4d, it %4d, count %3d  jl %3d trlen %3d\n",
                  k, it, count, jl, trlen);
        }
        /*-------------------- enough good candidates 1st time -> break */
        /* if ((count*evFrac >= nev-lock) && (prtrlen==-1)) */
        if (count*evFrac >= nev-lock) {
          fprintf(fstats, "count * evFrac >= nev-lock: %d * %d >= %d - %d\n", count, evFrac, nev, lock);
          break;
        }
        /*-------------------- count & jl unchanged since last test --> break */
        if ((count<=last_count) && (jl<=last_jl)) {
          fprintf(fstats, "inner: count <= last_count && jl <= last_jl: %d <= %d && %d <= %d\n", count, last_count, jl, last_jl);
          break;
        }
        last_count = count;
        last_jl = jl;
#endif
      } /* if [Restarting test] block */
      /*-------------------- end of inner (Lanczos) loop - Next: restart*/
    } /* while (k<mlan) loop */


    //SymEigenSolver(k, T, lanm1, EvecT, lanm1, Rval);
    //savedensemat(T, lanm1, k, k, "T.txt");
    //savedensemat(EvecT, lanm1, k, k, "Evec.txt");
    //save_vec(k, Rval, "eval.txt");
    /*--------------------   TWO passes to select good candidates */
    /*                       Pass-1: based on if ``p(Ritzvalue) > bar'' */
    jl = 0;
    for (i=0; i<k; i++) {
      //printf("resi[%d] = %.15e\n", i, fabs(beta*EvecT[i*lanm1+(k-1)]));
      /*--------------------   if this Ritz value is higher than ``bar'' */
      if (Rval[i] + EVSL_DBL_EPS_MULT * DBL_EPSILON >= bar) {
        /* move good eigenvectors/vals to front */
        if (i != jl) {
          evsl_dcopy(&k, EvecT+i*lanm1_l, &one, EvecT+jl*lanm1_l, &one);
          Rval[jl] = Rval[i];
          //resi[jl] = resi[i];
        }
        resi[jl] = fabs(beta*EvecT[i*lanm1_l+(k-1)]);
        //printf("beta = %.15e, resi[%d] = %.15e\n", beta, i, resi[jl]);
        jl++;
      }
    }
    //exit(0);
    //fprintf(fstats, "beta = %.1e\n", beta);
#ifdef EVSL_USING_CUDA_GPU
    evsl_memcpy_host_to_device(EvecT_device, EvecT, lanm1_l*lanm1_l*sizeof(double));
#endif
    /*---------------------- Compute the Ritz vectors:
     *                       Rvec(:,1:jl) = V(:,1:k) * EvecT(:,1:jl) */
    evsl_dgemm_device(&cN, &cN, &n, &jl, &k, &done, V, &n, EvecT_device, &lanm1, &dzero, Rvec, &n);
    if (ifGenEv) {
      evsl_dgemm_device(&cN, &cN, &n, &jl, &k, &done, Z, &n, EvecT_device, &lanm1, &dzero, BRvec, &n);
    }
    /*-------------------- Pass-2: check if Ritz vals of A are in [a,b] */
    /*                     number of Ritz values in [a,b] */
    ll = 0;
    /*-------------------- trlen = # Ritz vals that will go to TR set */
    prtrlen = trlen;
    trlen = 0;
    for (i=0; i<jl; i++) {
      double *y = Rvec + i*n_l;
      double *By = NULL;
      if (ifGenEv) {
        By = BRvec + i*n_l;
      }
      double *w = work;
      double *w2 = w + n;
      /*------------------ normalize just in case. */
      if (ifGenEv) {
        /* B-norm, w2 = B*y */
        matvec_B(y, w2);
        t = sqrt(evsl_ddot_device(&n, y, &one, w2, &one));
      } else {
        /* 2-norm */
        t = evsl_dnrm2_device(&n, y, &one);
      }
      /*-------------------- return code 2 --> zero eigenvector found */
      if (t == 0.0) {
        return 2;
      }
      /*-------------------- scal y */
      t = 1.0 / t;
      evsl_dscal_device(&n, &t, y, &one);
      /*-------------------- scal B*y */
      if (ifGenEv) {
        evsl_dscal_device(&n, &t, w2, &one);
      }
      /*-------------------- w = A*y */
      matvec_A(y, w);
      /*-------------------- Ritzval: t3 = (y'*w)/(y'*y) or
       *                              t3 = (y'*w)/(y'*B*y) */
      /*-------------------- Rayleigh quotient */
      double t3 = evsl_ddot_device(&n, y, &one, w, &one);
      /*--------------------  if lambda (==t3) is in [a,b] */
      if (t3 >= aa - EVSL_DBL_EPS_MULT * DBL_EPSILON && t3 <= bb + EVSL_DBL_EPS_MULT * DBL_EPSILON) {
        ll++;
        /*-------------------- compute residual wrt A for this pair */
        double nt3 = -t3;
        if (ifGenEv) {
          /* w = w - t3*w2, w2 = B*y,  (w=A*y-t3*B*y) */
          evsl_daxpy_device(&n, &nt3, w2, &one, w, &one);
        } else {
          /*-------------------- w = w - t3*y, (w=A*y-t3*y) */
          evsl_daxpy_device(&n, &nt3, y, &one, w, &one);
        }
        /*-------------------- res0 = 2-norm of w */
        res0 = evsl_dnrm2_device(&n, w, &one);
        /*-------------------- test res. of this Ritz pair against tol */
        /* r = resi[i];*/
        r = res0;
        if (r < tol) {
          //-------------------- check if need to realloc
          if (lock >= nev) {
            int old_nev = nev;
            nev += 1 + (int) (nev*nevInc);
            if (do_print) {
              fprintf(fstats, "-- More eigval found: realloc space for %d evs\n", nev);
            }
            Y = evsl_Realloc_device(Y, old_nev*n_l, double, nev*n_l, double);
            if (ifGenEv) {
              BY = evsl_Realloc_device(BY, old_nev*n_l, double, nev*n_l, double);
            }
            Lam = evsl_Realloc(Lam, nev, double);
            res = evsl_Realloc(res, nev, double);
          }
          /*--------------------   accept (t3, y) */
          evsl_dcopy_device(&n, y, &one, Y+lock*n_l, &one);
          if (ifGenEv) {
            evsl_dcopy_device(&n, By, &one, BY+lock*n_l, &one);
          }
          Lam[lock] = t3;
          res[lock] = res0;
          lock++;
        } else {
          /*-------------------- restart; move Ritz pair for TR to front */
          Rval[trlen] = Rval[i];
          evsl_dcopy_device(&n, y, &one, V+trlen*n_l, &one);
          if (ifGenEv) {
            evsl_dcopy_device(&n, By, &one, Z+trlen*n_l, &one);
          }
          /* special vector for TR that is the bottom row of
           * eigenvectors of Tm */
          s[trlen] = beta * EvecT[i*lanm1_l+(k-1)];
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
    if (do_print) {
      fprintf(fstats,"it %4d:  k %3d, jl %3d, ll %3d, lock %3d, trlen %3d\n",
              it, k, jl, ll, lock, trlen);
    }
    /*-------------------- TESTs for stopping */
    if ((prtrlen == 0) && (ll==0)) {
      /*--------------------  It looks like: Nothing left to compute  */
      if (do_print) {
        fprintf(fstats, "--------------------------------------------------\n");
        fprintf(fstats, " --> No ev.s left to be computed\n");
      }
      break;
    }
    if (it == maxit) {
      if (do_print) {
        fprintf(fstats, "--------------------------------------------------\n");
        fprintf(fstats, " --> Max its reached without convergence\n");
      }
      break;
    }
    /*-------------------- prepare to restart.  First zero out all T */
    memset(T, 0, lanm1_l*lanm1_l*sizeof(double));
    /*-------------------- move starting vector vector V(:,k+1);  V(:,trlen+1) = V(:,k+1) */
    evsl_dcopy_device(&n, V+k*n_l, &one, V+trlen*n_l, &one);
    if (ifGenEv) {
      evsl_dcopy_device(&n, Z+k*n_l, &one, Z+trlen*n_l, &one);
    }
  } /* outer loop (it) */

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
  evsl_Free_device(V);
  evsl_Free(T);
  evsl_Free(Rval);
  evsl_Free(resi);
  evsl_Free(EvecT);
#ifdef EVSL_USING_CUDA_GPU
  evsl_Free_device(EvecT_device);
#endif
  evsl_Free_device(Rvec);
  evsl_Free(s);
  evsl_Free_device(work);
  if (vrand) {
    evsl_Free_device(vrand);
  }
  if (ifGenEv) {
    evsl_Free_device(Z);
    evsl_Free_device(BY);
    evsl_Free_device(BRvec);
  }
  /*-------------------- record stats */
  tall = evsl_timer() - tall;
  /*-------------------- print stat */
  evslstat.t_iter = tall;

  return 0;
}

