#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

#define COMP_RES 0
/**
 * @brief Lanczos process for eigenvalue bounds [Thick restart version]
 *
 * @param lanm      Dimension of Krylov subspace [restart dimension]
 * 
 * @param maxit  max Num of outer Lanczos iterations (restarts) allowed -- 
 *         Each restart may or use the full lanm lanczos steps or fewer.
 * 
 * @param tol       tolerance for convergence
 * @param vinit     initial  vector for Lanczos -- [optional]
 *
 * @b Modifies:
 * @param[out] lammin   Lower bound of the spectrum
 * @param[out] lammax   Upper bound of the spectrum
 * @param[out] fstats File stream which stats are printed to
 *
 * @return Returns 0 on success 
 *
 **/
int LanTrbounds(int lanm, int maxit, double tol, double *vinit,
                int bndtype,
                double *lammin, double *lammax, FILE *fstats) {
  double lmin=0.0, lmax=0.0, t, t1, t2;
  int do_print = 1;
  /* handle case where fstats is NULL. Then no output. Needed for openMP. */
  if (fstats == NULL) {
    do_print = 0;
  }
  /* size of the matrix */
  int n = evsldata.n;
  /*--------------------- adjust lanm and maxit */
  lanm = min(lanm, n);
  int lanm1=lanm+1;
  /*  if use full lanczos, should not do more than n iterations */
  if (lanm == n) {
    maxit = min(maxit, n);
  }
  /*--------------------   some constants frequently used */
  /* char cT='T'; */
  char cN = 'N';
  int one = 1;
  double done=1.0, dmone=-1.0, dzero=0.0;
  int i;
  /*-----------------------------------------------------------------------* 
   * *thick restarted* Lanczos step 
   *-----------------------------------------------------------------------*/
  if (do_print) {
    fprintf(fstats, " LanTR for bounds: dim %d, maxits %d\n", lanm, maxit);
  }
  /*--------------------- the min number of steps to be performed for 
   *                      each innter loop, must be >= 1 - if not, 
   *                      it means that the Krylov dim is too small */
  int min_inner_step = 5;
  /*-------------------- it = number of Lanczos steps */
  int it = 0;
  /*-------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *T;
  Malloc(V, n*lanm1, double);
  /*-------------------- for gen eig prob, storage for Z = B * V */
  double *Z;
  if (evsldata.ifGenEv) {
    Malloc(Z, n*lanm1, double);
  } else {
    Z = V;
  }
  /*-------------------- T must be zeroed out initially */
  Calloc(T, lanm1*lanm1, double);
  /*-------------------- trlen = dim. of thick restart set */
  int trlen = 0;
  /*-------------------- Ritz values and vectors of p(A) */
  double *Rval, *Rvec, *BRvec=NULL;
  Malloc(Rval, lanm, double);
  /*-------------------- Only compute 2 Ritz vectors */
  Malloc(Rvec, n*2, double);
  if (evsldata.ifGenEv) {
    Malloc(BRvec, n*2, double);
  }
  /*-------------------- Eigen vectors of T */
  double *EvecT;
  Malloc(EvecT, lanm1*lanm1, double);
  /*-------------------- s used by TR (the ``spike'' of 1st block in Tm)*/
  double s[3];
  /*-------------------- alloc some work space */
  double *work;
  int work_size = 2*n;
  Malloc(work, work_size, double);  
  /*-------------------- copy initial vector to V(:,1)   */
  DCOPY(&n, vinit, &one, V, &one);
  /*-------------------- normalize it */
  if (evsldata.ifGenEv) {
    /* B norm */
    matvec_B(V, Z);
    t = 1.0 / sqrt(DDOT(&n, V, &one, Z, &one));
    /* z = B*v */
    DSCAL(&n, &t, Z, &one);
  } else {
    /* 2-norm */
    t = 1.0 / DNRM2(&n, V, &one);
  }
  /* unit B-norm or 2-norm */
  DSCAL(&n, &t, V, &one);
  /*-------------------- main (restarted Lan) outer loop */
  while (it < maxit) {
    /*-------------------- for ortho test */
    double wn = 0.0;
    int nwn = 0;
    /*  beta */
    double beta = 0.0;
    /*  start with V(:,k) */
    int k = trlen > 0 ? trlen + 1 : 0;
    /* ! add a test if dimension exceeds (m+1) 
     * (trlen + 1) + min_inner_step <= lanm + 1 */
    if (k+min_inner_step > lanm1) {
      fprintf(stderr, "Krylov dim too small for this problem. Try a larger dim\n");
      exit(1);
    }
    /*-------------------- thick restart special step */
    if (trlen > 0) {
      int k1 = k-1;
      /*------------------ a quick reference to V(:,k) */
      double *v = &V[k1*n];
      double *z = &Z[k1*n];
      /*------------------ next Lanczos vector */
      double *vnew = v + n;
      double *znew = z + n;
      /*------------------ znew = A * v */
      matvec_A(v, znew);
      /*-------------------- restart with 'trlen' Ritz values/vectors
                             T = diag(Rval(1:trlen)) */
      for (i=0; i<trlen; i++) {
        T[i*lanm1+i] = Rval[i];
        wn += fabs(Rval[i]);
      }
      /*--------------------- s(k) = V(:,k)'* znew */
      s[k1] = DDOT(&n, v, &one, znew, &one);
      /*--------------------- znew = znew - Z(:,1:k)*s(1:k) */
      DGEMV(&cN, &n, &k, &dmone, Z, &n, s, &one, &done, znew, &one);
      /*-------------------- expand T matrix to k-by-k, arrow-head shape
                             T = [T, s(1:k-1)] then T = [T; s(1:k)'] */
      for (i=0; i<k1; i++) {
        T[trlen*lanm1+i] = s[i];
        T[i*lanm1+trlen] = s[i];
        wn += 2.0 * fabs(s[i]);
      }
      T[trlen*lanm1+trlen] = s[k1];
      wn += fabs(s[k1]);
      if (evsldata.ifGenEv) {
        /*-------------------- vnew = B \ znew */
        solve_B(znew, vnew);
        /*-------------------- beta = (vnew, znew)^{1/2} */
        beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
      } else {
        /*-------------------- beta = norm(w) */
        beta = DNRM2(&n, vnew, &one);
      }
      wn += 2.0 * beta;
      nwn += 3*k;
      /*   beta ~ 0 */
      if (beta*nwn < orthTol*wn) {
        rand_double(n, vnew);
        if (evsldata.ifGenEv) {
          /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
          CGS_DGKS2(n, k, NGS_MAX, V, Z, vnew, work);          
          matvec_B(vnew, znew);
          beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
          double ibeta = 1.0 / beta;
          DSCAL(&n, &ibeta, vnew, &one);        
          DSCAL(&n, &ibeta, znew, &one);
          beta = 0.0;            
        } else {
          /*   vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
          /*   beta = norm(w) */
          CGS_DGKS(n, k, NGS_MAX, V, vnew, &beta, work);
          double ibeta = 1.0 / beta;
          DSCAL(&n, &ibeta, vnew, &one);
          beta = 0.0;
        }
      } else {
        /*------------------- w = w / beta */
        double ibeta = 1.0 / beta;
        DSCAL(&n, &ibeta, vnew, &one);
        if (evsldata.ifGenEv) {
          DSCAL(&n, &ibeta, znew, &one);
        }
      }
      /*------------------- T(k+1,k) = beta; T(k,k+1) = beta; */
      T[k1*lanm1+k] = beta;
      T[k*lanm1+k1] = beta;
    } /* if (trlen > 0) */
    /*-------------------- Done with TR step. Rest of Lanczos step */
    /*-------------------- regardless of trlen, *(k+1)* is the current 
     *                     number of Lanczos vectors in V */
    /*-------------------- pointer to the previous Lanczos vector */
    double *zold = k > 0 ? Z+(k-1)*n : NULL;
    /*------------------------------------------------------*/
    /*------------------ Lanczos inner loop ----------------*/
    /*------------------------------------------------------*/
    while (k < lanm && it < maxit) {
      k++;
      /*---------------- a quick reference to V(:,k) */
      double *v = &V[(k-1)*n];
      double *z = &Z[(k-1)*n];
      /*---------------- next Lanczos vector */
      double *vnew = v + n;
      double *znew = z + n;
      /*------------------ znew = A * v */
      matvec_A(v, znew);
      it++;
      /*-------------------- znew = znew - beta*zold */
      if (zold) {
        double nbeta = -beta;
        DAXPY(&n, &nbeta, zold, &one, znew, &one);
      }
      /*-------------------- alpha = znew'*v */
      double alpha = DDOT(&n, v, &one, znew, &one);
      /*-------------------- T(k,k) = alpha */
      T[(k-1)*lanm1+(k-1)] = alpha;
      wn += fabs(alpha);
      /*-------------------- znew = znew - alpha*z */
      double nalpha = -alpha;
      DAXPY(&n, &nalpha, z, &one, znew, &one);
      /*-------------------- FULL reortho to all previous Lan vectors */
      if (evsldata.ifGenEv) {
        /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
        CGS_DGKS2(n, k, NGS_MAX, Z, V, znew, work);
        /* vnew = B \ znew */
        solve_B(znew, vnew);
        /*-------------------- beta = (vnew, znew)^{1/2} */
        beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
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
        /* generate a new init vector*/
        rand_double(n, vnew);
        if (evsldata.ifGenEv) {
          /* vnew = vnew - V(:,1:k)*Z(:,1:k)'*vnew */
          CGS_DGKS2(n, k, NGS_MAX, V, Z, vnew, work);          
          matvec_B(vnew, znew);
          beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
          double ibeta = 1.0 / beta;
          DSCAL(&n, &ibeta, vnew, &one);          
          DSCAL(&n, &ibeta, znew, &one);
          beta = 0.0;            
        } else {
          /*   vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
          /*   beta = norm(w) */
          CGS_DGKS(n, k, NGS_MAX, V, vnew, &beta, work);          
          double ibeta = 1.0 / beta;
          DSCAL(&n, &ibeta, vnew, &one);
          beta = 0.0;
        }
      } else {
        /*---------------------- vnew = vnew / beta */
        double ibeta = 1.0 / beta;
        DSCAL(&n, &ibeta, vnew, &one);
        if (evsldata.ifGenEv) {
          /*-------------------- znew = znew / beta */
          DSCAL(&n, &ibeta, znew, &one);
        }
      }
      /*-------------------- T(k,k+1) = T(k+1,k) = beta */
      T[k*lanm1+(k-1)] = beta;
      T[(k-1)*lanm1+k] = beta;
    } /* while (k<mlan) loop */

    /*-------------------- solve eigen-problem for T(1:k,1:k)
                           vals in Rval, vecs in EvecT */
    SymEigenSolver(k, T, lanm1, EvecT, lanm1, Rval);

    /*-------------------- Rval is in ascending order */
    /*-------------------- Rval[0] is smallest, Rval[k-1] is largest */
    /*-------------------- special vector for TR that is the bottom row of 
                           eigenvectors of Tm */
    s[0] = beta * EvecT[k-1];
    s[1] = beta * EvecT[(k-1)*lanm1+(k-1)];
    /*---------------------- bounds */
    if (bndtype <= 1) {
      /*-------------------- BOUNDS type 1 (simple) */
      t1 = fabs(s[0]);
      t2 = fabs(s[1]);
    } else {
      /*-------------------- BOUNDS type 2 (Kato-Temple) */
      t1 = 2.0*s[0]*s[0] / (Rval[1] - Rval[0]);
      t2 = 2.0*s[1]*s[1] / (Rval[k-1] - Rval[k-2]);
    }
    lmin = Rval[0]   - t1;
    lmax = Rval[k-1] + t2;
    /*---------------------- Compute two Ritz vectors: 
     *                       Rvec(:,1) = V(:,1:k) * EvecT(:,1) 
     *                       Rvec(:,end) = V(:,1:k) * EvecT(:,end) */
    DGEMV(&cN, &n, &k, &done, V, &n, EvecT, &one, &dzero, Rvec, &one);
    DGEMV(&cN, &n, &k, &done, V, &n, EvecT+(k-1)*lanm1, &one, &dzero, Rvec+n, &one);
    if (evsldata.ifGenEv) {
      DGEMV(&cN, &n, &k, &done, Z, &n, EvecT, &one, &dzero, BRvec, &one);
      DGEMV(&cN, &n, &k, &done, Z, &n, EvecT+(k-1)*lanm1, &one, &dzero, BRvec+n, &one);
    }
    /*---------------------- Copy two Rval and Rvec to TR set */
    trlen = 2;
    for (i=0; i<2; i++) {
      double *y = Rvec + i*n;
      DCOPY(&n, y, &one, V+i*n, &one);
      if (evsldata.ifGenEv) {
        double *By = BRvec + i*n;
        DCOPY(&n, By, &one, Z+i*n, &one);
      }
    }
    Rval[1] = Rval[k-1];
    /*-------------------- recompute residual norm for debug only */
#if COMP_RES
    /*-------------------- compute residual */
    double *w1 = work;
    double *w2 = work + n;
    double rr[2];
    for (i=0; i<2; i++) {
      double *y = Rvec + i*n;
      double nt = -Rval[i];
      matvec_A(y, w1);
      if (evsldata.ifGenEv) {
        matvec_B(y, w2);
        DAXPY(&n, &nt, w2, &one, w1, &one);
        solve_B(w1, w2);
        rr[i] = sqrt(DDOT(&n, w1, &one, w2, &one));
      } else {
        DAXPY(&n, &nt, y, &one, w1, &one);
        rr[i] = DNRM2(&n, w1, &one);
      }
    } 
    if (do_print) {
      fprintf(fstats,"it %4d, k %3d: ritz %.15e %.15e, t1,t2 %e %e, res %.15e %.15e, comp res %.15e %.15e\n", 
              it, k, Rval[0], Rval[1], t1, t2, fabs(s[0]), fabs(s[1]), rr[0], rr[1]);
    }
#else
    if (do_print) {
      fprintf(fstats,"it %4d, k %3d: ritz %.15e %.15e, t1,t2 %e %e, res %.15e %.15e\n", 
              it, k, Rval[0], Rval[1], t1, t2, fabs(s[0]), fabs(s[1]));
    }
#endif
    /*---------------------- test convergence */
    if (t1+t2 < tol*(fabs(lmin)+fabs(lmax))) {
      break;
    }
    /*-------------------- prepare to restart.  First zero out all T */
    memset(T, 0, lanm1*lanm1*sizeof(double));
    /*-------------------- move starting vector vector V(:,k+1);  V(:,trlen+1) = V(:,k+1) */
    DCOPY(&n, V+k*n, &one, V+trlen*n, &one);
    if (evsldata.ifGenEv) {
      DCOPY(&n, Z+k*n, &one, Z+trlen*n, &one);
    }
  } /* outer loop (it) */

  /*-------------------- Done.  output : */
  *lammin = lmin;
  *lammax = lmax;
  /*-------------------- free arrays */
  free(V);
  free(T);
  free(Rval);
  free(EvecT);
  free(Rvec);
  free(work);
  if (evsldata.ifGenEv) {
    free(Z);
    free(BRvec);
  }

  return 0;
}

