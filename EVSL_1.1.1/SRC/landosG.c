#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "internal_header.h"

/**
 * @file SRC/landosG.c
 * @brief Function to use Lanczos method for approximating DOS for the
 * generalized eigenvalue problem.
 */

/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density) using Lanczos
 *    algorithm for the generalized eigenvalue problem.
 *
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] npts number of sample points used for the DOS curve
 *    @param[in] intv Stores the intervals of interest
 *      intv[0:1] = [a b] = interval where DOS is to be computed
 *      intv[2:3] = [lambda_min, lambda_max] \\
 *
 *    @param[out] xdos Length-npts long vector, x-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] ydos Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] neig  estimated number of eigenvalues
 *
 *
 *
 *    @note This works for both the standard and generalized eigenvalue
 *    problems.
 *    landos.c/LanDos is only for the standard eigenvalue problem.
 *----------------------------------------------------------------------*/

int LanDosG(const int nvec, const int msteps, int npts, double *xdos, double *ydos,
            double *neig, const double *const intv) {

  int i, j, k;
  int maxit = msteps, m;  /* Max number of iterations */
  /* size of the matrix */
  int n = evsldata.n;
  size_t n_l = n;
  const int ifGenEv = evsldata.ifGenEv;

  /* some RanGen requires n is multiple of 2, so alloc 1 more element for vinit */
  int n2 = ((n & 1) == 1) ? n+1 : n;

  double *vinit;
  vinit = evsl_Malloc_device(n2, double);

  int *ind;
  ind = evsl_Malloc(npts, int);
  double *y;
  y = evsl_Calloc(npts, double);

  /*-------------------- frequently used constants  */
  int one = 1;
  maxit = evsl_min(n, maxit);
  size_t maxit_l = maxit;
  double *gamma2;
  gamma2 = evsl_Calloc(maxit, double);
  /*-----------------------------------------------------------------------*
   * *Non-restarted* Lanczos iteration
   *-----------------------------------------------------------------------
   -------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *dT, *eT, *Z;
  V = evsl_Malloc_device(n_l * (maxit + 1), double);
  if (ifGenEv) {
    /* storage for Z = B * V */
    Z = evsl_Malloc_device(n_l * (maxit + 1), double);
  } else {
    /* Z and V are the same */
    Z = V;
  }
  /*-------------------- diag. subdiag of Tridiagional matrix */
  dT = evsl_Malloc(maxit, double);
  eT = evsl_Malloc(maxit, double);
  double *EvalT, *EvecT;
  EvalT = evsl_Malloc(maxit, double);          /* eigenvalues of tridia. matrix  T */
  EvecT = evsl_Malloc(maxit_l * maxit_l, double);  /* Eigen vectors of T */
  const double lm = intv[2];
  const double lM = intv[3];
  const double aa = evsl_max(intv[0], intv[2]);
  const double bb = evsl_min(intv[1], intv[3]);
  const double kappa = 1.25;
  const int M = evsl_min(msteps, 30);
  const double H = (lM - lm) / (M - 1);
  const double sigma = H / sqrt(8 * log(kappa));
  const double sigma2 = 2 * sigma * sigma;
  /*-------------------- If gaussian small than tol ignore point. */
  const double tol = 1e-08;
  double width = sigma * sqrt(-2.0 * log(tol));
  linspace(aa, bb, npts, xdos);  // xdos = linspace(lm,lM, npts);
  /*-------------------- u  is just a pointer. wk == work space */
  double *wk;
  const size_t wk_size = ifGenEv ? 6 * n_l : 4 * n_l;
  wk = evsl_Malloc_device(wk_size, double);

  for (m = 0; m < nvec; m++) {
    randn_double_device(n2, vinit);
    /*-------------------- copy initial vector to Z(:,1) */
    /* Filter the initial vector */
    if (ifGenEv) {
      solve_LT(vinit, V);
    }
    /*--------------------  normalize it */
    double t;
    if (ifGenEv) {
      /* B norm */
      matvec_B(V, Z);
      t = 1.0 / sqrt(evsl_ddot_device(&n, V, &one, Z, &one));
      evsl_dscal_device(&n, &t, Z, &one);
    } else {
      /* 2-norm */
      t = 1.0 / evsl_dnrm2_device(&n, vinit, &one);  // TODO: add a test here.
      evsl_dcopy_device(&n, vinit, &one, V, &one);
    }
    /* unit B^{-1}-norm or 2-norm */
    evsl_dscal_device(&n, &t, V, &one);
    /*-------------------- for ortho test */
    double wn = 0.0;
    int nwn = 0;
    /*-------------------- lanczos vectors updated by rotating pointer*/
    /*-------------------- pointers to Lanczos vectors */
    double *zold, *z, *znew;
    double *v, *vnew;
    /*--------------------  Lanczos recurrence coefficients */
    double alpha, nalpha, beta = 0.0, nbeta;
    /* ---------------- main Lanczos loop */
    for (k = 0; k < maxit; k++) {
      /*-------------------- quick reference to Z(:,k-1) when k>0*/
      zold = k > 0 ? Z + (k - 1) * n_l : NULL;
      /*-------------------- a quick reference to V(:,k) */
      v = &V[k * n_l];
      /*-------------------- a quick reference to Z(:,k) */
      z = &Z[k * n_l];
      /*-------------------- next Lanczos vector V(:,k+1)*/
      vnew = v + n;
      /*-------------------- next Lanczos vector Z(:,k+1)*/
      znew = z + n;
      matvec_A(v, znew);
      /*------------------ znew = znew - beta*zold */
      if (zold) {
        nbeta = -beta;
        evsl_daxpy_device(&n, &nbeta, zold, &one, znew, &one);
      }
      /*-------------------- alpha = znew'*v */
      alpha = evsl_ddot_device(&n, v, &one, znew, &one);
      /*-------------------- T(k,k) = alpha */
      dT[k] = alpha;
      wn += fabs(alpha);
      /*-------------------- znew = znew - alpha*z */
      nalpha = -alpha;
      evsl_daxpy_device(&n, &nalpha, z, &one, znew, &one);
      /*-------------------- FULL reortho to all previous Lan vectors */
      if (ifGenEv) {
        /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
        CGS_DGKS2(n, k + 1, NGS_MAX, Z, V, znew, wk);
        /* vnew = B \ znew */
        solve_B(znew, vnew);
        /*-------------------- beta = (vnew, znew)^{1/2} */
        beta = sqrt(evsl_ddot_device(&n, vnew, &one, znew, &one));
      } else {
        /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
        /* beta = norm(vnew) */
        CGS_DGKS(n, k + 1, NGS_MAX, V, vnew, &beta, wk);
      }
      wn += 2.0 * beta;
      nwn += 3;
      /*-------------------- lucky breakdown test */
      if (beta * nwn < orthTol * wn) {
        randn_double_device(n, vnew);
        if (ifGenEv) {
          /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
          CGS_DGKS2(n, k + 1, NGS_MAX, V, Z, vnew, wk);
          /* -------------- NOTE: B-matvec */
          matvec_B(vnew, znew);
          beta = sqrt(evsl_ddot_device(&n, vnew, &one, znew, &one));
          /*-------------------- vnew = vnew / beta */
          t = 1.0 / beta;
          evsl_dscal_device(&n, &t, vnew, &one);
          /*-------------------- znew = znew / beta */
          evsl_dscal_device(&n, &t, znew, &one);
          beta = 0.0;
        } else {
          /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
          /* beta = norm(vnew) */
          CGS_DGKS(n, k + 1, NGS_MAX, V, vnew, &beta, wk);
          /*-------------------- vnew = vnew / beta */
          t = 1.0 / beta;
          evsl_dscal_device(&n, &t, vnew, &one);
          beta = 0.0;
        }
      } else {
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        evsl_dscal_device(&n, &t, vnew, &one);
        if (ifGenEv) {
          /*-------------------- znew = znew / beta */
          evsl_dscal_device(&n, &t, znew, &one);
        }
      }
      /*-------------------- T(k+1,k) = beta */
      eT[k] = beta;
    }
    SymmTridEig(EvalT, EvecT, maxit, dT, eT);
    for (i = 0; i < maxit; i++) {
      /*-------------------- weights for Lanczos quadrature */
      /* Gamma2(i) = elementwise square of top entry of i-th eginvector
       * */
      gamma2[i] = EvecT[i * maxit_l] * EvecT[i * maxit_l];
    }
    /*-------------------- dos curve parameters
       Generate DOS from small gaussians centered at the ritz values */
    for (i = 0; i < msteps; i++) {
      // As msteps is width of ritzVal -> we get msteps eigenvectors
      const double t = EvalT[i];
      int numPlaced = 0;
      /*-------------------- Place elements close to t in ind */
      for (j = 0; j < npts; j++) {
        if (fabs(xdos[j] - t) < width) ind[numPlaced++] = j;
      }
      for (j = 0; j < numPlaced; j++) {
        y[ind[j]] += gamma2[i] *
                     exp(-((xdos[ind[j]] - t) * (xdos[ind[j]] - t)) / sigma2);
      }
    }
  } /* end of main loop */

  double scaling = 1.0 / (nvec * sqrt(sigma2 * EVSL_PI));
  /* y = ydos * scaling */
  evsl_dscal(&npts, &scaling, y, &one);
  evsl_dcopy(&npts, y, &one, ydos, &one);
  simpson(xdos, y, npts);
  *neig = y[npts - 1] * n;
  evsl_Free(gamma2);
  /*-------------------- free arrays */
  evsl_Free_device(vinit);
  evsl_Free_device(V);
  evsl_Free(dT);
  evsl_Free(eT);
  evsl_Free(EvalT);
  evsl_Free(EvecT);
  evsl_Free_device(wk);
  evsl_Free(y);
  evsl_Free(ind);
  if (ifGenEv) {
    evsl_Free_device(Z);
  }

  return 0;
}

