#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blaslapack.h"
#include "def.h"
#include "evsl.h"
#include "internal_proto.h"
#include "string.h"  //for memset
#include "struct.h"

/**
 * @file landosG.c
 * @brief Function to use Lanczos method for approximating DOS for the
 * general eigenvalue problem.
 */

/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density) using Lanczos
 *    algorithm for the general eigenvalue problem.
 *
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] degB Degree with which B should be approximated by \\
 *                Must be positive. Planned feature: if non-positive use
 *                Cholsky factorization.
 *    @param[in] npts number of sample points used for the DOS curve
 *    @param[in] *intv Stores the the intervals of interest
 *      intv[0:1] = [a b] = interval where DOS is to be computed
 *      intv[2:3] = [lambda_min, lambda_max] \\
 *    @param[in] tau Tolerence used
 *
 *    @param[out] *xdos Length-npts long vector, x-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] *ydos Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] neig  estimated number of eigenvalues
 *
 *    @note This works for both the standard and generalized eigenvalue
 *    problems.
 *    landos.c/LanDos is only for the standard eigenvalue problem.
 *----------------------------------------------------------------------*/

int LanDosG(const int nvec, const int msteps, int npts, double *xdos, double *ydos, 
            double *neig, const double *const intv) {

  int i, j, k;
  int maxit = msteps, m;  /* Max number of iterations */
  int n = evsldata.n; /* Number of elements in matrix */
  const int ifGenEv = evsldata.ifGenEv;
 
  double *vinit;
  Malloc(vinit, n, double);

  int *ind;
  Malloc(ind, npts, int);
  double *y;
  Calloc(y, npts, double);

  /*-------------------- frequently used constants  */
  int one = 1;
  /* size of the matrix */
  maxit = min(n, maxit);
  double *gamma2;
  Calloc(gamma2, maxit, double);
  /*-----------------------------------------------------------------------*
   * *Non-restarted* Lanczos iteration
   *-----------------------------------------------------------------------
   -------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
  double *V, *dT, *eT, *Z;
  Calloc(V, n * (maxit + 1), double);
  if (ifGenEv) {
    /* storage for Z = B * V */
    Calloc(Z, n * (maxit + 1), double);
  } else {
    /* Z and V are the same */
    Z = V;
  }
  /*-------------------- diag. subdiag of Tridiagional matrix */
  Malloc(dT, maxit, double);
  Malloc(eT, maxit, double);
  double *EvalT, *EvecT;
  Malloc(EvalT, maxit, double);          /* eigenvalues of tridia. matrix  T */
  Malloc(EvecT, maxit * maxit, double);  /* Eigen vectors of T */
  const double lm = intv[2];
  const double lM = intv[3];
  const double aa = max(intv[0], intv[2]);
  const double bb = min(intv[1], intv[3]);
  const double kappa = 1.25;
  const int M = min(msteps, 30);
  const double H = (lM - lm) / (M - 1);
  const double sigma = H / sqrt(8 * log(kappa));
  const double sigma2 = 2 * sigma * sigma;
  /*-------------------- If gaussian small than tol ignore point. */
  const double tol = 1e-08;
  double width = sigma * sqrt(-2.0 * log(tol));
  linspace(aa, bb, npts, xdos);  // xdos = linspace(lm,lM, npts);
  /*-------------------- u  is just a pointer. wk == work space */
  double *wk;
  const int wk_size = ifGenEv ? 6 * n : 4 * n;
  Malloc(wk, wk_size, double);
  for (m = 0; m < nvec; m++) {
    randn_double(n, vinit);
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
      t = 1.0 / sqrt(DDOT(&n, V, &one, Z, &one));
      DSCAL(&n, &t, Z, &one);
    } else {
      /* 2-norm */
      t = 1.0 / DNRM2(&n, vinit, &one);  // add a test here.
      DCOPY(&n, vinit, &one, V, &one);
    }
    /* unit B^{-1}-norm or 2-norm */
    DSCAL(&n, &t, V, &one);
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
      zold = k > 0 ? Z + (k - 1) * n : NULL;
      /*-------------------- a quick reference to V(:,k) */
      v = &V[k * n];
      /*-------------------- a quick reference to Z(:,k) */
      z = &Z[k * n];
      /*-------------------- next Lanczos vector V(:,k+1)*/
      vnew = v + n;
      /*-------------------- next Lanczos vector Z(:,k+1)*/
      znew = z + n;
      matvec_A(v, znew);
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
      if (ifGenEv) {
        /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
        CGS_DGKS2(n, k, NGS_MAX, Z, V, znew, wk);
        /* vnew = B \ znew */
        solve_B(znew, vnew);
        /*-------------------- beta = (vnew, znew)^{1/2} */
        beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
      } else {
        /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
        /* beta = norm(vnew) */
        CGS_DGKS(n, k + 1, NGS_MAX, V, vnew, &beta, wk);
      }
      wn += 2.0 * beta;
      nwn += 3;
      /*-------------------- lucky breakdown test */
      if (beta * nwn < orthTol * wn) {
        rand_double(n, vnew);
        if (ifGenEv) {
          /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
          CGS_DGKS2(n, k + 1, NGS_MAX, V, Z, vnew, wk);
          /* -------------- NOTE: B-matvec */
          matvec_B(vnew, znew);
          beta = sqrt(DDOT(&n, vnew, &one, znew, &one));
          /*-------------------- vnew = vnew / beta */
          t = 1.0 / beta;
          DSCAL(&n, &t, vnew, &one);
          /*-------------------- znew = znew / beta */
          DSCAL(&n, &t, znew, &one);
          beta = 0.0;
        } else {
          /* vnew = vnew - V(:,1:k)*V(:,1:k)'*vnew */
          /* beta = norm(vnew) */
          CGS_DGKS(n, k + 1, NGS_MAX, V, vnew, &beta, wk);
          /*-------------------- vnew = vnew / beta */
          t = 1.0 / beta;
          DSCAL(&n, &t, vnew, &one);
          beta = 0.0;
        }
      } else {
        /*-------------------- vnew = vnew / beta */
        t = 1.0 / beta;
        DSCAL(&n, &t, vnew, &one);
        if (ifGenEv) {
          /*-------------------- znew = znew / beta */
          DSCAL(&n, &t, znew, &one);
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
      gamma2[i] = EvecT[i * maxit] * EvecT[i * maxit];
    }
    /*-------------------- dos curve parameters
       Generate DOS from small gaussians centered at the ritz values */
    for (i = 0; i < msteps; i++) {
      // As msteps is width of ritzVal -> we get msteps eigenvectors
      const double t = EvalT[i];
      int numPlaced = 0;
      /*-------------------- Place elements close to t in ind */
      for (j = 0; j < npts; j++) {
        if (abs(xdos[j] - t) < width) ind[numPlaced++] = j;
      }
      for (j = 0; j < numPlaced; j++) {
        y[ind[j]] += gamma2[i] *
                     exp(-((xdos[ind[j]] - t) * (xdos[ind[j]] - t)) / sigma2);
      }
    }
  }
  double scaling = 1.0 / (nvec * sqrt(sigma2 * PI));
  /* y = ydos * scaling */
  DSCAL(&npts, &scaling, y, &one);
  DCOPY(&npts, y, &one, ydos, &one);
  simpson(xdos, y, npts);
  *neig = y[npts - 1] * n;
  free(gamma2);
  /*-------------------- free arrays */
  free(vinit);
  free(V);
  free(dT);
  free(eT);
  free(EvalT);
  free(EvecT);
  free(wk);
  free(y);
  free(ind);
  if (ifGenEv) {
    free(Z);
  }

  return 0;
}
