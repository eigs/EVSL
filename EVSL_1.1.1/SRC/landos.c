#include <stdio.h>
#include <stdlib.h>
#include "internal_header.h"

/**
 * @file SRC/landos.c
 * @brief Function to use Lanczos method for approximating DOS for the
 * standard eigenvalue problem.
 */

/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density)
 *
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] npts number of sample points used for the DOS curve
 *    @param[in] *intv Stores the two intervals of interest \\
 *      intv[0:1] = [a b] = interval where DOS is to be computed
 *      intv[2:3] = [lambda_min, lambda_max]\\
 *
 *    @param[out] *xdos Length-npts long vector, x-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] *ydos Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated before calling LanDos
 *
 *    @param[out] neig  estimated number of eigenvalues
 *
 *
 *
 *    @note This only works with the standard eigenvalue problem. Use landosG.c
 *    /LanDosG for the generalized eigenvalue problem.
 *
 *----------------------------------------------------------------------*/

int LanDos(const int nvec, int msteps, int npts, double *xdos, double *ydos,
           double *neig, const double *const intv) {
  //--------------------
  const int ifGenEv = evsldata.ifGenEv;
  if (ifGenEv) {
    fprintf(stderr,
            "ERROR: landos.c/LanDos only works with the standard eigenvalue "
            "problem. Please use landosG.c/LanDos for the generalized problem. "
            "\n");
    exit(-1);
  }
  double *alp, *bet, nbet, nalp, t, *V;
  int one = 1;
  int n, m, i, j;
  size_t n_l;

  n = evsldata.n;
  n_l = n;
  //-------------------- Variables that persist through iterations
  double *v, *y;  // v=Vector for current iteration; y Stores y values
  int *ind;
  v = evsl_Malloc_device(n, double);
  y = evsl_Calloc(npts, double);
  ind = evsl_Malloc(npts, int);
  /*-------------------- for tridiag. eigenvalue problem + lanczos
                         quadrature */
  double *S, *ritzVal, *gamma2;
  // S will contain a matrix compressed into a single array.
  size_t msteps_l = msteps;
  S = evsl_Malloc(msteps_l * msteps_l, double);
  gamma2 = evsl_Malloc(msteps, double);
  ritzVal = evsl_Malloc(msteps, double);
  const double lm = intv[2];
  const double lM = intv[3];
  const double tolBdwn = 1.e-13 * (fabs(lM) + fabs(lm));
  const double aa = evsl_max(intv[0], intv[2]);
  const double bb = evsl_min(intv[1], intv[3]);
  const double kappa = 1.25;
  const int M = evsl_min(msteps, 30);
  const double H = (lM - lm) / (M - 1);
  const double sigma = H / sqrt(8 * log(kappa));
  double sigma2 = 2 * sigma * sigma;
  //-------------------- If gaussian small than tol ignore point.
  const double tol = 1e-08;
  double width = sigma * sqrt(-2.0 * log(tol));
  linspace(aa, bb, npts, xdos);  // xdos = linspace(lm,lM, npts);

  alp = evsl_Malloc(msteps, double);
  bet = evsl_Malloc(msteps, double);
  V = evsl_Malloc_device((msteps + 1) * n_l, double);
  //-------------------- Lanczos loop for this vector
  for (m = 0; m < nvec; m++) {
    randn_double_device(n, v);  // w = randn(size(A,1),1);
    //--------------------Start of bulk of lanbound.c code
    t = evsl_ddot_device(&n, v, &one, v, &one);
    //-------------------- normalize vector
    //                     v = can also  use DNRM2 instead.
    t = 1.0 / sqrt(t);
    evsl_dscal_device(&n, &t, v, &one);
    evsl_dcopy_device(&n, v, &one, V, &one);
    double wn = 0.0;
    /*-------------------- main Lanczos loop */
    for (j = 0; j < msteps; j++) {
      // w = A*v
      matvec_A(&V[j * n_l], &V[(j + 1) * n_l]);
      // w = w - bet * vold
      if (j) {
        nbet = -bet[j - 1];
        evsl_daxpy_device(&n, &nbet, &V[(j - 1) * n_l], &one, &V[(j + 1) * n_l], &one);
      }
      /*-------------------- alp = w' * v */
      alp[j] = evsl_ddot_device(&n, &V[(j + 1) * n_l], &one, &V[j * n_l], &one);
      wn += alp[j] * alp[j];
      //-------------------- w = w - alp * v
      nalp = -alp[j];
      evsl_daxpy_device(&n, &nalp, &V[j * n_l], &one, &V[(j + 1) * n_l], &one);
      //-------------------- full reortho
      for (i = 0; i <= j; i++) {
        t = evsl_ddot_device(&n, &V[(j + 1) * n_l], &one, &V[i * n_l], &one);
        double mt = -t;
        evsl_daxpy_device(&n, &mt, &V[i * n_l], &one, &V[(j + 1) * n_l], &one);
      }
      bet[j] = evsl_ddot_device(&n, &V[(j + 1) * n_l], &one, &V[(j + 1) * n_l], &one);
      if (bet[j] * (j + 1) < orthTol * wn) {
        fprintf(stdout, "lanbounds: lucky break, j=%d, beta=%e, break\n", j, bet[j]);
        msteps = j + 1;
        break;
      }
      if (bet[j] > tolBdwn) {  // If it's not zero, continue as normal
        wn += 2.0 * bet[j];
        bet[j] = sqrt(bet[j]);
        t = 1.0 / bet[j];
        evsl_dscal_device(&n, &t, &V[(j + 1) * n_l], &one);
      } else {  // Otherwise generate a new vector and redo the previous
                // calculations on it
        randn_double_device(n, v);  // w = randn(size(A,1),1);
        for (i = 0; i <= j; i++) {
          t = evsl_ddot_device(&n, &V[(j + 1) * n_l], &one, &V[i * n_l], &one);
          double mt = -t;
          evsl_daxpy_device(&n, &mt, &V[i * n_l], &one, &V[(j + 1) * n_l], &one);
        }
        bet[j] = evsl_ddot_device(&n, &V[(j + 1) * n_l], &one, &V[(j + 1) * n_l], &one);
        wn += 2.0 * bet[j];
        bet[j] = sqrt(bet[j]);
        t = 1.0 / bet[j];
        evsl_dscal_device(&n, &t, &V[(j + 1) * n_l], &one);
        bet[j] = 0;
      }
    }
    /*-------------------- end Lanczos loop for this vector
      -------------------- diagonalize tridiagonal matrix */
    SymmTridEig(ritzVal, S, msteps, alp, bet);
    /* S = -eigvec
     ritzVal = diags of D */

    /* theta = ritzVal = sorted eigenvalues IN ASCENDING ORDER */
    for (i = 0; i < msteps; i++) {
      /*-------------------- weights for Lanczos quadrature
       Gamma2(i) = elementwise square of top entry of i-th eginvector*/
      gamma2[i] = S[i * msteps_l] * S[i * msteps_l];
    }
    /*-------------------- dos curve parameters
     Generate DOS from small gaussians centered at the ritz values */
    for (i = 0; i < msteps; i++) {
      /* As msteps is width of ritzVal -> we get msteps eigenvectors */
      const double t = ritzVal[i];
      int numPlaced = 0;
      /*-------------------- Place elements close to t in ind */
      for (j = 0; j < npts; j++) {
        if (fabs(xdos[j] - t) < width) ind[numPlaced++] = j;
      }

      for (j = 0; j < numPlaced; j++)
        y[ind[j]] += gamma2[i] *
                     exp(-((xdos[ind[j]] - t) * (xdos[ind[j]] - t)) / sigma2);
    }
    /*-------------------- end vector loop */
  }

  double scaling = 1.0 / (nvec * sqrt(sigma2 * EVSL_PI));
  /* y = ydos * scaling */
  evsl_dscal(&npts, &scaling, y, &one);
  evsl_dcopy(&npts, y, &one, ydos, &one);
  simpson(xdos, y, npts);

  *neig = y[npts - 1] * n;
  evsl_Free(gamma2);
  evsl_Free(S);
  evsl_Free(ritzVal);

  evsl_Free(alp);
  evsl_Free(bet);
  evsl_Free_device(V);

  evsl_Free_device(v);
  evsl_Free(y);
  evsl_Free(ind);

  return 0;
}
