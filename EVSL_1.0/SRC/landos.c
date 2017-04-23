#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blaslapack.h"
#include "def.h"
#include "internal_proto.h"
#include "string.h"  //for memset
#include "struct.h"
//#include "vector.h"

/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density)
 *
 *    @param[in] *A    matrix A
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] npts number of sample points used ofr the curve
 *
 *    @param[out] xdos Length-npts long vector, x-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *    @param[out] ydos Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *
 *----------------------------------------------------------------------*/

int LanDos(csrMat *A, int nvec, int msteps, int npts, double *xdos,
           double *ydos, int *intv) {
  // Allocations from lanbounds.c
  double *alp, *bet, nbet, nalp, t, *V;
  int one = 1;
  int n;

  n = A->nrows;
  const double lm = intv[0];
  const double lM = intv[1];
  const double tolBdwn = 1.e-13 * (abs(lM) + abs(lm));
  const double aa = max(intv[0], intv[2]);
  const double bb = max(intv[1], intv[3]);
  const double kappa = 1.25;
  const int M = min(msteps, 30);
  const double H = (lM - lm) / (M - 1);
  const double sigma = H / sqrt(8 * log(kappa));
  sigma2 = 2 * sigma * sigma;
  // If gaussian small than tol ignore point.
  const double tol = 1e-08;
  width = sigma * sqrt(-2 * log(tol));
  linspace(aa, bb, npts, xdos);       // xdos = linspace(lm,lM, npts);
  memset(y, 0, npts * sizeof(y[0]));  // y = zeros(size(xdos));

  Malloc(alp, msteps, double);
  Malloc(bet, msteps, double);
  Malloc(V, (msteps + 1) * n, double);

  // Variables that persist through iterations
  double *v;  // Vector for current iteration
  Malloc(v, n, double);
  double sigma2;  // A dos curve parameter
  double width;   // A dos curve parameter
  double *y;      // Stores y values
  Calloc(y, npts, double);

  for (int m = 0; m < nvec; m++) {
    randn_double(n, v);  // w = randn(size(A,1),1);
    //---------------------------------------
    // Start of bulk of lanbound.c code
    //---------------------------------------

    t = DDOT(&n, v, &one, v, &one);
    t = 1.0 / sqrt(t);  // t is now the number required to normalize vector.
    DSCAL(&n, &t, v, &one);
    DCOPY(&n, v, &one, V,
          &one);  // v = w/norm(w); Might be able to use DNRM2 instead.
    double wn = 0.0;
    /*-------------------- main Lanczos loop */
    int j;
    for (j = 0; j < msteps; j++) {
      // w = A*v
      matvec_genev(A, &V[j * n], &V[(j + 1) * n]);
      // w = w - bet * vold
      if (j) {
        nbet = -bet[j - 1];
        DAXPY(&n, &nbet, &V[(j - 1) * n], &one, &V[(j + 1) * n], &one);
      }
      /* alp = w' * v */
      alp[j] = DDOT(&n, &V[(j + 1) * n], &one, &V[j * n], &one);
      wn += alp[j] * alp[j];
      // w = w - alp * v
      nalp = -alp[j];
      DAXPY(&n, &nalp, &V[j * n], &one, &V[(j + 1) * n], &one);
      // full reortho
      int i;
      for (i = 0; i <= j; i++) {
        t = DDOT(&n, &V[(j + 1) * n], &one, &V[i * n], &one);
        double mt = -t;
        DAXPY(&n, &mt, &V[i * n], &one, &V[(j + 1) * n], &one);
      }
      bet[j] = DDOT(&n, &V[(j + 1) * n], &one, &V[(j + 1) * n], &one);
      if (bet[j] * (j + 1) < orthTol * wn) {
        fprintf(stdout, "lanbounds: lucky break, j=%d, beta=%e, break\n", j,
                bet[j]);
        msteps = j + 1;
        break;
      }
      wn += 2.0 * bet[j];
      bet[j] = sqrt(bet[j]);
      t = 1.0 / bet[j];
      DSCAL(&n, &t, &V[(j + 1) * n], &one);
    }

    double *S, *ritzVal;
    Malloc(S, msteps * msteps, double);
    // Note that S is a matrix compressed into a single array.
    Malloc(ritzVal, msteps, double);
    //-------------------- diagonalize tridiagonal matrix
    SymmTridEig(ritzVal, S, msteps, alp, bet);
    // S = -eigvec
    // ritzVal = diags of D

    //---------------------------------------
    // End of bulk of lanbound.c code
    //---------------------------------------

    // theta = ritzVal = sorted eigenvalues IN ASCENDING ORDER
    double *gamma2;
    Malloc(gamma2, msteps, double);
    for (int i = 0; i < msteps; i++) {
      gamma2[i] =
          S[i * msteps] *
          S[i * msteps];  // Note the difference due to row/column major order
    }
    printf("\n");

    // Gamma^2 is now elementwise square of smallest eginvector

    // dos curve parameters

    // Generate DOS from small gaussians centered at the ritz values
    for (int i = 0; i < msteps;
         i++) {  // As msteps is width of ritzVal -> we get msteps eigenvectors
      const double t = ritzVal[i];
      int *ind;
      int numind = 0;

      for (int j = 0; j < npts;
           j++) {  // Calculate number of elements matching pattern
        if (abs(xdos[j] - t) < width) {
          numind++;
        }
      }
      Calloc(ind, numind, int);
      int numPlaced = 0;
      for (int j = 0; j < npts;
           j++) {  // Place the elements matching said pattern
        if (abs(xdos[j] - t) < width) {
          ind[numPlaced++] = j;
        }
      }
      // ind now is = find(abs(xdos - t) < width);

      // This replaces y(ind) = y(ind) +
      // gamma2(i)*exp(-(xdos(ind)-t).^2/sigma2);
      for (int j = 0; j < numind; j++) {
        y[ind[j]] = y[ind[j]] +
                    gamma2[i] * exp(-((xdos[ind[j]] - t) * (xdos[ind[j]] - t)) /
                                    sigma2);
      }
      free(ind);
    }
    free(gamma2);
    free(S);
    free(ritzVal);
  }

  double sum = 0;
  memcpy(ydos, y, sizeof(double) * npts);
  for (int i = 0; i < npts; i++) {
    sum += ydos[i];
  }
  for (int i = 0; i < npts; i++) {
    ydos[i] /= (sum * (xdos[1] - xdos[0]));
  }

  free(alp);
  free(bet);
  free(V);
  // free(S);
  // free(ritzVal);

  free(v);
  free(y);

  return 0;
}
