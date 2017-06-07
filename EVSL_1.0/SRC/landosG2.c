#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blaslapack.h"
#include "def.h"
#include "evsl.h"
#include "internal_proto.h"
#include "string.h"  //for memset
#include "struct.h"
/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density)
 *
 *    @param[in] *A  -- used through calls to matvec_A
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] npts number of sample points used for the DOS curve
 *    @param[in] *intv Stores the two intervals of interest \\
 *      intv[0:1] = [lambda_min, lambda_max]\\
 *      intv[2:3] = [a b] = interval where DOS is to be computed
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

double fsqrt(const double a) { return sqrt(1 / a); }

double rec(const double a) { return 1 / a; }

void ones(int n, double* v) { int i = 0; for(i = 0; i < n; i++) { v[i] = 1; }}


/**
 * @brief @b Computes y=P(A) y, where pn is a Cheb. polynomial expansion 
 * 
 * This explicitly calls matvec, so it can be useful for implementing
 * user-specific matrix-vector multiplication.
 *
 * @param pol Struct containing the paramenters and expansion coefficient of
 * the polynomail.
 * @param v input vector
 *
 * @param[out] y p(A)v
 *
 * @b Workspace
 * @param w Work vector of length 3*n [allocate before call]
 * @param v is untouched
 **/
int pnav(polparams *pol, double *v, double *y, double *w, evslData data) {//Really just ChebAv
  printf("eln: %i", evsldata.n);
  int n = evsldata.n;
  printf("N: %i \n", n);
  /*-------------------- unpack pol */
  double *mu = pol->mu;
  double dd = pol->dd;
  double cc = pol->cc;
  int m = pol->deg;
  printf("mu[0]: %f", mu[0]);
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
  /*-------------------- degree loop. k IS the degree */
  for (k=1; k<=m; k++) {
    /*-------------------- y = mu[k]*Vk + y */    
    t = k == 1 ? t1 : t2; 
    /*-------------------- */    
    s = mu[k];
    if (evsldata.ifGenEv) {
      /*-------------------- Vkp1 = A*B\Vk - cc*Vk */    
      solve_B(vk, w2);
      matvec_A(w2, vkp1);
    } else {
      /*-------------------- Vkp1 = A*Vk - cc*Vk */    
      matvec_B(vk, vkp1);
    printf("vk0: %f \n", vk[0]);
    printf("vkp0: %f \n", vkp1[0]);
    }
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

  return 0;
}





int LanDosG2(const int nvec, int msteps, const int degB, int npts, double *xdos, double *ydos,
           double *neig, const double *const intv, const double tau) {
  //--------------------
  double *alp, *bet, nbet, nalp, t, *V, *Z;
  int one = 1;
  int n, m, i, j;
  const int mdeg = 200;

  n = evsldata.n;

  polparams pol_sqr;
  Malloc(pol_sqr.mu, mdeg, double);

  polparams pol_sol;
  Malloc(pol_sol.mu, mdeg, double);



  if (degB <= 0) {
	  printf("LanDos with degB <= not yet implemented!");
	  exit(-1);
  }
  else {
	  lsPol2(&intv[4], mdeg, fsqrt, tau, &pol_sqr);
	  lsPol2(&intv[4], mdeg, rec, tau, &pol_sol);
  }


  double lm = INFINITY;
  double lM = -INFINITY;
  //-------------------- Variables that persist through iterations
  //
  double *v, *y, *z, *vtmp;  // v=Vector for current iteration; y Stores y values
  int *ind;
  Malloc(v, n, double);
  Malloc(z, n, double);
  Calloc(y, npts, double);
  Malloc(ind, npts, int);
  /*-------------------- for tridiag. eigenvalue problem + lanczos
                         quadrature */
  double *S, *ritzVal, *gamma2;

  // S will contain a matrix compressed into a single array.
  Malloc(S, msteps * msteps, double);
  Malloc(gamma2, msteps, double);
  Malloc(ritzVal, msteps, double);
  // const double lm = intv[0];
  // const double lM = intv[1];
  const double tolBdwn = 1.e-13 * (abs(lM) + abs(lm));
  const double aa = max(intv[0], intv[2]);
  const double bb = min(intv[1], intv[3]);
  const double kappa = 1.25;
  const int M = min(msteps, 30);
  const double H = (lM - lm) / (M - 1);
  const double sigma = H / sqrt(8 * log(kappa));
  double sigma2 = 2 * sigma * sigma;
  //-------------------- If gaussian small than tol ignore point.
  const double tol = 1e-08;
  double width = sigma * sqrt(-2.0 * log(tol));
  linspace(aa, bb, npts, xdos);  // xdos = linspace(lm,lM, npts);

  double *ww; //Workspace for ChebAv
  Malloc(ww, 3 * n, double); //Should be length 3*(length pol)


  Malloc(alp, msteps, double);
  Malloc(bet, msteps, double);
  Malloc(V, (msteps + 1) * n, double);
  Malloc(Z, (msteps + 1) * n, double);
  Calloc(vtmp, n, double);
  //-------------------- Lanczos loop for this vector
  //Memory errors
  for (m = 0; m < nvec; m++) {
    ones(n, v);  // w = randn(size(A,1),1);
    printf("v: \n");
    for(i = 0; i < 5; i++) {
      printf("%f \n", v[i]);
    }
    printf("n: %i \n", n);
    printf("evsl.n: %i \n", evsldata.n);
    if(degB) {
      pnav(&pol_sqr, v, vtmp, ww, evsldata); //v = pnav(B, cm, hm mu_sqr, randn(n,1));
      //ChebAv(&pol_sqr, v, vtmp, ww); //v = pnav(B, cm, hm mu_sqr, randn(n,1));
      printf("vtmp1: %f \n", vtmp[0]);
      memcpy(v, vtmp, n * sizeof(double));
    }
    else {
      fprintf(stderr, "LanDos has not yet implemented non-positive degB\n");
      exit(-1);
    }
    // z = B*v
    matvec_B(v,z);


    //--------------------Start of bulk of lanbound.c code
    t = DDOT(&n, z, &one, v, &one); // t =z' * v;
    printf("t: %f \n", t);
    //-------------------- normalize vector
    //                     v = can also  use DNRM2 instead.
    t = 1.0 / sqrt(t); 
    DSCAL(&n, &t, v, &one); // v = v / t;
    DSCAL(&n, &t, z, &one); // z = z / t
    DCOPY(&n, v, &one, V, &one); // V(:, 1 ) = v;
    DCOPY(&n, z, &one, Z, &one); // V(:, 1 ) = v;
    printf("V[0]: %f \n", V[0]);
    printf("Z[0]: %f \n", Z[0]);
    printf("V[n]: %f \n", V[n-1]);
    printf("Z[n]: %f \n", Z[n-1]);
    double wn = 0.0;
    /*-------------------- main Lanczos loop */
    for (j = 0; j < msteps; j++) {
      // w = A*v
      matvec_A(&V[j * n], &V[(j + 1) * n]);
      // w = w - bet * vold
      if (j) {
        nbet = -bet[j - 1];
        DAXPY(&n, &nbet, &V[(j - 1) * n], &one, &V[(j + 1) * n], &one);
      }
      /*-------------------- alp = w' * v */
      alp[j] = DDOT(&n, &V[(j + 1) * n], &one, &V[j * n], &one);
      wn += alp[j] * alp[j];
      //-------------------- w = w - alp * v
      nalp = -alp[j];
      DAXPY(&n, &nalp, &V[j * n], &one, &V[(j + 1) * n], &one);
    printf("V[0]: %f \n", V[0]);
    printf("Z[0]: %f \n", Z[0]);
    printf("V[n]: %f \n", V[n-1]);
    printf("Z[n]: %f \n", Z[n-1]);
    printf("nalp: %f \n", nalp);
    printf("V[%i]: %f \n", (j+1)*n,V[(j+1)*n]);
    printf("Z[%i]: %f \n", (j+1)*n,Z[(j+1)*n]);
    printf("wn: %f \n", wn);
      //-------------------- full reortho
      for (i = 0; i <= j; i++) {
        t = DDOT(&n, &V[(j + 1) * n], &one, &V[i * n], &one);
        double mt = -t;
        DAXPY(&n, &mt, &Z[i * n], &one, &Z[(j + 1) * n], &one);
      }
    printf("Z[%i]: %f \n", (j+1)*n,Z[(j+1)*n]);
      if(degB) {
        pnav(&pol_sol,&V[(j+1) * n] , v, ww, evsldata); //v = pnav(B, cm, hm mu_sqr, randn(n,1));
        memcpy(&V[(j+1) * n], vtmp, n * sizeof(double));
      }
      else {
        fprintf(stderr, "LanDos does not yet support degB <= 0");
        exit(-1);
      }
    printf("V[0]: %f \n", V[0]);
    printf("Z[0]: %f \n", Z[0]);
    printf("V[n]: %f \n", V[n-1]);
    printf("Z[n]: %f \n", Z[n-1]);
    printf("nalp: %f \n", nalp);
    printf("V[%i]: %f \n", (j+1)*n,V[(j+1)*n]);
    printf("Z[%i]: %f \n", (j+1)*n,Z[(j+1)*n]);
    printf("wn: %f \n", wn);

      bet[j] = DDOT(&n, &V[(j + 1) * n], &one, &V[(j + 1) * n], &one);
      printf("Bet[j]: %f", bet[j]);
      exit(-1);
      if (bet[j] * (j + 1) < orthTol * wn) {
        fprintf(stdout, "lanbounds: lucky break, j=%d, beta=%e, break\n", j,
                bet[j]);
        msteps = j + 1;
        break;
      }
      if (bet[j] > tolBdwn) {  // If it's not zero, continue as normal
        wn += 2.0 * bet[j];
        bet[j] = sqrt(bet[j]);
        t = 1.0 / bet[j];
        DSCAL(&n, &t, &V[(j + 1) * n], &one);
      } else {  // Otherwise generate a new vector and redo the previous
                // calculations on it
        fprintf(stderr, "LanDos does not yet have berakdown support \n");
        exit(-1);
        // randn_double(n, v);  // w = randn(size(A,1),1);
        // for (i = 0; i <= j; i++) {
        //   t = DDOT(&n, &V[(j + 1) * n], &one, &V[i * n], &one);
        //   double mt = -t;
        //   DAXPY(&n, &mt, &V[i * n], &one, &V[(j + 1) * n], &one);
        // }
        // bet[j] = DDOT(&n, &V[(j + 1) * n], &one, &V[(j + 1) * n], &one);
        // wn += 2.0 * bet[j];
        // bet[j] = sqrt(bet[j]);
        // t = 1.0 / bet[j];
        // DSCAL(&n, &t, &V[(j + 1) * n], &one);
        // bet[j] = 0;
      }
    }
    //-------------------- end Lanczos loop for this vector
    //-------------------- diagonalize tridiagonal matrix
    SymmTridEig(ritzVal, S, msteps, alp, bet);
    // S = -eigvec
    // ritzVal = diags of D
    //---------------------------------------
    // End of bulk of lanbound.c code
    //---------------------------------------

    // theta = ritzVal = sorted eigenvalues IN ASCENDING ORDER
    for (i = 0; i < msteps; i++) {
      //-------------------- weights for Lanczos quadrature
      // Gamma2(i) = elementwise square of top entry of i-th eginvector
      gamma2[i] = S[i * msteps] * S[i * msteps];
    }
    //-------------------- dos curve parameters
    // Generate DOS from small gaussians centered at the ritz values
    for (i = 0; i < msteps; i++) {
      // As msteps is width of ritzVal -> we get msteps eigenvectors
      const double t = ritzVal[i];
      int numPlaced = 0;
      //-------------------- Place elements close to t in ind
      for (j = 0; j < npts; j++) {
        if (abs(xdos[j] - t) < width) ind[numPlaced++] = j;
      }

      for (j = 0; j < numPlaced; j++)
        y[ind[j]] += gamma2[i] *
                     exp(-((xdos[ind[j]] - t) * (xdos[ind[j]] - t)) / sigma2);
    }
    //-------------------- end vector loop
  }

  double scaling = 1.0 / (nvec * sqrt(sigma2 * PI));
  // y = ydos * scaling
  DSCAL(&npts, &scaling, y, &one);
  DCOPY(&npts, y, &one, ydos, &one);
  simpson2(xdos, y, npts);

  *neig = y[npts - 1] * n;
  free(gamma2);
  free(S);
  free(ritzVal);

  free(alp);
  free(bet);
  free(V);

  free(v);
  free(y);
  free(ind);

  return 0;
}
