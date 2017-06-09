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


double rec(const double a) { return 1.0 / a; } //Reciprocal 

double isqrt(const double a) { return 1.0 / sqrt(a); } //Inverse square root


/**
 * @brief @b Computes y=P(A) v, where pn is a Cheb. polynomial expansion 
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
int pnav(double* mu, const int m, const double cc, const double dd, double *v, double *y, double *w) {//Really just ChebAv
  int n = evsldata.n;
  /*-------------------- pointers to v_[k-1],v_[k], v_[k+1] from w */
  double *vk   = w;
  double *vkp1 = w+n;
  double *vkm1 = vkp1+n;
  /*-------------------- */
  int k, i;
  double t, s, *tmp, t1= 1.0 / dd, t2 = 2.0 / dd; 
  /*-------------------- vk <- v; vkm1 <- zeros(n,1) */
  memcpy(vk, v, n*sizeof(double));
  memset(vkm1, 0, n*sizeof(double));
  /*-------------------- special case: k == 0 */
  s = mu[0];
  for (i=0; i<n; i++) {
    y[i] = s*vk[i];
  }
  /*-------------------- degree loop. k IS the degree */
  for (k=1; k<=m; k++) {
    /*-------------------- y = mu[k]*Vk + y */    
    t = k == 1 ? t1 : t2; 
    /*-------------------- */    
    s = mu[k];
    matvec_B(vk, vkp1);

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




/**----------------------------------------------------------------------
 *
 *    Computes the density of states (DOS, or spectral density) using Lanczos
 *    algorithm for the general eigenvalue problem.
 *
 *
 *    @param[in] nvec  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] degB Degree with which B should be approximated by
 *    @param[in] npts number of sample points used for the DOS curve
 *    @param[in] *intv Stores the three intervals of interest
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

int LanDosG(const int nvec, const int msteps, const int degB, int npts, double *xdos, double *ydos,
	     double *neig, const double *const intv, const double tau) {
  //--------------------

  int maxit = msteps,m; //Max number of iteratios
  int n = evsldata.n;

  const int mdeg = 200; 

  polparams pol_sqr;
  polparams pol_sol;
  set_pol_def(&pol_sqr);
  set_pol_def(&pol_sol);

  double tall;

  double* mu_sqr = NULL;
  double* mu_sol = NULL;
  mu_sqr = (double*) malloc(mdeg * sizeof(double));
  mu_sol = (double*) malloc(mdeg * sizeof(double));
  pol_sqr.mu = mu_sqr;
  pol_sol.mu = mu_sol;

  double* vinit;
  Malloc(vinit, n, double);

  if (degB <= 0) {
    printf("LanDos with degB <= 0 not yet implemented");
    exit(-1);
  }
  else {
    lsPol(&intv[4], mdeg, isqrt, tau, &pol_sqr);
    lsPol(&intv[4], mdeg, rec, tau, &pol_sol);
  }
  int *ind;
  Malloc(ind, npts, int);
  double tm,  tmv=0.0;
  double *y; 
  Calloc(y, npts, double);
  //-------------------- to report timings/
  tall = cheblan_timer();
  int i, j, k;

  /*--------------------   frequently used constants  */
  int one = 1;
  /* size of the matrix */
  maxit = min(n, maxit);
  double *gamma2;
  Malloc(gamma2, maxit, double);
  /*-----------------------------------------------------------------------* 
   * *Non-restarted* Lanczos iteration 
   *-----------------------------------------------------------------------
   -------------------- Lanczos vectors V_m and tridiagonal matrix T_m */
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
  double *EvalT, *EvecT;
  Malloc(EvalT, maxit, double);       // eigenvalues of tridia. matrix  T
  Malloc(EvecT, maxit*maxit, double); // Eigen vectors of T
  const double lm = intv[0];
  const double lM = intv[1];
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
  /*-------------------- nmv counts  matvecs */
  int nmv = 0;
  /*-------------------- u  is just a pointer. wk == work space */
  double *wk, *vrand = NULL;
  int wk_size = evsldata.ifGenEv ? 6*n : 4*n;
  Malloc(wk, wk_size, double);
  for (m = 0; m < nvec; m++) {
    randn_double(n,vinit);
    /*-------------------- copy initial vector to Z(:,1) */
    /* Filter the initial vector */
    tm = cheblan_timer();
    pnav(pol_sqr.mu, pol_sqr.deg, pol_sqr.cc, pol_sqr.dd, vinit, V, wk); //Ish
    tmv += cheblan_timer() - tm;

    /*--------------------  normalize it */
    double t;
    if (evsldata.ifGenEv) {
      /* B norm */
      matvec_B(V, Z);
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
    /*-------------------- lanczos vectors updated by rotating pointer*/
    /*-------------------- pointers to Lanczos vectors */
    double *zold, *z, *znew;
    double *v, *vnew;
    /*--------------------  Lanczos recurrence coefficients */
    double alpha, nalpha, beta=0.0, nbeta;
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
      if (evsldata.ifGenEv) {
	/* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
	CGS_DGKS2(n, k, NGS_MAX, Z, V, znew, wk);
	/* -------------- NOTE: B-matvec 
	 *                vnew = B * znew */
	if (degB <= 0) {
	  printf("LanDos with degB <= 0 not yet implemented");
	  exit(-1);
	}else {
	  tm = cheblan_timer();
	  pnav(pol_sol.mu, pol_sol.deg, pol_sol.cc, pol_sol.dd, znew, vnew, wk);  //Ish
	  tmv += cheblan_timer() - tm;
	}
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
        rand_double(n, vnew);
	if (evsldata.ifGenEv) {
	  /* znew = znew - Z(:,1:k)*V(:,1:k)'*znew */
	  CGS_DGKS2(n, k+1, NGS_MAX, V, Z, vnew, wk);
	  /* -------------- NOTE: B-matvec */
	  matvec_B(vnew, znew);
	  nmv ++;
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
    }
    SymmTridEig(EvalT, EvecT, maxit, dT, eT);
    for (i = 0; i < maxit; i++) {
      //-------------------- weights for Lanczos quadrature
      // Gamma2(i) = elementwise square of top entry of i-th eginvector
      gamma2[i] = EvecT[i * maxit] * EvecT[i * maxit];
    }
    //-------------------- dos curve parameters
    // Generate DOS from small gaussians centered at the ritz values
    for (i = 0; i < msteps; i++) {
      // As msteps is width of ritzVal -> we get msteps eigenvectors
      const double t = EvalT[i];
      int numPlaced = 0;
      //-------------------- Place elements close to t in ind
      for (j = 0; j < npts; j++) {
        if (abs(xdos[j] - t) < width) ind[numPlaced++] = j;
      }
      for (j = 0; j < numPlaced; j++)
        y[ind[j]] += gamma2[i] *
	  exp(-((xdos[ind[j]] - t) * (xdos[ind[j]] - t)) / sigma2);
    }
  }
  double scaling = 1.0 / (nvec * sqrt(sigma2 * PI));
  // y = ydos * scaling
  DSCAL(&npts, &scaling, y, &one);
  DCOPY(&npts, y, &one, ydos, &one);
  simpson2(xdos, y, npts);
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
  free_pol(&pol_sqr);
  free_pol(&pol_sol);
  if (vrand) {
    free(vrand);
  }
  if (evsldata.ifGenEv) {
    free(Z);
  }
  /*-------------------- record stats */
  tall = cheblan_timer() - tall;
  /*-------------------- print stat */

  return 0;
}
