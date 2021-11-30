#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include "internal_header.h"

/**
 * @file ratfilter.c
 * @brief Computing and applying rational filters
 */
/**------------------------- Cauchy integration-based filter --------------
 * @brief Compute the locations of the poles
 *
 * @param[in] method    0 for Guass Legendre; 1 for midpoint
 * @param[in] n         Number of poles in the upper half plane
 * @param[out] zk   Vector of pole locations
 *
 *----------------------------------------------------------------------*/
void contQuad(int method, int n, EVSL_Complex* zk) {
  int i, m, INFO;
  double *beta, *D, *Z, *WORK;
  char JOBZ = 'V';
  EVSL_Complex tmp2;
  if (method == 0) {
    m = n-1;
    D = evsl_Malloc(n, double);
    Z = evsl_Malloc(n*n, double);
    WORK = evsl_Malloc(2*n-2, double);
    for (i=0; i<n; i++) {
      D[i] = 0.0;
    }
    beta = evsl_Malloc(m, double);
    for (i=0; i<m; i++) {
      beta[i] = 0.5/(sqrt(1-pow(2*(i+1),-2)));
    }
    evsl_dstev(&JOBZ, &n, D, beta, Z, &n, WORK, &INFO);
    for (i=0; i<n; i++) {
      tmp2 = I*EVSL_PI/2.0*(1.0-D[i]);
      zk[i] = cexp(tmp2);
    }
    evsl_Free(beta);
    evsl_Free(D);
    evsl_Free(Z);
    evsl_Free(WORK);
  } else if (method == 1) {
    for (i=0; i<n; i++) {
      tmp2 = EVSL_PI*I*(2*(i+1)-1)/(2.0*n);
      zk[i] = cexp(tmp2);
    }
  }
}

  /**------------------Multiple pole rational filter evaluation --------------
   * @brief Compute the function value of the multiple pole rational filter at real locations
   * @param[in] n      number of the pole
   * @param[in] mulp   multiplicity of the pole
   * @param[in] zk     array containing the poles.
   * @param[in] alp    fractional expansion coefficients
   * @param[in] m      number of locations to be evaluated
   * @param[in] z      real locations to be evaluated
   *
   * @param[out] xx    : function values at real locations z
   *
   *----------------------------------------------------------------------*/
void ratf2p2(int n, int *mulp, EVSL_Complex *zk, EVSL_Complex* alp, int m,
             double *z, double *xx) {
  EVSL_Complex y, x, t;
  int i, j, k, k1, k2;
  for (k2=0; k2<m; k2++) {
    k = 0;
    y = 0.0 + 0.0*I;
    for (j=0; j<n; j++) {
      x = 1.0 / (z[k2]-zk[j]);
      k1 = k + mulp[j];
      t = 0.0+0.0*I;
      for (i=k1-1;i>=k;i--) {
        t = x*(alp[i]+t);
      }
      y = y+t;
      k = k1;
    }
    xx[k2] = 2.0*creal(y);
  }
}



/**
 * @brief Get the fraction expansion of 1/[(z-s1)^k1 (z-s2)^k2]
 * */
void pfe2(EVSL_Complex s1, EVSL_Complex s2, int k1, int k2,
          EVSL_Complex* alp, EVSL_Complex* bet) {
  int i;
  EVSL_Complex d, xp;
  if (cabs(s1-s2) < 1.0e-12 * (cabs(s1)+cabs(s2))) {
    for (i=0; i<k1+k2; i++) {
      alp[i] = 0.0;
    }
    alp[k1+k2-1] = 1;
  } else if ((k1 == 1) && (k2 == 1)) {
    d = s1-s2;
    alp[k1-1] = 1.0 / d;
    bet[k1-1] = -alp[k1-1];
  } else {
    d = 1.0 + 0.0*I;
    xp = 1.0 / cpow((s1-s2),k2);
    for (i=0; i<k1; i++) {
      alp[k1-i-1] = d * xp;
      xp = xp / (s1-s2);
      d = -d * (k2+i) / (i+1.0);
    }
    d = 1.0 + 0.0*I;
    xp = 1.0 / cpow((s2-s1),k1);
    for (i=0; i<k2; i++) {
      bet[k2-i-1] = d * xp;
      xp = xp / (s2-s1);
      d = -d * (k1+i) / (i+1.0);
    }
  }
}

  /**
   * @brief Integration of 1/[(z-s1)^k1 (z-s2)^k2] from a to b
   */
EVSL_Complex integg2(EVSL_Complex s1, EVSL_Complex s2,
                       EVSL_Complex* alp, int k1, EVSL_Complex* bet,
                       int k2, double a, double b) {
  EVSL_Complex t, t1, t0, scal;
  int k;
  t = 0.0 + 0.0*I;
  t1 = 0.0 + 0.0*I;
  t0 = 0.0 +0.0*I;
  for (k=0; k<k1; k++) {
    scal = alp[k];
    if (k==0) {
      t0 = scal*clog((b-s1)/(a-s1));
    } else {
      t = t - (scal*1.0/k) * (1.0/cpow((b-s1),k)-1.0/cpow((a-s1),k));
    }
  }
  for (k=0; k<k2; k++) {
    scal = bet[k];
    if (k==0) {
      t1 = scal*clog((b-s2)/(a-s2));
    } else {
      t = t - (scal*1.0/k)*(1.0/cpow((b-s2),k)-1.0/cpow((a-s2),k));
    }
  }
  t = t + (t0+t1);
  return t;
}


/**
 *------------------multiple pole LS rational filter weights--------------
 * @brief Compute the LS weight for each multiple pole
 * @param[in] n      number of poles in the upper half plane
 * @param[in] zk     pole locations
 * @param[in] mulp    multiplicity of each pole
 * @param[in] lambda LS integration weight for [-1, 1]
 *
 * @param[out] omega LS weight for each pole
 *
 *----------------------------------------------------------------------*/
void weights(int n, EVSL_Complex* zk, int* mulp, double lambda,
             EVSL_Complex* omega) {
  int INFO;
  int nrhs = 1;
  int *ipiv;
  int m;
  double mu = 10.0;
  int i, j, ii, jj, ki, kj, n1, n2, nf=0, k1, k2;
  EVSL_Complex s1, s2, s3, t;
  EVSL_Complex *rhs, *A, *B, *mat, *alp, *bet;
  double scaling;
  for (i=0; i<n; i++) {
    nf += mulp[i];
  }
  m = 2*nf;
  ipiv = evsl_Malloc(m, int);
  rhs = evsl_Malloc(m, EVSL_Complex);
  A = evsl_Malloc(nf*nf, EVSL_Complex);
  B = evsl_Malloc(nf*nf, EVSL_Complex);
  mat = evsl_Malloc(4*nf*nf, EVSL_Complex);
  for(i=0; i<nf; i++) {
    for(j=0; j<nf; j++) {
      A[i*nf+j] = 0.0 + 0.0*I;
      B[i*nf+j] = 0.0 + 0.0*I;
    }
  }
  if (fabs(lambda) < 1.0e-12) {
    lambda = 1.0e-5;
  }
  ki = 0;
  for (ii=0; ii<n; ii++) {
    s1 = zk[ii];
    n1 = mulp[ii];
    kj = 0;
    s3 = conj(s1);
    for (i=0; i<n1; i++) {
      if (i==0) {
        rhs[ki+i] = lambda*clog((s3-1.0)/(s3+1.0));
      } else {
        rhs[ki+i] = -lambda*(1.0/(i))*(1/cpow((1.0-s3),i)-1.0/cpow((-1.0-s3),i));
      }
    }
    for (jj=0; jj<n; jj++) {
      s2 = zk[jj];
      n2 = mulp[jj];
      for (i=0; i<n1; i++) {
        for (j=0; j<n2; j++) {
          s3 = conj(s2);
          if (cabs(s1-s3) < 1.0e-12*(cabs(s1)+cabs(s3))) {
            alp = evsl_Malloc(i+j+2, EVSL_Complex);
            bet = evsl_Malloc(1, EVSL_Complex);
            k1 = i+1+j+1;
            k2 = 0;
          } else {
            alp = evsl_Malloc(i+1, EVSL_Complex);
            bet = evsl_Malloc(j+1, EVSL_Complex);
            k1 = i+1;
            k2 = j+1;
          }
          pfe2(s1, s3, k1, k2, alp, bet);
          t = integg2(s1, s3, alp, k1, bet, k2, -mu, mu);
          t += (lambda-1)*integg2(s1, s3, alp, k1, bet, k2, -1.0, 1.0);
          A[(ki+i)*nf+kj+j] = t;
          evsl_Free(bet);
          evsl_Free(alp);
          if (cabs(s1-s2) < 1.0e-12*(cabs(s1)+cabs(s2))) {
            alp = evsl_Malloc(i+j+2, EVSL_Complex);
            bet = evsl_Malloc(1, EVSL_Complex);
            k1 = i+1+j+1;
            k2 = 0;
          } else {
            alp = evsl_Malloc(i+1, EVSL_Complex);
            bet = evsl_Malloc(j+1, EVSL_Complex);
            k1 = i+1;
            k2 = j+1;
          }
          pfe2(s1, s2, k1, k2, alp, bet);
          t = integg2(s1, s2, alp, k1, bet, k2, -mu, mu);
          t += (lambda-1)*integg2(s1, s2, alp, k1, bet, k2, -1.0, 1.0);
          B[(ki+i)*nf+kj+j] = t;
          evsl_Free(alp);
          evsl_Free(bet);
        }
      }
      kj = kj+n2;
    }
    ki = ki+n1;
  }
  for (i=nf; i<2*nf; i++) {
    rhs[i] = conj(rhs[i-nf]);
  }
  /*---form mat = [A,B;conj(B),conj(A)]---*/
  /* Note that mat is stored column-wise for lapack routine */
  for (i=0; i<nf; i++) {
    for(j=0; j<nf; j++) {
      mat[i+j*m] = conj(A[i*nf+j]);
    }
  }
  for (i=0; i<nf; i++) {
    for (j=nf; j<m; j++) {
      mat[i+j*m] = conj(B[i*nf+j-nf]);
    }
  }
  for (i=nf; i<m; i++) {
    for (j=0; j<nf; j++) {
      mat[i+j*m] = B[(i-nf)*nf+j];
    }
  }
  for (i=nf; i<m; i++) {
    for (j=nf; j<m; j++) {
      mat[i+j*m] = A[(i-nf)*nf+j-nf];
    }
  }
  evsl_zgesv(&m, &nrhs, mat, &m, ipiv, rhs, &m, &INFO);
  for(i=0;i<nf;i++) {
    omega[i] = rhs[i];
  }

  /* Scale coefs to let the filter pass through [-1, 0.5] */
  double aa = 1.0;
  ratf2p2(n, mulp, zk, omega, 1, &aa, &scaling);
  scaling = 0.5 / scaling;
  for (i=0; i<nf; i++) {
    omega[i] *= scaling;
  }

  evsl_Free(A);
  evsl_Free(B);
  evsl_Free(rhs);
  evsl_Free(mat);
  evsl_Free(ipiv);
}


/** Transform poles and weights computed on [-1, 1] to [a, b]
   * @brief  Compute the weights and pole locations on [a, b]
   * @param[in] n     number of poles used in the upper half plane
   * @param[in] a,b   [a, b] is the interval of desired eigenvalues
   * @param[in,out] zk    location of the poles
   * @param[in] mulp   multiplicity of the poles
   *
   * @param[out] omegaM: multiple LS weights
   *
   *----------------------------------------------------------------------*/
int scaleweigthts(int n, double a, double b, EVSL_Complex *zk, int* mulp,
                  EVSL_Complex* omegaM) {
  int i, j, k, nf=0;
  double c, h;
  c = 0.5 * (a + b);
  h = 0.5 * (b - a);
  for (i=0; i<n; i++) {
    nf += mulp[i];
    zk[i] = h*zk[i]+c;
  }
  /* Transform the coefs for multiple poles */
  double tmp;
  k = -1;
  for (i=0; i<n; i++) {
    for (j=0; j<mulp[i]; j++) {
      k = k+1;
      omegaM[k] = omegaM[k]*cpow(h,j+1);
    }
  }
  /* Scale ration function to let it pass through [a, 1/2] */
  ratf2p2(n, mulp, zk, omegaM, 1, &a, &tmp);
  tmp = 0.5 / tmp;
  for (i=0; i<nf; i++) {
    omegaM[i] = omegaM[i] * tmp;
  }
  return 0;
}


/**
 * @brief Sets default values for ratparams struct
 * */
void set_ratf_def(ratparams *rat) {
  // -------------------- this sets default values for ratparams struct.
  rat->num    = 1;           // number of the poles
  rat->pw     = 2;           // default multplicity of each pole
  rat->method = 1;           // using poles from mid-point rule
  rat->beta   = 0.01;        // beta in LS approximation
  rat->bar    = 0.5;         // this is fixed for rational filter
  rat->aa     = -1.0;        // left endpoint of interval
  rat->bb     = 1.0;         // right endpoint of interval
  //rat->cc = 0.0;           // center of interval
  //rat->dd = 1.0;           // width of interval
}

/**----------------------------------------------------------------------
 * @param[in] intv  = an array of length 4
 *         [intv[0], intv[1]] is the interval of desired eigenvalues
 *         [intv[2], intv[3]] is the global interval of all eigenvalues
 *         it must contain all eigenvalues of A
 *
 * @param[out] rat
 * these are set in rat struct:\n
 *   omega : expansion coefficients of rational filter \n
 *    zk   : location of the poles used\n
 *    aa  : adjusted left endpoint\n
 *    bb  : adjusted right endpoint\n
 *    dd  : half-width and.. \n
 *    cc  : ..center of interval\n
 *
 *--------------------------------------------------------------------*/
int find_ratf(double *intv, ratparams *rat) {
  EVSL_Complex *omega; // weights of the poles
  EVSL_Complex *zk;    // location of the poles
  int *mulp;             // multiplicity of the each pole
  int n = rat->num, i, pow = 0, pw = rat->pw, method = rat->method;
  double beta = rat->beta;
  /*-------------------- A few parameters to be set or reset */
  mulp = evsl_Malloc(n, int);
  zk = evsl_Malloc(n, EVSL_Complex);
  for (i=0; i<n; i++) { // set the multiplicity of each pole
    mulp[i] = pw;
    pow += mulp[i];
  }
  rat->zk = zk;
  rat->mulp = mulp;
  rat->pow = pow; // total multiplicity of the poles
  omega = evsl_Malloc(pow, EVSL_Complex);
  rat->omega = omega;
  //-------------------- intervals related
  if (check_intv(intv, stdout) < 0) {
    return -1;
  }
  double aa, bb;
  aa = evsl_max(intv[0], intv[2]);  bb = evsl_min(intv[1], intv[3]);
  if (intv[0] < intv[2] || intv[1] > intv[3]) {
    fprintf(stdout, " warning [%s (%d)]: interval (%e, %e) is adjusted to (%e, %e)\n",
            __FILE__, __LINE__, intv[0], intv[1], aa, bb);
  }
  //double lmin = intv[2], lmax = intv[3];
  /*-------------------- */
  rat->aa = aa;
  rat->bb = bb;
  /*-------------------- cc, rr: center and half-width of [aa, bb] */
  //double cc = 0.5 * (aa + bb);
  //double dd = 0.5 * (bb - aa);
  //rat->cc = cc;
  //rat->dd = dd;
  /*------------ compute the location of the poles */
  contQuad(method, n, zk);
  /*------------ compute expansion coefficients of rational filter on [-1, 1] */
  weights(n, zk, mulp, beta, omega);
  /*-------------------- compute expansion coefficients on [aa, bb]*/
  scaleweigthts(n, aa, bb, zk, mulp, omega);

  rat->ASIGBsol = NULL;

  return 0;
}

void free_rat(ratparams *rat) {
  evsl_Free(rat->mulp);
  evsl_Free(rat->omega);
  evsl_Free(rat->zk);
  evsl_Free(rat->ASIGBsol);
}

/**
 * @brief Apply rational filter R to a vetor b
 *
 * For generalized e.v problem x = L' * (A-SB) \ L*b [w:=work]
 * x = L * b
 * w = (A-sB) \ x
 * x = L' * w
 *
 * @param[in] rat ratparams struct
 * @param[in] n Length of array
 * @param[in] b x = L * b
 * @param w6 Work array of size 4*n for standard ev problem,
 *                         size 6*n for generalized ev problem
 *
 * @param[out] x Becomes R(A)b
 *
 */
void RatFiltApply(int n, ratparams *rat, double *b, double *x, double *w6) {
#if EVSL_TIMING_LEVEL > 0
  double tt = evsl_timer();
#endif
  const int ifGenEv = evsldata.ifGenEv;
  int jj, kk, k=0, kf;
  int *mulp = rat->mulp;
  int num = rat->num;
  EVSL_Complex *omega = rat->omega;
  double dtwo = 2.0;
  double done = 1.0;
  int one = 1;

  double *xr, *xz, *bz, *br, *yr=NULL, *yz=NULL;
  double zkr, zkc;
  xr = w6;
  xz = xr + n;
  bz = xz + n;
  br = bz + n;
  if (ifGenEv) {
    yr = br + n;
    yz = yr + n;
  }
  /*------------------ loop through each pole */
  for (kk=0; kk<num; kk++) {
    /*---------------- solver for A-s[kk]*B */
    EVSLASIGMABSol *sol = &rat->ASIGBsol[kk];
    kf = k + mulp[kk];
    /*------------------ power loop */
    for (jj=kf-1; jj>=k; jj--) {
      /*---------------- weight */
      zkr = creal(omega[jj]);
      zkc = cimag(omega[jj]);
      /*---------------- initilize the right hand side */
      evsl_memcpy_device_to_device(br, b, n*sizeof(double));
      evsl_memcpy_device_to_device(bz, b, n*sizeof(double));
      evsl_dscal_device(&n, &zkr, br, &one);
      evsl_dscal_device(&n, &zkc, bz, &one);
      if (jj != kf-1) {
        evsl_daxpy_device(&n, &done, xr, &one, br, &one);
        evsl_daxpy_device(&n, &done, xz, &one, bz, &one);
      }
      /*---------------- solve shifted system */
      if (ifGenEv) {
        if (jj > k) {
          //(sol->func)(n, br, bz, yr, yz, sol->data);
          solve_ASigB(sol, n, br, bz, yr, yz);
          matvec_B(yr, xr);
          matvec_B(yz, xz);
        } else {
          /*------------- jj == k */
          //(sol->func)(n, br, bz, xr, xz, sol->data);
          solve_ASigB(sol, n, br, bz, xr, xz);
        }
      } else {
        //(sol->func)(n, br, bz, xr, xz, sol->data);
        solve_ASigB(sol, n, br, bz, xr, xz);
      }
    }
    /*------------------ solution (real part) */
    if (kk) {
      evsl_daxpy_device(&n, &dtwo, xr, &one, x, &one);
    } else {
      evsl_memcpy_device_to_device(x, xr, n*sizeof(double));
      evsl_dscal_device(&n, &dtwo, x, &one);
    }
    k = kf;
  }

#if EVSL_TIMING_LEVEL > 0
  evslstat.n_ratAv ++;
  evslstat.t_ratAv += evsl_timer() - tt;
#endif
}

