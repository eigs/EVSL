#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

/**----------------------------------------------------------------------
 *
 * @brief This function  computes the  coefficients of the  density of
 * states  in  the  chebyshev   basis.   It  also  returns  the
 * estimated number of eigenvalues in the interval given by intv.
 * @param Mdeg     degree of polynomial to be used. 
 * @param damping  type of damping to be used [0=none,1=jackson,2=sigma]
 * @param nvec     number of random vectors to use for sampling
 * @param intv   an array of length 4  \n
 *                 [intv[0] intv[1]] is the interval of desired eigenvalues 
 *                 that must be cut (sliced) into n_int  sub-intervals \n
 *                 [intv[2],intv[3]] is the global interval of eigenvalues 
 *                 it must contain all eigenvalues of A \n
 * @param[out] mu   array of Chebyshev coefficients 
 * @param[out] ecnt estimated num of eigenvalues in the interval of interest
 *
 *----------------------------------------------------------------------*/
int kpmdos(int Mdeg, int damping, int nvec, double *intv,
    double *mu, double *ecnt) {
  /*-------------------- initialize variables */
  int n;
  csrMat *A;
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    A = evsldata.A;
    n = A->nrows;
  }
  double *vkp1, *v, *vkm1, *vk, *jac;
  double *w = NULL;
  /*-------------------- workspace for generalized eigenvalue prob */
  if (evsldata.ifGenEv) {
    Malloc(w, n, double);
  }
  Malloc(vkp1, n, double);
  Malloc(v, n, double);
  Malloc(vkm1, n, double);
  Malloc(vk, n, double);
  Malloc(jac, Mdeg+1, double);
  double *tmp,  ctr, wid; 
  double scal, t, tcnt, beta1, beta2, aa, bb;
  int k, k1, i, m, mdegp1, one=1;
  //-------------------- check if the interval is valid
  if (check_intv(intv, stdout) < 0) {
    return -1;
  }

  aa = max(intv[0], intv[2]);  bb = min(intv[1], intv[3]);
  if (intv[0] < intv[2] || intv[1] > intv[3]) {
    fprintf(stdout, " warning [%s (%d)]: interval (%e, %e) is adjusted to (%e, %e)\n",
            __FILE__, __LINE__, intv[0], intv[1], aa, bb);
  }
  /*-------------------- some needed constants */
  ctr  = (intv[3]+intv[2])/2.0;
  wid  = (intv[3]-intv[2])/2.0;
  t = max(-1.0+DBL_EPSILON, (aa-ctr)/wid);
  beta1 = acos(t);
  t = min(1.0-DBL_EPSILON, (bb-ctr)/wid);
  beta2 = acos(t);
  /*-------------------- compute damping coefs. */
  dampcf(Mdeg, damping, jac);
  /*-------------------- readjust jac[0] it was divided by 2 */
  jac[0] = 1.0;
  memset(mu,0,(Mdeg+1)*sizeof(double));
  /*-------------------- seed the random generator */
  i = time_seeder();
  srand(i);
  tcnt = 0.0;
  /*-------------------- for random loop */
  for (m=0; m<nvec; m++) {
    if (evsldata.ifGenEv) {
      /* unit 2-norm v */
      rand_double(n, v);
      t = 1.0 / DNRM2(&n, v, &one);
      DSCAL(&n, &t, v, &one);  
      /*  w = L^{-T}v */
      evsldata.LTsol->func(v, w, evsldata.Bsol->data);
      /* v = B*w */
      matvec_B(w, v);
      t = DDOT(&n, v, &one, w, &one);
      memcpy(vk, w, n*sizeof(double));
    } else {
      /* unit 2-norm */
      rand_double(n, v);
      t = 1.0 / DNRM2(&n, v, &one);
      DSCAL(&n, &t, v, &one);
      memcpy(vk, v, n*sizeof(double));
    }
    
    memset(vkm1, 0, n*sizeof(double));
    mu[0] += jac[0];
    //-------------------- for eigCount
    tcnt -= jac[0]*(beta2-beta1);  
    /*-------------------- Chebyshev (degree) loop */
    for (k=0; k<Mdeg; k++){
      /*-------------------- Cheb. recurrence */
      if (evsldata.ifGenEv) {
        /* v_{k+1} := B \ A * v_k (partial result) */
        matvec_A(vk, w);
        evsldata.Bsol->func(w, vkp1, evsldata.Bsol->data);
      } else {
        matvec_A(vk, vkp1);
      }
      scal = k==0 ? 1.0 : 2.0;
      scal /= wid;
      for (i=0; i<n; i++) {
        vkp1[i] = scal*(vkp1[i]-ctr*vk[i]) - vkm1[i];
      }
      //-------------------- rotate pointers to exchange vectors
      tmp = vkm1;
      vkm1 = vk;
      vk = vkp1;
      vkp1 = tmp;
      /*-------------------- accumulate dot products for DOS expansion */
      k1 = k+1;
      t = 2*jac[k1] * DDOT(&n, vk, &one, v, &one);
      mu[k1] += t;
      /*-------------------- for eig. counts */
      tcnt -= t*(sin(k1*beta2)-sin(k1*beta1))/k1;  
    }
  }
  //--------------------change of interval + scaling in formula
  t = 1.0 /(((double)nvec)*PI) ;
  mdegp1 = Mdeg+1;
  DSCAL(&mdegp1, &t, mu, &one) ;
  tcnt *= t*((double) n);
  *ecnt = tcnt;
  /*-------------------- free memory  */
  free(vkp1);
  free(v);
  free(vkm1);
  free(vk);
  free(jac);
  if (evsldata.ifGenEv) {
    free(w);
  }

  return 0;
}

  /**  
  * @brief Computes the integrals \f$\int_{xi[0]}^{xi[j]} p(t) dt\f$
  *  where p(t) is the approximate DOS as given in the KPM method
  *  in the expanded form:  \f$\sum mu_i C_i /\sqrt{1-t^2}\f$
  **/
void intChx(int Mdeg, double *mu, int npts, double *xi, double *yi) {
  //
  int ndp1, j, k;
  double val0, theta0, *thetas;
  Malloc(thetas, npts, double);
  ndp1   = Mdeg+1; 
  //  if (xi[0]<-1.0) xi[0] = -1; 
  //if (xi[npts-1]> 1.0) xi[npts-1]  = 1; 

  for (j=0; j<npts; j++)
    thetas[j] = acos(xi[j]);
  theta0 = thetas[0];
  for (j=0; j<npts; j++) 
    yi[j] = mu[0]*(theta0 - thetas[j]);
  //-------------------- degree loop  
  for (k=1; k<ndp1; k++){
    val0 = sin(k*theta0)/k;
    //-------------------- points loop
    for (j=0; j<npts; j++)
      yi[j] += mu[k]*(val0 - sin(k*thetas[j])/k);
  }
  free (thetas);
}

/**----------------------------------------------------------------------- 
 * @brief given the dos function defined by mu find a partitioning
 * of sub-interval [a,b] of the spectrum so each 
 * subinterval has about the same number of eigenvalues
 * Mdeg = degree.. mu is of length Mdeg+1  [0---> Mdeg]
 * on return [ sli[i],sli[i+1] ] is a subinterval (slice).
 *
 * @param *sli  see above (output)
 * @param *mu   coeffs of polynomial (input)
 * @param Mdeg     degree of polynomial (input)
 * @param *intv  an array of length 4 
 *                [intv[0] intv[1]] is the interval of desired eigenvalues
 *                that must be cut (sliced) into n_int  sub-intervals
 *                [intv[2],intv[3]] is the global interval of eigenvalues
 *                it must contain all eigenvalues of A
 * @param n_int   number of slices wanted (input)
 * @param npts      number of points to use for discretizing the interval
 *                [a b]. The more points the more accurate the intervals. 
 *                it is recommended to set npts to a few times the number 
 *                of eigenvalues in the interval [a b] (input). 
 *
 *----------------------------------------------------------------------*/
int spslicer(double *sli, double *mu, int Mdeg, double *intv, int n_int, int npts) {
  int ls, ii, err=0;
  double  ctr, wid, aL, bL, target, aa, bb;

  if (check_intv(intv, stdout) < 0) {
    return -1;
  }

  // adjust a, b: intv[0], intv[1]
  aa = max(intv[0], intv[2]);  bb = min(intv[1], intv[3]);
  if (intv[0] < intv[2] || intv[1] > intv[3]) {
    fprintf(stdout, " warning [%s (%d)]:  interval (%e, %e) is adjusted to (%e, %e)\n",
        __FILE__, __LINE__, intv[0], intv[1], aa, bb);
  }

  //-------------------- 
  memset(sli,0,(n_int+1)*sizeof(double));
  //-------------------- transform to ref interval [-1 1]
  //-------------------- take care of trivial case n_int==1
  if (n_int == 1){
    sli[0] = intv[0];
    sli[1] = intv[1];
    return 0;
  }
  //-------------------- general case 
  ctr = (intv[3] + intv[2])/2;
  wid = (intv[3] - intv[2])/2;
  aL  = (aa - ctr)/wid;   // (a - ctr)/wid 
  bL  = (bb - ctr)/wid;   // (b - ctr)/wid
  aL = max(aL,-1.0);
  bL = min(bL,1.0);
  npts = max(npts,2*n_int+1);
  double *xi, *yi;
  Malloc(xi, npts, double);
  Malloc(yi, npts, double);
  linspace(aL, bL, npts, xi);
  //printf(" aL %15.3e bL %15.3e \n",aL,bL);
  //-------------------- get all integrals at the xi's 
  //-------------------- exact integrals used.
  intChx(Mdeg, mu, npts, xi, yi) ; 
  //-------------------- goal: equal share of integral per slice
  target = yi[npts-1] / (double)n_int;
  ls = 0;
  ii = 0;
  // use the unadjust left boundary
  sli[ls] = intv[0];
  //-------------------- main loop 
  while (++ls < n_int) {
    while (ii < npts && yi[ii] < target) {
      ii++;
    }
    if (ii == npts) {
      break;
    }
    //-------------------- take best of 2 points in interval
    if ( (target-yi[ii-1]) < (yi[ii]-target)) {
      ii--;
    }
    sli[ls] = ctr + wid*xi[ii];
    //-------------------- update target.. Slice size adjusted
    target = yi[ii] + (yi[npts-1] - yi[ii])/(n_int-ls);
    //printf("ls %d, n_int %d, target %e\n", ls, n_int, target);
  }

  // use the unadjust left boundary
  sli[n_int] = intv[1];

  //-------------------- check errors
  if (ls != n_int) {
    err = 1;
  }
  for (ii=1; ii<=n_int; ii++) {
    if (sli[ii] <= sli[ii-1]) {
      err += 2;
      break;
    }
  }

  if (err) {
    save_vec(npts, xi, "OUT/xi.out");
    save_vec(npts, yi, "OUT/yi.out");
  }

  /*-------------------- free arrays */
  free(xi);
  free(yi);
  return err;
}



