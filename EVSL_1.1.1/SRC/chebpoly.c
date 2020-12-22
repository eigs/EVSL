#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "internal_header.h"
/**
 * @file chebpoly.c
 * @brief Computing and applying polynomial filters
 */

/**
 * @brief set default values for polparams struct.
 **/
void set_pol_def(polparams *pol){
  pol->max_deg = 10000;    // max degree allowed
  pol->min_deg = 2;        // min allowed degree
  pol->damping = 2;        // damping. 0 = no damping, 1 = Jackson, 2 = Lanczos
  pol->thresh_ext = 0.20;  // threshold for accepting polynomial for end intervals
  pol->thresh_int = 0.8;   // threshold for accepting polynomial for interior
  pol->mu = NULL;
  pol->deg = 0;            // degree =0 means determine optimal degree.
  pol->intvtol = 1e-9;     // cut-off point of middle interval
}

/**
 * @brief Computes damping coefficient for cheb. expansions.
 *
 * @param[in] m         degree of the polynomial
 * @param[in] damping == 0 --> no damping \n
 *                == 1 --> Jackson \n
 *                == 2 --> Lanczos sigma damping
 * @param[out] jac  output array of dampened coefficients
 *
 * @warning jac must be pre allocated
 *
 * @return 0
 **/
int dampcf(int m, int damping, double *jac) {
  double thetJ = 0.0, thetL = 0.0, a1 = 0.0, a2 = 0.0, dm = (double) m;
  int k, k1;
  if (damping==1){
    thetJ = EVSL_PI/(dm+2);
    a1 = 1/(dm+2);
    a2 = sin(thetJ);
  }
  if (damping == 2)
    thetL = EVSL_PI/(dm+1);
  jac[0] = 0.5;    // <<-- Note this is instead of 1 - in order
  // to reflect  the 1/2 factor in zeroth term
  // of Chebyshev expansion

  for (k=1; k<=m; k++) {
    if (damping == 0) {
      jac[k] = 1.0;
    } else if (damping == 1){
      //-------------------- Note: slightly simpler formulas for jackson:
      k1 = k+1;
      jac[k] = a1*sin(k1*thetJ)/a2 + (1-k1*a1)*cos(k*thetJ) ;
      //-------------------- Lanczos sigma-damping:
    } else {
      jac[k] = sin(k*thetL)/(k*thetL);
    }
  }
  return (0);
}

/**
 *@brief function dif_eval(v, m, thc, jac)
 * evaluates the difference between the right and left
 * values of the polynomial expressed in chebyshev expansion
 * @param[in] v  vector of coefficients [see paper]
 * @param[in] m  degree
 * @param[in] thc  angle theta corresponding to peak of polynomial
 * @param[in] jac vector of damping coefficients
**/

double dif_eval(double *v, int m, double thc, double *jac){
  double fval = 0.0;
  int j;
  for (j=0; j<=m; j++)
    fval += v[j]*cos(j*thc)*jac[j];
  return (fval);
}

/**
 *
 *
 *  @brief function yi = chebxPltd computes yi = p_mu (xi),
 *
 *  where xi is a vectors of values. This can used for plotting
 *  the filter given by mu for example --
 *  Jackson (or other) dampings is not explicitly used here but
 *  is assumed to be multiplied by mu outside this routine.

 *  @param[in] m         degree of the polynomial = length(mu)-1
 *  @param[in] mu        Chev. expansion coefficients in KPM method
 *  @param[in] npts     = number of points in xi, yi
 *  @param[in] xi       = a vector of values where p(xi) is to be computed.
 *  @warning Note xi's must be in [-1 1]
 *
 *  @param[out] yi       = pn(xi(:) )
 *  @return 0
 **/
int chebxPltd(int m, double *mu, int npts, double *xi, double *yi) {
  int k, j, one = 1, n = npts;
  double scal, *vkm1, *vkp1, *vk;
  vkm1 = evsl_Malloc(n, double);
  vkp1 = evsl_Malloc(n, double);
  vk = evsl_Malloc(n, double);
  double *tmp;
  //-------------------- compute p(xi)
  memset(vkm1, 0, n*sizeof(double));
  vecset (n, 1.0, vk) ;
  //-------------------- start yi = mu[0]*vk = mu[0]*ones
  vecset (n, mu[0], yi) ;
  //-------------------- POLYNOMIAL loop
  for (k=0; k<m; k++) {
    //-------------------- generates term of degree k+1
    scal = (k == 0 ? 1.0 : 2.0);
    for (j=0; j<n; j++)
      vkp1[j] = scal*xi[j]*vk[j] - vkm1[j];
    tmp  = vkm1;
    vkm1 = vk;
    vk   = vkp1;
    vkp1 = tmp;
    //-------------------- accumulation of vector.
    evsl_daxpy(&n, &mu[k+1], vk, &one, yi, &one) ;
  }
  //-------------------- done
  evsl_Free(vkm1);
  evsl_Free(vkp1);
  evsl_Free(vk);
  return 0;
}

/**
 * @brief Determines polynomial for end interval cases.
 *
 *
 *  In these cases, polynomial is just a scaled Chebyshev polynomial. However
 *  we need to express it in the same basis as in the other (middle interval)
 *  cases. This function determines this expansion
 *
 * @param[in,out] pol A struct containing the parameters of polynomial.
 * @param[in] aIn The start index of the transformed interval
 * @param[in] bIn The end index of the transformed interval
 * @note [aIn, bIn] is the transformed interval
 *
 *
 * @b Modifies
 * mu Expansion coefficients of best polynomial found.
 * deg: Degree of polynomial
 * gam: Site of delta function that is expanded.
 *      Accurate 'balancing' is done:
 *             If p(t) is best approximation to delta function at gam
 *              then gam is selected so that p(a)=p(b) - within the
 *              tolerance tolBal (set in this function to 1.e-10)
 * bar: If \f$P(\lambda_{i}) \geq \f$ bar, accept eigenvalue as belonging to
 * interval; else reject.
 *
 * @note [a b] is now "interval for eigenvalues to damp",
 * [aIn, bIn] is the interval of interest
 * del is used to expand interval of interest slightly so that pn(t)<bar by
 * some margin.
 *
 *
 *
 * * */
void chext(polparams *pol, double aIn, double bIn){
  int max_deg = pol->max_deg;
  // int min_deg = pol->min_deg;   NOT used
  double thresh = pol->thresh_ext;
  double *mu = pol->mu;
  //-------------------- local variables
  double del = 0.1*sqrt((bIn-aIn)*0.5);
  //double del = 0.0;
  double eps = 1e-13;
  //double eps = 0.0;
  double a, b, x, e, c, sigma, sigma1, sigma_new, g0, g1, gnew, s1, s2, s3;
  double *t0, *t1, *tnew; // coef. of three consecutive cheby. expansions
  double bar, gam;
  int mbest=0;
  int m1 = max_deg+1, j, k;
  //-------------------- local work space
  int work_size = 3*m1;
  //-------------------- this is for the forced degree case
  if (pol->deg > 0){
    thresh = -1.0;
    max_deg = pol->deg;
  }

  if (bIn >= (1.0-eps)) {
    /*-------------------- right interval case */
    x = aIn;
    b = x - del;
    a = -1.0;
    gam = 1.0;
  } else {
    /*-------------------- left interval case */
    x = bIn;
    a = x+del;
    b = 1.0;
    gam = -1.0;
  }
  e = (b-a)*0.5;
  c = (b+a)*0.5;
  sigma1 = e/(gam-c);
  sigma  = sigma1;
  /*-------------------- alloc some work space for coef. of cheby. expansions*/
  t0 = evsl_Calloc(work_size, double);
  t1 = t0+m1;
  tnew = t1+m1;
  // t0 = 1
  t0[0] = 1.0;
  // t1(t) = (t-c)*sig1/e
  t1[0] = -c*(sigma1/e);
  t1[1] = (sigma1/e);
  //-------------------- for evaluating polyn. at x
  g0 = 1.0;
  g1 = (x-c)*sigma1/e;
  if (g1 < thresh){
    mbest = 1;
  }else{
    //-------------------- degree loop : [! make sure to start at k=2]
    for(k=2;k<=max_deg;k++){
      sigma_new = 1.0/(2.0/sigma1-sigma);
      s1 = sigma_new/e;
      s2 = 2*s1*c;
      s3 = sigma*sigma_new;
      for(j=2;j<=k;j++)
        tnew[j-1] = s1*(t1[j]+t1[j-2]);
      tnew[1] = tnew[1] + s1*t1[0];
      tnew[0] = s1*t1[1];
      tnew[k] = s1*t1[k-1];
      for(j=0;j<=k;j++)
        tnew[j] = tnew[j] - s2*t1[j] - s3*t0[j];
      for(j=0;j<=k;j++){
        t0[j] = t1[j];
        t1[j] = tnew[j];
      }
      //-------------------- recurrence to evaluate pn(x)
      gnew = 2*(x-c)*s1*g1 - s3*g0;
      g0 = g1;
      g1 = gnew;
      //-------------------- best degree
      mbest = k;
      if (g1<thresh)
        break;
      sigma = sigma_new;
    }
  }
  memcpy(mu,t1,(mbest+1)*sizeof(double));
  bar = g1;
  pol->deg = mbest;
  pol->bar = bar;
  pol->gam = gam;
  evsl_Free(t0);
}

/**
 * @brief Find the indexofSmallestElement in array
 * @param[in] array Array to find the smallest index of
 * @param[in] size size of the arry
 * @return index of smallest entry in array
**/
int indexofSmallestElement(double *array, int size){
  int index = 0, i;
  for(i = 1; i < size; i++){
    if(array[i] < array[index])
      index = i;
  }
  return index;
}

/**
 * @brief Finds the roots of linear combination of chebyshev polynomials
 * @param[in] m   degree of polynomial
 * @param[in] v difference between cosines on left and right [(3.12) in paper]
 * @param[in] jac   damping coefficients
 * @param[in] tha    theta_a [refer to paper]
 * @param[in] thb    theta_b [refer to paper]
 * @param[out] mu     expansion coefficients.
 * @param[out] thcOut  value of theta_c
 *
 * @warning mu must be preallocated
**/
int rootchb(int m, double *v, double* jac, double tha, double thb, double *mu,
            double *thcOut){
  int MaxIterBalan = 30;     // max steps in Newton to balance interval
  double tolBal;
  // do 2 newton steps -- if OK exit otherwise
  // continue to get root by solving eigv. pb
  int j, it;
  double fval = 0.0, d;
  double fa, fb, thN, thc;
  tolBal = fabs(tha-thb)*1.e-13;
  thc = 0.5*(tha+thb);
  /*-------------------- check whether or not this will work */
  fb = dif_eval(v, m, thb, jac);
  fa = dif_eval(v, m, tha, jac);
  //--------------------this will not work -exit + use higher deg.
  if ((fa>0) || (fb<0))
    return 1 ;
  /*-------------------- Newton iteration to balance the interval*/
  for (it=0; it<=MaxIterBalan; it++) {
    fval = dif_eval(v, m, thc, jac);
    /*-------------------- do one newton step - d= derivative
      thN = thetaNewton*/
    d = 0.0;
    for (j=1; j<=m; j++)
      d += jac[j]*j*sin(j*(thc))*v[j];
    thN = thc + fval/d;
    /*-------------------- test for stopping */
   if ((fabs(fval) < tolBal) || fabs(thc - thN) < DBL_EPSILON * fabs(thc)) {
     break;
   }
   /*-------------------- test for doing a form of bisection */
    if (fval >0){
      if((thN < thb) || (thN > tha))
         thN = 0.5*(thc+tha);
      thb = thc;
      thc = thN;
    } else {
      if((thN < thb) || (thN > tha) )
         thN = 0.5*(thc+thb);
      tha = thc;
      thc = thN;
    }
  }
  /*-------------------- done - return mu and thethac */
  for(j=0; j<=m; j++)
    mu[j] = cos(j*(thc))*jac[j];
  *thcOut = thc;
  return 0;
}

/**
 *
 * @brief Sets the values in pol
 *
 * @param[in] intv An array of length 4, [intv[0, intv[1]] is the interval of
 * desired eigenvalues. [intv[2], intv[3]] is the global interval of all
 * eigenvalues it must contain all eigenvalues of A.
 *
 * @param[in,out] pol The polynomial struct to set the values of.
 *
 * @warning Set the following values of pol \n
 *  mu Expansion coefficients of best polynomial found  \n
 *  gam  Site of delta function that is expanded.
 *          accurate 'balancing' is done:
 *          if p(t) is best approximation to delta function at gam
 *           then gam is selected so that p(a)=p(b) - within the
 *          tolerance tolBal (set in this function to 1.e-10) \n
 *  bar If \f$p(\lambda_i)) \geq \f$ bar, accept eignevalue as belonging to
 *  interval; else reject.\n
 *  deg degree of polynomial \n
 *  dd half-width and... \n
 *  cc ... Center of interval containing all eigenvalues [these are used by
 *  ChebAv]
 *
**/
int find_pol(double *intv, polparams *pol) {
  double *mu, *v, *jac, t=0.0, itv[2],  vals[2];
  int max_deg=pol->max_deg, min_deg=pol->min_deg, damping=pol->damping;
  double tha=0.0, thb=0.0, thc=0.0;
  double gam,  thresh;
  int m, j, nitv, mbest;
  //-------------------- intervals related
  if (check_intv(intv, stdout) < 0) {
    return -1;
  }
  /*-------------------- A few parameters to be set or reset */
  mu = evsl_Malloc(max_deg+1, double);
  pol->mu = mu;
  v = evsl_Malloc(max_deg+1, double);
  jac = evsl_Malloc(max_deg+1, double);
  /*-------------------- A parameter for interval check */
  // double IntTol = 2*DBL_EPSILON; // If b*IntTol>1 accept [a b] extreme
  //double IntTol = 0.0005;
  double IntTol = pol->intvtol;
  double aa, bb;
  aa = evsl_max(intv[0], intv[2]);  bb = evsl_min(intv[1], intv[3]);
  if (intv[0] < intv[2] || intv[1] > intv[3]) {
    fprintf(stdout, " warning [%s (%d)]: interval (%e, %e) is adjusted to (%e, %e)\n",
        __FILE__, __LINE__, intv[0], intv[1], aa, bb);
  }
  double lmin = intv[2], lmax = intv[3];
  /*-------------------- cc, rr: center and half-width of [lmin, lmax] */
  double cc = 0.5 * (lmax + lmin);
  double dd = 0.5 * (lmax - lmin);
  pol->cc = cc;
  pol->dd = dd;
  /*-------------------- adjust intervals just in case. */
  //a = evsl_max(a, lmin);
  //b = evsl_min(b, lmax);
  /*   transform [lmin, lmax] to [-1,1] by y = (x-cc) / dd
   * transform [a, b] to [aT, bT] accordingly */
  aa  = (aa - cc) / dd;
  bb  = (bb - cc) / dd;
  aa  = evsl_max(aa, -1.0);
  bb  = evsl_min(bb,  1.0);
  //printf("transformed interval [%.15e %.15e]\n", a,b);
  thb = acos(bb);
  tha = acos(aa);
  //printf(" a = %e   b = %e  tha = %e thb = %e\n",a,b,tha,thb);
  // right interval. Same thing on other end
  /*-------------------- deal with extremal interval cases */
  if (aa-IntTol <= -1.0) {
    /*-------------------- left interval */
    thc = tha;
    //--------------------  note: p(itv) == where to evaluate p to
    //                      obtain bar values.
    nitv = 1;
    aa = -1.0;
    gam = -1.0;            // set center for a
    pol->type = 1;
  } else if (bb + IntTol >= 1.0) {
    /*-------------------- right interval */
    thc = thb;
    nitv   = 1;
    bb = 1;
    gam = 1;               // set center for b
    pol->type = 2;
  } else {
    /*-------------------- middle interval */
    itv[0] = aa;
    itv[1] = bb;
    nitv = 2;
    pol->type = 0;
  }
  /*-------------------- threshold for finding the best poly */
  if (nitv == 1) {
    /*-------------------- end interval case done separately */
    chext(pol,aa,bb);
  } else {
    /*-------------------- give a starting degree - around 1/(half gap) */
    min_deg = 2 + 0.5/(bb-aa);
    min_deg = evsl_max(min_deg, 2);
    min_deg = evsl_min(min_deg, (max_deg+1)/2);
    // min_deg = 2;
    thresh = pol->thresh_int;
    //-------------------- this is a short-circuit for the
    //                     case of a forced degree
    if (pol->deg > 0){
      thresh = -1;
      max_deg = pol->deg;
    }
    /*-------------------- initialize v vector */
    for (j=0; j<min_deg; j++)
      v[j] = cos(j*thb) - cos(j*tha);
    /*-------------------- DEGREE LOOP -------------------- */
    for (m=min_deg; m < max_deg; m++){
      dampcf(m, damping, jac);
      //-------------------- update v: add one more entry
      v[m] = cos(m*thb)-cos(m*tha);
      //---------------------Balacing the interval + get new mu
      if (rootchb(m, v, jac, tha, thb, mu, &thc)) {
      //-------------------- if m<0 degree is too low - skip
        //printf("rootchb == 1, m = %d, [%.15e, %.15e], tha %e thb %e\n",
        //       m, aa, bb, tha, thb);
        continue;
      }
      //----------------------New center
      gam = cos(thc);
      //-------------------- for scaling
      chebxPltd(m, mu, 1, &gam, &t);
      chebxPltd(m, mu, nitv, itv, vals);
      //-------------------- test for acceptance of this pol.
      //printf("mindeg=%d, m = %d, val %e %e, t %e, thresh %e\n",
      //       min_deg, m, vals[0], vals[1], t, thresh);
      if (vals[0] <= t*thresh && vals[1] <= t*thresh) {
        m++;
        break;
      }
    }
    mbest = m - 1;
    //-------------------- scale the polynomial
    for (j=0; j<=m; j++)
      mu[j] /= t;
    pol->bar = evsl_min(vals[0], vals[1])/t;
    pol->gam = gam;
    pol->deg = mbest;
  }
  /* printf("min_deg %d, max_deg %d, deg %d\n", min_deg, max_deg, mbest); */
  //save_vec(pol->deg+1, mu, "OUT/mu.mtx");
  evsl_Free(v);
  evsl_Free(jac);
  return 0;
}

/**
 * Frees a polparams struct's mu
 *
 * @param[in,out] pol struct to free
 * */
void free_pol(polparams *pol) {
  if (pol->mu) {
    evsl_Free(pol->mu);
  }
}

/**
 * @brief @b Computes y=P(A) v, where pn is a Cheb. polynomial expansion
 *
 * This explicitly calls matvec, so it can be useful for implementing
 * user-specific matrix-vector multiplication.
 *
 * @param[in] pol Struct containing the paramenters and expansion coefficient of
 * the polynomail.
 * @param[in] v input vector
 *
 * @param[out] y p(A)v
 *
 * @b Workspace
 * @param[in, out] w Work vector of length 3*n [allocate before call]
 * @param[in, out] v is untouched
 **/
int ChebAv(polparams *pol, double *v, double *y, double *w) {
#if EVSL_TIMING_LEVEL > 0
  double tt = evsl_timer();
#endif
  const int ifGenEv = evsldata.ifGenEv;
  int n = evsldata.n;
  /*-------------------- unpack pol */
  double *mu = pol->mu;
  double dd = pol->dd;
  double cc = pol->cc;
  int m = pol->deg;
  int one = 1;
  /*-------------------- pointers to v_[k-1],v_[k], v_[k+1] from w */
  double *vk   = w;
  double *vkp1 = w+n;
  double *vkm1 = vkp1+n;
  double *w2 = ifGenEv ? vkm1 + n : NULL;
  /*-------------------- */
  int k;
  double t, s, *tmp, t1= 1.0 / dd, t2 = 2.0 / dd;
  /*-------------------- vk <- v; vkm1 <- zeros(n,1) */
  /* we have to do this copy, because we don't want to alter the
   * elements in v */
  evsl_memcpy_device_to_device(vk, v, n*sizeof(double));
  //memset(vkm1, 0, n*sizeof(double));
  /*-------------------- special case: k == 0 */
  s = mu[0];
  evsl_memcpy_device_to_device(y, v, n*sizeof(double));
  evsl_dscal_device(&n, &s, y, &one);

  /*-------------------- degree loop. k IS the degree */
  for (k=1; k<=m; k++) {
    /*-------------------- y = mu[k]*Vk + y */
    t = k == 1 ? t1 : t2;
    /*-------------------- */
    s = mu[k];
    if (ifGenEv) {
      /*-------------------- Vkp1 = A*B\Vk - cc*Vk */
      solve_B(vk, w2);
      matvec_A(w2, vkp1);
    } else {
      /*-------------------- Vkp1 = A*Vk - cc*Vk */
      matvec_A(vk, vkp1);
    }

    /*
    for (i=0; i<n; i++) {
      vkp1[i] = t*(vkp1[i]-cc*vk[i]) - vkm1[i];
      y[i] += s*vkp1[i];
    }
    */

#ifdef EVSL_USING_CUDA_GPU
    /* fuse 3 blas-1 into 1 kernel */
    evsl_chebAv_device(n, k, t, s, cc, vkp1, vk, vkm1, y);
#else
    double ncc = -cc;
    double dmone = -1.0;
    evsl_daxpy_device(&n, &ncc, vk, &one, vkp1, &one);
    evsl_dscal_device(&n, &t, vkp1, &one);
    if (k > 1) {
      evsl_daxpy_device(&n, &dmone, vkm1, &one, vkp1, &one);
    }
    evsl_daxpy_device(&n, &s, vkp1, &one, y, &one);
#endif

    /*-------------------- next: rotate vectors via pointer exchange */
    tmp = vkm1;
    vkm1 = vk;
    vk = vkp1;
    vkp1 = tmp;
  }

#if EVSL_TIMING_LEVEL > 0
  evslstat.n_polAv ++;
  evslstat.t_polAv += evsl_timer() - tt;
#endif

  return 0;
}

