#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

void set_pol_def(polparams *pol){
  // -------------------- this sets default values for polparams struct. 
  pol->max_deg = 300;      // max degree allowed
  pol->min_deg = 2;        // min allowed degree
  pol->damping = 2;        // damping. 0 = no damping, 1 = Jackson, 2 = Lanczos
  pol->thresh_ext = 0.50;  // threshold for accepting polynomial for end intervals 
  pol->thresh_int = 0.8;   // threshold for accepting polynomial for interior
  //intervals 
  pol-> deg = 0;           // degree =0 means determine optimal degree.
}

int dampcf(int m, int damping, double *jac){
  //------------------------------------------------------------
  // computes damping coefficient for cheb. expansions
  // damping =  0  -->  no damping
  // damping =  1  -->  Jackson
  // damping =  2  -->  Lanczos sigma damping
  //------------------------------------------------------------
  double thetJ = 0.0, thetL = 0.0, a1 = 0.0, a2 = 0.0, dm = (double) m;
  int k, k1;
  if (damping==1){
    thetJ = PI/(dm+2);
    a1 = 1/(dm+2);
    a2 = sin(thetJ);
  } 
  if (damping == 2) 
    thetL = PI/(dm+1);   
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

int chebxPltd(int m, double *mu, int npts, double *xi, double *yi) {
  /*--------------------------------------------------------------  
    function yi = chebxPltd -- computes yi = p_mu (xi),
    where xi is a vectors of values. This can used for plotting 
    the filter given by mu for example --  
    Jackson (or other) dampings is not explicitly used here but
    is assumed to be multiplied by mu outside this routine.

    m        = degree of the polynomial = length(mu)-1
    mu       = Chev. expansion coefficients in KPM method
    npts     = number of points in xi, yi
    xi       = a vector of values where p(xi) is to be computed.
    Note xi's must be in [-1 1] 
    OUT:	  
    yi       = pn(xi(:) )
--------------------------------------------------------------*/
  int k, j, one = 1, n = npts;
  double scal;
  double *vkm1 = malloc(n*sizeof(double));
  double *vkp1 = malloc(n*sizeof(double));
  double *vk   = malloc(n*sizeof(double));
  double *tmp;
  //-------------------- compute p(xi)
  memset(vkm1, 0, n*sizeof(double));
  vecset (n, 1.0, vk) ; 
  //-------------------- start yi = mu[0]*vk = mu[0]*ones
  vecset (n, mu[0], yi) ; 
  //-------------------- POLYNOMIAL loop
  for (k=0; k<m; k++) {
    //-------------------- generates term of degree k+1
    scal = (k==0? 1.0 : 2.0);
    for (j=0; j<n; j++)
      vkp1[j] = scal*xi[j]*vk[j] - vkm1[j];  
    tmp  = vkm1;
    vkm1 = vk;
    vk   = vkp1;
    vkp1 = tmp;
    //-------------------- accumulation of vector.
    DAXPY(&n, &mu[k+1], vk, &one, yi, &one) ;
  }
  //-------------------- done
  free(vkm1);
  free(vkp1);
  free(vk);
  return 0;
}

void chext(polparams *pol, double aIn, double bIn){
  //	   double *bar, int *mbest)
  /* determines  polynomial for  end interval  cases.  In these cases,
     polynomial is just a scaled Chebyshev polynomial. However we need
     to express it in the same basis as in the other (middle interval)
     cases. This function determines this expansion 
     [aIn, bIn] = transformed interval 
     pol        = struct containing the parameters of polynomial.
     On return these are set in pol:
     /   mu   : expansion coefficients of best polynomial found 
     /   deg  : degree of polymomial
     /   gam  : site of delta function that is expanded.
     /          accurate 'balancing' is done:
     /          if p(t) is best approximation to delta function at gam 
     /           then gam is selected so that p(a)=p(b) - within the 
     /          tolerance tolBal (set in this function to 1.e-10)
     /   bar  : if (p(\lambda_i)) >= bar, accept eigenvalue as belonging
     /          to interval; else reject. 
     / Note : [a b] is now "interval for eigenvalues to damp", 
     [aIn, bIn] is the interval of interest 
     del is used to expand interval of interest slightly so 
     that pn(t)<bar by some margin ------------------------------*/
  int max_deg = pol->max_deg;
  // int min_deg = pol->min_deg;   NOT used 
  double thresh = pol->thresh_ext;
  double *mu = pol->mu;
  //-------------------- local variables 
  double del = 0.1*sqrt((bIn-aIn)*0.5); 
  double eps = 1e-13;
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
  Calloc(t0, work_size, double);
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
  free(t0);
}

// find the indexofSmallestElement in array
int indexofSmallestElement(double *array, int size){
  int index = 0, i;
  for(i = 1; i < size; i++){
    if(array[i] < array[index])
      index = i;              
  }
  return index;
}

int rootchb(int m, double *v, double* jac, double thcIn, double tha, 
    double thb, double *mu, double *thc){
  int MaxIterBalan= 2;    // max steps in Newton to balance interval
  double tolBal=1.e-6;  // tolerance for Newton - based on f val only
  // do 2 newton steps -- of OK exit otherwise
  // continue to get root by solving eigv. pb
  int j, it, mm, one = 1;
  double fval = 0.0, d, *T, *g;
  *thc = thcIn;
  /*-------------------- Newton iteration */
  for (it=0; it<=MaxIterBalan; it++) {
    //-------------------- Balacing the interval
    for (j=0; j<=m; j++) 
      mu[j] = cos(j* (*thc))*jac[j];
    mm = m+1;

    fval = DDOT(&mm,mu,&one,v,&one);
    // printf("     --1 thc %f fval %e \n",*thc,fval);
    if (fabs(fval)<tolBal){
      return(0) ;
    } 
    /*-------------------- do one newton step */
    d = 0.0;
    for (j=1; j<=m; j++)
      d += jac[j]*j*sin(j*(*thc))*v[j];
    *thc += fval/d;
    //--------------------value not good - exit immediatly
    if (( (*thc) < thb) || ( (*thc) > tha)) {
      break;
    }
  }
  //--------------------if acceptable exit. Otherwise get *thc 
  //                    using sol. of eigenvalue problem  
  mm = m+1;
  /*-------------------- T must be zeroed out initially*/
  Calloc(T, m*m, double);
  Malloc(g, mm, double);
  // I store T in a upper hess. form in order to use dhseqr.f 
  //printf("entering non-newton\n");
  //printf(" tha %f thc %f thb %f\n",tha,*thc,thb); 
  for (j=0; j<=m; j++)
    g[j] = v[j]*jac[j];
  for (j=0; j<=m; j++)
    g[j] /= g[m];
  for (j=1; j<m; j++) {
    T[j*m+j-1] = -1;
    T[j*m-m+j] = -1;
  } 
  T[1] = -2;
  for (j=0; j<m; j++) {
    T[(m-1)*m+j] += g[j];
  }   
  /*-------------------- compute eigenvalues of T */
  char jobz='E';/* want eigenvalues only */
  char compz = 'N';/* no Schur vectors */
  int info, ilo = 1, ihi = m, ldz = 1, lwork = m;
  double *evr, *evi, *z=NULL, *work;
  Malloc(evr, m, double);
  Malloc(evi, m, double);
  Malloc(work, lwork, double);
  DHSEQR(&jobz, &compz, &m, &ilo, &ihi, T, &m, evr, evi, 
      z, &ldz, work, &lwork,&info);
  if (info != 0) {
    fprintf(stdout, "DHSEQR error: %d\n", info);
    return(1);
  }
  for (j=0; j<m; j++) 
    evr[j]*= -0.5; 
  //-------------------- root closest to thc
  double target = cos(*thc);
  for (j=0; j<m; j++)
    evi[j] = fabs(evr[j]-target);

  int index = indexofSmallestElement(evi, m);
  double y = evr[index];
  *thc = acos(y);
  if (( (*thc) < thb) || ( (*thc) > tha)) 
    *thc = thcIn;

  for(j=0; j<=m; j++)
    mu[j] = cos(j*(*thc))*jac[j];
  free(evr);
  free(evi);
  free(work);
  free(g);
  free(T);
  return 0;
}

int find_pol(double *intv, polparams *pol){ 
  /*----------------------------------------------------------------------
    / IN:
    / intv  = an array of length 4 
    /         [intv[0], intv[1]] is the interval of desired eigenvalues
    /         [intv[2], intv[3]] is the global interval of all eigenvalues
    /         it must contain all eigenvalues of A
    /  
    / OUT:
    / these are set in pol struct:
    /   mu   : expansion coefficients of best polynomial found 
    /   gam  : site of delta function that is expanded.
    /          accurate 'balancing' is done:
    /          if p(t) is best approximation to delta function at gam 
    /           then gam is selected so that p(a)=p(b) - within the 
    /          tolerance tolBal (set in this function to 1.e-10)
    /   bar  : if (p(\lambda_i)) >= bar, accept eigenvalue as belonging
    /          to interval; else reject. 
    /   deg  : degree of polymomial
    /    dd  : half-width and.. 
    /    cc  : ..center of interval containing all eigenvalues 
    /         [these are used by ChebAv] 
   * --------------------------------------------------------------------*/
  double *mu, *v, *jac, t=0.0, itv[2],  vals[2];
  int max_deg=pol->max_deg, min_deg=pol->min_deg, damping=pol->damping;
  double tha=0.0, thb=0.0, thc=0.0, thcIn=0.0;
  double gam,  thresh;
  int m, j, nitv,  mbest; 
  /*-------------------- A few parameters to be set or reset */
  Malloc(mu, max_deg+1, double);
  pol->mu = mu;
  Malloc(v, max_deg+1, double);
  Malloc(jac, max_deg+1, double);
  /*-------------------- A parameter for interval check */
  // double IntTol = 2*DBL_EPSILON; // If b*IntTol>1 accept [a b] extreme
  double IntTol = 0.0005;
  //-------------------- intervals related
  if (check_intv(intv, stdout) < 0) {
    return -1;
  }

  double aa, bb;
  aa = max(intv[0], intv[2]);  bb = min(intv[1], intv[3]);
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
  //a = max(a, lmin);
  //b = min(b, lmax);
  /*   transform [lmin, lmax] to [-1,1] by y = (x-cc) / dd
   * transform [a, b] to [aT, bT] accordingly */
  aa  = (aa - cc) / dd;
  bb  = (bb - cc) / dd;
  aa  = max(aa, -1.0);
  bb  = min(bb,  1.0);
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
    gam = -1.0;               // set center for a
  } else if (bb + IntTol >= 1.0) {
    /*-------------------- right interval */
    thc = thb;
    nitv   = 1;
    bb = 1;
    gam = 1; // set center for b 
  } else {
    /*-------------------- middle interval */
    thc = 0.5*(tha+thb); 
    itv[0] = aa;
    itv[1] = bb;
    nitv = 2;
  }
  /*-------------------- threshold for finding the best poly */
  if (nitv == 1) { 
    /*-------------------- end interval case done separately */
    chext(pol,aa,bb);
  } else {
    /*-----------------------------------------------*/ 
    min_deg = max(min_deg,2);
    thresh = pol->thresh_int;
    thcIn = thc;
    //-------------------- this is a short-circuit for the 
    // case of a forced degree 
    if (pol->deg > 0){
      thresh = -1;
      max_deg = pol->deg;
    }
    /*-------------------- initialize v vector */
    for (j=0; j<min_deg; j++)
      v[j] = cos(j*thb) - cos(j*tha);
    /*-------------------- DEGREE LOOP -------------------- */
    for (m=min_deg; m < max_deg;m++){
      dampcf(m, damping, jac);
      //-------------------- update v: add one more entry
      v[m] = cos(m*thb)-cos(m*tha);
      //---------------------Balacing the interval + get new mu
      rootchb(m, v, jac, thcIn, tha, thb, mu, &thc);
      //----------------------New center 
      gam = cos(thc);
      //-------------------- for scaling
      chebxPltd(m, mu, 1, &gam, &t);
      chebxPltd(m, mu, nitv, itv, vals);
      //-------------------- test for acceptance of this pol. 
      if (vals[0] <= t*thresh && vals[1] <= t*thresh) {
        m++;
        break;
      }
    }
    mbest = m - 1;
    //-------------------- scale the polynomial
    for (j=0; j<=m; j++) 
      mu[j] /= t;
    pol->bar = min(vals[0], vals[1])/t;
    pol->gam = gam;
    pol->deg = mbest;
  }
  free(v);
  free(jac);
  return 0;
}


int ChebAv(csrMat *A, polparams *pol, double *v, double *y, double *w) {
  /*----------------------------------------------------------------------
    Computes  y = p(A) y, where pn is a Cheb. polynomial expansion 
    [this does not call matvec - but does the sparse matrix vectror product
    internally] 

    INPUT
    csrMat *A  = matrix A
    pol        = struct containing the parameters and expansion
    coefficients of the polynomial. 
    double *v  = input vector

    OUT: 
    double *y,= p(A)v

    WORK SPACE: 
    double *w = work vector of length 3*n [allocate before call] 
    v is untouched. 
----------------------------------------------------------------------*/
  //-------------------- unpack A 
  int n = A->nrows;
  int  *ia = A->ia;
  int  *ja = A->ja;
  double  *a = A->a;
  //-------------------- unpack pol
  double *mu = pol->mu;
  double dd = pol->dd;
  double cc = pol->cc;
  int m = pol->deg;
  double zer = 0.0;
  //-------------------- pointers to v_[k-1],v_[k], v_[k+1]  from w
  double *vk   = w;
  double *vkp1 = w+n;
  double *vkm1 = vkp1+n;
  //-------------------- 
  int k, i, j;
  double r, s, t, *tmp;

  double t1 = 1./dd, t2 = 2.0/dd;
  //-------------------- vk <- v; vkm1 <- zeros(n,1)
  memcpy(vk, v, n*sizeof(double));
  memset(vkm1,zer,n*sizeof(double));
  /*-------------------- special case: k == 0 */
  s = mu[0];
  for (i=0; i<n; i++) 
    y[i] = s*vk[i];
  //-------------------- degree loop. k IS the degree.  
  for (k=1; k<=m; k++) {
    /*-------------------- y = mu[k]*Vk + y */
    t = (k==1 ? t1 : t2);
    /*-------------------- Vkp1 = A*Vk - cc*Vk; */
    s = mu[k]; 
    for (i=0; i<n; i++){
      r = -cc*vk[i];
      for (j=ia[i]; j<ia[i+1];j++) 	
        r += vk[ja[j]]*a[j];
      //-------------------- t* ( A*Vk - cc*Vk) - vkm1
      r = r*t - vkm1[i] ;
      //-------------------- save v_{k+1} and update y
      vkp1[i] = r;
      y[i]+=s*r;
    }   
    //-------------------- next step: rotate vectors via pointer exchange
    tmp = vkm1;
    vkm1 = vk;
    vk = vkp1;
    vkp1 = tmp;
  }
  return 0;
}

int ChebAv0(csrMat *A, polparams *pol, double *v, double *y, double *w) {
  /*----------------------------------------------------------------------
    / Computes  y = p(A) y, where pn is a Cheb. polynomial expansion 
    / This is unused but left here on purpose -- simpler but a bit slower 
    / routine.. This calls explicitly matvec -- so it can be useful for
    / implementing user-specific matrix-vector multiplication 

    INPUT
    csrMat *A  = matrix A
    pol        = struct containing the parameters and expansion
    coefficients of the polynomial. 
    double *v  = input vector

    OUT: 
    double *y,= p(A)v

    WORK SPACE: 
    double *w = work vector of length 3*n [allocate before call] 
    v is untouched. 
  ----------------------------------------------------------------------*/
  //-------------------- unpack n from A  
  int n = A->nrows;
  //-------------------- unpack pol
  double *mu = pol->mu;
  double dd = pol->dd;
  double cc = pol->cc;
  int m = pol->deg;
  //-------------------- pointers to v_[k-1],v_[k], v_[k+1]  from w
  double *vk   = w;
  double *vkp1 = w+n;
  double *vkm1 = vkp1+n;
  //-------------------- 
  int k, i; 
  double t,  s, *tmp, t1= 1./dd, t2 = 2.0/dd; 
  double zer = 0.0;
  //-------------------- vk <- v; vkm1 <- zeros(n,1)
  memcpy(vk, v, n*sizeof(double));
  memset(vkm1,zer,n*sizeof(double));
  /*-------------------- special case: k == 0 */
  s = mu[0] ;
  for (i=0; i<n; i++) 
    y[i] = s*vk[i]; 
  //-------------------- degree loop. k IS the degree.     
  for (k=1; k<=m; k++) {
    /*-------------------- y = mu[k]*Vk + y */    
    t = (k==1 ? t1:t2); 
    /*-------------------- Vkp1 = A*Vk - cc*Vk; */    
    s = mu[k];
    matvec (A, vk, vkp1);
    for (i=0; i<n; i++){
      vkp1[i] = t*(vkp1[i]-cc*vk[i]) - vkm1[i];
      //-------------------- for degree 2 and up: 
      y[i] += s*vkp1[i];
    }
    //-------------------- next: rotate vectors via pointer exchange
    tmp = vkm1;
    vkm1 = vk;
    vk = vkp1;
    vkp1 = tmp;
  }
  return 0;
}

