#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "internal_header.h"

#define NBUF 2
/**
 * @file chebsi.c
 * @brief Polynomial Filtered Subspace Iteration
 */

/**
 * @brief Chebyshev polynomial filtering Subspace Iteration
 *
 *   @param[in] nev         Estimate of number of eigenvalues in the interval --
 *           ideally nev == exact number or a little larger.
 *           ChebSI stops when  at least nev eigenvalues are
 *           found or when no more candidates are left in interval.
 *   @param[in] intv        An array of length 4
 *           [intv[0], intv[1]] is the interval of desired eigenvalues
 *           [intv[2], intv[3]] is the global interval of all eigenvalues
 *           it must contain all eigenvalues of A
 *   @param[in] maxit       Max Num of outer subspace iterations allowed
 *   @param[in] tol         Tolerance for convergence. stop when ||res||< tol
 *   @param[in] vinit       Nev initial vectors [to be made optional]
 *   @param[in] pol         A struct containing the parameters of the polynomial..
 *
 *   @b Modifies:
 * @param[out] nevo     Number of eigenvalues/vectors computed
 * @param[out] Yo       A set of eigenvectors  [n x nev2 matrix]
 * @param[out] lamo     Associated eigenvalues [nev2 x 1 vector]
 * @param[out] reso     Associated residual norms [nev x 1 vector]
 * @param[out] fstats   File stream which stats are printed to
 *
 *   @warning Memory allocation for Yo/lamo/reso within this function
 */

int ChebSI(int nev, double *intv, int maxit,
           double tol, double *vinit, polparams *pol, int *nevo,
           double **lamo, double **Yo, double **reso, FILE *fstats) {
  /*-------------------- for stats */
  double tm, tall = 0.0;
  int icol, nmv = 0;
  double tolP = tol*0.01;
  tall = evsl_timer();
  //    int max_deg = pol->max_deg,   min_deg = pol->min_deg;
  /*-------------------   size of A */
  int n = evsldata.n;
  /*--------------------   some constants frequently used */
  char cT = 'T';
  char cN = 'N';
  int one = 1;
  int nev2 = nev*nev;
  int nnev = n*nev;
  double done=1.0,dmone=-1.0,dzero=0.0;
  if (check_intv(intv, fstats) < 0) {
    *nevo = 0;
    return 0;
  }
  double aa = intv[0];
  double bb = intv[1];
  /*--------------------   frequency of convergence checks (check every cvcheck iterations) */
  int cvcheck = 5;
  /*-------------------   misc */
  int i;
  // this code(step 1)is the same as in cheblan; to be moved
  // to a separate file for reuse
  /*-------------------- unpack some values from pol */
  int deg =pol->deg;
  double gamB=pol->gam, bar=pol->bar;
  fprintf(fstats, "Cheb Poly- deg = %d, gam = %.15e, bar: %.15e\n",
      deg, gamB, bar);
  /*-------------------- gamB must be within [-1, 1] */
  if (gamB > 1.0 || gamB < -1.0) {
    fprintf(stdout, "gamB error %.15e\n", gamB);
    return -1;
  }
  /*-------------------- filtered subspace iteration */
  fprintf(fstats, "Step 2: ChebSI, block size = %d\n", nev);
  /*-------------------- it = number of the current iteration */
  int it = 0;
  /*-------------------- memory for projected matrix T=V'p(A)V and its eigenpairs (evecT,evalT) */
  double *T, *T_host, *evecT, *evecT_device, *evalT;
  T = evsl_Malloc_device(nev2, double);
  evecT = evsl_Malloc(nev2, double);
#ifdef EVSL_USING_CUDA_GPU
  T_host = evsl_Malloc(nev2, double);
  evecT_device = evsl_Malloc_device(nev2, double);
#else
  T_host = T;
  evecT_device = evecT;
#endif
  evalT = evsl_Malloc(nev, double);
  ///*-------------------- memory for converged Ritz pairs (Y, Lam) */
  //    double *Y, *Lam;
  //    Y = evsl_Malloc(nev2, double);
  //    Lam = evsl_Malloc(nev2, double);

  /*-------------------- memory for block of approx. eigenvectors/evals and their products with p(A)*/
  double *V, *Lam, *PV, *R, *res;
  double *V_out, *Lam_out, *res_out;
  V = evsl_Malloc_device(nnev, double);
  PV = evsl_Malloc_device(nnev, double);
  Lam = evsl_Malloc(nev, double);
  R = evsl_Malloc_device(nnev, double);   // block of residuals
  res = evsl_Malloc(nev, double);  // residual norms w.r.t. A

  /*-------------------- number of locked, and unconverged (active) pairs */
  int nlock, nact, nlock_ab;
  int *act_idx, *lock_idx;  // arrays used to store idices of active and locked columns of V
  act_idx = evsl_Malloc(nev, int);
  lock_idx = evsl_Malloc(nev, int);
  nlock = 0;
  nact  = nev;
  nlock_ab = 0;  // number of locked eigenvalues that are in [a,b]

  /*-------------------- alloc some work space */
  double *work, *buf;
  int work_size = 3*n;
  work = evsl_Malloc_device(work_size, double);
  buf = evsl_Malloc_device(nnev, double);  // buffer for evsl_dgemm calls

  /*-------------------- orthonormalize initial V and compute PV = P(A)V */
  orth(vinit,n,nev,V,work);

  for (icol = 0; icol<nev; icol++) {
    ChebAv(pol, V+icol*n, PV+icol*n, work);
  }
  nmv += deg*nev;

  int find_more = 1;
  /*-------------------- Filtered Subspace itertion: MAIN LOOP */
  while ( (it < maxit) && (find_more) ) {
    /*---  V <- orth(PV)  */
    int nnlock = n*nlock;                                      // pointer to the first active column
    int nnact  = n*nact;
    evsl_dcopy_device(&nnact, PV+nnlock, &one, V+nnlock, &one);
#if EVSL_TIMING_LEVEL > 0
    tm = evsl_timer();
#endif
    /*---  Orthogonalize active columns V(:,nlock:nev-1) against locked columns V(:,0:nlocked-1)*/
    if ( nlock>0 ) {
      evsl_dgemm_device(&cT, &cN, &nlock, &nact, &n, &done, V, &n, V+nnlock, &n, &dzero, T, &nev);
      evsl_dgemm_device(&cN, &cN, &n, &nact, &nlock, &dmone, V, &n, T, &nev, &done, V+nnlock, &n);
    }
#if EVSL_TIMING_LEVEL > 0
    evslstat.t_siorth = evsl_timer() - tm;
#endif
    /*--- Orthogonormalize columns of V(:,nlock:nev-1) */
    orth(V+nnlock, n, nact, buf, work);
    evsl_dcopy_device(&nnact, buf, &one, V+nnlock, &one);
    /*---  PV <- P(A)*V */
    //polyfilt(A, deg, mu, dd, cc, V+nnlock, nact, PV+nnlock, work);
    for (icol = nlock; icol<nlock+nact; icol++) {
      ChebAv(pol, V+icol*n, PV+icol*n, work);
    }
    nmv += deg*nact;
#if EVSL_TIMING_LEVEL > 0
    tm = evsl_timer();
#endif
    // orthogonalize PVact against Vlocked
    if ( nlock>0 ) {
      evsl_dgemm_device(&cT, &cN, &nlock, &nact, &n, &done, V, &n, PV+nnlock, &n, &dzero, T, &nev);
      evsl_dgemm_device(&cN, &cN, &n, &nact, &nlock, &dmone, V, &n, T, &nev, &done, PV+nnlock, &n);
    }
#if EVSL_TIMING_LEVEL > 0
    evslstat.t_siorth = evsl_timer() - tm;
#endif
    // end orthogonalize PVact against Vlocked
    /*---  Lock converged pairs */
    if ( ((it+1)%cvcheck == 0) || ((it+1)==maxit)  ) {
      //// Rayleigh--Ritz with unconverged columns V
      /*---  T = V'p(A)V; [evecT,evalT]=eig(T);    */
#if EVSL_TIMING_LEVEL > 0
      tm = evsl_timer();
#endif
      evsl_dgemm_device(&cT,&cN,&nact,&nact,&n,&done,V+nnlock,&n,PV+nnlock,&n,&dzero,T,&nev);

#if EVSL_TIMING_LEVEL > 0
      evslstat.t_ritz = evsl_timer() - tm;
#endif

#ifdef EVSL_USING_CUDA_GPU
      evsl_memcpy_device_to_host(T_host, T, nev2*sizeof(double));
#endif

      SymEigenSolver(nact, T_host, nev, evecT, nev, evalT+nlock);

#ifdef EVSL_USING_CUDA_GPU
      evsl_memcpy_host_to_device(evecT_device, evecT, nev2*sizeof(double));
#endif

#if EVSL_TIMING_LEVEL > 0
      tm = evsl_timer();
#endif
      /*---  V <- V*evecT; p(A)V <- PV*evecT (update only active columns)    */
      evsl_dgemm_device(&cN,&cN,&n,&nact,&nact,&done,V+nnlock,&n,evecT_device,&nev,&dzero,buf,&n);
      evsl_dcopy_device(&nnact, buf, &one, V+nnlock, &one);
      evsl_dgemm_device(&cN,&cN,&n,&nact,&nact,&done,PV+nnlock,&n,evecT_device,&nev,&dzero,buf,&n);
      evsl_dcopy_device(&nnact, buf, &one, PV+nnlock, &one);

#if EVSL_TIMING_LEVEL > 0
      evslstat.t_ritz = evsl_timer() - tm;
#endif

      /*---  Compute active residuals R = PV - V*diag(evalT)    */
      for (i=nlock; i<nev; i++) {
        double t = -evalT[i];
        evsl_dcopy_device(&n, PV+i*n, &one, R+i*n, &one);
        evsl_daxpy_device(&n, &t, V+i*n, &one, R+i*n, &one);
        //for (j=0; j<n; j++) {
        //  R[i*n+j] = PV[i*n+j]+t*V[i*n+j];
        //}
      }
      /*---  Detect converged eigenpairs w.r.t. p(A) and A (candidates). Move them to first nlock columns    */
      int nlock_new = 0;    // number of newly locked pairs
      int nlock_ab_new = 0;
      nact = 0;
      for (i=nlock; i<nev; i++) {
        /*---  Compute norm of R(:,i)   */
        double resP = sqrt(evsl_ddot_device(&n, R+i*n, &one, R+i*n, &one));
        if (resP < tolP) {
          /*---  Compute norm of AV(:,i) - V(:,i)*Lambda(i)   */
          matvec_A(V+i*n, buf);
          nmv++;
          double rq = evsl_ddot_device(&n, V+i*n, &one, buf, &one);  // Rayleigh Quotient for A
          double nrq = -rq;
          evsl_dcopy_device(&n, buf, &one, R+i*n, &one);
          evsl_daxpy_device(&n, &nrq, V+i*n, &one, R+i*n, &one);
          //for (j=0; j < n; j++) {
          //  R[i*n+j] = buf[j] - rq*V[i*n+j];
          //}
          double resA = sqrt(evsl_ddot_device(&n, R+i*n, &one, R+i*n, &one));
          if (resA < tol) {
            lock_idx[nlock_new] = i;
            res[nlock+nlock_new] = resA;
            Lam[nlock+nlock_new] = rq;
            nlock_new++;
            /*---  Determine if the newly locked eigenvalue is in [a,b] */
            if ( rq >= aa - EVSL_DBL_EPS_MULT * DBL_EPSILON && rq <= bb + EVSL_DBL_EPS_MULT * DBL_EPSILON ) {
              nlock_ab_new++;
            }
          } else {
            act_idx[nact] = i;
            nact++;
          }
        } else {
          act_idx[nact] = i;
          nact++;
        }
      }

      /*---  Move newly locked eigenvectors to  columns nlock:nlock+nlock_new-1 of V */
      for (i = 0; i<nlock_new; i++) {
        evsl_dcopy_device(&n, V+lock_idx[i]*n, &one, buf+nnlock+i*n, &one);
        evsl_dcopy_device(&n, PV+lock_idx[i]*n, &one, R+nnlock+i*n, &one);
      }
      /*---  Move active columns to V(:, nlock:nev)*/
      for (i = 0; i<nact; i++) {
        evsl_dcopy_device(&n, V+act_idx[i]*n, &one, buf+nnlock+(nlock_new+i)*n, &one);
        evsl_dcopy_device(&n, PV+act_idx[i]*n, &one, R+nnlock+(nlock_new+i)*n, &one);
      }
      evsl_dcopy_device(&nnact, buf+nnlock, &one, V+nnlock, &one);
      evsl_dcopy_device(&nnact, R+nnlock, &one, PV+nnlock, &one);

      nlock += nlock_new;
      nlock_ab += nlock_ab_new;

      fprintf(fstats, "it %4d:   nMV %7d,  nlock %3d, nlock_ab  %3d\n",
          it+1, nmv, nlock, nlock_ab);

      /*---  Decide if iteration should be stopped */
      /*     Stop if number of locked pairs that are in [a,b] are by NBUF smaller than total
             number of locked pairs or if all nev pairs in the block have been locked    */
      if ( ((nlock-nlock_ab)>=NBUF) || (nlock == nev) ) {
        find_more =0;
        fprintf(fstats, "-------------------------------------\n");
        fprintf(fstats, " No ev.s left to be computed\n");
        fprintf(fstats, " Number of evals found = %d\n", nlock_ab);
        fprintf(fstats, "-------------------------------------------------------------------\n");
      }
    }
    it++;
  }

  /*-------------------- Postprocessing : remove eigenpairs not associated with [a,b]*/
  Lam_out = evsl_Malloc(nev, double);
  res_out = evsl_Malloc(nev, double);
  V_out = evsl_Malloc_device(nnev, double);
  int idx=0;
  for (i=0; i<nlock;i++) {
    double t = Lam[i];
    if ( t >= aa - EVSL_DBL_EPS_MULT * DBL_EPSILON && t <= bb + EVSL_DBL_EPS_MULT * DBL_EPSILON ) {
      Lam_out[idx] = t;
      res_out[idx] = res[i];
      evsl_dcopy_device(&n, V+i*n, &one, V_out+idx*n, &one);
      idx++;
    }
  }

  /*-------------------- Done.  output : */
  *nevo = idx;
  *lamo = Lam_out;
  *Yo = V_out;
  *reso = res_out;
  /*-------------------- free arrays */
  evsl_Free_device(T);
  evsl_Free(evalT);
  evsl_Free(evecT);
#ifdef EVSL_USING_CUDA_GPU
  evsl_Free_device(evecT_device);
  evsl_Free(T_host);
#endif
  evsl_Free_device(R);
  evsl_Free_device(V);
  evsl_Free_device(PV);
  evsl_Free(Lam);
  evsl_Free(res);
  evsl_Free(act_idx);
  evsl_Free(lock_idx);
  evsl_Free_device(buf);
  evsl_Free_device(work);

  /*-------------------- record stats */
  tall = evsl_timer() - tall;
  /*-------------------- print stat */
  fprintf(fstats, "------This slice consumed: \n");
  fprintf(fstats, "Matvecs :        %d\n", nmv);
  fprintf(fstats, "total  time :    %.2f\n", tall);

  evslstat.t_iter = tall;

  return 0;
}

