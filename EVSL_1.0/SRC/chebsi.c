#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

#define NBUF 2

/**
 * @brief Chebyshev polynomial filtering Subspace Iteration
 *
 *   @param nev         Estimate of number of eigenvalues in the interval --
 *           ideally nev == exact number or a little larger.
 *           ChebSI stops when  at least nev eigenvalues are              
 *           found or when no more candidates are left in interval.
 *   @param intv        An array of length 4 
 *           [intv[0], intv[1]] is the interval of desired eigenvalues
 *           [intv[2], intv[3]] is the global interval of all eigenvalues
 *           it must contain all eigenvalues of A
 *   @param maxit       Max Num of outer subspace iterations allowed 
 *   @param tol         Tolerance for convergence. stop when ||res||< tol
 *   @param vinit       Nev initial vectors [to be made optional]
 *   @param pol         A struct containing the parameters of the polynomial..
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
  double tm, tall=0.0, tmv=0.0;
  int icol, nmv = 0;
  double tolP = tol*0.01;
  tall = cheblan_timer();
  //    int max_deg = pol->max_deg,   min_deg = pol->min_deg;
  /*-------------------   size of A */
  int n;
  /* if users provided their own matvec function, input matrix A will be ignored */
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
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
  int i,j;
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
  double *T, *evecT, *evalT;
  Malloc(T, nev2, double);
  Malloc(evecT, nev2, double);
  Malloc(evalT, nev, double);
  ///*-------------------- memory for converged Ritz pairs (Y, Lam) */
  //    double *Y, *Lam;
  //    Malloc(Y, nev2, double);
  //    Malloc(Lam, nev2, double);

  /*-------------------- memory for block of approx. eigenvectors/evals and their products with p(A)*/
  double *V, *Lam, *PV, *R, *res; 
  double *V_out, *Lam_out, *res_out; 
  Malloc(V, nnev, double);
  Malloc(PV, nnev, double);
  Malloc(Lam, nev, double);
  Malloc(R, nnev, double);   // block of residuals 
  Malloc(res, nev, double);  // residual norms w.r.t. A

  /*-------------------- number of locked, and unconverged (active) pairs */
  int nlock, nact, nlock_ab;  
  int *act_idx, *lock_idx;  // arrays used to store idices of active and locked columns of V 
  Malloc(act_idx, nev, int);
  Malloc(lock_idx, nev, int);
  nlock = 0;
  nact  = nev;    
  nlock_ab = 0;  // number of locked eigenvalues that are in [a,b]

  /*-------------------- alloc some work space */
  double *work, *buf;
  int work_size = 3*n;
  Malloc(work, work_size, double);
  Malloc(buf, nnev, double);  // buffer for DGEMM calls

  /*-------------------- orthonormalize initial V and compute PV = P(A)V */
  orth(vinit,n,nev,V,work);

  tm = cheblan_timer();
  for (icol = 0; icol<nev; icol++) {
    ChebAv(pol, V+icol*n, PV+icol*n, work);
  }     
  tmv += cheblan_timer() - tm;
  nmv += deg*nev;

  int find_more = 1;
  /*-------------------- Filtered Subspace itertion: MAIN LOOP */
  while ( (it < maxit) && (find_more) ) {

    /*---  V <- orth(PV)  */ 
    int nnlock = n*nlock;                                      // pointer to the first active column
    int nnact  = n*nact;
    DCOPY(&nnact, PV+nnlock, &one, V+nnlock, &one);
    /*---  Orthogonalize active columns V(:,nlock:nev-1) against locked columns V(:,0:nlocked-1)*/ 
    if ( nlock>0 ) { 
      DGEMM(&cT, &cN, &nlock, &nact, &n, &done, V, &n, V+nnlock, &n, &dzero, T, &nev);
      DGEMM(&cN, &cN, &n, &nact, &nlock, &dmone, V, &n, T, &nev, &done, V+nnlock, &n);
    }
    /*--- Orthogonormalize columns of V(:,nlock:nev-1) */ 
    orth(V+nnlock, n, nact, buf, work);
    DCOPY(&nnact, buf, &one, V+nnlock, &one);
    /*---  PV <- P(A)*V */
    tm = cheblan_timer();
    //polyfilt(A, deg, mu, dd, cc, V+nnlock, nact, PV+nnlock, work); 
    for (icol = nlock; icol<nlock+nact; icol++)
      ChebAv(pol, V+icol*n, PV+icol*n, work);
    tmv += cheblan_timer() - tm;
    nmv += deg*nact;
    // orthogonalize PVact against Vlocked
    if ( nlock>0 ) { 
      DGEMM(&cT, &cN, &nlock, &nact, &n, &done, V, &n, PV+nnlock, &n, &dzero, T, &nev);
      DGEMM(&cN, &cN, &n, &nact, &nlock, &dmone, V, &n, T, &nev, &done, PV+nnlock, &n);
    }
    // end orthogonalize PVact against Vlocked
    /*---  Lock converged pairs */
    if ( ((it+1)%cvcheck == 0) || ((it+1)==maxit)  ) {
      //// Rayleigh--Ritz with unconverged columns V
      /*---  T = V'p(A)V; [evecT,evalT]=eig(T);    */
      DGEMM(&cT,&cN,&nact,&nact,&n,&done,V+nnlock,&n,PV+nnlock,&n,&dzero,T,&nev);     
      SymEigenSolver(nact, T, nev, evecT, nev, evalT+nlock);
      /*---  V <- V*evecT; p(A)V <- PV*evecT (update only active columns)    */ 
      DGEMM(&cN,&cN,&n,&nact,&nact,&done,V+nnlock,&n,evecT,&nev,&dzero,buf,&n);
      DCOPY(&nnact, buf, &one, V+nnlock, &one);
      DGEMM(&cN,&cN,&n,&nact,&nact,&done,PV+nnlock,&n,evecT,&nev,&dzero,buf,&n);
      DCOPY(&nnact, buf, &one, PV+nnlock, &one);
      /*---  Compute active residuals R = PV - V*diag(evalT)    */
      for (i=nlock; i<nev; i++) {
        double t = -evalT[i]; 
        //DSCAL(&n, &t, R+i*n, &one); 
        for (j=0; j<n; j++) {
          R[i*n+j] = PV[i*n+j]+t*V[i*n+j];    
        }
      }
      /*---  Detect converged eigenpairs w.r.t. p(A) and A (candidates). Move them to first nlock columns    */
      int nlock_new = 0;    // number of newly locked pairs 
      int nlock_ab_new = 0; 
      nact = 0;
      for (i=nlock; i<nev; i++) {
        /*---  Compute norm of R(:,i)   */
        double resP = sqrt(DDOT(&n, R+i*n, &one, R+i*n, &one));
        if (resP < tolP) {
          /*---  Compute norm of AV(:,i) - V(:,i)*Lambda(i)   */
          tm = cheblan_timer();
          matvec_A(V+i*n, buf);
          tmv += cheblan_timer() - tm;
          nmv++;
          double rq = DDOT(&n, V+i*n, &one, buf, &one);  // Rayleigh Quotient for A
          for (j=0; j < n; j++) {
            R[i*n+j] = buf[j] - rq*V[i*n+j]; 
          }
          double resA = sqrt(DDOT(&n, R+i*n, &one, R+i*n, &one));
          if (resA < tol) {
            lock_idx[nlock_new] = i;
            res[nlock+nlock_new] = resA;
            Lam[nlock+nlock_new] = rq;
            nlock_new++;
            /*---  Determine if the newly locked eigenvalue is in [a,b] */
            if ( rq >= aa - DBL_EPSILON && rq <= bb + DBL_EPSILON ) {
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
        DCOPY(&n, V+lock_idx[i]*n, &one, buf+nnlock+i*n, &one);
        DCOPY(&n, PV+lock_idx[i]*n, &one, R+nnlock+i*n, &one);
      }    
      /*---  Move active columns to V(:, nlock:nev)*/
      for (i = 0; i<nact; i++) {
        DCOPY(&n, V+act_idx[i]*n, &one, buf+nnlock+(nlock_new+i)*n, &one);
        DCOPY(&n, PV+act_idx[i]*n, &one, R+nnlock+(nlock_new+i)*n, &one);
      }
      DCOPY(&nnact, buf+nnlock, &one, V+nnlock, &one);
      DCOPY(&nnact, R+nnlock, &one, PV+nnlock, &one);

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
  Malloc(Lam_out, nev, double); 
  Malloc(res_out, nev, double); 
  Malloc(V_out, nnev, double); 
  int idx=0;
  for (i=0; i<nlock;i++) {
    double t = Lam[i];
    if ( t >= aa - DBL_EPSILON && t <= bb + DBL_EPSILON ) {
      Lam_out[idx] = t;
      res_out[idx] = res[i];
      DCOPY(&n, V+i*n, &one, V_out+idx*n, &one);
      idx++;  
    }   
  } 
 
  /*-------------------- Done.  output : */
  *nevo = idx;
  *lamo = Lam_out;
  *Yo = V_out;
  *reso = res_out;  
  /*-------------------- free arrays */
  free(T);
  free(evalT);
  free(evecT);
  free(R);
  free(V);
  free(PV);
  free(Lam);
  free(res);
  free(act_idx);
  free(lock_idx);
  free(buf);
  free(work);
  /*-------------------- record stats */
  tall = cheblan_timer() - tall;
  /*-------------------- print stat */
  fprintf(fstats, "------This slice consumed: \n");
  fprintf(fstats, "Matvecs :        %d\n", nmv);
  fprintf(fstats, "total  time :    %.2f\n", tall);
  fprintf(fstats, "matvec time :    %.2f\n", tmv);
  fprintf(fstats,"======================================================\n");
  return 0;

}

