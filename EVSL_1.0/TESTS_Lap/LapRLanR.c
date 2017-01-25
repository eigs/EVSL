#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "evsl.h"
#include "io.h"
#include "cholmod.h"
#include "umfpack.h"

int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv);
void solvefunc(int n, double *br, double *bz, double *xr, double *xz, void *data);

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/*------------ data needed by umfpack for solving a linear system */
typedef struct _umfdata {
  SuiteSparse_long *Ap;
  SuiteSparse_long *Ai;
  double *Ax;
  double *Az;
  void *Numeric;
} umfdata;

int main(int argc, char *argv[]) {
  /*------------------------------------------------------------
    generates a laplacean matrix on an nx x ny x nz mesh 
    and computes all eigenvalues in a given interval [a  b]
    The default set values are
    nx = 41; ny = 53; nz = 1;
    a = 0.4; b = 0.8;
    nslices = 1 [one slice only] 
    other parameters 
    tol [tolerance for stopping - based on residual]
    Mdeg = pol. degree used for DOS
    nvec  = number of sample vectors used for DOS 
    This uses:
    Thick-restart Lanczos with rational filtering
    ------------------------------------------------------------*/
  int n, nnz,nx, ny, nz, i, j, k, npts, nslices, nvec, Mdeg, nev, 
      status, mlan, max_its, ev_int, sl, flg, ierr;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol,   *sli, *mu;
  double xintv[4];
  /* initial vector: random */
  double *vinit;
  /* Input for umfpack*/
  SuiteSparse_long *Ap, *Ai;
  int *diag; 
  /* parameters for rational filter */
  int num = 1; // number of poles used for each slice
  int pw = 2; // multiplicity of each pole
  double beta = 0.01; // beta in the LS approximation
  FILE *fstats = NULL;
  if (!(fstats = fopen("OUT/LapRLanR","w"))) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
  /*-------------------- default values */
  nx   = 41;
  ny   = 53;
  nz   = 1;
  a    = 0.4;
  b    = 0.8;
  nslices = 4;
  //-----------------------------------------------------------------------
  //-------------------- reset some default values from command line [Yuanzhe/]
  /* user input from command line */
  flg = findarg("help", NA, NULL, argc, argv);
  if (flg) {
    printf("Usage: ./testL.ex -nx [int] -ny [int] -nz [int] -a [double] -b [double] -nslices [int]\n");
    return 0;
  }
  findarg("nx", INT, &nx, argc, argv);
  findarg("ny", INT, &ny, argc, argv);
  findarg("nz", INT, &nz, argc, argv);
  findarg("a", DOUBLE, &a, argc, argv);
  findarg("b", DOUBLE, &b, argc, argv);
  findarg("nslices", INT, &nslices, argc, argv);
  fprintf(fstats,"used nx = %3d ny = %3d nz = %3d",nx,ny,nz);
  fprintf(fstats," [ a = %4.2f  b= %4.2f],  nslices=%2d \n",a,b,nslices);
  //----------------------------------------------------------------------- DONE
  //-------------------- eigenvalue bounds set by hand.
  lmin = 0.0;  
  lmax =  ((nz == 1)? 8.0 :12) ;
  xintv[0] = a;
  xintv[1] = b;
  xintv[2] = lmin;
  xintv[3] = lmax;
  tol = 1e-8;
  /*-------------------- generate 2D/3D Laplacian matrix 
   *                     saved in coo format */
  n = nx * ny * nz;
  ierr = lapgen(nx, ny, nz, &Acoo);  
  /*-------------------- convert coo to csr */
  ierr = cooMat_to_csrMat(0, &Acoo, &Acsr); 
  nnz = Acoo.nnz;
  /* output the problem settings */
  fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(fstats, "Laplacian: %d x %d x %d, n = %d, nnz = %d\n", nx, ny, nz, n, Acoo.nnz);
  fprintf(fstats, "Interval: [%20.15f, %20.15f]  -- %d slices \n", a, b, nslices);
  fprintf(fstats, "Step 0: Eigenvalue bound s for A: [%.15e, %.15e]\n", lmin, lmax);
  
  // format convertion for umfpack
  solveShift solshift;
  diag = (int *)malloc(n*sizeof(int));
  Ap = (SuiteSparse_long *)malloc((n+1)*sizeof(SuiteSparse_long));    // pointer to starting point of each row
  Ai = (SuiteSparse_long *)malloc(nnz*sizeof(SuiteSparse_long));  // column indices
  // each shifted matrix shares the same Ap and Ai
  for (i=0; i<n+1; i++) {
    Ap[i] = Acsr.ia[i];
  }
  for (i=0; i<nnz; i++) {
    Ai[i] = Acsr.ja[i];
  }
  //-------------------- find the location of the diagonal entries in CSR format n_intv 
  //                     shifted matrices share the same diag, Ax and Az
  for (i = 0 ; i < n ; i++){
    for (j = Ap[i] ; j < Ap[i+1] ; j++){      
      k = Ai[j];
      if(i == k){
        diag[i] = j;
        break;
      }
    }
  }

  /*-------------------- call kpmdos to get the DOS for dividing the spectrum*/
  /*-------------------- define kpmdos parameters */
  Mdeg = 100;
  nvec = 60;
  mu = malloc((Mdeg+1)*sizeof(double));
  //-------------------- call kpmdos 
  ierr = kpmdos(&Acsr, Mdeg, 1, nvec, xintv, mu, &ecount);
  if (ierr) {
    printf("kpmdos error %d\n", ierr);
    return 1;
  }
  fprintf(fstats, " estimated eig count in interval: %10.2e \n",ecount);
  //-------------------- call splicer to slice the spectrum
  npts = 10 * ecount; 
  sli = malloc((nslices+1)*sizeof(double));
  ierr = spslicer(sli, mu, Mdeg, xintv, nslices,  npts);  
  if (ierr) {
    printf("spslicer error %d\n", ierr);
    return 1;
  }
  printf("====================  SLICES FOUND  ====================\n");
  for (j=0; j<nslices;j++)
    printf(" %2d: [% .15e , % .15e]\n", j+1, sli[j],sli[j+1]);
  //-------------------- # eigs per slice
  ev_int = (int) (1 + ecount / ((double) nslices));

  //-------------------- initial vector  
  vinit = (double*) malloc(n*sizeof(double));
  rand_double(n, vinit);
  //-------------------- debug only :
  //  save_vec(n, vinit, "OUT/vinit.mtx");
  
  //-------------------- For each slice call RatLanrTr
  for (sl=0; sl<nslices; sl++){
    double **Ax, **Az;
    Ax = (double **)malloc(num*sizeof(double *));//Ax[i] = real part of the ith shifted matrix
    Az = (double **)malloc(num*sizeof(double *));//Az[i]: imaginary part of the ith shfited matrix
    // Allocate memory for the ith shifted matrix  A - zk[i]I
    for (i=0; i< num; i++){
      Ax[i] = (double *)malloc(nnz*sizeof(double));//real part of each shifted matrix entries
      Az[i] = (double *)malloc(nnz*sizeof(double));//imaginary part of each shifted matrix entries
    }
    for(i=0; i<num; i++){
      for(j=0; j<nnz; j++){
	Az[i][j] = 0.0;
      }
    }
    printf("======================================================\n");
    int nev2;
    double *lam, *Y, *res;
    int *ind;
    int nev_ex;
    double *lam_ex;
    //-------------------- 
    a = sli[sl];
    b = sli[sl+1];
    printf(" subinterval: [%.4e , %.4e]\n", a, b); 
    double intv[4];
    intv[0] = a;
    intv[1] = b;
    intv[2] = lmin;
    intv[3] = lmax;
    // Find the rational filter on this slice
    ratparams rat;
    // Set up default parameters for rat
    set_ratf_def(&rat);
    // change some default values here:
    rat.num = num;
    rat.beta = beta;
    rat.pw = pw;
    // Now determine rational filter
    find_ratf(intv, &rat);    
    void **Symbolic, **Numeric;
    Symbolic= (void**)malloc(num*sizeof(void *));
    Numeric = (void**)malloc(num*sizeof(void *));
    umfdata *umf = (umfdata *) malloc(num*sizeof(umfdata));
    
    for (i= 0; i< num; i++){
      //------------ Form the ith shifted matrix A_shifted = A - zk[i]*I
      memcpy(Ax[i],Acsr.a,nnz*sizeof(double));
      complex double *zk = rat.zk;
      double zkr = creal(zk[i]);
      double zkc = cimag(zk[i]);
      for(j = 0 ; j <n ; j++){
	Ax[i][diag[j]] = Ax[i][diag[j]] - zkr;
	Az[i][diag[j]] = - zkc;
      }
      //Symbolic Factorization
      status = umfpack_zl_symbolic (n, n, Ap, Ai, Ax[i], Az[i], &Symbolic[i], NULL, NULL);
      if (status < 0) {
	printf("umfpack_zl_symbolic failed, %d\n", status);
	return 1;
      }
      //Numerical Factorization
      status = umfpack_zl_numeric(Ap, Ai, Ax[i], Az[i], Symbolic[i], &Numeric[i], NULL, NULL);
      if (status < 0) {
	printf("umfpack_zl_numeric failed and exit, %d\n",status);
	return 1;
      }
      //Release the memory
      umfpack_zl_free_symbolic(&Symbolic[i]);
        
      umf[i].Ap = Ap;  umf[i].Ai = Ai;  umf[i].Ax = Ax[i];  umf[i].Az = Az[i];  umf[i].Numeric = Numeric[i];
    }
    //-------------------- approximate number of eigenvalues wanted
    nev = ev_int+2;
    //-------------------- Dimension of Krylov subspace and maximal iterations
    mlan = max(4*nev, 100);  mlan = min(mlan, n); 
    max_its = 3*mlan;
    // structure for direct solver passed to RatLanNr
    // func [function pointer]: function to solve complex shifted systems
    // it must be a function of the following type:
    // void (*solveShiftFunc)(int n, double *br, double *bz,
    //                        double *xr, double *xz, void *data);
    // void *data: encapsulates all data needed by the direct solver
    // An example is given here for using UMFPACK by T.Davis
    solshift.func = &solvefunc;
    solshift.data = malloc(num*sizeof(void *));
    for (i=0; i<num; i++) {
      solshift.data[i] = &umf[i];
    }
    //-------------------- RationalLanTr
    ierr = RatLanTr(&Acsr, &solshift, mlan, nev, intv, &rat, max_its, tol, vinit, &nev2, &lam, &Y, &res, fstats);
    if (ierr) {
      printf("RatLanTr error %d\n", ierr);
      return 1;
    }
    /* [compute residual] already computed in res */
    /* sort the eigenvals: ascending order
     * ind: keep the orginal indices */
    ind = (int *) malloc(nev2*sizeof(int));
    sort_double(nev2, lam, ind);
    /* compute exact eigenvalues */
    exeiglap3(nx, ny, nz, a, b, &nev_ex, &lam_ex);
    printf(" number of eigenvalues: %d, found: %d\n", nev_ex, nev2);
    /* print eigenvalues */
    fprintf(fstats, "                                   Eigenvalues in [a, b]\n");
    fprintf(fstats, "    Computed [%d]       ||Res||              Exact [%d]", nev2, nev_ex);
    if (nev2 == nev_ex) {
      fprintf(fstats, "                 Err");
    }
    fprintf(fstats, "\n");
    for (i=0; i<max(nev2, nev_ex); i++) {
      if (i < nev2) {
        fprintf(fstats, "% .15e  %.1e", lam[i], res[ind[i]]);
      } else {
        fprintf(fstats, "                               ");
      }
      if (i < nev_ex) { 
        fprintf(fstats, "        % .15e", lam_ex[i]);
      }
      if (nev2 == nev_ex) {
        fprintf(fstats, "        % .1e", lam[i]-lam_ex[i]);
      }
      fprintf(fstats,"\n");
      if (i>50) {
        fprintf(fstats,"                        -- More not shown --\n");
        break;
      } 
    }
    //-------------------- free allocated space withing this scope
    if (lam) free(lam);
    if (Y) free(Y);
    if (res) free(res);
    for(i=0; i<num; i++){
      free(Ax[i]);
      free(Az[i]);
      umfpack_zl_free_numeric (&Numeric[i]);
    }
    free(Ax);
    free(Az);
    free(Numeric);
    free(Symbolic);
    free(ind);
    free(lam_ex);
    free(rat.omega);
    free(rat.mulp);
    free(rat.zk);
    free(umf);
    free(solshift.data);    
  }
  //-------------------- free other allocated space 
  free(vinit);
  free(sli);
  free(Ap);
  free(diag);
  free(Ai);
  free_coo(&Acoo);
  free_csr(&Acsr);
  free(mu);
  fclose(fstats);
  //
  return 0;
}

void solvefunc(int n, double *br, double *bz, double *xr, double *xz, void *data) {
  /*-------------------------------------------------------------------------------
   * complex linear solver routine passed to evsl
   * NOTE: This function MUST be of this prototype
   * INPUT:
   *   n: size of the system
   *   br, bz: vectors of length n, complex right-hand side (real and imaginary)
   *   data: all data that are needed for solving the system
   * OUTPUT:
   *   xr, xz: vectors of length n, complex solution (real and imaginary)
   *------------------------------------------------------------------------------*/
  umfdata *umf = (umfdata *) data;
  SuiteSparse_long* Ap = umf->Ap;
  SuiteSparse_long* Ai = umf->Ai;
  double* Ax = umf->Ax;
  double* Az = umf->Az;  
  void* Numeric = umf->Numeric;
  double Control[UMFPACK_CONTROL]; 
  umfpack_zl_defaults(Control);
  Control[UMFPACK_IRSTEP] = 0; // no iterative refinement for umfpack 
  umfpack_zl_solve(UMFPACK_A, Ap, Ai, Ax, Az, xr, xz, br, bz, Numeric, Control, NULL); 
}

