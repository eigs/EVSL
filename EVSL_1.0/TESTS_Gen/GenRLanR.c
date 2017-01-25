#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <complex.h>
#include "evsl.h"
#include "io.h"
#include "cholmod.h"
#include "umfpack.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define TRIV_SLICER 0 

/*------------ data needed by umfpack for solving a linear system */
typedef struct _umfdata {
  SuiteSparse_long *Ap;
  SuiteSparse_long *Ai;
  double *Ax;
  double *Az;
  void *Numeric;
} umfdata;

/*-------------------- Protos */
int read_coo_MM(const char *matfile, int idxin, int idxout,   cooMat *Acoo); 
int get_matrix_info(FILE *fmat, io_t *pio);
void solvefunc(int n, double *br, double *bz, double *xr, double *xz, void *data);
/*-------------------- End Protos */

int main () { 
  int ierr = 0;
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix
   * read in matrix format -- using
   * non-restart Lanczos with rational filtering.
   *-------------------------------------------------------------*/
  int n=0, nnz=0, sl, i, j, k, mlan, nev, totcnt; 
  double a, b, ecount, xintv[4];
  double lmin, lmax; 
  double *alleigs; 
  int n_intv;      // number of subintervals (slices)  
  int npnts;       // number of integration points for eigenvalue count 
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr; 
  /*-------------------- Pass pointers from Acsr to Ap, Ai, Ax, Az for umfpack*/
  SuiteSparse_long *Ap, *Ai;
  int *diag; 
  double **Ax, **Az;   
  /* tolerance to accept ev */
  double tol;
  /* total #ev computed; the computed eig val/vec */
  int nevOut;
  /*-------------------- max Lanczos steps */
  int max_its;
  /*-------------------- poly. deg. and # of vectors for kpmdos */
  int Mdeg, nvec;
  /*-------------------- IO */
  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE]; 
  /* initial vector: random */
  double *vinit;
  tol = 1e-8;
  /* parameters for rational filter */
  int pow = 2; // multiplicity of the pole 
  double beta = 0.01; // beta in the LS approximation
  // slicer parameters 
  npnts = 1000;
  Mdeg = 100;  nvec = 60;      
  // interior eigensolver parameters  
  double *mu = malloc((Mdeg+1)*sizeof(double)); // coeff. for kpmdos
  int *counts; // #ev computed in each slice
  double *sli; // endpoints of partitioned slices
  /*------------------ file "matfile" contains paths to matrices */
  if( NULL == ( fmat = fopen( "matfile", "r" ) ) ) {
    fprintf( flog, "Can't open matfile...\n" );
    exit(2);
  }
  /*-------------------- read number of matrices ..*/  
  memset( line, 0, MAX_LINE );
  fgets( line, MAX_LINE, fmat );
  if( ( numat = atoi( line ) ) <= 0 ) {
    fprintf( flog, "Invalid count of matrices...\n" );
    exit(3);
  }
  /*-------------------- LOOP through matrices -*/
  for(mat = 1; mat <= numat; mat++) {
    if(get_matrix_info(fmat, &io) != 0) {
      fprintf(flog, "Invalid format in matfile ...\n");
      exit(5);
    }
    /*----------------input matrix and interval information -*/
    fprintf(flog, "MATRIX: %s...\n", io.MatNam);
    a = io.a; // left endpoint of input interval
    b = io.b; // right endpoint of input interval
    n_intv = io.n_intv;
    char path[1024] ;   // path to write the output files
    strcpy( path, "OUT/GenRLanR_");
    strcat( path, io.MatNam);
    fstats = fopen(path,"w"); // write all the output to the file io.MatNam 
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX: %s...\n", io.MatNam);
    fprintf(fstats,"Partition the interval of interest [%f,%f] into %d slices\n", a,b,n_intv);
    counts = malloc(n_intv*sizeof(int)); 
    sli = malloc( (n_intv+1)*sizeof(double));
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB){
      ierr =read_coo_MM(io.Fname, 1, 0, &Acoo); 
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
        nnz = Acoo.nnz; 
        n = Acoo.nrows;
        printf("size of matrix is %d\n", n);
        printf("nnz of matrix is %d\n", nnz);
      }
      else {
        fprintf(flog, "read_coo error = %d\n", ierr);
        exit(6);
      }
      /*-------------------- conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
    }
    if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT  not supported (yet) * \n");
      exit(7);
    }
    // format convertion for umfpack 
    solveShift solshift;
    diag = (int *)malloc(n*sizeof(int)); // location of diagonal entries 
    Ap = (SuiteSparse_long *)malloc((n+1)*sizeof(SuiteSparse_long));    // pointer to starting point of each row
    Ai = (SuiteSparse_long *)malloc(nnz*sizeof(SuiteSparse_long));  // column indices
    Ax = (double **)malloc(n_intv*sizeof(double *));//Ax[i] = real part of the ith shifted matrix 
    Az = (double **)malloc(n_intv*sizeof(double *));//Az[i]: imaginary part of the ith shfited matrix
    // each shifted matrix shares the same Ap and Ai
    for (i=0; i<n+1; i++) {
      Ap[i] = Acsr.ia[i];
    }
    for (i=0; i<nnz; i++) {
      Ai[i] = Acsr.ja[i];
    }
    //--------------------find the location of the diagonal entries in CSR format n_intv shifted 
    //--------------------matrices share the same diag, Ax and Az
    for (i = 0 ; i < n ; i++){
      for (j = Ap[i] ; j < Ap[i+1] ; j++){      
        k = Ai[j];
        if(i == k){
          diag[i] = j;
          break;
        }
      }
    }
    /*-------------------- define parameters for DOS*/
    alleigs = malloc(n*sizeof(double)); 
    vinit = (double*) malloc(n*sizeof(double));
    rand_double(n, vinit);
    /*-------------------- get lambda_min lambda_max estimates */
    ierr = LanBounds(&Acsr, 60, vinit, &lmin, &lmax);
    fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
    /*-------------------- define [a b] now so we can get estimates now    
      on number of eigenvalues in [a b] from kpmdos */
    fprintf(fstats," --> interval: a  %9.3e  b %9.3e \n",a, b);
    /*-------------------- define kpmdos parameters */
    xintv[0] = a;  xintv[1] = b;
    xintv[2] = lmin; xintv[3] = lmax;
    //-------------------- trivial slicer or kpmdos? 
    if (TRIV_SLICER) {
      linspace(a, b, n_intv+1,  sli);
    } else {
      double t = cheblan_timer();
      ierr = kpmdos(&Acsr, Mdeg, 0, nvec, xintv, mu, &ecount);
      t = cheblan_timer() - t;
      if (ierr) {
	printf("kpmdos error %d\n", ierr);
	return 1;
      }
      fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);      
      fprintf(fstats, "Step 1a: Estimated eig count in interval - %10.2e \n",ecount);
      if ((ecount <0) | (ecount > n)) {
        printf(" e-count estimate is incorrect \n ");
        exit(7);
      }
      /*-------------------- error in ecount */
      /*
      if ((ecount <1) | (ecount>n))
        exit(6);
      */
      /*-------------------- define slicer parameters */
      npnts = 10 * ecount;
      ierr = spslicer(sli, mu, Mdeg, xintv, n_intv, npnts);
      if (ierr) {
        printf("spslicer error %d\n", ierr);
        return 1;
      }
    }
    fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",Mdeg, nvec, npnts);
    fprintf(fstats, "Step 1b: Slices found: \n");
    for (j=0; j<n_intv;j++)
      fprintf(fstats,"[% 12.4e , % 12.4e]\n", sli[j],sli[j+1]);
    //-------------------- # eigs per slice
    //-------------------- approximate number of eigenvalues wanted
    double fac = 1.2;   // this factor insures that # of eigs per slice is slightly overestimated 
    printf("Total number of eigenvalues estimated = %d \n", (int)(ecount));
    nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
    nev = (int)(fac*nev);                        // want an overestimate of ev_int 
    fprintf(fstats,"Step 2: In each slice compute %d eigenvalues ... \n", nev);
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    totcnt = 0;
    for (sl =0; sl<n_intv; sl++){
      // Allocate memory for the slth shifted matrix  A - (mid+width*1i)I    
      Ax[sl] = (double *)malloc(nnz*sizeof(double));//real part of each shifted matrix entries
      Az[sl] = (double *)malloc(nnz*sizeof(double));//imaginary part of each shifted matrix entries    
      for(i=0; i<nnz; i++){
        Az[sl][i] = 0.0;
      } 
      fprintf(fstats,"======================================================\n");
      double *lam, *Y, *res;
      int *ind;
      //-------------------- 
      a = sli[sl];
      b = sli[sl+1];
      fprintf(fstats, " subinterval: [%12.4e , %12.4e]\n", a, b); 
      //-------------------- Parameters for RatLanTr
      mlan = max(4*nev,100); // maximal Krylov subspace dim.
      mlan = min(mlan, n); 
      max_its = 3*mlan;
      fprintf(fstats, "Thick Restarted Lanczos with dimension %d\n", mlan);
      fprintf(fstats, "Max Lanczos steps %d\n", max_its);
      double intv[4];
      intv[0] = a; 
      intv[1] = b;
      intv[2] = lmin;
      intv[3] = lmax;
      printf("==============Compute the %dth subinterval: [%.4e, %.4e] out of %d==============\n",
          sl+1, a, b, n_intv); 
      // Find the rational filter on this slice
      ratparams rat;
      // Set up default parameters for rat
      set_ratf_def(&rat);
      // change some default values here:
      rat.beta = beta;
      rat.pow = pow;
      // Now determine rational filter
      find_ratf(intv, &rat);    
      double mid = rat.cc;
      double width = rat.dd;       
      // Form the shifted matrix A - (mid+width*1i)*I
      memcpy(Ax[sl],Acsr.a,nnz*sizeof(double));
      for(i = 0 ; i <n ; i++) {
        Ax[sl][diag[i]] = Ax[sl][diag[i]] - mid;
        Az[sl][diag[i]] = -width;
      }
      //------------Call Umfpack to factor A - (mid+width*1i)*I
      //double Control [UMFPACK_CONTROL]; umfpack_zl_defaults(Control); Control[UMFPACK_PRL] = 3;
      void *Symbolic, *Numeric;
      //------------Umfpack symbolic factorization
      int status;
      status = umfpack_zl_symbolic (n, n, Ap, Ai, Ax[sl], Az[sl], &Symbolic, NULL, NULL);
      if (status < 0) {
        printf("umfpack_zl_symbolic failed and exit\n") ;
        return 1;
      }
      //------------Umfpack numerical factorization 
      double tall = 0.0;
      tall = cheblan_timer();
      status = umfpack_zl_numeric(Ap, Ai, Ax[sl], Az[sl], Symbolic, &Numeric, NULL, NULL);
      if (status < 0) {
        printf("umfpack_zl_numeric failed and exit\n") ;
        return 0;
      }
      //Release the memory
      umfpack_zl_free_symbolic(&Symbolic);
      tall = cheblan_timer()-tall;
      printf("factorization  time :    %.2f\n", tall);

      // special case for 1 pole
      umfdata *umf = (umfdata *) malloc(1*sizeof(umfdata));
      umf[0].Ap = Ap;  umf[0].Ai = Ai;  umf[0].Ax = Ax[sl];  umf[0].Az = Az[sl];  umf[0].Numeric = Numeric;
      // structure for linear solver passed to RatLanNr
      // func [function pointer]: function to solve complex shifted system
      // it must be a function of the following type:
      // void (*solveShiftFunc)(int n, double *br, double *bz,
      //                        double *xr, double *xz, void *data);
      // void **data: encapsulates all data needed by the direct solver
      // An example is given here for using UMFPACK by T.Davis
      solshift.func = &solvefunc;  
      solshift.data = malloc(1*sizeof(void *));
      solshift.data[0] = &umf[0];

      //-------------------- RationalLanTr
      ierr = RatLanTr(&Acsr, &solshift, mlan, nev, intv, &rat, max_its, tol, vinit, &nevOut, &lam, &Y, &res, fstats);

      if (ierr) {
        printf("RatLanTr error %d\n", ierr);
        return 1;
      }
      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = (int *) malloc(nevOut*sizeof(int));
      sort_double(nevOut, lam, ind);
      /* print eigenvalues */
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      fprintf(fstats, "                                   Eigenvalues in [%f, %f]\n",a,b);
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      fprintf(fstats, "    Computed [%d out of %d estimated]           ||Res||     ", nevOut, nev);
      fprintf(fstats, "\n");
      for (i=0; i<nevOut; i++) {
        fprintf(fstats, "        % .15e                 %.1e", lam[i], res[ind[i]]);
        fprintf(fstats,"\n");
      }
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

      memcpy(&alleigs[totcnt],lam,nevOut*sizeof(double));
      totcnt += nevOut;
      counts[sl] = nevOut;
      //-------------------- free memory allocated within loop
      free(Ax[sl]);
      free(Az[sl]);
      if (lam)  free(lam);
      if (Y)  free(Y);
      if (res)  free(res);
      free(ind);
      free(rat.omega);
      umfpack_zl_free_numeric (&Numeric);
      free(umf);
      free(solshift.data);
      /*-------------------- end slice loop */
    }
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);
    sprintf(path, "OUT/EigsOut_RLanR_%s", io.MatNam);
    FILE *fmtout = fopen(path,"w");
    if (fmtout) {
      for (j=0; j<totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }
    /*-------------------- free memory for matrix loop */
    free(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free(sli);
    free(counts);
    free(Ap);
    free(diag);
    free(Ai);
    free(Ax);
    free(Az);
    free(alleigs);
    if (fstats != stdout) fclose(fstats);
    /*-------------------- end matrix loop */
  }
  //-------------------- free rest of memory 
  free(mu);
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
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

