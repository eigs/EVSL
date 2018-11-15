#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "evsl.h"
#include "io.h"
#include "evsl_direct.h"

#define TRIV_SLICER 0

int main () {
  int ierr = 0;
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix
   * read in matrix format
   * compute lmin and lmax internally. needs interval of eigenvalues
   * to compute and number of slices - both of which are entered via
   * the file matfile/
   *-------------------------------------------------------------*/
  int n=0, sl, i, j,  nev, totcnt;
  //int nnz;
  double a, b, ecount, xintv[4];
  double lmin, lmax;
  double *alleigs;
  int n_intv;      // number of subintervals (slices)
  int npnts;       // number of integration points for eigenvalue count
  //--------------------related to kpmdos
  //int Mdeg = 180,  nvec = 100;    // benzene
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
  /* tolerance to accept ev */
  double tol;
  /* stats */
  /* #ev computed; the computed eig val/vec,
   * residuals of the computed eig pairs */
  int nevOut;
  //  double *lam=NULL, *Y=NULL, *res=NULL ;
  /*-------------------- mim/max degree of  polynomial, max Lanczos steps */
  int max_its, Mdeg, nvec;
  /*-------------------- IO */
  FILE *flog = stdout, *fmat = NULL, *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  /* initial vector: random */
  double *vinit;

  /* setup solver parameters, all in one place */
  /* ------- C60 --------------*/
  //  a = -1.134037080628220E+01; // = lmin
  //  /* ------- Si5H12 --------------*/
  //  a = -1;  b = 5.896;

  //  /* ------- Ga10As10H30 --------------*/
  //  a = -1.215;  b = 1.6877208;
  tol = 1e-12;
  // set parameters
  npnts = 1000;
  Mdeg = 300;
  nvec = 60;
  /*-------------------- start EVSL */
  EVSLStart();

  // interior eigensolver parameters
  max_its = 1000; //1000;  // max number of iterations

  double *mu = evsl_Malloc(Mdeg+1, double);
  int *counts;
  double *sli;
  /* -------------------------------- */

  /*------------------ file "matfile" contains paths to matrices */
  if( NULL == ( fmat = fopen( "matfile", "r" ) ) ) {
    fprintf( flog, "Can't open matfile...\n" );
    exit(2);
  }
  /*-------------------- read number of matrices ..*/
  memset( line, 0, MAX_LINE );
  if (NULL == fgets( line, MAX_LINE, fmat )) {
    fprintf( flog, "error in reading matfile...\n" );
    exit(2);
  }
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
    /*--------------------------input matrix and interval information -*/
    fprintf(flog, "MATRIX: %s...\n", io.MatNam1);
    a = io.a;
    b = io.b;
    n_intv = io.n_intv;

    struct stat st = {0}; /* Make sure OUT directory exists */
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    char path[1024] ;   // path to write the output files
    strcpy( path, "OUT/MMRSI_");
    strcat( path, io.MatNam1);
    //-------------------- where to write output
    fstats = fopen(path,"w"); // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX: %s...\n", io.MatNam1);
    fprintf(fstats,"Partition the interval of interest [%f,%f] into %d slices\n", a,b,n_intv);
    counts = evsl_Malloc(n_intv, int);
    sli = evsl_Malloc(n_intv+1, double);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB){
      ierr =read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
        //nnz = Acoo.nnz;
        n = Acoo.nrows;
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
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- define ChebLanTr parameters */
    alleigs = evsl_Malloc(n, double);
    vinit = evsl_Malloc(n, double);
    rand_double(n, vinit);
    /*-------------------- get lambda_min lambda_max estimates */
    ierr = LanTrbounds(50, 200, 1e-8, vinit, 1, &lmin, &lmax, fstats);
    //    lmin = lmin-0.1; lmax = lmax+0.1;
    fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
    fprintf(fstats," --> interval: a  %9.3e  b %9.3e \n",a, b);
    /*-------------------- define kpmdos parameters */
    xintv[0] = a;    xintv[1] = b;
    xintv[2] = lmin; xintv[3] = lmax;

    ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
    if (ierr) {
      printf("kpmdos error %d\n", ierr);
      return 1;
    }
    if (TRIV_SLICER) {
      linspace(a, b, n_intv+1,  sli);
    } else {
      //    fprintf(flog," -->  ecount = %e\n",ecount);
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
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    totcnt = 0;
    for (sl =0; sl<n_intv; sl++){
      fprintf(fstats,"======================================================\n");
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      a = sli[sl];
      b = sli[sl+1];
      fprintf(fstats, " subinterval: [% 12.4e , % 12.4e]\n", a, b);

      double fac = 1.2;
      nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
      nev = evsl_max( (int)(fac*nev), nev+5 );                        // want an overestimate of ev_int
      max_its = 1000;

      fprintf(fstats, "Filtered subspace iteration with block size %d\n", nev);
      fprintf(fstats, "Max steps %d\n", max_its);
   
      // rat filter
      double intv[4];
      intv[0] = a;
      intv[1] = b;
      intv[2] = lmin;
      intv[3] = lmax;
      // find the rational filter on this slice
      ratparams rat;
      // setup default parameters for rat
      set_ratf_def(&rat);
      // now determine rational filter
      find_ratf(intv, &rat);
      
      // set up direct solver
      void **solshiftdata = evsl_Malloc(rat.num, void *);
      SetupASIGMABSolDirect(&Acsr, NULL, rat.num, nev, rat.zk, solshiftdata);
      for (i=0; i<rat.num; i++) {
         SetASigmaBSol(&rat, i, ASIGMABSolDirect, solshiftdata[i]);
      }

      ierr = RatSI(nev, intv, max_its, tol, &rat, &nevOut, &lam, &Y, &res);

      if (ierr) {
        printf("RatSI error %d\n", ierr);
      }

      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = evsl_Malloc(nevOut, int);
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

      fprintf(fstats,"======================================================\n");
      /*-------------------- print results [eigenvalues] */
      memcpy(&alleigs[totcnt],lam,nevOut*sizeof(double));
      totcnt += nevOut;
      counts[sl] = nevOut;
      //-------------------- free memory allocated in this loop
      if (lam)  evsl_Free(lam);
      if (Y)  evsl_Free(Y);
      if (res)  evsl_Free(res);
      FreeASIGMABSolDirect(rat.num, solshiftdata);
      evsl_Free(solshiftdata);
      free_rat(&rat);
      /*-------------------- end slice loop */
    }
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);
    sprintf(path, "OUT/EigsOut_PSI_%s",io.MatNam1);
    FILE *fmtout = fopen(path,"w");
    if (fmtout) {
      for (j=0; j<totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }
    /*-------------------- free rest of allocated memory */
    evsl_Free(counts);
    evsl_Free(sli);
    evsl_Free(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
    evsl_Free(alleigs);
    if (fstats != stdout) fclose(fstats);
    /*-------------------- end matrix loop */
  }
  evsl_Free(mu);
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

