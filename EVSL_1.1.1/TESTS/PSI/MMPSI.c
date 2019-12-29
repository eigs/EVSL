#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "evsl.h"
#include "io.h"

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
  polparams pol;
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
  tol = 1e-10;
  // set parameters
  npnts = 1000;
  Mdeg = 300;
  nvec = 60;
  /*-------------------- start EVSL */
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
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
    strcpy( path, "OUT/MMPSI_");
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
#ifdef EVSL_USING_CUDA_GPU
    /*-------------------- set matrix A (csr on GPU) */
    csrMat Acsr_gpu;
    evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
    SetAMatrix_device_csr(&Acsr_gpu);
#else
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
#endif
    /*-------------------- define ChebLanTr parameters */
    alleigs = evsl_Malloc(n, double);
    vinit = evsl_Malloc_device(n, double);
    rand_double_device(n, vinit);
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
    //-------------------- # eigs per slice
    //-------------------- approximate number of eigenvalues wanted
    //    nev = 2 + (int) (1 + ecount / ((double) n_intv));
    double fac = 1.2;   // this factor insures that # of eigs per slice is slightly overestimated
    nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
    nev = (int)(fac*nev);                        // want an overestimate of ev_int
    fprintf(fstats,"Step 2: In each slice compute %d eigenvalues ... \n", nev);
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    totcnt = 0;
    for (sl =0; sl<n_intv; sl++){
      fprintf(fstats,"======================================================\n");
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      StatsReset();
      a = sli[sl];
      b = sli[sl+1];
      fprintf(fstats, " subinterval: [% 12.4e , % 12.4e]\n", a, b);
      //-------------------- Parameters for ChebLanTr
      fprintf(fstats, "Filtered subspace iteration with block size %d\n", nev);
      fprintf(fstats, "Max steps %d\n", max_its);
      xintv[0] = a;
      xintv[1] = b;
      //-------------------- random initial guess
      double *V0;
      V0 = evsl_Malloc_device(n*nev, double);
      rand_double_device(n*nev, V0);
      //-------------------- filtered subspace iteration
      set_pol_def(&pol);
      // can change default values here e.g.
      pol.damping = 0;  pol.thresh_int = 0.5;
      find_pol(xintv, &pol);
      fprintf(fstats, " polynomial deg %d, bar %e gam %e\n",pol.deg,pol.bar, pol.gam);
      ierr = ChebSI(nev, xintv, max_its, tol, V0,
          &pol, &nevOut, &lam, &Y, &res, fstats);

      if (ierr) {
        printf("ChebSI error %d\n", ierr);
        return 1;
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
      if (Y)  evsl_Free_device(Y);
      if (res)  evsl_Free(res);
      evsl_Free(ind);
      free_pol(&pol);
      evsl_Free_device(V0);
      StatsPrint(fstats);
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
    evsl_Free_device(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
#ifdef EVSL_USING_CUDA_GPU
    evsl_free_csr_gpu(&Acsr_gpu);
#endif
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

