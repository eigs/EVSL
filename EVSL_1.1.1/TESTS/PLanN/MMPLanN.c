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
   * read in sparse matrix format
   *-------------------------------------------------------------*/
  int n=0, sl, i, j, mlan, nev, totcnt;
  double a, b, ecount, xintv[4];
  double lmin, lmax;
  double *alleigs;
  int n_intv;      // number of subintervals (slices)
  int npts;       // number of integration points for eigenvalue count
  //--------------------related to kpmdos
  //int Mdeg = 180,  nvec = 100;
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
  /*-------------------- tolerance to accept ev */
  double tol;
  /* stats */
  /* #ev computed; the computed eig val/vec,
   * residuals of the computed eig pairs */
  int nevOut;
  //  double *lam=NULL, *Y=NULL, *res=NULL ;
  /*-------------------- Max degree of polynomial, or  max Lanczos steps for DOS; numbr smpling vectors */
  int Mdeg=80, nvec=30;
  /*-------------------- IO */
  FILE *flog = stdout, *fmat = NULL, *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  /* initial vector: random */
  double *vinit, *xdos, *ydos;
  polparams pol;
  //-------------------- tolerance for stopping criterion
  tol = 1e-6;
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
  /*-------------------- start EVSL */
  EVSLStart();
  //-------------------- interior eigensolver parameters
  int *counts;
  double *sli;
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
    /*----------------input matrix and interval information -*/
    fprintf(flog, "MATRIX: %s...\n", io.MatNam1);
    a = io.a;
    b = io.b;
    n_intv = io.n_intv;

    struct stat st = {0}; /* Make sure OUT directory exists */
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    /*-------------------- path to write the output files*/
    char path[1024] ;
    strcpy( path, "OUT/MMPLanN_");
    strcat( path, io.MatNam1);
    fstats = fopen(path,"w"); // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX: %s...\n", io.MatNam1);
    fprintf(fstats,"Partition the interval of interest [%f,%f] into %d slice (s)\n", a,b,n_intv);
    counts = evsl_Malloc(n_intv, int);
    sli = evsl_Malloc(n_intv+1, double);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr =read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
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
    fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
    /*-------------------- define [a b] now so we can get estimates now
                           on number of eigenvalues in [a b] from kpmdos */
    fprintf(fstats," --> interval: a  %9.3e  b %9.3e \n",a, b);
    /*-------------------- define LanczosDOS parameters */
    xintv[0] = a;    xintv[1] = b;
    xintv[2] = lmin; xintv[3] = lmax;
    //-------------------- landos or triv_slicer
    if (TRIV_SLICER) {
      linspace(a, b, n_intv+1,  sli);
    } else {
      double t = evsl_timer();
      npts = 100*n_intv;
      xdos = evsl_Malloc(npts, double);
      ydos = evsl_Malloc(npts, double);

      //fprintf(stdout," %d %d \n",npts,n_intv);
      LanDos(nvec, Mdeg, npts, xdos, ydos, &ecount, xintv);
      //fprintf(stdout," %f \n",ecount);
      t = evsl_timer() - t;

      fprintf(fstats, " Time to build DOS (Landos) was : %10.2f  \n",t);
      fprintf(fstats, " estimated eig count in interval: %.15e \n",ecount);
      //-------------------- call splicer to slice the spectrum
      if ((ecount <0) | (ecount > n)) {
        printf(" e-count estimate is incorrect \n ");
        exit(7);
      }
      /*-------------------- error in ecount */
      if ((ecount <1) | (ecount>n))
        exit(6);
      /*-------------------- define slicer parameters */

      fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",
              Mdeg, nvec, npts);
      spslicer2(xdos, ydos, n_intv,  npts, sli);
      evsl_Free(xdos);
      evsl_Free(ydos);
    }
    /*-------------------- Dos determined+ slicing done */
    fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npts = %d\n",
            Mdeg, nvec, npts);
    fprintf(fstats, "Step 1b: Slices found: \n");
    for (j=0; j<n_intv;j++) {
      fprintf(fstats, " %2d: [%.15e , %.15e]\n", j+1, sli[j],sli[j+1]);
    }
    //-------------------- # eigs per slice
    //-------------------- approximate number of eigenvalues wanted
    double fac = 1.2;   // this factor ensures that # of eigs per slice is slightly overestimated
    nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
    nev = (int)(fac*nev);                        // want an overestimate of ev_int
    fprintf(fstats,"Step 2: In each slice compute %d eigenvalues ... \n", nev);
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    totcnt = 0;
    for (sl = 0; sl < n_intv; sl++) {
      fprintf(fstats,"======================================================\n");
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      a = sli[sl];
      b = sli[sl+1];
      StatsReset();
      fprintf(fstats, " subinterval: [% 12.4e , % 12.4e]\n", a, b);
      //-------------------- Parameters for ChebLanTr
      mlan = evsl_max(5*nev,300);  mlan = evsl_min(mlan, n);
      //mlan  = 2000;  // max number of Lanczos iterations
      fprintf(fstats, "Non-Restarted Lanczos with dimension %d\n", mlan);
      xintv[0] = a;
      xintv[1] = b;
      //-------------------- set up default parameters for pol.
      set_pol_def(&pol);
      // can change default values here e.g.
      pol.damping = 2;         // use lanczos damping
      pol.thresh_int = 0.8;    // change thresholds for getting a sharper
      pol.thresh_ext = 0.2;   // [higher deg] polynomial
      //-------------------- Now determine polymomial
      find_pol(xintv, &pol);

      fprintf(fstats, " polynomial [type = %d], deg %d, bar %e gam %e\n",
              pol.type, pol.deg, pol.bar, pol.gam);
      //-------------------- Call ChebLanNr
      ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nevOut, &lam, &Y,
                       &res, fstats);
      if (ierr) {
        printf("ChebLanNr error %d\n", ierr);
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
      memcpy(&alleigs[totcnt],lam,nevOut*sizeof(double));
      totcnt += nevOut;
      counts[sl] = nevOut;
      StatsPrint(stdout);
      //-------------------- free memory allocated within this loop
      if (lam) {
         evsl_Free(lam);
      }
      if (Y) {
         evsl_Free_device(Y);
      }
      if (res) {
         evsl_Free(res);
      }
      evsl_Free(ind);
      free_pol(&pol);
      /*-------------------- end slice loop */
    }
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);
    sprintf(path, "OUT/EigsOut_PLanN_%s", io.MatNam1);
    FILE *fmtout = fopen(path,"w");
    if (fmtout) {
      for (j=0; j<totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }
    /*-------------------- free memory */
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
  //-------------------- free rest of memory
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  /*-------------------- finalize EVSL */
  EVSLFinish();

  return 0;
}

