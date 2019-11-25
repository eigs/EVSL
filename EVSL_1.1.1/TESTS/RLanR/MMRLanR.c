#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include <complex.h>
#include "evsl.h"
#include "io.h"
#include "evsl_direct.h"

#define TRIV_SLICER 0

int main () {
  int ierr = 0;
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix
   * read in matrix format -- using
   * non-restart Lanczos with rational filtering.
   *-------------------------------------------------------------*/
  int n=0, sl, i, j, mlan, nev, totcnt;
  double a, b, ecount, xintv[4];
  double lmin, lmax;
  double *alleigs;
  int n_intv;      // number of subintervals (slices)
  int npnts;       // number of integration points for eigenvalue count
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
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
  /* slicer parameters */
  npnts = 1000;
  Mdeg = 300;
  nvec = 100;
  /*-------------------- start EVSL */
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
  EVSLStart();
  /* interior eigensolver parameters */
  double *mu = evsl_Malloc(Mdeg+1, double); // coeff. for kpmdos
  int *counts; // #ev computed in each slice
  double *sli; // endpoints of partitioned slices
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
    a = io.a; // left endpoint of input interval
    b = io.b; // right endpoint of input interval
    n_intv = io.n_intv;

    struct stat st = {0}; /* Make sure OUT directory exists */
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    char path[1024];   // path to write the output files
    strcpy( path, "OUT/MMRLanR_");
    strcat( path, io.MatNam1);
    fstats = fopen(path,"w"); // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX: %s...\n", io.MatNam1);
    fprintf(fstats,"Partition the interval of interest [%f,%f] into %d slices\n",
            a,b,n_intv);
    counts = evsl_Malloc(n_intv, int);
    sli = evsl_Malloc(n_intv+1, double);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB){
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
    /*-------------------- define parameters for DOS */
    alleigs = evsl_Malloc(n, double);
    vinit = evsl_Malloc_device(n, double);
    rand_double_device(n, vinit);
    /*-------------------- get lambda_min lambda_max estimates */
    ierr = LanTrbounds(50, 200, 1e-8, vinit, 1, &lmin, &lmax, fstats);
    fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n",
            lmin, lmax);
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
      double t = evsl_timer();
      ierr = kpmdos(Mdeg, 0, nvec, xintv, mu, &ecount);
      t = evsl_timer() - t;
      if (ierr) {
        printf("kpmdos error %d\n", ierr);
        return 1;
      }
      fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);
      fprintf(fstats, "Step 1a: Estimated eig count in interval - %10.2e \n",
              ecount);
      if (ecount < 0 || ecount > n) {
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
    fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",
            Mdeg, nvec, npnts);
    fprintf(fstats, "Step 1b: Slices found: \n");
    for (j=0; j<n_intv;j++)
      fprintf(fstats,"[% 12.4e , % 12.4e]\n", sli[j],sli[j+1]);
    //-------------------- # eigs per slice
    //-------------------- approximate number of eigenvalues wanted
    double fac = 1.2;   // this factor insures that # of eigs per slice is slightly overestimated
    printf("Total number of eigenvalues estimated = %d \n", (int)(ecount));
    nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
    nev = (int)(fac*nev);                        // want an overestimate of ev_int
    fprintf(fstats, "Step 2: In each slice compute %d eigenvalues ... \n", nev);
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    totcnt = 0;
    for (sl=0; sl<n_intv; sl++) {
      fprintf(fstats,"======================================================\n");
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      a = sli[sl];
      b = sli[sl+1];
      fprintf(fstats, " subinterval: [% 12.4e , % 12.4e]\n", a, b);
      //-------------------- Parameters for RatLanTr
      mlan = evsl_max(4*nev,300); // maximal Krylov subspace dim.
      mlan = evsl_min(mlan, n);
      max_its = 3*mlan;
      fprintf(fstats, "Thick Restarted Lanczos with dimension %d\n", mlan);
      fprintf(fstats, " Max Lanczos steps %d\n", max_its);
      double intv[4]; // endpoints of this slice and spectrum
      intv[0] = a;
      intv[1] = b;
      intv[2] = lmin;
      intv[3] = lmax;
      printf("=== Compute subinterval %d: [%.4e, %.4e] out of %d ===\n",
             sl+1, a, b, n_intv);
      //------------Find the rational filter on this slice
      ratparams rat;
      //------------Set up default parameters for rat
      set_ratf_def(&rat);
      //------------Change some default values here:
      rat.beta = beta;
      rat.pow = pow;
      //-------------Now determine rational filter
      find_ratf(intv, &rat);
      /*------------ use direct solver function  */
      void **solshiftdata = evsl_Malloc(rat.num, void *);
      /*------------ factoring the shifted matrices and store the factors */
      SetupASIGMABSolDirect(&Acsr, NULL, rat.num, rat.zk, solshiftdata);
      /*------------ set the solver for A-sI in rat */
      for (i=0; i<rat.num; i++) {
        SetASigmaBSol(&rat, i, ASIGMABSolDirect, solshiftdata[i]);
      }
      //-------------------- RationalLanTr
      ierr = RatLanTr(mlan, nev, intv, max_its, tol, vinit, &rat, &nevOut, &lam,
                      &Y, &res, fstats);
      if (ierr) {
        printf("RatLanTr error %d\n", ierr);
        return 1;
      }
      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = evsl_Malloc(nevOut, int);
      sort_double(nevOut, lam, ind);

      /* print eigenvalues */
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      fprintf(fstats, "                                   Eigenvalues in [%f, %f]\n",
              a,b);
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      fprintf(fstats, "    Computed [%d out of %d estimated]           ||Res||     ",
              nevOut, nev);
      fprintf(fstats, "\n");
      for (i=0; i<nevOut; i++) {
        fprintf(fstats, "        % .15e                 %.1e", lam[i], res[ind[i]]);
        fprintf(fstats,"\n");
      }
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      memcpy(&alleigs[totcnt],lam,nevOut*sizeof(double));
      totcnt += nevOut;
      counts[sl] = nevOut;
      /*-------------------- free memory within  this loop */
      if (lam)  evsl_Free(lam);
      if (Y)  evsl_Free_device(Y);
      if (res)  evsl_Free(res);
      evsl_Free(ind);
      FreeASIGMABSolDirect(rat.num, solshiftdata);
      evsl_Free(solshiftdata);
      free_rat(&rat);
      /*-------------------- end slice loop */
    }
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);
    sprintf(path, "OUT/EigsOut_RLanR_%s", io.MatNam1);
    FILE *fmtout = fopen(path,"w");
    if (fmtout) {
      for (j=0; j<totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }
    /*-------------------- free memory within  this loop */
    evsl_Free_device(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
#ifdef EVSL_USING_CUDA_GPU
    evsl_free_csr_gpu(&Acsr_gpu);
#endif
    evsl_Free(sli);
    evsl_Free(counts);
    evsl_Free(alleigs);
    if (fstats != stdout) fclose(fstats);
    /*-------------------- end matrix loop */
  }
  //-------------------- free rest of memory
  evsl_Free(mu);
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

