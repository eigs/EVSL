#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <omp.h>
#include "evsl.h"
#include "io.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define TRIV_SLICER 0

int main () { 
  int ierr = 0;
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix
   * read in sparse matrix format
   * Uses:
   * Thick-restart Lanczos with pol. filtering
   * this driver implements thread parallelism across slices/ 
   *-------------------------------------------------------------*/
  int n=0, sl, i, j, mlan, nev, totcnt;
  //int nnz;
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
  /* tolerance to accept ev */
  double tol;
  /* stats */
  /* #ev computed; the computed eig val/vec, 
   * residuals of the computed eig pairs */
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
  tol = 1e-8;
  //-------------------- slicer parameters
  Mdeg = 300;
  nvec = 60;
  /*-------------------- start EVSL */
  EVSLStart();
  //-------------------- interior eigensolver parameters  
  double *mu = malloc((Mdeg+1)*sizeof(double));
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
    fprintf(flog, "MATRIX: %s...\n", io.MatNam);
    a = io.a;
    b = io.b;
    n_intv = io.n_intv;
    /*-------------------- path to write the output files*/
    char path[1024];
    strcpy( path, "OUT/MMPLanR_OMP_");
    strcat( path, io.MatNam);
    fstats = fopen(path,"w"); // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX: %s...\n", io.MatNam);
    fprintf(fstats,"Partition the interval of interest [%f,%f] into %d slices\n", a,b,n_intv);
    counts = malloc(n_intv*sizeof(int)); 
    sli = malloc((n_intv+1)*sizeof(double));
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr =read_coo_MM(io.Fname, 1, 0, &Acoo); 
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
    alleigs = (double *) malloc(n*sizeof(double)); 
    vinit = (double *) malloc(n*sizeof(double));
    rand_double(n, vinit);
    /*-------------------- get lambda_min lambda_max estimates */
    ierr = LanBounds(60, vinit, &lmin, &lmax);
    fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
    /*-------------------- define [a b] now so we can get estimates now    
                           on number of eigenvalues in [a b] from kpmdos */
    fprintf(fstats," --> interval: a  %9.3e  b %9.3e \n",a, b);
    /*-------------------- define kpmdos parameters */
    xintv[0] = a;  xintv[1] = b;
    xintv[2] = lmin; xintv[3] = lmax;
    //-------------------- kpmdos or triv_slicer
    if (TRIV_SLICER) {
      linspace(a, b, n_intv+1,  sli);
    } else {
      double t = cheblan_timer();
      ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
      t = cheblan_timer() - t;
      if (ierr) {
        printf("kpmdos error %d\n", ierr);
        return 1;
      }
      fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);
      fprintf(fstats, "Step 1a: Estimated eig count in interval - %.15e \n",ecount);
      if ((ecount <0) | (ecount > n)) {
        printf(" e-count estimate is incorrect \n ");
        exit(7);
      }
      /*-------------------- error in ecount */
      /*
      if (ecount < 1 || ecount > n)
        exit(6);
      */
      /*-------------------- define slicer parameters */
      npts = 10 * ecount;
      ierr = spslicer(sli, mu, Mdeg, xintv, n_intv, npts);
      if (ierr) {
        printf("spslicer error %d\n", ierr);
        return 1;
      }
    }
    fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npts = %d\n",
            Mdeg, nvec, npts);
    fprintf(fstats, "Step 1b: Slices found: \n");
    for (j=0; j<n_intv;j++) {
      fprintf(fstats, " %2d: [%.15e , %.15e]\n", j+1, sli[j],sli[j+1]);
    }
    //-------------------- # eigs per slice
    //-------------------- approximate number of eigenvalues wanted
    double fac = 1.2;   // this factor insures that # of eigs per slice is slightly overestimated 
    nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
    nev = (int)(fac*nev);                        // want an overestimate of ev_int 
    fprintf(fstats,"Step 2: In each slice compute %d eigenvalues ... \n", nev);
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    double tsolve = cheblan_timer();
    totcnt = 0;
    mlan = max(4*nev,100);   mlan = min(n, mlan);
    max_its = 3*mlan; 
    // shared array for storing eigenvalues and residuals for each slice 
    double *lam_global, *res_global;
    lam_global = (double *) malloc(nev*n_intv*sizeof(double));
    res_global = (double *) malloc(nev*n_intv*sizeof(double));

#pragma omp parallel for private(sl)
    for (sl=0; sl<n_intv; sl++) {
      int nevOut, j;
      double *lam, *Y, *res;
      double *res_sorted;
      int *ind;
      double intv[4];
      polparams pol;
      /*
         int tmp1 = omp_get_thread_num()+1; 
         int tmp2 = omp_get_num_threads(); 
         #pragma omp critical
         fprintf(flog,"Thread %d out of %d:  \n", tmp1, tmp2);
      */
      //-------------------- ChebLanTr
      intv[0] = sli[sl];
      intv[1] = sli[sl+1];
      intv[2] = lmin; 
      intv[3] = lmax;

      set_pol_def(&pol);
      find_pol(intv, &pol);
      ierr = ChebLanTr(mlan, nev, intv, max_its, tol, vinit,
                       &pol, &nevOut, &lam, &Y, &res, NULL);
      if (ierr) {
        printf("ChebLanTr error %d\n", ierr);
      }
      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = (int *) malloc(nevOut*sizeof(int));
      sort_double(nevOut, lam, ind);
      // sort residuals accordingly 
      res_sorted = (double *) malloc(nevOut*sizeof(double));
      for (j=0; j<nevOut; j++){
        res_sorted[j] = res[ind[j]];
      }       

      /* copy result to global arrays */
      memcpy(&lam_global[sl*nev],lam,nevOut*sizeof(double));
      memcpy(&res_global[sl*nev],res_sorted,nevOut*sizeof(double));
      counts[sl] = nevOut;
      //-------------------- free memory within this loop
      if (lam)  free(lam);
      if (Y)  free(Y);
      if (res)  free(res);
      free(res_sorted);
      free(ind);
      free_pol(&pol);
      /*-------------------- end slice loop */
    } // end of parallel for-loop

    tsolve = cheblan_timer() - tsolve;

    // Now output result and move computed eigenvalues to alleigs
    for (sl=0; sl<n_intv; sl++) {
      //-------------------- 
      a = sli[sl];
      b = sli[sl+1];

      fprintf(fstats,"======================================================\n");
      fprintf(fstats, " subinterval: [% 12.4e , % 12.4e ]\n", a, b); 
      fprintf(fstats, "Thick Restarted Lanczos with dimension %d\n", mlan);
      fprintf(fstats, "Max Lanczos steps %d\n", max_its);

      /* print eigenvalues */
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      fprintf(fstats, "                                   Eigenvalues in [%f, %f]\n",a,b);
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      fprintf(fstats, "    Computed [%d out of %d estimated]           ||Res||     ", counts[sl], nev);
      fprintf(fstats, "\n");
      for (i=0; i<counts[sl]; i++) {
        fprintf(fstats, "        % .15e                 %.1e",
                lam_global[sl*nev+i], res_global[sl*nev+i]);
        fprintf(fstats,"\n");
      }
      fprintf(fstats, "- - - -  - - - - - - - - - - - - - - - - - - - - - -\n");

      memcpy(&alleigs[totcnt],&lam_global[sl*nev],counts[sl]*sizeof(double));
      totcnt += counts[sl];
    }
    fprintf(fstats, "Solution  time :    %.2f\n", tsolve);
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);
    //FILE *fmtout = fopen("EigsOut","w");
    //for (j=0; j<totcnt; j++)
    //  fprintf(fmtout, "%22.15e\n", alleigs[j]);
    //fclose(fmtout);
    /*-------------------- free memory */
    free(counts);
    free(sli);
    free(lam_global);
    free(res_global);
    free(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free(alleigs);
    if (fstats != stdout) fclose(fstats);
    /*-------------------- end matrix loop */
  }
  //-------------------- free rest of memory 
  free(mu);
  if( flog != stdout ) fclose ( flog );
  fclose( fmat );
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}
