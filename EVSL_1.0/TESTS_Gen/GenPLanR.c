#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include "evsl.h"
#include "io.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define TRIV_SLICER 0

/*-------------------- Protos */
int read_coo_MM(const char *matfile, int idxin, int idxout,   cooMat *Acoo); 
int get_matrix_info( FILE *fmat, io_t *pio );
/*-------------------- End Protos */

int main () { 
  int ierr = 0;
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix
   * read in sparse matrix format
   * !! needed: fix lanbounds 
   *-------------------------------------------------------------*/
  int n=0, sl, i, j, mlan, nev, totcnt; 
  // int nnz;
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
  polparams pol;
  //-------------------- tolerance for stopping criterion
  tol = 1e-10;
  
  //-------------------- slicer parameters 
  npnts = 1000;
  Mdeg = 100;  nvec = 100;      // C60
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
    a = io.a;
    b = io.b;
    n_intv = io.n_intv;
    /*-------------------- path to write the output files*/
    char path[1024] ;   
    strcpy( path, "OUT/GenPLanR_");
    strcat( path, io.MatNam);
    // write all the output to the file io.MatNam 
    fstats = fopen(path,"w"); 
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
    /*-------------------- define ChebLanTr parameters */
    alleigs = malloc(n*sizeof(double)); 
    vinit = (double *) malloc(n*sizeof(double));
    rand_double(n, vinit);
    /*-------------------- get lambda_min lambda_max estimates */
    ierr = LanBounds(&Acsr, 60, vinit, &lmin, &lmax);
    //    lmin = lmin-0.1; lmax = lmax+0.1;
    //fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
    fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
    /*-------------------- define [a b] now so we can get estimates now    
      on number of eigenvalues in [a b] from kpmdos */
    fprintf(fstats," --> interval: a  %9.3e  b %9.3e \n",a, b);
    /*-------------------- define kpmdos parameters */
    xintv[0] = a;  xintv[1] = b;
    xintv[2] = lmin; xintv[3] = lmax;

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
      if ((ecount <1) | (ecount>n))
        exit(6);
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
      fprintf(fstats, " %2d: [%.15e , %.15e]\n", j+1, sli[j],sli[j+1]);
    //-------------------- # eigs per slice
    //-------------------- approximate number of eigenvalues wanted
    double fac = 1.2;   // this factor insures that # of eigs per slice is slightly overestimated 
    nev = (int) (1 + ecount / ((double) n_intv));  // # eigs per slice
    nev = (int)(fac*nev);                        // want an overestimate of ev_int 
    fprintf(fstats,"Step 2: In each slice compute %d eigenvalues ... \n", nev);
    /*-------------------- MAIN intv LOOP: for each sclice Do: */
    totcnt = 0;
    for (sl =0; sl<n_intv; sl++) {
      fprintf(fstats,"======================================================\n");
      double *lam, *Y, *res;
      int *ind;
      //-------------------- 
      a = sli[sl];
      b = sli[sl+1];
      fprintf(fstats, " subinterval: [% 12.4e , % 12.4e]\n", a, b); 
      //-------------------- Parameters for ChebLanTr
      mlan = max(4*nev,100);  mlan = min(mlan, n); 
      max_its = 3*mlan;  // max number of Lanczos iterations
      fprintf(fstats, "Thick Restarted Lanczos with dimension %d\n", mlan);
      fprintf(fstats, "Max Lanczos steps %d\n", max_its);
      xintv[0] = a; 
      xintv[1] = b;
      //-------------------- set up default parameters for pol.
      set_pol_def(&pol);
      // can change default values here e.g.
      // pol.max_deg = 500; pol.damping = 1; 
      //-------------------- Now determine polymomial
      find_pol(xintv, &pol); 
      //-------------------- Call ChebLanTr        
      ierr = ChebLanTr(&Acsr, mlan, nev, xintv, max_its, tol, vinit,
          &pol, &nevOut, &lam, &Y, &res, fstats);

      if (ierr) {
        printf("ChebLanTr error %d\n", ierr);
        return 1;
      }

      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = (int*) malloc(nevOut*sizeof(int));
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
      //-------------------- free memory allocated within this loop
      if (lam)  free(lam);
      if (Y)  free(Y);
      if (res)  free(res);
      free(ind);
      free(pol.mu);
      /*-------------------- end slice loop */
    }
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);

    sprintf(path, "OUT/EigsOut_PLanR_%s", io.MatNam);
    FILE *fmtout = fopen(path,"w");
    if (fmtout) {
      for (j=0; j<totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }
    /*-------------------- free memory */
    free(counts);
    free(sli);
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
  return 0;
}
