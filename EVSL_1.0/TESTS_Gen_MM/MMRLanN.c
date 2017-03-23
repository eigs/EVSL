#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "evsl.h"
#include "io.h"
#include "evsl_suitesparse.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/*-------------------- Protos */
int read_coo_MM(const char *matfile, int idxin, int idxout,   cooMat *Acoo); 
int get_matrix_info(FILE *fmat, io_t *pio);
/*-------------------- End Protos */

int main() {
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix pair
   * read in matrix format -- using
   * Thick-Restarted Lanczos with rational filtering.
   *-------------------------------------------------------------*/
  int n=0, nnz =0, i, j, npts, nslices, nvec, Mdeg, nev, 
    max_its, ev_int, sl, ierr, totcnt;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli, *mu;
  double xintv[4];
  double *alleigs; 
  int *counts; // #ev computed in each slice  
  /* initial vector: random */
  double *vinit;
  /* parameters for rational filter */
  int num = 1; // number of poles used for each slice
  int pow = 2; // multiplicity of each pole
  double beta = 0.01; // beta in the LS approximation
  /*-------------------- matrices A, B: coo format and csr format */
  cooMat Acoo, Bcoo;
  csrMat Acsr, Bcsr;
  /* slicer parameters */  
  Mdeg = 300;
  nvec = 60;
  mu = malloc((Mdeg+1)*sizeof(double));
  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  /*-------------------- Bsol */
  BSolDataSuiteSparse Bsol;
  /*-------------------- stopping tol */
  tol = 1e-8;
  /*-------------------- start EVSL */
  EVSLStart();
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
  for(mat = 1; mat <= numat; mat++) {
    if(get_matrix_info(fmat, &io) != 0) {
      fprintf(flog, "Invalid format in matfile ...\n");
      exit(5);
    }
    /*----------------input matrix and interval information -*/
    fprintf(flog, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(flog, "MATRIX B: %s...\n", io.MatNam2);    
    a = io.a; // left endpoint of input interval
    b = io.b; // right endpoint of input interval
    nslices = io.n_intv;
    char path[1024];   // path to write the output files
    strcpy( path, "OUT/MMRLanN_");
    strcat( path, io.MatNam1);
    fstats = fopen(path,"w"); // write all the output to the file io.MatNam 
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(fstats, "MATRIX B: %s...\n", io.MatNam2);
    fprintf(fstats,"Partition the interval of interest [%f,%f] into %d slices\n",
            a,b,nslices);
    counts = malloc(nslices*sizeof(int));             
    sli = malloc( (nslices+1)*sizeof(double));
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB){
      ierr =read_coo_MM(io.Fname1, 1, 0, &Acoo); 
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
        nnz = Acoo.nnz; 
        n = Acoo.nrows;
        printf("size of A is %d\n", n);
        printf("nnz of  A is %d\n", nnz);
      }
      else {
        fprintf(flog, "read_coo error for A = %d\n", ierr);
        exit(6);
      }
      ierr =read_coo_MM(io.Fname2, 1, 0, &Bcoo); 
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
        nnz = Bcoo.nnz; 
        n = Bcoo.nrows;
        printf("size of B is %d\n", n);
        printf("nnz of  B is %d\n", nnz);
      }
      else {
        fprintf(flog, "read_coo error for B = %d\n", ierr);
        exit(6);
      }
      /*-------------------- conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr); 
      ierr = cooMat_to_csrMat(0, &Bcoo, &Bcsr);
    }
    if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT  not supported (yet) * \n");
      exit(7);
    }
    alleigs = malloc(n*sizeof(double)); 
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- set the right-hand side matrix B */
    SetBMatrix(&Bcsr);
    /*-------------------- use SuiteSparse as the solver for B */
    SetupBSolSuiteSparse(&Bcsr, &Bsol);
    /*-------------------- set the solver for B and LT */
    SetBSol(BSolSuiteSparse, (void *) &Bsol);
    SetLTSol(LTSolSuiteSparse);
    /*-------------------- for generalized eigenvalue problem */
    SetGenEig();
    /*-------------------- step 0: get eigenvalue bounds */
    //-------------------- initial vector  
    vinit = (double *) malloc(n*sizeof(double));
    rand_double(n, vinit);
    ierr = LanBounds(60, vinit, &lmin, &lmax);
    fprintf(fstats, "Step 0: Eigenvalue bound s for B^{-1}*A: [%.15e, %.15e]\n", 
	    lmin, lmax);
    /*-------------------- interval and eig bounds */
    xintv[0] = a;
    xintv[1] = b;
    xintv[2] = lmin;
    xintv[3] = lmax;
    /*-------------------- call kpmdos to get the DOS for dividing the spectrum*/
    /*-------------------- define kpmdos parameters */
    //-------------------- call kpmdos 
    double t = cheblan_timer();
    ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
    t = cheblan_timer() - t;
    if (ierr) {
      printf("kpmdos error %d\n", ierr);
      return 1;
    }
    fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);
    fprintf(fstats, " estimated eig count in interval: %.15e \n",ecount);
    //-------------------- call splicer to slice the spectrum
    npts = 10 * ecount; 
    sli = malloc((nslices+1)*sizeof(double));
    fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",
	    Mdeg, nvec, npts);
    ierr = spslicer(sli, mu, Mdeg, xintv, nslices,  npts);
    if (ierr) {
      printf("spslicer error %d\n", ierr);
      return 1;
    }
    printf("====================  SLICES FOUND  ====================\n");
    for (j=0; j<nslices; j++) {
      printf(" %2d: [% .15e , % .15e]\n", j+1, sli[j],sli[j+1]);
    }
    //-------------------- # eigs per slice
    ev_int = (int) (1 + ecount / ((double) nslices));
    totcnt = 0;
    //-------------------- For each slice call RatLanrNr
    for (sl=0; sl<nslices; sl++) {
      printf("======================================================\n");
      int nev2;
      double *lam, *Y, *res;
      int *ind;
      //-------------------- 
      a = sli[sl];
      b = sli[sl+1];
      printf(" subinterval: [%.4e , %.4e]\n", a, b); 
      double intv[4];
      intv[0] = a;
      intv[1] = b;
      intv[2] = lmin;
      intv[3] = lmax;
      // find the rational filter on this slice
      ratparams rat;
      // setup default parameters for rat
      set_ratf_def(&rat);
      // change some default parameters here:
      rat.pw = pow;
      rat.num = num;
      rat.beta = beta;
      // now determine rational filter
      find_ratf(intv, &rat);
      // use the solver function from UMFPACK
      void **solshiftdata = (void **) malloc(num*sizeof(void *));
      /*------------ factoring the shifted matrices and store the factors */
      SetupASIGMABSolSuiteSparse(&Acsr, &Bcsr, num, rat.zk, solshiftdata);
      /*------------ give the data to rat */
      SetASigmaBSol(&rat, NULL, ASIGMABSolSuiteSparse, solshiftdata);
      //-------------------- approximate number of eigenvalues wanted
      nev = ev_int+2;
      //-------------------- maximal Lanczos iterations   
      max_its = max(3*nev,100);  max_its = min(max_its, n);
      //-------------------- RationalLanNr
      ierr = RatLanNr(intv, max_its, tol, vinit, &rat, &nev2, &lam, 
		      &Y, &res, fstats);
      if (ierr) {
	printf("RatLanNr error %d\n", ierr);
	return 1;
      }

      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = (int *) malloc(nev2*sizeof(int));
      sort_double(nev2, lam, ind);
      printf(" number of eigenvalues found: %d\n", nev2);
      /* print eigenvalues */
      fprintf(fstats, "    Eigenvalues in [a, b]\n");
      fprintf(fstats, "    Computed [%d]        ||Res||\n", nev2);
      for (i=0; i<nev2; i++) {
	fprintf(fstats, "% .15e  %.1e\n", lam[i], res[ind[i]]);
      }
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
      memcpy(&alleigs[totcnt],lam,nev2*sizeof(double));
      totcnt += nev2;
      counts[sl] = nev2;      
      //-------------------- free allocated space withing this scope
      if (lam) free(lam);
      if (Y) free(Y);
      if (res) free(res);
      FreeASIGMABSolSuiteSparse(rat.num, solshiftdata);
      free(solshiftdata);
      free_rat(&rat);
      free(ind);
    } //for (sl=0; sl<nslices; sl++)
    //-------------------- free other allocated space 
    fprintf(fstats," --> Total eigenvalues found = %d\n",totcnt);
    sprintf(path, "OUT/EigsOut_RLanN_(%s, %s)", io.MatNam1, io.MatNam2);
    FILE *fmtout = fopen(path,"w");
    if (fmtout) {
      for (j=0; j<totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }    
    free(vinit);
    free(sli);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free_coo(&Bcoo);
    free_csr(&Bcsr);
    FreeBSolSuiteSparseData(&Bsol);
    free(alleigs);
    free(counts);
    if (fstats != stdout) fclose(fstats);
    /*-------------------- end matrix loop */
  }
  free(mu);
  if( flog != stdout) fclose( flog );
  fclose( fmat ); 
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

