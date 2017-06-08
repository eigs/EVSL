#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "evsl.h"
#include "io.h"

/*-------------------- protos */
int exDOS(double *vals, int n, int npts, 
	  double *x, double *y, double *intv);
int read_coo_MM(const char *matfile, int idxin, int idxout, cooMat *Acoo);
int get_matrix_info(FILE *fmat, io_t *pio);


#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
/*
 * Reads in a diag matrix.
 *
 * @parm[in,out] mat UNallocated cooMat to read into.
 * @parm[in] filename file to read from, where the first line contains number
 * of elements/width/height of matrix, and the rest of the lines contain the
 * values. */
int readDiagMat(const char* filename, cooMat* mat) {
  int width, i; 
  FILE* ifp = fopen(filename, "r");
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file \n");
    exit(1);
  }
  fscanf(ifp, "%i", &width);
  // Setup cooMat
  mat->ncols = width;
  mat->nrows = width;
  mat->nnz = width;
  int* ir = (int*)malloc(sizeof(int) * width);
  int* jc = (int*)malloc(sizeof(int) * width);
  mat->ir = ir;
  mat->jc = jc;
  double* vv = (double*)malloc(sizeof(double) * width);
  mat->vv = vv;

  for (i = 0; i < width; i++) {
    mat->ir[i] = i;
    mat->jc[i] = i;
    fscanf(ifp, "%lf", &mat->vv[i]);
  }
  fclose(ifp);
  return 0;
}

/*
 *-----------------------------------------------------------------------
 * Tests landos.c -- Only the  following variable are the outputs. The
 * noted  values are  the values  when the  randn_double in  landos is
 * replaced with a vector of ones.
 *-----------------------------------------------------------------------
 */
int main() {
  const int msteps = 30;
  const int degB = 20;
  const int npts = 200;
  const int nvec = 50;
  const double tau = 1e-4;
  double intv[6] = {-2.739543872224533e-13, 0.0325, -2.739543872224533e-13, 0.0325, 0.5479, 2.5000};
  int n=0, i, j,  nslices, Mdeg, nev, 
      mlan, max_its, ev_int, sl, ierr, totcnt;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli, *mu;
  double xintv[4];
  double *alleigs; 
  int *counts; // #ev computed in each slice  
  /* initial vector: random */
  double *vinit;
  polparams pol;
  cooMat Acoo, Bcoo, cooMat;
  csrMat Acsr, Bcsr, csrMat;

  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  /*-------------------- Bsol */
  /*-------------------- stopping tol */
  tol = 1e-8;
  /*-------------------- start EVSL */
  EVSLStart();
  /*------------------ file "matfile" contains paths to matrices */
  if( NULL == ( fmat = fopen( "matfile", "r" ) ) ) {
    fprintf( flog, "Can't open matfile...\n" );
    exit(2);
  }
  else {
    printf("Read matfile \n");
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
    strcpy( path, "OUT/MMPLanR_");
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
      ierr = read_coo_MM(io.Fname1, 1, 0, &Acoo); 
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
        //nnz = Acoo.nnz; 
        n = Acoo.nrows;
        //printf("size of A is %d\n", n);
        //printf("nnz of  A is %d\n", nnz);
      }
      else {
        fprintf(flog, "read_coo error for A = %d\n", ierr);
        exit(6);
      }
      ierr = read_coo_MM(io.Fname2, 1, 0, &Bcoo); 
      if (ierr == 0) {
        fprintf(fstats,"matrix read successfully\n");
        if (Bcoo.nrows != n) {
          return 1;
        }
        //nnz = Bcoo.nnz; 
        //n = Bcoo.nrows;
        //printf("size of B is %d\n", n);
        //printf("nnz of  B is %d\n", nnz);
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

    SetGenEig();
  }

  //-------------------- Read in a test matrix
  readDiagMat("testmat.dat", &cooMat);
  cooMat_to_csrMat(0, &cooMat, &csrMat);

  //-------------------- Define some constants to test with
  //-------------------- reset to whole interval
  int ret;
  double neig;
  //-------------------- exact histogram and computed DOS
  double* xHist = (double*)calloc(npts, sizeof(double));
  double* yHist = (double*)calloc(npts, sizeof(double));
  double* xdos = (double*)calloc(npts, sizeof(double));
  double* ydos = (double*)calloc(npts, sizeof(double));

  ret = LanDosG2(nvec, msteps, degB,  npts, xdos, ydos, &neig, intv, tau);
  fprintf(stdout, " LanDos ret %d \n",ret) ;
  
  ret = exDOS(cooMat.vv, cooMat.ncols, npts, xHist, yHist, intv) ; 
  free_coo(&cooMat);
  fprintf(stdout, " exDOS ret %d \n",ret) ;

  //--------------------Make OUT dir if it does'nt exist
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
	  mkdir("OUT", 0700);
  }

  //-------------------- Write to  output files
  FILE* ofp = fopen("OUT/myydos.txt", "w");
  for (i = 0; i < npts; i++) 
    fprintf(ofp," %10.4f  %10.4f\n",xdos[i],ydos[i]);
  fclose(ofp);

  //-------------------- save exact DOS 
  ofp = fopen("OUT/Exydos.txt", "w");
  for (i = 0; i < npts; i++) 
    fprintf(ofp," %10.4f  %10.4f\n",xHist[i],yHist[i]);
  fclose(ofp);
  //-------------------- invoke gnuplot for plotting ...
  system("gnuplot < tester.gnuplot");
  //-------------------- and gv for visualizing /
  system("gv tester.eps");
  //-------------------- done 
  free(xHist);
  free(yHist);
  free(xdos);
  free(ydos);
  free_csr(&csrMat);
  return 0;
}