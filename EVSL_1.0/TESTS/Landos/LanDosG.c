#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "evsl.h"
#include "io.h"

/*-------------------- protos */
int exDOS(double* vals, int n, int npts, double* x, double* y, double* intv);
int read_coo_MM(const char* matfile, int idxin, int idxout, cooMat* Acoo);
int get_matrix_info(FILE* fmat, io_t* pio);

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
/*
 * Reads in a vector as an nx1 matrix.
 *
 * @parm[in,out] mat UNallocated cooMat to read into.
 * @parm[in] filename file to read from, where the first line contains number
 * of elements/width/height of matrix, and the rest of the lines contain the
 * values. */
int readVec(const char* filename, cooMat* mat) {
  int numEles, i;
  FILE* ifp = fopen(filename, "r");
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file \n");
    exit(1);
  }
  fscanf(ifp, "%i", &numEles);
  // Setup cooMat
  mat->ncols = numEles;
  mat->nrows = 1;
  mat->nnz = numEles;
  int* ir = (int*)malloc(sizeof(int) * numEles);
  int* jc = (int*)malloc(sizeof(int) * numEles);
  mat->ir = ir;
  mat->jc = jc;
  double* vv = (double*)malloc(sizeof(double) * numEles);
  mat->vv = vv;
  for (i = 0; i < numEles; i++) {
    mat->ir[i] = 1;
    mat->jc[i] = i;
    fscanf(ifp, "%lf", &mat->vv[i]);
  }
  fclose(ifp);
  return 0;
}

/*
 *-----------------------------------------------------------------------
 * Tests landosG.c , the Lanczos DOS approximate for the general eigenvalue
 * problem. Includes graphical comparison of calculated vs exact DOS
 *-----------------------------------------------------------------------
 */
int main() {
  const int msteps = 30;    // Number of steps
  const int degB = 20;      // Degree to aproximate B with
  const int npts = 200;     // Number of points
  const int nvec = 50;      // Number of random vectors to use
  const double tau = 1e-4;  // Tolerance in polynomial approximation
  // ---------------- Intervals of interest
  double intv[6] = {-2.739543872224533e-13,
                    0.0325,
                    -2.739543872224533e-13,
                    0.0325,
                    0.5479,
                    2.5000};
  int n = 0, i, nslices, ierr;
  double a, b;

  cooMat Acoo, Bcoo, cooMat;  // A, B, and eigenvalue matrices in COO form
  csrMat Acsr, Bcsr, csrMat;  // A, B, and eigenvalue matrices in CSR form

  FILE *flog = stdout, *fmat = NULL;
  FILE* fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  /*-------------------- start EVSL */
  EVSLStart();
  /*------------------ file "matfile" contains paths to matrices */
  if (NULL == (fmat = fopen("matfile", "r"))) {
    fprintf(flog, "Can't open matfile...\n");
    exit(2);
  } else {
    printf("Read matfile \n");
  }
  /*-------------------- read number of matrices ..*/
  memset(line, 0, MAX_LINE);
  if (NULL == fgets(line, MAX_LINE, fmat)) {
    fprintf(flog, "error in reading matfile...\n");
    exit(2);
  }
  if ((numat = atoi(line)) <= 0) {
    fprintf(flog, "Invalid count of matrices...\n");
    exit(3);
  }
  for (mat = 1; mat <= numat; mat++) {
    if (get_matrix_info(fmat, &io) != 0) {
      fprintf(flog, "Invalid format in matfile ...\n");
      exit(5);
    }
    /*----------------input matrix and interval information -*/
    fprintf(flog, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(flog, "MATRIX B: %s...\n", io.MatNam2);
    a = io.a;  // left endpoint of input interval
    b = io.b;  // right endpoint of input interval
    nslices = io.n_intv;
    char path[1024];  // path to write the output files
    strcpy(path, "OUT/MMPLanR_");
    strcat(path, io.MatNam1);
    fstats = fopen(path, "w");  // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(fstats, "MATRIX B: %s...\n", io.MatNam2);
    fprintf(fstats,
            "Partition the interval of interest [%f,%f] into %d slices\n", a, b,
            nslices);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr = read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        // nnz = Acoo.nnz;
        n = Acoo.nrows;
        // printf("size of A is %d\n", n);
        // printf("nnz of  A is %d\n", nnz);
      } else {
        fprintf(flog, "read_coo error for A = %d\n", ierr);
        exit(6);
      }
      ierr = read_coo_MM(io.Fname2, 1, 0, &Bcoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        if (Bcoo.nrows != n) {
          return 1;
        }
        // nnz = Bcoo.nnz;
        // n = Bcoo.nrows;
        // printf("size of B is %d\n", n);
        // printf("nnz of  B is %d\n", nnz);
      } else {
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
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- set the right-hand side matrix B */
    SetBMatrix(&Bcsr);

    SetGenEig();
  }

  //-------------------- Read in the eigenvalues
  readVec("ev.dat", &cooMat);
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

  // ------------------- Calculate the approximate DOS
  ret = LanDosG(nvec, msteps, degB, npts, xdos, ydos, &neig, intv, tau);
  fprintf(stdout, " LanDos ret %d \n", ret);

  // -------------------- Calculate the exact DOS
  ret = exDOS(cooMat.vv, cooMat.ncols, npts, xHist, yHist, intv);

  free_coo(&cooMat);
  free_coo(&Acoo);
  free_coo(&Bcoo);
  fprintf(stdout, " exDOS ret %d \n", ret);

  //--------------------Make OUT dir if it doesn't exist
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0700);
  }

  //-------------------- Write to  output files
  FILE* ofp = fopen("OUT/myydosG.txt", "w");
  for (i = 0; i < npts; i++)
    fprintf(ofp, " %10.4f  %10.4f\n", xdos[i], ydos[i]);
  fclose(ofp);

  //-------------------- save exact DOS
  ofp = fopen("OUT/ExydosG.txt", "w");
  for (i = 0; i < npts; i++)
    fprintf(ofp, " %10.4f  %10.4f\n", xHist[i], yHist[i]);
  fclose(ofp);
  //-------------------- invoke gnuplot for plotting ...
  system("gnuplot < testerG.gnuplot");
  //-------------------- and gv for visualizing /
  system("gv testerG.eps");
  //-------------------- done
  free(xHist);
  free(yHist);
  free(xdos);
  free(ydos);
  free_csr(&csrMat);
  free_csr(&Acsr);
  free_csr(&Bcsr);
  return 0;
}
