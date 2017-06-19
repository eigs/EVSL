#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "evsl.h"
#include "io.h"

#define BsolPol 1

/*-------------------- protos */
int exDOS(double* vals, int n, int npts, double* x, double* y, double* intv);
int read_coo_MM(const char* matfile, int idxin, int idxout, cooMat* Acoo);
int get_matrix_info(FILE* fmat, io_t* pio);

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
/*
 * Reads in a vector as an nx1 matrix.
 *
 * @parm[out] npts pointer to an int to store # of points
 * @parm[out] vec UNallocated space to read vector to
 * @parm[in] filename file to read from, where the first line contains number
 * of elements/width/height of matrix, and the rest of the lines contain the
 * values. */
int readVec(const char* filename, int* npts, double** vec) {
  int i;
  FILE* ifp = fopen(filename, "r");
  fscanf(ifp, "%i", npts);
  *vec = (double*)malloc(sizeof(double) * *npts);
  for (i = 0; i < (*npts); i++) {
    fscanf(ifp, "%lf", (&(*vec)[i]));
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
  srand(time(NULL));
  const int msteps = 30;    // Number of steps
  const int degB = 40;      // Degree to aproximate B with
  const int npts = 200;     // Number of points
  const int nvec = 30;      // Number of random vectors to use
  const double tau = 1e-4;  // Tolerance in polynomial approximation
  // ---------------- Intervals of interest
  // intv[0] and intv[1] are the smallest and largest eigenvalues of (A,B)
  // intv[2] and intv[3] are the input interval of interest [a. b]
  // intv[4] and intv[5] are the smallest and largest eigenvalues of B after
  // diagonal scaling
  double intv[6] = {-2.739543872224533e-13,
                     0.0325,
                    -2.739543872224533e-13,
                     0.0325,
                     0.5479,
                     2.5000};
  int n = 0, i, nslices, ierr;
  double a, b;

  cooMat Acoo, Bcoo;  // A, B
  csrMat Acsr, Bcsr;  // A, B
  double* sqrtdiag = NULL;

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
    fprintf(fstats, "Partition the interval of interest [%f,%f] into %d slices\n", a, b,
            nslices);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr = read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        // nnz = Acoo.nnz;
        n = Acoo.nrows;
        if(n <= 0) {
          fprintf(stderr, "non-positive number of rows");
          exit(7);
        }


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
      /*------------------ diagonal scaling for Acoo and Bcoo */
      sqrtdiag = (double*)calloc(n, sizeof(double));
      extractDiag(&Bcoo, sqrtdiag);
      diagScaling(&Acoo, &Bcoo, sqrtdiag);
      // save_vec(n, diag, "OUT/diag.txt");
      /*-------------------- conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
      if(ierr) {
        fprintf(flog, "Could not convert matrix to csr: %i", ierr);
        exit(8);
      }
      ierr = cooMat_to_csrMat(0, &Bcoo, &Bcsr);
      if(ierr) {
        fprintf(flog, "Could not convert matrix to csr: %i", ierr);
        exit(9);
      }
    }
    if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT  not supported (yet) * \n");
      exit(7);
    }
    if (1) {
      /*----------------  compute the range of the spectrum of B */
      SetStdEig();
      SetAMatrix(&Bcsr);
      double* vinit = (double*)malloc(n * sizeof(double));
      rand_double(n, vinit);
      double lmin = 0.0, lmax = 0.0;
      ierr = LanTrbounds(50, 200, 1e-8, vinit, 1, &lmin, &lmax, fstats);
      SetGenEig();
      if(ierr) {
        fprintf(flog, "Could not run LanTrBounds: %i", ierr);
        exit(10);
      }
      /*------------- get the bounds for B ------*/
      intv[4] = lmin;
      intv[5] = lmax;
      /*-------------------- set the left-hand side matrix A */
      SetAMatrix(&Acsr);
      /*-------------------- set the right-hand side matrix B */
      SetBMatrix(&Bcsr);
#if BsolPol
      /*-------------------- Use polynomial to solve B */
      BSolDataPol Bsol;
      Bsol.intv[0] = lmin;
      Bsol.intv[1] = lmax;    
      SetupBSolPol(&Bcsr, &Bsol);
      SetBSol(BSolPol, (void *) &Bsol);
#else
      /*-------------------- Use Choleksy factorization to solve B */
      BSolDataSuiteSparse Bsol;
      /*-------------------- use SuiteSparse as the solver for B */
      SetupBSolSuiteSparse(&Bcsr, &Bsol);
      /*-------------------- set the solver for B */
      SetBSol(BSolSuiteSparse, (void *) &Bsol);
#endif
      SetGenEig();
      rand_double(n, vinit);
      ierr = LanTrbounds(40, 200, 1e-10, vinit, 1, &lmin, &lmax, fstats);
      if(ierr) {
        fprintf(flog, "Could not run LanTrBounds: %i", ierr);
        exit(10);
      }
      free(vinit);
#if BsolPol
      FreeBSolPolData(&Bsol);
#else
      FreeBSolSuiteSparseData(&Bsol);
#endif
      /*----------------- get the bounds for (A, B) ---------*/
      intv[0] = lmin;
      intv[1] = lmax;
      /*----------------- plotting the DOS on [a, b] ---------*/
      intv[2] = lmin;
      intv[3] = lmax;
    } else {
      /*-------------------- set the left-hand side matrix A */
      SetAMatrix(&Acsr);
      /*-------------------- set the right-hand side matrix B */
      SetBMatrix(&Bcsr);
      SetGenEig();
    }
    // printf("%lf, %lf, %lf, %lf \n", lmin, lmax, intv[0], intv[1]);
    //-------------------- Read in the eigenvalues
    double* ev;
    int numev;
    readVec("NM1AB_eigenvalues.dat", &numev, &ev);

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

    double t0 = cheblan_timer();
    ret = LanDosG(nvec, msteps, degB, npts, xdos, ydos, &neig, intv, tau);
    double t1 = cheblan_timer();
    fprintf(stdout, " LanDos ret %d  in %0.04fs\n", ret, t1 - t0);

    // -------------------- Calculate the exact DOS
    ret = exDOS(ev, numev, npts, xHist, yHist, intv);


    free_coo(&Acoo);
    free_coo(&Bcoo);
    free_csr(&Acsr);
    free_csr(&Bcsr);
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
    ierr = system("gnuplot < testerG.gnuplot");
    if(ierr) {
      printf("Could not run gnuplot: %i", ierr);
    }
    //-------------------- and gv for visualizing /
    ierr = system("gv testerG.eps");
    if(ierr) {
      printf("Could not run gv: %i", ierr);
    }
    free(xHist);
    free(yHist);
    free(xdos);
    free(ydos);
    free(ev);
    fclose(fstats);
    if (sqrtdiag) free(sqrtdiag);
  }
  fclose(fmat);
  EVSLFinish();
  //-------------------- done
  return 0;
}
