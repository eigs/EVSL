#include "evsl.h"
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

/*-------------------- protos */
int exDOS(double *vals, int n, int npts, double *x, double *y, double *intv);

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
/*
 * Reads in a diag matrix.
 *
 * @parm[in,out] mat UNallocated cooMat to read into.
 * @parm[in] filename file to read from, where the first line contains number
 * of elements/width/height of matrix, and the rest of the lines contain the
 * values. */
int readDiagMat(const char *filename, cooMat *mat) {
  int width, i;
  FILE *ifp = fopen(filename, "r");
  if (ifp == NULL) {
    fprintf(stderr, "Can't open input file \n");
    exit(1);
  }
  fscanf(ifp, "%i", &width);
  // Setup cooMat
  mat->ncols = width;
  mat->nrows = width;
  mat->nnz = width;
  int *ir = (int *)malloc(sizeof(int) * width);
  int *jc = (int *)malloc(sizeof(int) * width);
  mat->ir = ir;
  mat->jc = jc;
  double *vv = (double *)malloc(sizeof(double) * width);
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
 * Tests landos.c -- Includes graphical comparison of exact DOS and calculated.
 *-----------------------------------------------------------------------
 */
int main() {
  srand(time(NULL));
  cooMat cooMat;
  csrMat csrMat;
  //-------------------- Read in a test matrix
  readDiagMat("testmat.dat", &cooMat);
  cooMat_to_csrMat(0, &cooMat, &csrMat);

  //-------------------- Define some constants to test with
  const int msteps = 40; // Steps to perform
  const int npts = 200;  // Number of points
  const int nvec = 100;  // Number of random vectors to generate
  double intv[6] = {-2.448170338612495,
                    11.868902203167497,
                    0,
                    0,
                    5,
                    8};    // Interval of interest
  const int degB = 1;      // Positive value to work with landosG
  const double tau = 1e-4; // Tolerance
  //-------------------- reset to whole interval
  intv[2] = intv[0];
  intv[3] = intv[1];
  int i, ret;
  double neig, neig2; // Number eigenvalues
  //-------------------- exact histogram and computed DOS
  double *xHist = (double *)malloc(npts * sizeof(double)); // Exact DOS x values
  double *yHist = (double *)malloc(npts * sizeof(double)); // Exact DOS y values
  double *xdos = (double *)malloc(npts * sizeof(double));  // Calculated DOS x
                                                           // vals
  double *ydos = (double *)malloc(npts * sizeof(double));  // Calculated DOS y
  double *xdos2 = (double *)malloc(npts * sizeof(double)); // Calculated DOS x
                                                           // vals
  double *ydos2 = (double *)malloc(npts * sizeof(double)); // Calculated DOS y
                                                           // vals

  SetStdEig();
  EVSLStart();
  SetAMatrix(&csrMat);
  double t0 = cheblan_timer();
  ret = LanDosG(nvec, msteps, degB, npts, xdos, ydos, &neig, intv,
                tau); // Calculate DOS
  double t1 = cheblan_timer();
  fprintf(stdout, " LanDos ret %d \n", ret);
  double t2 = cheblan_timer();
  ret = LanDos(nvec, msteps, npts, xdos2, ydos2, &neig2, intv); // Calculate DOS
  double t3 = cheblan_timer();
  fprintf(stdout, " LanDos ret %d \n", ret);
  printf("neig1: %f, neig2: %f \n", neig, neig2);
  printf("1: %f, 2: %f \n", t1 - t0, t3 - t2);

  ret = exDOS(cooMat.vv, cooMat.ncols, npts, xHist, yHist, intv); // Exact DOS
  EVSLFinish();
  free_coo(&cooMat);
  fprintf(stdout, " exDOS ret %d \n", ret);
  //--------------------Make OUT dir if it does'nt exist
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0750);
  }

  //-------------------- Write to  output files
  FILE *ofp = fopen("OUT/myydos.txt", "w");
  FILE *ofp2 = fopen("OUT/myydos2.txt", "w");
  for (i = 0; i < npts; i++) {
    fprintf(ofp, " %10.4f  %10.4f\n", xdos[i], ydos[i]);
    fprintf(ofp2, " %10.4f  %10.4f\n", xdos2[i], ydos2[i]);
  }
  fclose(ofp);
  fclose(ofp2);

  //-------------------- save exact DOS
  ofp = fopen("OUT/Exydos.txt", "w");
  for (i = 0; i < npts; i++)
    fprintf(ofp, " %10.4f  %10.4f\n", xHist[i], yHist[i]);
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
  free(xdos2);
  free(ydos2);
  free_csr(&csrMat);
  return 0;
}
