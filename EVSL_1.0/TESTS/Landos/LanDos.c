#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "evsl.h"

/*-------------------- protos */
int exDOS(double *vals, int n, int npts, 
	  double *x, double *y, double *intv);

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
 * Tests landos.c -- onlu 
 * The following variable are the outputs. The noted values are the values
 * when the randn_double in landos is replaced with a vector of ones.
 *
 *
 *
 */
int main() {
  cooMat cooMat;
  csrMat csrMat;
  // Read in a test matrix
  readDiagMat("testmat.dat", &cooMat);
  cooMat_to_csrMat(0, &cooMat, &csrMat);
 
  // Define some constants to test with
  const int msteps = 30;
  const int npts = 200;
  const int nvec = 100;
  double intv[4] = {-2.448170338612495, 11.868902203167497, 5, 8};
  intv[2] = intv[0];
  intv[3] = intv[1];
  int i;

  double* xHist = (double*)calloc(npts, sizeof(double));
  double* yHist = (double*)calloc(npts, sizeof(double));
  double* xdos = (double*)calloc(npts, sizeof(double));
  double* ydos = (double*)calloc(npts, sizeof(double));

  double neig;  
  int ret;

  SetAMatrix(&csrMat);
  ret = LanDos(nvec, msteps, npts, xdos, ydos, &neig, intv);
  fprintf(stdout, " LanDos ret %d \n",ret) ;
  
  ret = exDOS(cooMat.vv, cooMat.ncols, npts, xHist, yHist, intv) ; 
  free_coo(&cooMat);
  fprintf(stdout, " exDOS ret %d \n",ret) ;

  //Make OUT dir if it does'nt exist
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
	  mkdir("OUT", 0700);
  }

  // Write to an output file
  FILE* ofp = fopen("OUT/myydos.txt", "w");
  for (i = 0; i < npts; i++) 
    fprintf(ofp," %10.4f  %10.4f\n",xdos[i],ydos[i]);
  fclose(ofp);


  // Write to an output file
  ofp = fopen("OUT/Exydos.txt", "w");
  for (i = 0; i < npts; i++) 
    fprintf(ofp," %10.4f  %10.4f\n",xHist[i],yHist[i]);
  fclose(ofp);

  system("gnuplot < tester.gnuplot");
  system("gv tester.eps");
  

  free(xHist);
  free(yHist);
  free(xdos);
  free(ydos);

  free_csr(&csrMat);

  return 0;
}
