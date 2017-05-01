#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "evsl.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define TRIV_SLICER 0
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
 * Tests landos.c, simposon.c and spslicer2.c.
 * The following variable are the outputs. The noted values are the values
 * when the randn_double in landos is replaced with a vector of ones.
 *
 *
 * xdos = [5, 5.015075, 5.030151, 5.045226 ... 7.964925, 8.000000 ]
 * ydos = [0, 0.062025, 0.062213, 0.062406 ... 0.100559, 0.100404 ]
 * si   = [0, 0.000936, 0.001876, 0.002818 ... 0.254688, 0.256203 ]
 * sli  = [5, 5.949749, 6.716593, 6.718593 ]
 * neig = 153.721918
 *
 *
 */
int main() {
  cooMat cooMat;
  csrMat csrMat;

  // Read in a test matrix
  readDiagMat("testmat.dat", &cooMat);
  cooMat_to_csrMat(0, &cooMat, &csrMat);
  free_coo(&cooMat);

  // Define some constants to test with
  const int msteps = 30;
  const int npts = 200;
  const int nvec = 100;
  const double intv[4] = {-2.448170338612495, 11.868902203167497, 5, 8};
  
  int i;

  double* xdos = (double*)calloc(npts, sizeof(double));
  double* ydos = (double*)calloc(npts, sizeof(double));

  double neig;  
  int ret;

  SetAMatrix(&csrMat);
  ret = LanDos(nvec, msteps, npts, xdos, ydos, &neig, intv);
  fprintf(stdout, " ret %d \n",ret) ;

  double* si = (double*)calloc(npts, sizeof(double));
  simpson(xdos, ydos, npts, si);

  const int n_int = 4;
  double* sli = (double*)calloc(n_int, sizeof(double));
  spslicer2(xdos, si, n_int, npts, sli);

  //Make OUT dir if it does'nt exist
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
	  mkdir("OUT", 0700);
  }

  // Write to an output file
  FILE* ofp = fopen("OUT/myydos.txt", "w");
  for (i = 0; i < npts; i++) {
    fprintf(ofp, "%lf\n", ydos[i]);
  }
  fclose(ofp);


  ofp = fopen("OUT/myxdos.txt", "w");
  for (i = 0; i < npts; i++) {
    fprintf(ofp, "%lf\n", xdos[i]);
  }
  fclose(ofp);

  ofp = fopen("OUT/mysi.txt", "w");
  for (i = 0; i < npts; i++) {
    fprintf(ofp, "%lf\n", si[i]);
  }
  fclose(ofp);

  ofp = fopen("OUT/mysli.txt", "w");
  for (i = 0; i < n_int; i++) {
    fprintf(ofp, "%lf\n", sli[i]);
  }
  fclose(ofp);

  ofp = fopen("OUT/myneig.txt", "w");
  fprintf(ofp, "%lf", neig);
  fclose(ofp);

  free(si);
  free(sli);
  free(xdos);
  free(ydos);
  free_csr(&csrMat);

  return 0;
}
