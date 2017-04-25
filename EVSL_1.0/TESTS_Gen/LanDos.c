#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "evsl.h"
#include "io.h"

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
  int width;
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

  for (int i = 0; i < width; i++) {
    mat->ir[i] = i;
    mat->jc[i] = i;
    fscanf(ifp, "%lf", &mat->vv[i]);
  }
  fclose(ifp);
  return 0;
}
int main() {
  cooMat cooMat;
  csrMat csrMat;
  readDiagMat("testmat.dat", &cooMat);
  cooMat_to_csrMat(0, &cooMat, &csrMat);
  free_coo(&cooMat);
  const int msteps = 30;
  const int npts = 200;
  const int nvec = 100;
  const double intv[4] = {-2.448170338612495, 11.868902203167497, 5, 8};

  double* xdos = (double*)calloc(npts, sizeof(double));
  double* ydos = (double*)calloc(npts, sizeof(double));

  int ret = LanDos(&csrMat, nvec, msteps, npts, xdos, ydos, intv);

  printf("xdos:\n");
  for (int i = 0; i < npts; i++) {
    printf("%lf\t", xdos[i]);
  }
  printf("\n");
  printf("ydos:\n");
  for (int i = 0; i < npts; i++) {
    printf("%lf\t", ydos[i]);
  }
  FILE* ofp = fopen("myydos.txt", "w");
  for (int i = 0; i < npts; i++) {
    fprintf(ofp, "%lf\n", ydos[i]);
  }
  double* si = (double*)calloc(npts, sizeof(double));
  simpson(xdos, ydos, npts, si);
  free(xdos);
  free(ydos);
  free_csr(&csrMat);

  return 0;
}
