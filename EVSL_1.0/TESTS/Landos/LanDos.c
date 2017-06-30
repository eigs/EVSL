#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <time.h>
#include <unistd.h>
#include "evsl.h"
#include "io.h"

/*-------------------- protos */
int exDOS(double *vals, int n, int npts, double *x, double *y, double *intv);
int findarg(const char *argname, ARG_TYPE type, void *val, int argc,
            char **argv);

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

/**
 * Reads in a diag matrix.
 *
 * @param[in,out] mat UNallocated cooMat to read into.
 * @param[in] filename file to read from, where the first line contains number
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
  /* Setup cooMat */
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

/**
 *-----------------------------------------------------------------------
 * Tests landos.c -- Includes graphical comparison of exact DOS and calculated.
 *
 * use -graph_exact_dos 1 to enable graphing the exact DOS
 *-----------------------------------------------------------------------
 */
int main(int argc, char *argv[]) {
  int ierr;
  srand(time(NULL));
  cooMat cooMat;
  csrMat csrMat;

  int graph_exact_dos_tmp = 0;
  findarg("graph_exact_dos", INT, &graph_exact_dos_tmp, argc, argv);
  const int graph_exact_dos = graph_exact_dos_tmp;
  /*-------------------- Read in a test matrix */
  readDiagMat("testmat.dat", &cooMat);
  cooMat_to_csrMat(0, &cooMat, &csrMat);

  /*-------------------- Define some constants to test with */
  const int msteps = 40; /* Steps to perform */
  const int npts = 200;  /* Number of points */
  const int nvec = 100;  /* Number of random vectors to generate */
  double intv[6] = {-2.448170338612495,
                    11.868902203167497,
                    0,
                    0,
                    5,
                    8}; /* Interval of interest */
  /*-------------------- reset to whole interval */
  intv[2] = intv[0];
  intv[3] = intv[1];
  int i, ret;
  double neig; /* Number eigenvalues */
  /*-------------------- exact histogram and computed DOS */
  double *xHist;
  double *yHist;
  if (graph_exact_dos) {
    xHist = (double *)malloc(npts * sizeof(double)); /*Exact DOS x values */

    yHist = (double *)malloc(npts * sizeof(double)); /* Exact DOS y values */
  }
  double *xdos =
      (double *)malloc(npts * sizeof(double)); /* Calculated DOS x vals */
  double *ydos = (double *)malloc(npts * sizeof(double)); /* Calculated DOS y */

  SetStdEig();
  EVSLStart();
  SetAMatrix(&csrMat);
  /* Calculate computed DOS */
  ret = LanDos(nvec, msteps, npts, xdos, ydos, &neig, intv);
  fprintf(stdout, " LanDos ret %d \n", ret);

  /* Calculate the exact DOS */
  if (graph_exact_dos) {
    ret = exDOS(cooMat.vv, cooMat.ncols, npts, xHist, yHist, intv);
    fprintf(stdout, " exDOS ret %d \n", ret);
  }
  EVSLFinish();
  free_coo(&cooMat);
  /* --------------------Make OUT dir if it doesn't exist */
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0750);
  }

  /*-------------------- Write to  output files */
  FILE *ofp = fopen("OUT/LanDos_Approx_DOS.txt", "w");
  for (i = 0; i < npts; i++) {
    fprintf(ofp, " %10.4f  %10.4f\n", xdos[i], ydos[i]);
  }
  fclose(ofp);

  if (graph_exact_dos) {
    /*-------------------- save exact DOS */
    ofp = fopen("OUT/LanDos_Exact_DOS.txt", "w");
    for (i = 0; i < npts; i++)
      fprintf(ofp, " %10.4f  %10.4f\n", xHist[i], yHist[i]);
    fclose(ofp);
  }
  printf("The data output is located in  OUT/ \n");
  struct utsname buffer;
  errno = 0;
  if (uname(&buffer) != 0) {
    perror("uname");
    exit(EXIT_FAILURE);
  }

  /*-------------------- invoke gnuplot for plotting ... */
  if (graph_exact_dos) {
    ierr = system("gnuplot < tester_ex.gnuplot");
  } else {
    ierr = system("gnuplot < tester.gnuplot");
  }

  if (ierr) {
    fprintf(stderr,
            "Error using 'gnuplot < tester.gnuplot', \n"
            "postscript plot could not be generated \n");
  } else {
    printf("A postscript graph has been placed in OUT/tester.eps \n");
    /*-------------------- and gv for visualizing */
    if (!strcmp(buffer.sysname, "Linux")) {
      ierr = system("gv OUT/tester.eps");
      if (ierr) {
        fprintf(stderr, "Error using 'gv OUT/tester.eps' \n");
      }

    } else {
      printf(
          "To view the postscript graph use a postcript viewer such as  "
          "gv \n");
    }
  }
  if (graph_exact_dos) {
    free(xHist);
    free(yHist);
  }
  free(xdos);
  free(ydos);
  return 0;
}
