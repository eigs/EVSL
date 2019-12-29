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
#include "lapl.h"

int exDOS(double *vals, int n, int npts, double *x, double *y, double *intv);

/**
 *-----------------------------------------------------------------------
 * Tests landos.c -- Includes graphical comparison of exact DOS and calculated.
 *
 * use -graph_exact_dos 0 to disable graphing the exact DOS
 *-----------------------------------------------------------------------
 */
int main(int argc, char *argv[]) {
  srand(time(NULL));
  const int msteps = 40; /* Steps to perform */
  const int npts = 200;  /* Number of points */
  const int nvec = 100;  /* Number of random vectors to generate */
  double intv[4];
  const int num = 20;

  cooMat Acoo;
  csrMat Acsr;
#ifdef EVSL_USING_CUDA_GPU
  csrMat Acsr_gpu;
#endif

  /*-------------------- IO */
  FILE *fstats = NULL;
  io_t io;
  int ierr, n = 0, graph_exact_dos = 1;

  /* initial vector: random */
  /*-------------------- start EVSL */
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
  EVSLStart();
  /*  Generate the Laplacian matrix (which we can trivially get eigenvalues for*/
  lapgen(num, num, num, &Acoo);
  cooMat_to_csrMat(0, &Acoo, &Acsr);
#ifdef EVSL_USING_CUDA_GPU
  evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
#endif

  // /*-------------------- path to write the output files*/
  char path[1024];
  strcpy(path, "OUT/LanDos_lap.log");
  fstats = fopen(path, "w");  // write all the output to the file io.MatNam
  if (!fstats) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  n = Acoo.nrows;
  /*-------------------- Define some constants to test with */
  /*-------------------- Intervals of interest
    intv[0] and intv[1] are the input interval of interest [a. b]
    intv[2] and intv[3] are the smallest and largest eigenvalues of (A,B)
    */
  int i, ret;
  double neig; /* Number eigenvalues */
  /*-------------------- exact histogram and computed DOS */
  double *xHist = NULL;
  double *yHist = NULL;
  if (graph_exact_dos) {
    xHist = evsl_Malloc(npts, double); /*Exact DOS x values */
    yHist = evsl_Malloc(npts, double); /* Exact DOS y values */
  }
  double *xdos = evsl_Malloc(npts, double); /* Calculated DOS x vals */
  double *ydos = evsl_Malloc(npts, double); /* Calculated DOS y */

  SetStdEig();
#ifdef EVSL_USING_CUDA_GPU
  SetAMatrix_device_csr(&Acsr_gpu);
#else
  SetAMatrix(&Acsr);
#endif

  double *vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);
  double lmin = 0.0, lmax = 0.0;
  ierr = LanTrbounds(50, 200, 1e-8, vinit, 1, &lmin, &lmax, fstats);

  intv[0] = lmin;
  intv[1] = lmax;
  /*----------------- get the bounds for (A, B) ---------*/
  intv[2] = lmin;
  intv[3] = lmax;

  double *ev = NULL;
  int numev;
  if (graph_exact_dos) {
    /*  Calculate the exact eigenvalues */
    exeiglap3(num, num, num, lmin, lmax, &numev, &ev);
  }

  /* Calculate computed DOS */
  ret = LanDos(nvec, msteps, npts, xdos, ydos, &neig, intv);
  fprintf(stdout, " LanDos ret %d \n", ret);

  /* Calculate the exact DOS */
  if (graph_exact_dos) {
    ret = exDOS(ev, numev, npts, xHist, yHist, intv);
    fprintf(stdout, " exDOS ret %d \n", ret);
  }
  free_coo(&Acoo);
  free_csr(&Acsr);
#ifdef EVSL_USING_CUDA_GPU
  evsl_free_csr_gpu(&Acsr_gpu);
#endif
  evsl_Free_device(vinit);
  /* --------------------Make OUT dir if it doesn't exist */
  struct stat st = {0};
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0750);
  }

  /*-------------------- Write to  output files */
  char computed_path[1024];
  strcpy(computed_path, "OUT/LanDos_Approx_DOS_");
  strcat(computed_path, io.MatNam1);
  FILE *ofp = fopen(computed_path, "w");
  for (i = 0; i < npts; i++) {
    fprintf(ofp, " %10.4f  %10.4f\n", xdos[i], ydos[i]);
  }
  fclose(ofp);

  if (graph_exact_dos) {
    /*-------------------- save exact DOS */
    strcpy(path, "OUT/LanDos_Exact_DOS_");
    strcat(path, io.MatNam1);
    ofp = fopen(path, "w");
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
  char command[1024];
  strcpy(command, "gnuplot ");
  strcat(command, " -e \"filename='");
  strcat(command, computed_path);

  if (graph_exact_dos) {
    strcat(command, "';exact_dos='");
    strcat(command, path);
    strcat(command, "'\" tester_ex.gnuplot");
    ierr = system(command);
  } else {
    strcat(command, "'\" tester.gnuplot");
    ierr = system(command);
  }

  if (ierr) {
    fprintf(stderr,
            "Error using 'gnuplot < tester.gnuplot', \n"
            "postscript plot could not be generated \n");
  } else {
    printf("A postscript graph has been placed in %s%s\n", computed_path,
           ".eps");
    /*-------------------- and gv for visualizing / */
    if (!strcmp(buffer.sysname, "Linux")) {
      strcpy(command, "gv ");
      strcat(command, computed_path);
      strcat(command, ".eps &");
      ierr = system(command);
      if (ierr) {
        fprintf(stderr, "Error using 'gv %s' \n", command);
        printf(
            "To view the postscript graph use a postcript viewer such as  "
            "gv \n");
      }
    } else {
      printf(
          "To view the postscript graph use a postcript viewer such as  "
          "gv \n");
    }
  }

  if (graph_exact_dos) {
    evsl_Free(xHist);
    evsl_Free(yHist);
  }
  evsl_Free(xdos);
  evsl_Free(ydos);
  if (ev) {
    evsl_Free(ev);
  }
  fclose(fstats);
  EVSLFinish();
  return 0;
}

