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
 * Reads in a vector as an nx1 matrix.
 *
 * @param[out] npts pointer to an int to store # of points
 * @param[out] vec UNallocated space to read vector to
 * @param[in] filename file to read from, where the first line contains number
 * of elements/width/height of matrix, and the rest of the lines contain the
 * values. */
int readVec(const char *filename, int *npts, double **vec) {
  int i;
  FILE *ifp = fopen(filename, "r");
  if (ifp) {
    fscanf(ifp, "%i", npts);
    *vec = (double *)malloc(sizeof(double) * *npts);
    for (i = 0; i < (*npts); i++) {
      fscanf(ifp, "%lf", (&(*vec)[i]));
    }
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
  srand(time(NULL));
  const int msteps = 40; /* Steps to perform */
  const int npts = 200;  /* Number of points */
  const int nvec = 100;  /* Number of random vectors to generate */
  double intv[4];

  cooMat Acoo;
  csrMat Acsr;

  /*-------------------- IO */
  FILE *flog = stdout, *fmat = NULL, *fstats = NULL;
  io_t io;
  int numat, mat, ierr, n, graph_exact_dos = 0;
  char line[MAX_LINE];

  findarg("graph_exact_dos", INT, &graph_exact_dos, argc, argv);
  /* initial vector: random */
  /*-------------------- start EVSL */
  EVSLStart();
  //-------------------- interior eigensolver parameters
  /*------------------ file "matfile" contains paths to matrices */
  if (NULL == (fmat = fopen("matfile", "r"))) {
    fprintf(flog, "Can't open matfile...\n");
    exit(2);
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
  /*-------------------- LOOP through matrices -*/
  for (mat = 1; mat <= numat; mat++) {
    if (get_matrix_info(fmat, &io) != 0) {
      fprintf(flog, "Invalid format in matfile ...\n");
      exit(5);
    }
    /*----------------input matrix and interval information -*/
    fprintf(flog, "MATRIX: %s...\n", io.MatNam);

    /*-------------------- path to write the output files*/
    char path[1024];
    strcpy(path, "OUT/LanDos");
    strcat(path, io.MatNam);
    strcat(path, ".log");
    fstats = fopen(path, "w");  // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX: %s...\n", io.MatNam);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr = read_coo_MM(io.Fname, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        // nnz = Acoo.nnz;
        n = Acoo.nrows;
      } else {
        fprintf(flog, "read_coo error = %d\n", ierr);
        exit(6);
      }
      /*-------------------- conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
    }
    if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT  not supported (yet) * \n");
      exit(7);
    }
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- Read in a test matrix */

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
      xHist = (double *)malloc(npts * sizeof(double)); /*Exact DOS x values */
      yHist = (double *)malloc(npts * sizeof(double)); /* Exact DOS y values */
    }
    double *xdos =
        (double *)malloc(npts * sizeof(double)); /* Calculated DOS x vals */
    double *ydos =
        (double *)malloc(npts * sizeof(double)); /* Calculated DOS y */

    SetStdEig();
    EVSLStart();
    SetAMatrix(&Acsr);

    double *vinit = (double *)malloc(n * sizeof(double));
    rand_double(n, vinit);
    double lmin = 0.0, lmax = 0.0;
    ierr = LanTrbounds(50, 200, 1e-8, vinit, 1, &lmin, &lmax, fstats);

    intv[0] = lmin;
    intv[1] = lmax;
    /*----------------- get the bounds for (A, B) ---------*/
    intv[2] = lmin;
    intv[3] = lmax;

    double *ev;
    int numev;
    if (graph_exact_dos) {
      readVec("NM1AB_eigenvalues.dat", &numev, &ev);
    }

    /* Calculate computed DOS */
    ret = LanDos(nvec, msteps, npts, xdos, ydos, &neig, intv);
    fprintf(stdout, " LanDos ret %d \n", ret);

    /* Calculate the exact DOS */
    if (graph_exact_dos) {
      ret = exDOS(Acoo.vv, Acoo.ncols, npts, xHist, yHist, intv);
      fprintf(stdout, " exDOS ret %d \n", ret);
    }
    free_coo(&Acoo);
    free_csr(&Acsr);
    free(vinit);
    /* --------------------Make OUT dir if it doesn't exist */
    struct stat st = {0};
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    /*-------------------- Write to  output files */
    char computed_path[1024];
    strcpy(computed_path, "OUT/LanDosG_Approx_DOS_");
    strcat(computed_path, io.MatNam);
    FILE *ofp = fopen(computed_path, "w");
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
    char command[1024];
    strcpy(command, "gnuplot ");
    strcat(command, " -e \"filename='");
    strcat(command, computed_path);

    if (graph_exact_dos) {
      strcat(command, "';exact_dos='");
      strcat(command, path);
      strcat(command, "'\" testerG_ex.gnuplot");
      ierr = system(command);
    } else {
      strcat(command, "'\" testerG.gnuplot");
      ierr = system(command);
    }

    if (ierr) {
      fprintf(stderr,
              "Error using 'gnuplot < testerG.gnuplot', \n"
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
      free(xHist);
      free(yHist);
    }
    free(xdos);
    free(ydos);
    if(graph_exact_dos) {
      free(ev);
    }
    fclose(fstats);
  }
  EVSLFinish();
  return 0;
}
