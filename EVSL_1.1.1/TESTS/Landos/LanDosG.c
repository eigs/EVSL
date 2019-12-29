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
  fscanf(ifp, "%i", npts);
  *vec = evsl_Malloc(*npts, double);
  for (i = 0; i < (*npts); i++) {
    fscanf(ifp, "%lf", (&(*vec)[i]));
  }
  fclose(ifp);
  return 0;
}

/**
 *-----------------------------------------------------------------------
 * Tests landosG.c , the Lanczos DOS computed for the generalized eigenvalue
 * problem. Includes graphical comparison of calculated vs exact DOS
 *
 * use -graph_exact_dos 1 to enable graphing the exact DOS (from file)
 *
 *-----------------------------------------------------------------------
 */
int main(int argc, char *argv[]) {
  srand(time(NULL));
  const int msteps = 30;   /* Number of steps */
  const int degB = 40;     /* Degree to aproximate B with */
  const int npts = 200;    /* Number of points */
  const int nvec = 30;     /* Number of random vectors to use */
  const double tau = 1e-4; /* Tolerance in polynomial approximation */
  /* ---------------- Intervals of interest
     intv[0] and intv[1] are the input interval of interest [a. b]
     intv[2] and intv[3] are the smallest and largest eigenvalues of (A,B)
     */
  double intv[4];
  int n = 0, i;

  cooMat Acoo, Bcoo; /* A, B */
  csrMat Acsr, Bcsr; /* A, B */
#ifdef EVSL_USING_CUDA_GPU
  csrMat Acsr_gpu;
  csrMat Bcsr_gpu;
#endif
  double *sqrtdiag = NULL;

  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  int graph_exact_dos = 0;
  findarg("graph_exact_dos", INT, &graph_exact_dos, argc, argv);
  /*-------------------- start EVSL */
  EVSLStart();
  /*------------------ file "matfile" contains paths to matrices */
  if (NULL == (fmat = fopen("matfileG", "r"))) {
    fprintf(flog, "Can't open matfileG...\n");
    exit(2);
  } else {
    printf(" Read matfileG \n");
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
    int ierr;
    if ( (ierr = get_matrix_info(fmat, &io)) != 0 ) {
      fprintf(flog, "Invalid format in matfile %d...\n", ierr);
      exit(5);
    }
    /*----------------input matrix and interval information -*/
    fprintf(flog, " MATRIX A: %s...\n", io.MatNam1);
    fprintf(flog, " MATRIX B: %s...\n", io.MatNam2);
    char path[1024]; /* path to write the output files */
    strcpy(path, "OUT/LanDosG_");
    strcat(path, io.MatNam1);
    strcat(path, ".log");
    fstats = fopen(path, "w"); /* write all the output to the file io.MatNam */
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(fstats, "MATRIX B: %s...\n", io.MatNam2);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr = read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        /* nnz = Acoo.nnz; */
        n = Acoo.nrows;
        if (n <= 0) {
          fprintf(stderr, "non-positive number of rows");
          exit(7);
        }

        /* printf("size of A is %d\n", n);
           printf("nnz of  A is %d\n", nnz); */
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
      } else {
        fprintf(flog, "read_coo error for B = %d\n", ierr);
        exit(6);
      }
      /*------------------ diagonal scaling for Acoo and Bcoo */
      sqrtdiag = evsl_Calloc(n, double);
      /*-------------------- conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
      if (ierr) {
        fprintf(flog, "Could not convert matrix to csr: %i", ierr);
        exit(8);
      }
      ierr = cooMat_to_csrMat(0, &Bcoo, &Bcsr);
      if (ierr) {
        fprintf(flog, "Could not convert matrix to csr: %i", ierr);
        exit(9);
      }
    } else if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT  not supported (yet) * \n");
      exit(7);
    }

    /*-------------------- diagonal scaling for L-S poly. approx.
     *                     of B^{-1} and B^{-1/2},
     *                     which will be used in the DOS */
    /*-------------------- sqrt of diag(B) */
    extrDiagCsr(&Bcsr, sqrtdiag);
    for (i = 0; i < n; i++) {
      sqrtdiag[i] = sqrt(sqrtdiag[i]);
    }
    diagScalCsr(&Acsr, sqrtdiag);
    diagScalCsr(&Bcsr, sqrtdiag);

#ifdef EVSL_USING_CUDA_GPU
    evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
    evsl_create_csr_gpu(&Bcsr, &Bcsr_gpu);
#endif

    /*---------------- Set EVSL to solve std eig problem to
     *                 compute the range of the spectrum of B */
    SetStdEig();
#ifdef EVSL_USING_CUDA_GPU
    SetAMatrix_device_csr(&Bcsr_gpu);
#else
    SetAMatrix(&Bcsr);
#endif
    double *vinit = evsl_Malloc_device(n, double);
    rand_double_device(n, vinit);
    double lmin = 0.0, lmax = 0.0;
    ierr = LanTrbounds(50, 200, 1e-8, vinit, 1, &lmin, &lmax, fstats);
    if (ierr) {
      fprintf(flog, "Could not run LanTrBounds: %i", ierr);
      exit(10);
    }
    /*-------------------- Use polynomial to solve B and sqrt(B) */
    BSolDataPol BInv, BSqrtInv;
    /*-------------------- Setup the Bsol and Bsqrtsol struct */
    SetupPolRec(n, degB, tau, lmin, lmax, &BInv);
    SetupPolSqrt(n, degB, tau, lmin, lmax, &BSqrtInv);
    SetBSol(BSolPol, (void *)&BInv);
    SetLTSol(BSolPol, (void *)&BSqrtInv);
    printf(" The degree for LS polynomial approximations to B^{-1} and B^{-1/2} "
           "are %d and %d\n", BInv.deg, BSqrtInv.deg);
#ifdef EVSL_USING_CUDA_GPU
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix_device_csr(&Acsr_gpu);
    /*-------------------- set the right-hand side matrix B */
    SetBMatrix_device_csr(&Bcsr_gpu);
#else
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- set the right-hand side matrix B */
    SetBMatrix(&Bcsr);
#endif
    /*-------------------- Solve Gen Eig problem */
    SetGenEig();
    /*-------------------- eig bounds for B\A */
    rand_double_device(n, vinit);
    ierr = LanTrbounds(40, 200, 1e-10, vinit, 1, &lmin, &lmax, fstats);
    if (ierr) {
      fprintf(flog, "Could not run LanTrBounds: %i", ierr);
      exit(10);
    }
    fprintf(stdout, " Eig bounds for B\\A [%e %e]\n", lmin, lmax);
    /*----------------- plotting the DOS on [a, b] ---------*/
    intv[0] = io.a; //lmin;
    intv[1] = io.b; //lmax;
    /*----------------- get the bounds for (A, B) ---------*/
    intv[2] = lmin;
    intv[3] = lmax;
    /*-------------------- Read in the eigenvalues */
    double *ev = NULL;
    int numev;
    if (graph_exact_dos) {
      /*  Read in the exact eigenvalues from the dat file */
      readVec("NM1AB_eigenvalues.dat", &numev, &ev);
    }

    int ret;
    double neig;
    double *xHist = NULL;
    double *yHist = NULL;
    /*-------------------- exact histogram and computed DOS */
    if (graph_exact_dos) {
      xHist = evsl_Malloc(npts, double);
      yHist = evsl_Malloc(npts, double);
    }
    double *xdos = evsl_Malloc(npts, double);
    double *ydos = evsl_Malloc(npts, double);

    double t0 = evsl_timer();
    /* ------------------- Calculate the computed DOS */
    fprintf(stdout, " Running LanDosG on [%e, %e]\n", intv[0], intv[1]);
    ret = LanDosG(nvec, msteps, npts, xdos, ydos, &neig, intv);
    double t1 = evsl_timer();
    fprintf(stdout, " Estimated Eig count %f\n", neig);
    if (ret) {
      fprintf(stdout, " LanDosG error %d\n", ret);
      exit(1);
    }
    fprintf(stdout, " LanDosG time %.2f s\n", t1 - t0);

    /* -------------------- Calculate the exact DOS */
    if (graph_exact_dos) {
      ret = exDOS(ev, numev, npts, xHist, yHist, intv);
      if (ret) {
        fprintf(stdout, " exDos error %d\n", ret);
        exit(1);
      }
    }

    free_coo(&Acoo);
    free_coo(&Bcoo);
    free_csr(&Acsr);
    free_csr(&Bcsr);
#ifdef EVSL_USING_CUDA_GPU
    evsl_free_csr_gpu(&Acsr_gpu);
    evsl_free_csr_gpu(&Bcsr_gpu);
#endif
    evsl_Free_device(vinit);
    FreeBSolPolData(&BInv);
    FreeBSolPolData(&BSqrtInv);

    /*--------------------Make OUT dir if it doesn't exist */
    struct stat st = {0};
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    /*-------------------- Write to  output files */
    char computed_path[1024];
    strcpy(computed_path, "OUT/LanDosG_Approx_DOS_");
    strcat(computed_path, io.MatNam1);
    FILE *ofp = fopen(computed_path, "w");
    for (i = 0; i < npts; i++) {
      fprintf(ofp, " %.15e  %.15e\n", xdos[i], ydos[i]);
    }
    fclose(ofp);
    printf(" Wrote to:%s \n", computed_path);

    if (graph_exact_dos) {
      /*-------------------- save exact DOS */
      strcpy(path, "OUT/LanDosG_Exact_DOS_");
      strcat(path, io.MatNam1);
      ofp = fopen(path, "w");
      for (i = 0; i < npts; i++)
        fprintf(ofp, " %.15e  %.15e\n", xHist[i], yHist[i]);
      fclose(ofp);
    }

    printf(" The data output is located in  OUT/ \n");
    struct utsname buffer;
    errno = 0;
    if (uname(&buffer) != 0) {
      perror("uname");
      exit(EXIT_FAILURE);
    }

    /*-------------------- invoke gnuplot for plotting ... */
    char command[1024];
    strcpy(command, "gnuplot ");
    strcat(command, " -e \"filename='");
    strcat(command, computed_path);

    if (graph_exact_dos) {
      strcat(command, "';exact_dos='");
      strcat(command, path);
      strcat(command, "'\" testerG_ex.gnuplot");
      printf("Run command %s\n", command);
      ierr = system(command);
    } else {
      strcat(command, "'\" testerG.gnuplot");
      printf("Run command %s\n", command);
      ierr = system(command);
    }

    if (ierr) {
      fprintf(stderr,
              "Error using 'gnuplot < testerG.gnuplot', \n"
              "postscript plot could not be generated \n");
    } else {
      printf(" A postscript graph has been placed in %s%s\n", computed_path,
             ".eps");
      /*-------------------- and gv for visualizing / */
      if (!strcmp(buffer.sysname, "Linux")) {
        strcpy(command, "gv ");
        strcat(command, computed_path);
        strcat(command, ".eps &");
        printf("Run command %s\n", command);
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
    if (ev)  evsl_Free(ev);
    fclose(fstats);
    if (sqrtdiag) {
      evsl_Free(sqrtdiag);
    }
  } /* matrix loop */

  fclose(fmat);

  EVSLFinish();
  return 0;
}

