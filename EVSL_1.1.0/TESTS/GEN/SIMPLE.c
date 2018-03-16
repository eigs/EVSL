
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include "io.h"
#include "evsl.h"
#include "evsl_direct.h"
#include "blaslapack.h"


#ifdef __cplusplus
extern "C" {
#define max(a, b) std::max(a, b)
#define min(a, b) std::min(a, b)
#else
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

#define NO_SOLVES /*  If defined uses no solves */
/*-------------------- Protos */
int read_coo_MM(const char *matfile, int idxin, int idxout, cooMat *Acoo);
int get_matrix_info(FILE *fmat, io_t *pio);
/*-------------------- End Protos */

int main(int argc, char** argv) {
  /**--------------------------------------------------------------
   * @file TESTS/GEN/SIMPLE.c
   * @brief
   * This is a simple driver for the generalized eigenvalue problem.
   * It uses no factorizations, and only requires the matrix and, interval
   * and an estimate of the number of eigenvalues in the interval (greater than
   * or equal to the exact number).
   *
   * While this is pretty long, the majority is boilerplate-esque code, and code
   * to handle I/O. 
   *
   * The matrix is specified in the matfile, while the number of eigenvalues
   * are specified as the only command argument.
   *-------------------------------------------------------------*/
  int n = 0, i, j, nslices, nev, mlan, ev_int, ierr, totcnt;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol;
  double xintv[4]; /**< Specifies the interval */
  double
      sli[2]; /**< Contains the endpoints of each sub-interval, in this case 
                   only 2*/
  double *alleigs;
  /* initial vector: random */
  double *vinit;
  polparams pol; /**< Parameters for polynomial*/
  /*-------------------- matrices A, B: coo format and csr format */
  cooMat Acoo, Bcoo; /**< COO matrices (for input)*/
  csrMat Acsr, Bcsr; /**< CSR matrices (which EVSL operates on)*/
  /* csrMat Acsr0, Bcsr0; */
  /* slicer parameters */
  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io; /** Used for input */
  int numat, mat;/**< Used for input */
  char line[MAX_LINE]; /**< Used to read Matrix Market format */
  /*-------------------- Bsol */
  int         Bsol_direct = 0; /* if using direct solver for B^{-1} */
  void       *Bsol = NULL;
  BSolDataPol BsolPol, BsqrtsolPol;
  int         BsolDeg = 200;  /* Max degree to aproximate B with */
  double      BsolTol = 1e-6; /* Tolerance in polynomial approximation */
  int         DiagScalB = 1;  /* Apply diagonal scaling to A and B */
  double     *Sqrtdiag = NULL;

  char* p;

  assert(argc == 2); /**< Confirm number of eigenvalues are provided */
  long conv= strtol(argv[1], &p, 10);
  if (errno != 0 || *p != '\0' || conv > INT_MAX) {
    fprintf(stderr, "Could not use eigenvalue count\n");
    exit(-1);
  } else {
    ecount = (double) conv;
  }

  /*-------------------- stopping tol */
  tol = 1e-6;
  /*-------------------- start EVSL */
  EVSLStart();
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
  numat = 1; /*Limit to a single matrix, as only a single eigenvalue count
  can be input */

  for (mat = 1; mat <= numat; mat++) {

    if (get_matrix_info(fmat, &io) != 0) {
      fprintf(flog, "Invalid format in matfile ...\n");
      exit(5);
    }
    /*----------------input matrix and interval information -*/
    fprintf(flog, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(flog, "MATRIX B: %s...\n", io.MatNam2);
    a = io.a; // left endpoint of input interval
    b = io.b; // right endpoint of input interval
    nslices = 1;

    struct stat st = {0}; /* Make sure OUT directory exists */
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    char path[1024]; // path to write the output files
    strcpy(path, "OUT/SIMPLE_");
    strcat(path, io.MatNam1);
    fstats = fopen(path, "w"); // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(fstats, "MATRIX B: %s...\n", io.MatNam2);
    fprintf(fstats, "Partition the interval of interest [%f,%f] into %d slices\n", 
            a, b, nslices);

    /* Because we are are not slicing, the slices are just the input */
    sli[0] = a;
    sli[1] = b;

    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr = read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        n = Acoo.nrows;
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
      /*-------------------- conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
      ierr = cooMat_to_csrMat(0, &Bcoo, &Bcsr);
    } else if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT  not supported (yet) * \n");
      exit(7);
    }


    /* Finish with input  */

    alleigs = malloc(n * sizeof(double));
    /*-------------------- initial vector */
    vinit = (double *)malloc(n * sizeof(double));
    rand_double(n, vinit);/*  Fill with random */
    /*-------------------- B^{-1} */
    {
      if (DiagScalB) {
        /* For FEM matrices, diagonal scaling is often found useful to 
         * improve condition */
        Sqrtdiag = (double *)calloc(n, sizeof(double));
        extrDiagCsr(&Bcsr, Sqrtdiag);
        for (i = 0; i < n; i++) {
          Sqrtdiag[i] = sqrt(Sqrtdiag[i]);
        }
        diagScalCsr(&Acsr, Sqrtdiag);
        diagScalCsr(&Bcsr, Sqrtdiag);
      }
      /*-------------------- use polynomial approx. to B^{-1} and B^{-1/2}
       *                     B is assumbed to be ``well-conditioned'' */
      /* compute eig bounds for B, set to std eig prob for now */
      SetStdEig();
      SetAMatrix(&Bcsr);
      double lminB;
      double lmaxB;
      /*  Calculate the bounds */
      LanTrbounds(50, 200, 1e-10, vinit, 1, &lminB, &lmaxB, fstats);
      SetupPolRec (n, BsolDeg, BsolTol, lminB, lmaxB, &BsolPol);
      fprintf(fstats, "B-INV  Pol deg %d\n", BsolPol.deg);
      SetupPolSqrt(n, BsolDeg, BsolTol, lminB, lmaxB, &BsqrtsolPol);
      fprintf(fstats, "B-SQRT Pol deg %d\n", BsqrtsolPol.deg);
      /*-------------------- set the solver for B and LT */
      SetBSol (BSolPol, &BsolPol);
      SetLTSol(BSolPol, &BsqrtsolPol);
    }
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- set the right-hand side matrix B */
    SetBMatrix(&Bcsr);
    /*-------------------- for generalized eigenvalue problem */
    SetGenEig();
    ierr = LanTrbounds(50, 200, 1e-12, vinit, 1, &lmin, &lmax, fstats);
    /*-------------------- step 0: get eigenvalue bounds */
    fprintf(fstats, "Step 0: Eigenvalue bound s for B^{-1}*A: [%.15e, %.15e]\n",
            lmin, lmax);
    /*-------------------- interval and eig bounds */
    xintv[0] = a;
    xintv[1] = b;
    xintv[2] = lmin;
    xintv[3] = lmax;
    //linspace(a, b, nslices+1, sli);

    //-------------------- # eigs per slice
    ev_int = (int)(1 + ecount / ((double)nslices));
    totcnt = 0;
    /* We only have one slice. View other example for multiple slices */
    {
      printf("======================================================\n");
      int nev2;
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      a = sli[0]; /* Specify the left end of this sub interval */
      b = sli[1]; /* Specify the right end of this sub interval */
      printf(" subinterval: [%.15e , %.15e]\n", a, b);
      xintv[0] = a; /*  Left end of the subinterval */
      xintv[1] = b; /* Right end of the subinterval */
      xintv[2] = lmin; /* Minimum eigenvalue */
      xintv[3] = lmax; /* Maximum eigenvalue */

      //-------------------- set up default parameters for pol.
      set_pol_def(&pol);
      // can change default values here e.g.
      pol.damping = 2;
      pol.thresh_int = 0.8;
      pol.thresh_ext = 0.2;
      // pol.max_deg  = 300;

      //-------------------- Now determine polymomial
      find_pol(xintv, &pol);
      fprintf(fstats, " polynomial [type = %d], deg %d, bar %e gam %e\n",
              pol.type, pol.deg, pol.bar, pol.gam);
      // save_vec(pol.deg+1, pol.mu, "OUT/mu.txt");
      //-------------------- approximate number of eigenvalues wanted
      nev = ev_int + 2;
      //-------------------- Dimension of Krylov subspace and maximal iterations
      mlan = max(5 * nev, 300);
      mlan = min(mlan, n);
      //-------------------- then call ChenLanNr
      ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, &lam, &Y, &res,
                       fstats);
      if (ierr) {
        printf("ChebLanTr error %d\n", ierr);
        return 1;
      }

      /* in the case of digonal scaling, recover Y and recompute residual */
      if (Sqrtdiag) {
        double *v1 = (double *) malloc(n*sizeof(double));
        double *v2 = (double *) malloc(n*sizeof(double));
        int one = 1;
        for (i = 0; i < nev2; i++) {
          double *yi = Y+i*n;
          double t = -lam[i];
          matvec_csr(yi, v1, &Acsr);
          matvec_csr(yi, v2, &Bcsr);
          DAXPY(&n, &t, v2, &one, v1, &one);
          for (j = 0; j < n; j++) {
            v1[j] *= Sqrtdiag[j];
          }
          res[i] = DNRM2(&n, v1, &one);
          
          for (j = 0; j < n; j++) {
            yi[j] /= Sqrtdiag[j];
          }
        }
        free(v1);
        free(v2);
      }

      /* sort the eigenvals: ascending order
       * ind: keep the orginal indices */
      ind = (int *)malloc(nev2 * sizeof(int));

      sort_double(nev2, lam, ind);

      printf(" number of eigenvalues found: %d\n", nev2);
      /* print eigenvalues */
      fprintf(fstats, "    Eigenvalues in [a, b]\n");
      fprintf(fstats, "    Computed [%d]        ||Res||\n", nev2);
      for (i = 0; i < nev2; i++) {
        fprintf(fstats, "% .15e  %.15e\n", lam[i], res[ind[i]]);
      }
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - "
                      "- - - - - - - - - - - - - - - - - -\n");
      /*  Copy to array of all eigenvalues */
      memcpy(&alleigs[totcnt], lam, nev2 * sizeof(double));

      totcnt += nev2;
      //-------------------- free allocated space withing this scope
      if (lam) {
        free(lam);
      }
      if (Y) {
        free(Y);
      }
      if (res) {
        free(res);
      }
      free_pol(&pol);
      free(ind);
    } 
    //-------------------- free other allocated space
    fprintf(fstats, " --> Total eigenvalues found = %d\n", totcnt);
    sprintf(path, "OUT/EigsOut_Lan_MMPLanN_(%s_%s)", io.MatNam1, io.MatNam2);
    FILE *fmtout = fopen(path, "w");
    if (fmtout) {
      for (j = 0; j < totcnt; j++) {
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      }
      fclose(fmtout);
    }
    free(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free_coo(&Bcoo);
    free_csr(&Bcsr);
    free(alleigs);
    if (fstats != stdout) {
      fclose(fstats);
    }
    if (Sqrtdiag) {
      free(Sqrtdiag);
    }
    if (Bsol_direct) {
      FreeBSolDirectData(Bsol);
    } else { 
      FreeBSolPolData(&BsolPol);
      FreeBSolPolData(&BsqrtsolPol);
    }
    /*-------------------- end matrix loop */
  }
  if (flog != stdout) {
    fclose(flog);
  }
  fclose(fmat);
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

#ifdef __cplusplus
}
#endif
