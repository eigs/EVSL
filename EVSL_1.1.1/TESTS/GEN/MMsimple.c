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
#include "../../SRC/blas/evsl_blas.h"

#define NO_SOLVES /*  If defined uses no solves */

int main(int argc, char** argv) {
  /**--------------------------------------------------------------
   * @file TESTS/GEN/MMsimple.c
   * @brief
   * This is a simple driver for the generalized eigenvalue problem.
   * It uses no factorizations, and only requires the matrix and, interval
   * and an estimate of the number of eigenvalues in the interval (greater than
   * or equal to the exact number). Instead of reading the matrix info from
   * `matfile' it enters it in the driver itself. Only one matrix is tested
   * and only one interval is used [no slicing]
   *
   * While this is still a little long, the majority is boilerplate-esque code, and code
   * to handle I/O.
   *
   * The matrix is specified in the matfile, while the number of eigenvalues
   * are specified as the only command argument.
   *-------------------------------------------------------------*/
  int n = 0, i, j,  mlan, ierr, totcnt;
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
#ifdef EVSL_USING_CUDA_GPU
  csrMat Acsr_gpu;
  csrMat Bcsr_gpu;
#endif
  /* csrMat Acsr0, Bcsr0; */
  /* slicer parameters */
  FILE *flog = stdout, *fstats = NULL;
  io_t io; /** Used for input */
  /*-------------------- Bsol */
  int         Bsol_direct = 0; /* if using direct solver for B^{-1} */
  void       *Bsol = NULL;
  BSolDataPol BsolPol, BsqrtsolPol;
  int         BsolDeg = 200;  /* Max degree to approximate inv(B) with */
  double      BsolTol = 1e-6; /* Tolerance in polynomial approximation */
  int         DiagScalB = 1;  /* Apply diagonal scaling to A and B */
  double     *Sqrtdiag = NULL, *Sqrtdiag_device = NULL;
  /*-------------------- stopping tol */
  tol = 1e-08;
  /*-------------------- start EVSL */
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
  EVSLStart();
  /*----------------input matrix and interval information -*/
  /*-------------------- In this driver we avoid reading matrix from file - input info here:" */
  strcpy(io.Fname1,"../../MATRICES/NM1A.mtx");   // path for matrixA
  strcpy(io.Fname2,"../../MATRICES/NM1B.mtx");   // path for matrixB
  strcpy(io.MatNam1,"nm1A");                     // short name for A
  strcpy(io.MatNam2,"nm1B");                     // short name for B
  a = 3.947842e-07;                              // a of interval [a b]
  b = 3.947842e-05;                              // b of interval [a b]
  io.Fmt = 3;                                    // MM1 format
  ecount = 60;                                   // estimate of eigenvalue count
  /*-------------------- print some info */
  fprintf(flog, "MATRIX A: %s...\n", io.MatNam1);
  fprintf(flog, "MATRIX B: %s...\n", io.MatNam2);
  /*-------------------- output files*/
  struct stat st = {0}; /* Make sure OUT directory exists */
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    char path[1024]; // path to write the output files
    strcpy(path, "OUT/MMsimple_");
    strcat(path, io.MatNam1);
    fstats = fopen(path, "w"); // write all the output to the file io.MatNam
    if (!fstats) {
      printf(" failed in opening output file in OUT/\n");
      fstats = stdout;
    }
    fprintf(fstats, "MATRIX A: %s...\n", io.MatNam1);
    fprintf(fstats, "MATRIX B: %s...\n", io.MatNam2);
     /*-------------------- No slicing: only one interval -- */
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
    /*-------------------- Finish with input  */
    alleigs = evsl_Malloc(n, double);
    /*-------------------- initial vector */
    vinit = evsl_Malloc_device(n,  double);
    rand_double_device(n, vinit);/*  Fill with random */
    /*-------------------- B^{-1} */
    if (DiagScalB) {
      /* For FEM matrices, diagonal scaling is often found useful to
       * improve condition */
      Sqrtdiag = evsl_Calloc(n, double);
      extrDiagCsr(&Bcsr, Sqrtdiag);
      for (i = 0; i < n; i++) {
        Sqrtdiag[i] = sqrt(Sqrtdiag[i]);
      }
      diagScalCsr(&Acsr, Sqrtdiag);
      diagScalCsr(&Bcsr, Sqrtdiag);
#ifdef EVSL_USING_CUDA_GPU
      Sqrtdiag_device = evsl_Malloc_device(n, double);
      evsl_memcpy_host_to_device(Sqrtdiag_device, Sqrtdiag, n*sizeof(double));
#else
      Sqrtdiag_device = Sqrtdiag;
#endif
    }
#ifdef EVSL_USING_CUDA_GPU
    evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
    evsl_create_csr_gpu(&Bcsr, &Bcsr_gpu);
#endif
    /*-------------------- use polynomial approx. to B^{-1} and B^{-1/2}
     *                     B is assumed to be ``well-conditioned'' */
    /*       compute eig bounds for B, set to std eig prob for now */
    SetStdEig();
#ifdef EVSL_USING_CUDA_GPU
      SetAMatrix_device_csr(&Bcsr_gpu);
#else
    SetAMatrix(&Bcsr);
#endif
    double lminB, lmaxB;
    /*--------------------  Calculate the bounds */
    LanTrbounds(50, 200, 1e-10, vinit, 1, &lminB, &lmaxB, fstats);
    SetupPolRec (n, BsolDeg, BsolTol, lminB, lmaxB, &BsolPol);
    fprintf(fstats, "B-INV  Pol deg %d\n", BsolPol.deg);
    SetupPolSqrt(n, BsolDeg, BsolTol, lminB, lmaxB, &BsqrtsolPol);
    fprintf(fstats, "B-SQRT Pol deg %d\n", BsqrtsolPol.deg);
    /*-------------------- set the solver for B and LT */
    SetBSol (BSolPol, &BsolPol);
    SetLTSol(BSolPol, &BsqrtsolPol);
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
    /*-------------------- back to the generalized eigenvalue problem */
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
    //-------------------- # eigs per slice
    totcnt = 0;
    /*-------------------- We only have one slice. View other example for multiple slices */
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
    //-------------------- Now determine polymomial
    find_pol(xintv, &pol);
    fprintf(fstats, " polynomial [type = %d], deg %d, bar %e gam %e\n",
            pol.type, pol.deg, pol.bar, pol.gam);
    // save_vec(pol.deg+1, pol.mu, "OUT/mu.txt");
    //-------------------- approximate number of eigenvalues wanted
    //-------------------- Dimension of Krylov subspace and maximal iterations
    mlan = evsl_max(5 * (ecount+2), 300);
    mlan = evsl_min(mlan, n);
    //-------------------- then call ChenLanNr
    ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, &lam, &Y, &res, fstats);
    if (ierr) {
      printf("ChebLanTr error %d\n", ierr);
      return 1;
    }
    /* in the case of diagonal scaling, recover Y and recompute residual */
    if (Sqrtdiag) {
      double *v1 = evsl_Malloc_device(n, double);
      double *v2 = evsl_Malloc_device(n, double);
      int one = 1;
      for (i = 0; i < nev2; i++) {
        double *yi = Y+i*n;
        double t = -lam[i];
        evsldata.Amv->func(yi, v1, evsldata.Amv->data);
        evsldata.Bmv->func(yi, v2, evsldata.Bmv->data);
        evsl_daxpy_device(&n, &t, v2, &one, v1, &one);
        /* v1[j] *= Sqrtdiag[j] */
        evsl_element_mult_device(n, Sqrtdiag_device, v1);
        res[i] = evsl_dnrm2_device(&n, v1, &one);
        /* yi[j] /= Sqrtdiag[j] */
        evsl_element_divide_device(n, Sqrtdiag_device, yi);
      }
      evsl_Free_device(v1);
      evsl_Free_device(v2);
    }
    /* sort the eigenvals: ascending order
     * ind: keep the original indices */
    ind = evsl_Malloc(nev2, int);
    sort_double(nev2, lam, ind);
    printf(" number of eigenvalues found: %d\n", nev2);
    /*-------------------- print eigenvalues */
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
      evsl_Free(lam);
    }
    if (Y) {
      evsl_Free_device(Y);
    }
    if (res) {
      evsl_Free(res);
    }
    free_pol(&pol);
    evsl_Free(ind);
    //-------------------- free other allocated space
    fprintf(fstats, " --> Total eigenvalues found = %d\n", totcnt);
    sprintf(path, "OUT/EigsOut_MMsimple_(%s_%s)", io.MatNam1, io.MatNam2);
    FILE *fmtout = fopen(path, "w");
    if (fmtout) {
      for (j = 0; j < totcnt; j++) {
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      }
      fclose(fmtout);
    }
    evsl_Free_device(vinit);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free_coo(&Bcoo);
    free_csr(&Bcsr);
#ifdef EVSL_USING_CUDA_GPU
    evsl_free_csr_gpu(&Acsr_gpu);
    evsl_free_csr_gpu(&Bcsr_gpu);
#endif
    evsl_Free(alleigs);
    if (fstats != stdout) {
      fclose(fstats);
    }
    if (Sqrtdiag) {
      evsl_Free(Sqrtdiag);
    }
#ifdef EVSL_USING_CUDA_GPU
    if (Sqrtdiag_device) {
      evsl_Free_device(Sqrtdiag_device);
    }
#endif
    if (Bsol_direct) {
      FreeBSolDirectData(Bsol);
    } else {
      FreeBSolPolData(&BsolPol);
      FreeBSolPolData(&BsqrtsolPol);
    }
    if (flog != stdout) {
      fclose(flog);
    }
    /*-------------------- finalize EVSL */
    EVSLFinish();
    return 0;
}

