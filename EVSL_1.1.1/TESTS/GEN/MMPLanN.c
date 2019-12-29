#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "io.h"
#include "evsl.h"
#include "evsl_direct.h"
#include "../../SRC/blas/evsl_blas.h"

int main() {
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix pair
   * read in matrix format -- using
   * Non-Restart Lanczos with polynomial filtering.
   *-------------------------------------------------------------*/
  int n = 0, i, j, npts, nslices, nvec, nev, mlan, ev_int, sl, ierr, totcnt;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli;
  double xintv[4];
  double *alleigs;
  int *counts; // #ev computed in each slice
  /* initial vector: random */
  double *vinit;
  polparams pol;
  /*-------------------- matrices A, B: coo format and csr format */
  cooMat Acoo, Bcoo;
  csrMat Acsr, Bcsr;
#ifdef EVSL_USING_CUDA_GPU
  csrMat Acsr_gpu;
  csrMat Bcsr_gpu;
#endif
  /* slicer parameters */
  int msteps = 40;
  nvec = 10;
  npts = 200;
  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io;
  int numat, mat;
  char line[MAX_LINE];
  /*-------------------- Bsol */
  int         Bsol_direct = 0; /* if using direct solver for B^{-1} */
  void       *Bsol = NULL;
  BSolDataPol BsolPol, BsqrtsolPol;
  int         BsolDeg = 200;  /* Max degree to aproximate B with */
  double      BsolTol = 1e-6; /* Tolerance in polynomial approximation */
  int         DiagScalB = 1;  /* Apply diagonal scaling to A and B */
  double     *Sqrtdiag = NULL, *Sqrtdiag_device = NULL;

  /*-------------------- stopping tol */
  tol = 1e-6;
  /*-------------------- start EVSL */
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
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
    nslices = io.n_intv;

    struct stat st = {0}; /* Make sure OUT directory exists */
    if (stat("OUT", &st) == -1) {
      mkdir("OUT", 0750);
    }

    char path[1024]; // path to write the output files
    strcpy(path, "OUT/Lan_MMPLanN_");
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
    counts = evsl_Malloc(nslices, int);
    sli = evsl_Malloc(nslices+1, double);
    /*-------------------- Read matrix - case: COO/MatrixMarket formats */
    if (io.Fmt > HB) {
      ierr = read_coo_MM(io.Fname1, 1, 0, &Acoo);
      if (ierr == 0) {
        fprintf(fstats, "matrix read successfully\n");
        // nnz = Acoo.nnz;
        n = Acoo.nrows;
        // printf("size of A is %d\n", n);
        // printf("nnz of  A is %d\n", nnz);
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
    alleigs = evsl_Malloc(n, double);
    /*-------------------- initial vector */
    vinit = evsl_Malloc_device(n, double);
    rand_double_device(n, vinit);
    /*-------------------- B^{-1} */
    if (Bsol_direct) {
      /*-------------------- use direct solver as the solver for B */
      SetupBSolDirect(&Bcsr, &Bsol);
      /*-------------------- set the solver for B and LT */
      SetBSol(BSolDirect, Bsol);
      SetLTSol(LTSolDirect, Bsol);
#ifdef EVSL_USING_CUDA_GPU
      evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
      evsl_create_csr_gpu(&Bcsr, &Bcsr_gpu);
#endif
    } else {
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
       *                     B is assumbed to be ``well-conditioned'' */
      /* compute eig bounds for B, set to std eig prob for now */
      SetStdEig();
#ifdef EVSL_USING_CUDA_GPU
      SetAMatrix_device_csr(&Bcsr_gpu);
#else
      SetAMatrix(&Bcsr);
#endif
      double lminB, lmaxB;
      LanTrbounds(50, 200, 1e-10, vinit, 1, &lminB, &lmaxB, fstats);
      fprintf(fstats, "lmin, lmax of B: %e %e (lmax/lmin = %e)\n", lminB, lmaxB, lmaxB/lminB);
      SetupPolRec (n, BsolDeg, BsolTol, lminB, lmaxB, &BsolPol);
      fprintf(fstats, "B-INV  Pol deg %d\n", BsolPol.deg);
      SetupPolSqrt(n, BsolDeg, BsolTol, lminB, lmaxB, &BsqrtsolPol);
      fprintf(fstats, "B-SQRT Pol deg %d\n", BsqrtsolPol.deg);
      /*-------------------- set the solver for B and LT */
      SetBSol (BSolPol, &BsolPol);
      SetLTSol(BSolPol, &BsqrtsolPol);
    }
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
    /*-------------------- for generalized eigenvalue problem */
    SetGenEig();
    /*-------------------- step 0: get eigenvalue bounds */
    ierr = LanTrbounds(50, 200, 1e-12, vinit, 1, &lmin, &lmax, fstats);
    fprintf(fstats, "Step 0: Eigenvalue bound s for B^{-1}*A: [%.15e, %.15e]\n",
            lmin, lmax);
    /*-------------------- interval and eig bounds */
    xintv[0] = a;
    xintv[1] = b;
    xintv[2] = lmin;
    xintv[3] = lmax;
    /*-------------------- call LanczosDOS for spectrum slicing */
    /*-------------------- define landos parameters */
    double t = evsl_timer();
    double *xdos = evsl_Calloc(npts, double);
    double *ydos = evsl_Calloc(npts, double);
    ierr = LanDosG(nvec, msteps, npts, xdos, ydos, &ecount, xintv);
    t = evsl_timer() - t;
    if (ierr) {
      printf("Landos error %d\n", ierr);
      return 1;
    }
    fprintf(fstats, " Time to build DOS (Landos) was : %10.2f  \n", t);
    fprintf(fstats, " estimated eig count in interval: %.15e \n", ecount);
    //-------------------- call splicer to slice the spectrum
    fprintf(fstats, "DOS parameters: msteps = %d, nvec = %d, npnts = %d\n",
            msteps, nvec, npts);
    spslicer2(xdos, ydos, nslices, npts, sli);
    printf("====================  SLICES FOUND  ====================\n");
    for (j = 0; j < nslices; j++) {
      printf(" %2d: [% .15e , % .15e]\n", j + 1, sli[j], sli[j + 1]);
    }
    //-------------------- # eigs per slice
    ev_int = (int)(1 + ecount / ((double)nslices));
    totcnt = 0;
    //-------------------- For each slice call RatLanrNr
    for (sl = 0; sl < nslices; sl++) {
      printf("======================================================\n");
      int nev2;
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      a = sli[sl];
      b = sli[sl + 1];
      printf(" subinterval: [%.15e , %.15e]\n", a, b);
      xintv[0] = a;
      xintv[1] = b;
      xintv[2] = lmin;
      xintv[3] = lmax;
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
      mlan = evsl_max(5 * nev, 300);
      mlan = evsl_min(mlan, n);
      //-------------------- then call ChenLanNr
      ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, &lam, &Y, &res,
                       fstats);
      if (ierr) {
        printf("ChebLanTr error %d\n", ierr);
        return 1;
      }

      /* in the case of digonal scaling, recover Y and recompute residual */
      if (Sqrtdiag) {
        double *v1 = evsl_Malloc_device(n, double);
        double *v2 = evsl_Malloc_device(n, double);
        for (i = 0; i < nev2; i++) {
          int one = 1;
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
       * ind: keep the orginal indices */
      ind = evsl_Malloc(nev2, int);
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
      memcpy(&alleigs[totcnt], lam, nev2 * sizeof(double));
      totcnt += nev2;
      counts[sl] = nev2;
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
      StatsPrint(fstats);
    } // for (sl=0; sl<nslices; sl++)

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
    evsl_Free_device(vinit);
    evsl_Free(sli);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free_coo(&Bcoo);
    free_csr(&Bcsr);
#ifdef EVSL_USING_CUDA_GPU
    evsl_free_csr_gpu(&Acsr_gpu);
    evsl_free_csr_gpu(&Bcsr_gpu);
#endif
    evsl_Free(alleigs);
    evsl_Free(counts);
    evsl_Free(xdos);
    evsl_Free(ydos);
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

