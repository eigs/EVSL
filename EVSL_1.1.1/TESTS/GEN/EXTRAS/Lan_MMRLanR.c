#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "io.h"
#include "evsl.h"
#include "evsl_direct.h"

int main() {
  /*--------------------------------------------------------------
   * this tests the spectrum slicing idea for a generic matrix pair
   * read in matrix format -- using
   * Thick-Restarted Lanczos with rational filtering.
   *-------------------------------------------------------------*/
  int n = 0, i, j, npts, nslices, nvec, nev, mlan, max_its, ev_int, sl, ierr,
      totcnt;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli;
  double xintv[4];
  double *alleigs;
  int *counts; // #ev computed in each slice
  /* initial vector: random */
  double *vinit;
  /* parameters for rational filter */
  int num = 1;        // number of poles used for each slice
  int pow = 2;        // multiplicity of each pole
  double beta = 0.01; // beta in the LS approximation
  /*-------------------- matrices A, B: coo format and csr format */
  cooMat Acoo, Bcoo;
  csrMat Acsr, Bcsr, Acsr0, Bcsr0;
  double *sqrtdiag = NULL;
  /* slicer parameters */
  int msteps = 40;
  nvec = 10;
  npts = 200;
  FILE *flog = stdout, *fmat = NULL;
  FILE *fstats = NULL;
  io_t io;
  const int degB = 200;    // Max degree to aproximate B with
  const double tau = 1e-4; // Tolerance in polynomial approximation
  int numat, mat;
  char line[MAX_LINE];
  /*-------------------- stopping tol */
  tol = 1e-6;
  /*-------------------- Polynomial approximation to B and sqrtB*/
  BSolDataPol Bsol, Bsqrtsol;
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
    strcpy(path, "OUT/Lan_MMRLanR_");
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
      /*------------------ diagonal scaling for Acoo and Bcoo */
      sqrtdiag = evsl_Calloc(n, double);
      /*------------------ conversion from COO to CSR format */
      ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
      ierr = cooMat_to_csrMat(0, &Bcoo, &Bcsr);
    } else if (io.Fmt == HB) {
      fprintf(flog, "HB FORMAT not supported (yet) * \n");
      exit(7);
    }

    /*-------------------- diagonal scaling for L-S poly. approx.
     *                     of B^{-1} and B^{-1/2},
     *                     which will be used in the DOS */
    /*-------------------- sqrt of diag(B) */
    extrDiagCsr(&Bcsr, sqrtdiag);
    for (i=0; i<n; i++) {
      sqrtdiag[i] = sqrt(sqrtdiag[i]);
    }
    /*-------------------- backup A and B */
    csr_copy(&Acsr, &Acsr0, 1); /* 1 stands for memory alloc */
    csr_copy(&Bcsr, &Bcsr0, 1);
    /*-------------------- Scale A and B */
    diagScalCsr(&Acsr, sqrtdiag);
    diagScalCsr(&Bcsr, sqrtdiag);
    if (sqrtdiag) {
      evsl_Free(sqrtdiag);
    }

    /*---------------- Set EVSL to solve std eig problem to
     *---------------- compute the range of the spectrum of B */
    SetStdEig();
    SetAMatrix(&Bcsr);
    vinit = evsl_Malloc(n, double);
    rand_double(n, vinit);
    ierr = LanTrbounds(50, 200, 1e-10, vinit, 1, &lmin, &lmax, fstats);
    /*-------------------- Use polynomial to solve B and sqrt(B) */
    /*-------------------- Setup the Bsol and Bsqrtsol struct */
    SetupPolRec (n, degB, tau, lmin, lmax, &Bsol);
    SetupPolSqrt(n, degB, tau, lmin, lmax, &Bsqrtsol);
    SetBSol (BSolPol, (void *)&Bsol);
    SetLTSol(BSolPol, (void *)&Bsqrtsol);
    printf("The degree for LS polynomial approximations to B^{-1} and B^{-1/2} "
           "are %d and %d\n", Bsol.deg, Bsqrtsol.deg);
    /*-------------------- set the left-hand side matrix A */
    SetAMatrix(&Acsr);
    /*-------------------- set the right-hand side matrix B */
    SetBMatrix(&Bcsr);
    /*-------------------- for generalized eigenvalue problem */
    SetGenEig();
    /*-------------------- step 0: get eigenvalue bounds */
    //-------------------- initial vector
    rand_double(n, vinit);
    ierr = LanTrbounds(50, 200, 1e-10, vinit, 1, &lmin, &lmax, fstats);
    fprintf(fstats, "Step 0: Eigenvalue bound s for B^{-1}*A: [%.15e, %.15e]\n",
            lmin, lmax);
    /*-------------------- interval and eig bounds */
    xintv[0] = a;
    xintv[1] = b;
    xintv[2] = lmin;
    xintv[3] = lmax;
    /*-------------------- call landos to get the DOS for dividing the
     *                     spectrum*/
    /*-------------------- define landos parameters */
    double t = evsl_timer();
    double *xdos = evsl_Malloc(npts, double);
    double *ydos = evsl_Malloc(npts, double);
    ierr = LanDosG(nvec, msteps, npts, xdos, ydos, &ecount, xintv);
    t = evsl_timer() - t;
    if (ierr) {
      printf("Landos error %d\n", ierr);
      return 1;
    }
    fprintf(fstats, " Time to build DOS (Landos) was : %10.2f  \n", t);
    fprintf(fstats, " estimated eig count in interval: %.15e \n", ecount);
    //-------------------- call splicer to slice the spectrum
    sli = evsl_Malloc(nslices+1, double);
    fprintf(fstats, "DOS parameters: msteps = %d, nvec = %d, npnts = %d\n",
            msteps, nvec, npts);
    spslicer2(xdos, ydos, nslices, npts, sli);
    printf("====================  SLICES FOUND  ====================\n");
    for (j = 0; j < nslices; j++) {
      printf(" %2d: [% .15e , % .15e]\n", j + 1, sli[j], sli[j + 1]);
    }
    /*-------------------- DONE WITH DOS */
    FreeBSolPolData(&Bsol);
    FreeBSolPolData(&Bsqrtsol);

    //-------------------- # eigs per slice
    ev_int = (int)(1 + ecount / ((double)nslices));
    totcnt = 0;
    alleigs = evsl_Malloc(n, double);

    /* recover the original matrices A and B before scaling
     * Note that B-sol and sqrt(B)-sol will not be needed in RatLan,
     * so we can recover them */
    csr_copy(&Acsr0, &Acsr, 0); /* 0 stands for no memory alloc */
    csr_copy(&Bcsr0, &Bcsr, 0);

    //-------------------- For each slice call RatLanrNr
    for (sl = 0; sl < nslices; sl++) {
      printf("======================================================\n");
      int nev2;
      double *lam, *Y, *res;
      int *ind;
      //--------------------
      a = sli[sl];
      b = sli[sl + 1];
      printf(" subinterval: [%.4e , %.4e]\n", a, b);
      double intv[4];
      intv[0] = a;
      intv[1] = b;
      intv[2] = lmin;
      intv[3] = lmax;
      // find the rational filter on this slice
      ratparams rat;
      // setup default parameters for rat
      set_ratf_def(&rat);
      // change some default parameters here:
      rat.pw = pow;
      rat.num = num;
      rat.beta = beta;
      // now determine rational filter
      find_ratf(intv, &rat);
      // use direct solver function
      void **solshiftdata = evsl_Malloc(num, void *);
      /*------------ factoring the shifted matrices and store the factors */
      SetupASIGMABSolDirect(&Acsr, &Bcsr, num, rat.zk, solshiftdata);
      /*------------ give the data to rat */
      SetASigmaBSol(&rat, NULL, ASIGMABSolDirect, solshiftdata);
      //-------------------- approximate number of eigenvalues wanted
      nev = ev_int + 2;
      //-------------------- Dimension of Krylov subspace and maximal iterations
      mlan = evsl_max(4 * nev, 100);
      mlan = evsl_min(mlan, n);
      max_its = 3 * mlan;
      //-------------------- RationalLanNr
      ierr = RatLanTr(mlan, nev, intv, max_its, tol, vinit, &rat, &nev2, &lam,
                      &Y, &res, fstats);
      if (ierr) {
        printf("RatLanNr error %d\n", ierr);
        return 1;
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
        fprintf(fstats, "% .15e  %.1e\n", lam[i], res[ind[i]]);
      }
      fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - "
                      "- - - - - - - - - - - - - - - - - -\n");
      memcpy(&alleigs[totcnt], lam, nev2 * sizeof(double));
      totcnt += nev2;
      counts[sl] = nev2;
      //-------------------- free allocated space withing this scope
      if (lam)
        evsl_Free(lam);
      if (Y)
        evsl_Free(Y);
      if (res)
        evsl_Free(res);
      FreeASIGMABSolDirect(rat.num, solshiftdata);
      evsl_Free(solshiftdata);
      free_rat(&rat);
      evsl_Free(ind);
    } // for (sl=0; sl<nslices; sl++)
    //-------------------- free other allocated space
    fprintf(fstats, " --> Total eigenvalues found = %d\n", totcnt);
    sprintf(path, "OUT/EigsOut_Lan_MMRLanR_(%s_%s)", io.MatNam1, io.MatNam2);
    FILE *fmtout = fopen(path, "w");
    if (fmtout) {
      for (j = 0; j < totcnt; j++)
        fprintf(fmtout, "%.15e\n", alleigs[j]);
      fclose(fmtout);
    }
    evsl_Free(vinit);
    evsl_Free(sli);
    free_coo(&Acoo);
    free_csr(&Acsr);
    free_coo(&Bcoo);
    free_csr(&Bcsr);
    free_csr(&Acsr0);
    free_csr(&Bcsr0);
    evsl_Free(alleigs);
    evsl_Free(counts);
    evsl_Free(xdos);
    evsl_Free(ydos);
    if (fstats != stdout) {
      fclose(fstats);
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
