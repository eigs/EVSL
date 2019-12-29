#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include "evsl.h"
#include "evsl_direct.h"
#include "io.h"
#include "lapl.h"

int main(int argc, char *argv[]) {
  /*------------------------------------------------------------
    generates a laplacean matrix on an nx x ny x nz mesh
    and computes all eigenvalues in a given interval [a  b]
    The default set values are
    nx = 41; ny = 53; nz = 1;
    a = 0.4; b = 0.8;
    nslices = 1 [one slice only]
    other parameters
    tol [tolerance for stopping - based on residual]
    Mdeg = pol. degree used for DOS
    nvec  = number of sample vectors used for DOS
    This uses:
    Thick-restart Lanczos with rational filtering
    ------------------------------------------------------------*/
  int n, nx, ny, nz, i, j, npts, nslices, nvec, Mdeg, nev, mlan, max_its,
      ev_int, sl, flg, ierr;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli, *mu;
  double xintv[4];
  /* initial vector: random */
  double *vinit;
  /* parameters for rational filter */
  int num = 1;         // number of poles used for each slice
  int pow = 2;         // multiplicity of each pole
  double beta = 0.01;  // beta in the LS approximation

  struct stat st = {0}; /* Make sure OUT directory exists */
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0750);
  }

  FILE *fstats = NULL;
  if (!(fstats = fopen("OUT/LapRLanR.out", "w"))) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
  /*-------------------- default values */
  nx = 41;
  ny = 53;
  nz = 8;
  a = 0.4;
  b = 0.8;
  nslices = 4;
  //-----------------------------------------------------------------------
  //-------------------- reset some default values from command line [Yuanzhe/]
  /* user input from command line */
  flg = findarg("help", NA, NULL, argc, argv);
  if (flg) {
    printf("Usage: ./testL.ex -nx [int] -ny [int] -nz [int] -a [double] -b "
           "[double] -nslices [int]\n");
    return 0;
  }
  findarg("nx", INT, &nx, argc, argv);
  findarg("ny", INT, &ny, argc, argv);
  findarg("nz", INT, &nz, argc, argv);
  findarg("a", DOUBLE, &a, argc, argv);
  findarg("b", DOUBLE, &b, argc, argv);
  findarg("nslices", INT, &nslices, argc, argv);
  fprintf(fstats, "used nx = %3d ny = %3d nz = %3d", nx, ny, nz);
  fprintf(fstats, " [a = %4.2f  b= %4.2f],  nslices=%2d \n", a, b, nslices);
  //-------------------- eigenvalue bounds set by hand.
  lmin = 0.0;
  lmax = nz == 1 ? 8.0 : 12.0;
  xintv[0] = a;
  xintv[1] = b;
  xintv[2] = lmin;
  xintv[3] = lmax;
  tol = 1e-8;
  n = nx * ny * nz;
  /*-------------------- generate 2D/3D Laplacian matrix
   *                     saved in coo format */
  ierr = lapgen(nx, ny, nz, &Acoo);
  /*-------------------- convert coo to csr */
  ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
  /* output the problem settings */
  fprintf(fstats,
          "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(fstats, "Laplacian: %d x %d x %d, n = %d, nnz = %d\n", nx, ny, nz, n,
          Acoo.nnz);
  fprintf(fstats, "Interval: [%20.15f, %20.15f]  -- %d slices \n", a, b,
          nslices);
  fprintf(fstats, "Step 0: Eigenvalue bound s for A: [%.15e, %.15e]\n", lmin,
          lmax);
  /*-------------------- call kpmdos to get the DOS for dividing the spectrum*/
  /*-------------------- define kpmdos parameters */
  Mdeg = 300;
  nvec = 60;
  /*-------------------- start EVSL */
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
  EVSLStart();
  /*-------------------- set the left-hand side matrix A */
#ifdef EVSL_USING_CUDA_GPU
  csrMat Acsr_gpu;
  evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
  SetAMatrix_device_csr(&Acsr_gpu);
#else
  SetAMatrix(&Acsr);
#endif
  //-------------------- call kpmdos
  mu = evsl_Malloc(Mdeg+1, double);
  double t = evsl_timer();
  ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
  t = evsl_timer() - t;
  if (ierr) {
    printf("kpmdos error %d\n", ierr);
    return 1;
  }
  fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f \n", t);
  fprintf(fstats, " estimated eig count in interval: %10.2e \n", ecount);
  //-------------------- call splicer to slice the spectrum
  npts = 10 * ecount;
  sli = evsl_Malloc(nslices+1, double);
  ierr = spslicer(sli, mu, Mdeg, xintv, nslices, npts);
  if (ierr) {
    printf("spslicer error %d\n", ierr);
    return 1;
  }
  printf("====================  SLICES FOUND  ====================\n");
  for (j = 0; j < nslices; j++) {
    printf(" %2d: [% .15e , % .15e]\n", j + 1, sli[j], sli[j + 1]);
  }
  //-------------------- # eigs per slice
  ev_int = (int)(1 + ecount / ((double)nslices));
  //-------------------- initial vector
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);
  //-------------------- For each slice call RatLanrTr
  for (sl = 0; sl < nslices; sl++) {
    printf("======================================================\n");
    int nev2;
    double *lam, *Y, *res;
    int *ind;
    int nev_ex;
    double *lam_ex;
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
    SetupASIGMABSolDirect(&Acsr, NULL, num, rat.zk, solshiftdata);
    /*------------ set the solver for A-sI in rat */
    for (i=0; i<rat.num; i++) {
       SetASigmaBSol(&rat, i, ASIGMABSolDirect, solshiftdata[i]);
    }
    //-------------------- approximate number of eigenvalues wanted
    nev = ev_int + 2;
    //-------------------- Dimension of Krylov subspace and maximal iterations
    mlan = evsl_max(4 * nev, 300);
    mlan = evsl_min(mlan, n);
    max_its = 3 * mlan;
    //-------------------- RationalLanTr: note that when using CUDA
    //                     vinit (input) and Y (output) are device memory
    ierr = RatLanTr(mlan, nev, intv, max_its, tol, vinit, &rat, &nev2, &lam, &Y,
                    &res, fstats);
    if (ierr) {
      printf("RatLanTr error %d\n", ierr);
      return 1;
    }
    /* [compute residual] already computed in res */
    /* sort the eigenvals: ascending order
     * ind: keep the orginal indices */
    ind = evsl_Malloc(nev2, int);
    sort_double(nev2, lam, ind);
    /* compute exact eigenvalues */
    exeiglap3(nx, ny, nz, a, b, &nev_ex, &lam_ex);
    printf(" number of eigenvalues: %d, found: %d\n", nev_ex, nev2);
    /* print eigenvalues */
    fprintf(fstats,
            "                                   Eigenvalues in [a, b]\n");
    fprintf(fstats, "    Computed [%d]       ||Res||              Exact [%d]",
            nev2, nev_ex);
    if (nev2 == nev_ex) {
      fprintf(fstats, "                 Err");
    }
    fprintf(fstats, "\n");
    for (i = 0; i < evsl_max(nev2, nev_ex); i++) {
      if (i < nev2) {
        fprintf(fstats, "% .15e  %.1e", lam[i], res[ind[i]]);
      } else {
        fprintf(fstats, "                               ");
      }
      if (i < nev_ex) {
        fprintf(fstats, "        % .15e", lam_ex[i]);
      }
      if (nev2 == nev_ex) {
        fprintf(fstats, "        % .1e", lam[i] - lam_ex[i]);
      }
      fprintf(fstats, "\n");
      if (i > 50) {
        fprintf(fstats, "                        -- More not shown --\n");
        break;
      }
    }
    //-------------------- free allocated space withing this scope
    if (lam) evsl_Free(lam);
    if (Y) evsl_Free_device(Y);
    if (res) evsl_Free(res);
    FreeASIGMABSolDirect(rat.num, solshiftdata);
    evsl_Free(solshiftdata);
    evsl_Free(ind);
    evsl_Free(lam_ex);
    free_rat(&rat);
  }  // for (sl=0; sl<nslices; sl++)
  //-------------------- free other allocated space
  evsl_Free_device(vinit);
  evsl_Free(sli);
  free_coo(&Acoo);
  free_csr(&Acsr);
#ifdef EVSL_USING_CUDA_GPU
  evsl_free_csr_gpu(&Acsr_gpu);
#endif
  evsl_Free(mu);
  fclose(fstats);
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

