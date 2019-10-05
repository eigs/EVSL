#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "evsl.h"
#include "io.h"
#include "evsl_direct.h"
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
    Non-restart Lanczos with rational filtering
    ------------------------------------------------------------*/
  int n, nx, ny, nz, i, j, npts, nslices, nvec, Mdeg, nev,
      max_its, ev_int, sl, flg, ierr;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli, *mu;
  double xintv[4];
  /* initial vector: random */
  double *vinit;
  /* parameters for rational filter */
  int num = 1; // number of poles used for each slice
  int pow = 2; // multiplicity of each pole
  double beta = 0.01; // beta in the LS approximation

  struct stat st = {0}; /* Make sure OUT directory exists */
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0750);
  }

  FILE *fstats = NULL;
  if (!(fstats = fopen("OUT/LapRLanN.out","w"))) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
  /*-------------------- default values */
  nx   = 41;
  ny   = 53;
  nz   = 1;
  a    = 0.4;
  b    = 0.8;
  nslices = 4;
  //-----------------------------------------------------------------------
  //-------------------- reset some default values from command line [Yuanzhe/]
  /* user input from command line */
  flg = findarg("help", NA, NULL, argc, argv);
  if (flg) {
    printf("Usage: ./testL.ex -nx [int] -ny [int] -nz [int] -a [double] -b [double] -nslices [int]\n");
    return 0;
  }
  findarg("nx", INT, &nx, argc, argv);
  findarg("ny", INT, &ny, argc, argv);
  findarg("nz", INT, &nz, argc, argv);
  findarg("a", DOUBLE, &a, argc, argv);
  findarg("b", DOUBLE, &b, argc, argv);
  findarg("nslices", INT, &nslices, argc, argv);
  fprintf(fstats,"used nx = %3d ny = %3d nz = %3d",nx,ny,nz);
  fprintf(fstats," [a = %4.2f  b= %4.2f],  nslices=%2d \n",a,b,nslices);
  //-------------------- eigenvalue bounds set by hand.
  lmin = 0.0;
  lmax =  ((nz == 1)? 8.0 :12) ;
  xintv[0] = a;
  xintv[1] = b;
  xintv[2] = lmin;
  xintv[3] = lmax;
  tol = 1e-6;
  n = nx * ny * nz;
  /*-------------------- generate 2D/3D Laplacian matrix
   *                     saved in coo format */
  ierr = lapgen(nx, ny, nz, &Acoo);
  /*-------------------- convert coo to csr */
  ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
  /* output the problem settings */
  fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(fstats, "Laplacian: %d x %d x %d, n = %d, nnz = %d\n", nx, ny, nz, n, Acoo.nnz);
  fprintf(fstats, "Interval: [%20.15f, %20.15f]  -- %d slices \n", a, b, nslices);
  fprintf(fstats, "Step 0: Eigenvalue bound s for A: [%.15e, %.15e]\n", lmin, lmax);
  /*-------------------- call kpmdos to get the DOS for dividing the spectrum*/
  /*-------------------- define kpmdos parameters */
  Mdeg = 300;
  nvec = 60;
  /*-------------------- start EVSL */
  EVSLStart();
  /*-------------------- set the left-hand side matrix A */
  SetAMatrix(&Acsr);
  //-------------------- call kpmdos
  mu = evsl_Malloc(Mdeg+1, double);
  double t = evsl_timer();
  ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
  t = evsl_timer() - t;
  if (ierr) {
    printf("kpmdos error %d\n", ierr);
    return 1;
  }
  fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);
  fprintf(fstats, " estimated eig count in interval: %10.2e \n",ecount);
  //-------------------- call splicer to slice the spectrum
  npts = 10 * ecount;
  sli = evsl_Malloc(nslices+1, double);
  ierr = spslicer(sli, mu, Mdeg, xintv, nslices,  npts);
  if (ierr) {
    printf("spslicer error %d\n", ierr);
    return 1;
  }
  printf("====================  SLICES FOUND  ====================\n");
  for (j=0; j<nslices;j++) {
    printf(" %2d: [% .15e , % .15e]\n", j+1, sli[j],sli[j+1]);
  }
  //-------------------- # eigs per slice
  ev_int = (int) (1 + ecount / ((double) nslices));
  //-------------------- initial vector
  vinit = evsl_Malloc(n, double);
  rand_double(n, vinit);
  //-------------------- For each slice call RatLanrNr
  for (sl=0; sl<nslices; sl++) {
    printf("======================================================\n");
    int nev2;
    double *lam, *Y, *res;
    int *ind;
    int nev_ex;
    double *lam_ex;
    //--------------------
    StatsReset();
    a = sli[sl];
    b = sli[sl+1];
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
    /*------------ use direct solver function  */
    void **solshiftdata = evsl_Malloc(rat.num, void *);
    /*------------ factoring the shifted matrices and store the factors */
    SetupASIGMABSolDirect(&Acsr, NULL, rat.num, rat.zk, solshiftdata);
    /*------------ set the solver for A-sI in rat */
    for (i=0; i<rat.num; i++) {
       SetASigmaBSol(&rat, i, ASIGMABSolDirect, solshiftdata[i]);
    }
    //-------------------- approximate number of eigenvalues wanted
    nev = ev_int+2;
    //-------------------- maximal Lanczos iterations
    max_its = evsl_max(4*nev,300);  max_its = evsl_min(max_its, n);
    //-------------------- RationalLanNr
    ierr = RatLanNr(intv, max_its, tol, vinit, &rat, &nev2, &lam,
                    &Y, &res, fstats);
    if (ierr) {
      printf("RatLanNr error %d\n", ierr);
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
    fprintf(fstats, "                                   Eigenvalues in [a, b]\n");
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
        fprintf(fstats, "        % .1e", lam[i]-lam_ex[i]);
      }
      fprintf(fstats,"\n");
      if (i>50) {
        fprintf(fstats,"                        -- More not shown --\n");
        break;
      }
    }
    //-------------------- free allocated space withing this scope
    if (lam) evsl_Free(lam);
    if (Y) evsl_Free(Y);
    if (res) evsl_Free(res);
    evsl_Free(ind);
    evsl_Free(lam_ex);
    FreeASIGMABSolDirect(rat.num, solshiftdata);
    evsl_Free(solshiftdata);
    free_rat(&rat);
    StatsPrint(fstats);
  } //for (sl=0; sl<nslices; sl++)
  //-------------------- free other allocated space
  evsl_Free(vinit);
  evsl_Free(sli);
  free_coo(&Acoo);
  free_csr(&Acsr);
  evsl_Free(mu);
  fclose(fstats);
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

