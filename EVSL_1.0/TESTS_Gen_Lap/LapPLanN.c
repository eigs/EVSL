#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "evsl.h"
#include "io.h"
#include "evsl_suitesparse.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv);
int lapgen(int nx, int ny, int nz, cooMat *Acoo);

int main(int argc, char *argv[]) {
  /*------------------------------------------------------------
    generates a laplacean matrix on an nx x ny x nz mesh and
    matrix B on an (nx x ny) x nz x 1 mesh,
    and computes all eigenvalues of A x = lambda B x
    in a given interval [a  b]
    The default set values are
    nx = 41; ny = 53; nz = 1;
    a = 0.4; b = 0.8;
    nslices = 1 [one slice only] 
    other parameters 
    tol [tolerance for stopping - based on residual]
    Mdeg = pol. degree used for DOS
    nvec  = number of sample vectors used for DOS
    This uses:
    Non-restart Lanczos with polynomial filtering
    ------------------------------------------------------------*/
  int n, nx, ny, nz, i, j, npts, nslices, nvec, Mdeg, nev, 
      mlan, ev_int, sl, flg, ierr;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli, *mu;
  double xintv[4];
  /* initial vector: random */
  double *vinit;
  polparams pol;
  FILE *fstats = NULL;
  if (!(fstats = fopen("OUT/LapPLanN.out","w"))) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo, Bcoo;
  csrMat Acsr, Bcsr;
  /*-------------------- Bsol */
  BSolDataSuiteSparse Bsol;
  /*-------------------- default values */
  nx   = 10;
  ny   = 10;
  nz   = 20;
  a    = 1.5;
  b    = 2.5;
  nslices = 4;
  /*-------------------- user input from command line */
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
  /*-------------------- matrix size */
  n = nx * ny * nz;
  /*-------------------- stopping tol */
  tol = 1e-8;
  /*-------------------- output the problem settings */
  fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(fstats, "Laplacian A : %d x %d x %d, n = %d\n", nx, ny, nz, n);
  fprintf(fstats, "Laplacian B : %d x %d, n = %d\n", nx*ny,  nz, n);
  fprintf(fstats, "Interval: [%20.15f, %20.15f]  -- %d slices \n", a, b, nslices);
  /*-------------------- generate 1D/3D Laplacian matrix 
   *                     saved in coo format */
  ierr = lapgen(nx, ny, nz, &Acoo);
  ierr = lapgen(nx*ny*nz, 1, 1, &Bcoo);
  /*-------------------- convert coo to csr */
  ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
  ierr = cooMat_to_csrMat(0, &Bcoo, &Bcsr);
  /*-------------------- start EVSL */
  EVSLStart();
  /*-------------------- set the left-hand side matrix A */
  SetAMatrix(&Acsr);
  /*-------------------- set the right-hand side matrix B */
  SetBMatrix(&Bcsr);
  /*-------------------- use SuiteSparse as the solver for B */
  SetupBSolSuiteSparse(&Bcsr, &Bsol);
  /*-------------------- set the solver for B  and L^{T}*/
  SetBSol(BSolSuiteSparse, (void *) &Bsol);
  SetLTSol(LTSolSuiteSparse); 
  /*-------------------- for generalized eigenvalue problem */
  SetGenEig();
  /*-------------------- step 0: get eigenvalue bounds */
  /*-------------------- initial vector */
  vinit = (double *) malloc(n*sizeof(double));
  rand_double(n, vinit);
  ierr = LanBounds(60, vinit, &lmin, &lmax);
  fprintf(fstats, "Step 0: Eigenvalue bound s for A: [%.15e, %.15e]\n",
          lmin, lmax);
  /*-------------------- interval and eig bounds */
  xintv[0] = a;
  xintv[1] = b;
  xintv[2] = lmin;
  xintv[3] = lmax;
  /*-------------------- call kpmdos to get the DOS for dividing the spectrum*/
  /*-------------------- define kpmdos parameters */
  Mdeg = 300;
  nvec = 60;
  mu = malloc((Mdeg+1)*sizeof(double));
  /*-------------------- call kpmdos */
  double t = cheblan_timer();
  ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
  t = cheblan_timer() - t;
  if (ierr) {
    printf("kpmdos error %d\n", ierr);
    return 1;
  }
  fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);
  fprintf(fstats, " estimated eig count in interval: %.15e \n",ecount);
  /*-------------------- call splicer to slice the spectrum */
  npts = 10 * ecount; 
  sli = malloc((nslices+1)*sizeof(double));

  fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n", 
          Mdeg, nvec, npts);
  ierr = spslicer(sli, mu, Mdeg, xintv, nslices,  npts);
  if (ierr) {
    printf("spslicer error %d\n", ierr);
    return 1;
  }
  printf("====================  SLICES FOUND  ====================\n");
  for (j=0; j<nslices;j++) {
    printf(" %2d: [% .15e , % .15e]\n", j+1, sli[j],sli[j+1]);
  }
  /*-------------------- # eigs per slice */
  ev_int = (int) (1 + ecount / ((double) nslices));
  /*-------------------- For each slice call ChebLanr */
  for (sl=0; sl<nslices; sl++) {
    printf("======================================================\n");
    int nev2;
    double *lam, *Y, *res;
    int *ind;
    //-------------------- 
    a = sli[sl];
    b = sli[sl+1];
    printf(" subinterval: [%.4e , %.4e]\n", a, b); 
    //-------------------- approximate number of eigenvalues wanted
    nev = ev_int+2;
    //-------------------- Dimension of Krylov subspace 
    mlan = max(4*nev, 100);
    mlan = min(mlan, n);

    //-------------------- ChebLanNr
    xintv[0] = a;     xintv[1] = b;
    xintv[2] = lmin;  xintv[3] = lmax;
    //-------------------- set up default parameters for pol.      
    set_pol_def(&pol);
    //-------------------- this is to show how you can reset some of the
    //                     parameters to determine the filter polynomial
    pol.damping = 0;
    //-------------------- use a stricter requirement for polynomial
    pol.thresh_int = 0.25;
    pol.thresh_ext = 0.15;
    pol.max_deg  = 300;
    // pol.deg = 20 //<< this will force this exact degree . not recommended
    //                   it is better to change the values of the thresholds
    //                   pol.thresh_ext and plot.thresh_int
    //-------------------- Now determine polymomial to use
    find_pol(xintv, &pol);       

    fprintf(fstats, " polynomial deg %d, bar %e gam %e\n",
            pol.deg,pol.bar, pol.gam);
    //-------------------- then call ChenLanNr
    ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, 
                     &lam, &Y, &res, fstats);
    if (ierr) {
      printf("ChebLanNr error %d\n", ierr);
      return 1;
    }

    /* sort the eigenvals: ascending order
     * ind: keep the orginal indices */
    ind = (int *) malloc(nev2*sizeof(int));
    sort_double(nev2, lam, ind);
    printf(" number of eigenvalues found: %d\n", nev2);
    /* print eigenvalues */
    fprintf(fstats, "    Eigenvalues in [a, b]\n");
    fprintf(fstats, "    Computed [%d]        ||Res||\n", nev2);
    for (i=0; i<nev2; i++) {
      fprintf(fstats, "% .15e  %.1e\n", lam[i], res[ind[i]]);
      if (i>50) {
        fprintf(fstats,"                        -- More not shown --\n");
        break;
      } 
    }
    //-------------------- free allocated space withing this scope
    if (lam)  free(lam);
    if (Y)  free(Y);
    if (res)  free(res);
    free_pol(&pol);
    free(ind);
  } //for (sl=0; sl<nslices; sl++)
  //-------------------- free other allocated space 
  free(vinit);
  free(sli);
  free_coo(&Acoo);
  free_csr(&Acsr);
  free_coo(&Bcoo);
  free_csr(&Bcsr);
  FreeBSolSuiteSparseData(&Bsol);
  free(mu);
  fclose(fstats);
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

