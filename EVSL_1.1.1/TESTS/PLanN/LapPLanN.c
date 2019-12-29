#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "evsl.h"
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
    Non-restart Lanczos with polynomial filtering
    ------------------------------------------------------------*/
  int n, nx, ny, nz, i, j, npts, nslices, nvec, Mdeg, nev,
      mlan, ev_int, sl, flg, ierr;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol, *sli, *mu;
  double xintv[4];
  double *vinit;
  polparams pol;

  struct stat st = {0}; /* Make sure OUT directory exists */
  if (stat("OUT", &st) == -1) {
    mkdir("OUT", 0750);
  }

  FILE *fstats = NULL;
  if (!(fstats = fopen("OUT/LapPLanN.out","w"))) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  /*-------------------- matrix A: coo format and csr format */
  cooMat Acoo;
  csrMat Acsr;
  /*-------------------- default values */
  nx   = 16;
  ny   = 16;
  nz   = 20;
  a    = 0.4;
  b    = 0.8;
  nslices = 4;
  //-----------------------------------------------------------------------
  //-------------------- reset some default values from command line
  /* user input from command line */
  flg = findarg("help", NA, NULL, argc, argv);
  if (flg) {
    printf("Usage: ./testL.ex -nx [int] -ny [int] -nz [int] -a [double] -b [double] -nslices [int]\n");
    return 0;
  }
  /* print output on screen */
  flg = findarg("stdout", NA, NULL, argc, argv);
  if (flg) {
     fstats = stdout;
  }
  findarg("nx", INT, &nx, argc, argv);
  findarg("ny", INT, &ny, argc, argv);
  findarg("nz", INT, &nz, argc, argv);
  findarg("a", DOUBLE, &a, argc, argv);
  findarg("b", DOUBLE, &b, argc, argv);
  findarg("nslices", INT, &nslices, argc, argv);
  fprintf(fstats,"used nx = %3d ny = %3d nz = %3d",nx,ny,nz);
  fprintf(fstats," [a = %f  b = %f],  nslices=%2d \n",a,b,nslices);
  //-------------------- eigenvalue bounds set by hand.
  lmin = 0.0;
  lmax = nz == 1 ? 8.0 : 12.0;
  xintv[0] = a;
  xintv[1] = b;
  xintv[2] = lmin;
  xintv[3] = lmax;
  tol  = 1e-6;
  n = nx * ny * nz;
  /*-------------------- output the problem settings */
  fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(fstats, "Laplacian: %d x %d x %d, n = %d\n", nx, ny, nz, n);
  fprintf(fstats, "Interval: [%20.15f, %20.15f]  -- %d slices \n", a, b, nslices);
  /*-------------------- generate 2D/3D Laplacian matrix
   *                     saved in coo format */
  ierr = lapgen(nx, ny, nz, &Acoo);
  /*-------------------- convert coo to csr */
  ierr = cooMat_to_csrMat(0, &Acoo, &Acsr);
  /*-------------------- step 0: get eigenvalue bounds */
  fprintf(fstats, "Step 0: Eigenvalue bounds for A: [%.15e, %.15e]\n", lmin, lmax);
  /*-------------------- call landos to get the DOS for dividing the spectrum*/
  /*-------------------- define landos parameters */
  Mdeg = 300;
  nvec = 60;
#ifdef EVSL_USING_CUDA_GPU
  evsl_device_query(0);
#endif
  /*-------------------- start EVSL */
  EVSLStart();
#ifdef EVSL_USING_CUDA_GPU
  /*-------------------- set matrix A (csr on GPU) */
  csrMat Acsr_gpu;
  evsl_create_csr_gpu(&Acsr, &Acsr_gpu);
  SetAMatrix_device_csr(&Acsr_gpu);
#else
  /*-------------------- set matrix A */
  SetAMatrix(&Acsr);
#endif
  /*-------------------- call LanDos */
  mu = evsl_Malloc(Mdeg+1, double);
  double t = evsl_timer();
#if 0
  //-------------------- number of points - determines fine-ness of slices
  npts = 100*nslices;
  double *xdos, *ydos;
  xdos = evsl_Malloc(npts, double);
  ydos = evsl_Malloc(npts, double);

  //fprintf(stdout," %d %d \n",npts,nslices);
  LanDos(nvec, Mdeg, npts, xdos, ydos, &ecount, xintv);
  //fprintf(stdout," %f \n",ecount);

  //for (j=0; j<npts;j++) {
  //  printf(" %10.4f %10.4f \n",xdos[j],ydos[j]);}

  t = evsl_timer() - t;

  fprintf(fstats, " Time to build DOS (Landos) was : %10.2f  \n",t);
  fprintf(fstats, " estimated eig count in interval: %.15e \n",ecount);
  //-------------------- call splicer to slice the spectrum

  sli = evsl_Malloc(nslices+1, double);

  fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",
          Mdeg, nvec, npts);
  spslicer2(xdos, ydos, nslices,  npts, sli);
  evsl_Free(xdos);
  evsl_Free(ydos);
#else
  ierr = kpmdos(Mdeg, 1, nvec, xintv, mu, &ecount);
  t = evsl_timer() - t;
  if (ierr) {
    printf("kpmdos error %d\n", ierr);
    return 1;
  }
  fprintf(fstats, " Time to build DOS (kpmdos) was : %10.2f  \n",t);
  fprintf(fstats, " estimated eig count in interval: %.15e \n",ecount);
  fprintf(stdout, " estimated eig count in interval: %.15e \n",ecount);
  //-------------------- call spslicer to slice the spectrum
  npts = 10 * ecount;
  sli = evsl_Malloc(nslices+1, double);

  fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",Mdeg, nvec, npts);
  ierr = spslicer(sli, mu, Mdeg, xintv, nslices,  npts);
#endif
  //-------------------- slicing done
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
  vinit = evsl_Malloc_device(n, double);
  rand_double_device(n, vinit);
  //-------------------- debug only :
  //  save_vec(n, vinit, "OUT/vinit.mtx");
  //-------------------- For each slice call ChebLanr
  for (sl=0; sl<nslices; sl++){
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
    //-------------------- approximate number of eigenvalues wanted
    nev = ev_int+2;
    //-------------------- Dimension of Krylov subspace
    mlan = evsl_max(5*nev, 300);
    mlan = evsl_min(mlan, n);
    //-------------------- Interval
    xintv[0] = a;     xintv[1] = b;
    xintv[2] = lmin;  xintv[3] = lmax;
    //-------------------- set up default parameters for pol.
    set_pol_def(&pol);
    //-------------------- this is to show how you can reset some of the
    //                     parameters to determine the filter polynomial
    pol.damping = 2;
    //-------------------- use a stricter requirement for polynomial
    pol.thresh_int = 0.8;
    pol.thresh_ext = 0.2;
    pol.max_deg  = 3000;
    // pol.deg = 20 //<< this will force this exact degree . not recommended
    //                   it is better to change the values of the thresholds
    //                   pol.thresh_ext and plot.thresh_int
    //-------------------- Now determine polymomial to use
    find_pol(xintv, &pol);

    fprintf(fstats, " polynomial [type = %d], deg %d, bar %e gam %e\n",
            pol.type, pol.deg, pol.bar, pol.gam);
    /*-------------------- then call ChenLanNr: note that when using CUDA
                           vinit (input) and Y (output) are device memory */
    ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, &lam, &Y, &res, fstats);
    if (ierr) {
      printf("ChebLanNr error %d\n", ierr);
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
    fprintf(fstats, "    Computed [%d]       ||Res||              Exact [%d]", nev2, nev_ex);
    if (nev2 == nev_ex) {
      fprintf(fstats, "                 Err");
    }
    fprintf(fstats, "\n");
    for (i=0; i < evsl_max(nev2, nev_ex); i++) {
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
    evsl_Free(lam_ex);
    StatsPrint(fstats);
  }
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

