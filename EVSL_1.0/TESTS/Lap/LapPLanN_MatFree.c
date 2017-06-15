#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "evsl.h"
#include "io.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv);
int lapgen(int nx, int ny, int nz, cooMat *Acoo);
int exeiglap3(int nx, int ny, int nz, double a, double b, int *m, double **vo);
/* Example for matrix-free solvers: 2D/3D (constant coefficient) Laplacian matrices
 * Ax = \lambda x
 * The matrix A is not formed.
 * We provide a matvec routine, which only needs the stencil and grid sizes 
 * The matvec routine and the associated data will need to be registered */
/* matvec routine [it must be of this prototype] */
void Lap2D3DMatvec(double *x, double *y, void *data);
/* datatype for performing matvec for Laplacians */
typedef struct _lapmv_t {
  int nx, ny, nz;
  double *stencil;
} lapmv_t;

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
      mlan,  ev_int, sl, flg, ierr;
  /* find the eigenvalues of A in the interval [a,b] */
  double a, b, lmax, lmin, ecount, tol,   *sli, *mu;
  double xintv[4];
  double *vinit;
  polparams pol;
  FILE *fstats = NULL;
  if (!(fstats = fopen("OUT/LapPLanN_MatFree.out","w"))) {
    printf(" failed in opening output file in OUT/\n");
    fstats = stdout;
  }
  /*-------------------- default values */
  nx   = 15;
  ny   = 16;
  nz   = 25;
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
  //------ 5pt or 7pt stencil: [self, x-1, x+1, y-1, y+1, z-1, z+1]
  // assume that grid points are ordered by x first then y and then z
  double stencil[7] = {0.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  stencil[0] = nz <= 1 ? 4.0 : 6.0;
  /* this struct will be passed to the matvec routine */
  lapmv_t lapmv;
  lapmv.nx = nx;  lapmv.ny = ny;  lapmv.nz = nz;  lapmv.stencil = stencil;
  //-------------------- eigenvalue bounds set by hand.
  lmin = 0.0;  
  lmax = nz == 1 ? 8.0 : 12.0;
  xintv[0] = a;
  xintv[1] = b;
  xintv[2] = lmin;
  xintv[3] = lmax;
  tol  = 1e-8;
  n = nx * ny * nz;
  /*-------------------- output the problem settings */
  fprintf(fstats, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  fprintf(fstats, "Laplacian: %d x %d x %d, n = %d\n", nx, ny, nz, n);
  fprintf(fstats, "Interval: [%20.15f, %20.15f]  -- %d slices \n", a, b, nslices);
  /*-------------------- step 0: get eigenvalue bounds */
  fprintf(fstats, "Step 0: Eigenvalue bound s for A: [%.15e, %.15e]\n", lmin, lmax);
  /*-------------------- call kpmdos to get the DOS for dividing the spectrum*/
  /*-------------------- define kpmdos parameters */
  Mdeg = 300;
  nvec = 60;
  /*-------------------- start EVSL */
  EVSLStart();
  /*-------------------- without forming the matrix, 
   *                     just setup the matvec function and data */
  SetAMatvec(n, &Lap2D3DMatvec, (void*) &lapmv);
  /*-------------------- call kpmdos */
  mu = malloc((Mdeg+1)*sizeof(double));
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

  fprintf(fstats,"DOS parameters: Mdeg = %d, nvec = %d, npnts = %d\n",Mdeg, nvec, npts);
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
  vinit = (double *) malloc(n*sizeof(double));
  rand_double(n, vinit);
  //-------------------- debug only :
  //  save_vec(n, vinit, "OUT/vinit.mtx");
  //-------------------- For each slice call ChebLanr
  for (sl =0; sl<nslices; sl++){
    printf("======================================================\n");
    int nev2;
    double *lam, *Y, *res;
    int *ind;
    int nev_ex;
    double *lam_ex;
    //-------------------- 
    a = sli[sl];
    b = sli[sl+1];
    printf(" subinterval: [%.4e , %.4e]\n", a, b); 
    //-------------------- approximate number of eigenvalues wanted
    nev = ev_int+2;
    //-------------------- Dimension of Krylov subspace 
    mlan = max(4*nev, 100);
    mlan = min(mlan, n);
    //-------------------- ChebLanTr
    xintv[0] = a;     xintv[1] = b;
    xintv[2] = lmin;  xintv[3] = lmax;
    //-------------------- set up default parameters for pol.      
    set_pol_def(&pol);
    //-------------------- this is to show how you can reset some of the
    //                     parameters to determine the filter polynomial
    pol.damping = 0;
    //-------------------- use a stricter requirement for polynomial
    pol.thresh_int = 0.5;
    pol.thresh_ext = 0.15;
    pol.max_deg  = 300;
    // pol.deg = 20 //<< this will force this exact degree . not recommended
    //                   it is better to change the values of the thresholds
    //                   pol.thresh_ext and plot.thresh_int
    //-------------------- Now determine polymomial to use
    find_pol(xintv, &pol);       

    fprintf(fstats, " polynomial deg %d, bar %e gam %e\n",pol.deg,pol.bar, pol.gam);
    //-------------------- then call ChenLanNr
    ierr = ChebLanNr(xintv, mlan, tol, vinit, &pol, &nev2, &lam, &Y, &res, fstats);
    if (ierr) {
      printf("ChebLanNr error %d\n", ierr);
      return 1;
    }

    /* [compute residual] already computed in res */
    /* sort the eigenvals: ascending order
     * ind: keep the orginal indices */
    ind = (int *) malloc(nev2*sizeof(int));
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
    for (i=0; i<max(nev2, nev_ex); i++) {
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
    if (lam)  free(lam);
    if (Y)  free(Y);
    if (res)  free(res);
    free_pol(&pol);
    free(ind);
    free(lam_ex);
  }
  //-------------------- free other allocated space 
  free(vinit);
  free(sli);
  free(mu);
  fclose(fstats);
  /*-------------------- finalize EVSL */
  EVSLFinish();
  return 0;
}

/*----------------- external matvec routine provided by users */
void Lap2D3DMatvec(double *x, double *y, void *data) {
  /* y = A * x
   * data: pointer to a struct that contains all needed data
   */ 
  lapmv_t *lapmv = (lapmv_t *) data;
  int nx = lapmv->nx;
  int ny = lapmv->ny;
  int nz = lapmv->nz;
  double *stencil = lapmv->stencil;
  int i,j,k,p;

  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
        p = k*nx*ny + j*nx + i;
        y[p] = stencil[0] * x[p];
        // x-1, x+1
        if (i>0)    { y[p] += stencil[1] * x[p-1]; }
        if (i<nx-1) { y[p] += stencil[2] * x[p+1]; }
        // y-1, y+1
        if (j>0)    { y[p] += stencil[3] * x[p-nx]; }
        if (j<ny-1) { y[p] += stencil[4] * x[p+nx]; }
        // z-1, z+1
        if (k>0)    { y[p] += stencil[5] * x[p-nx*ny]; }
        if (k<nz-1) { y[p] += stencil[6] * x[p+nx*ny]; }
      }
    }
  }
}

