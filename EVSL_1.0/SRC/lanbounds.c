#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

int LanBounds(csrMat *A, int msteps, double *v, double *lmin, double *lmax) {
/*----------------------------------------------------------------------
/    INPUT:
/    csrMat *A,   = matrix A
/    int msteps,  = number of Lanczos steps
/    double *v,   = initial vector
/    RETURN:
/    double *lmin,= [lmin lmax] is the desired interval containing  
/    double *lmax   all eigenvalues of A
*----------------------------------------------------------------------*/   
  double *alp, *bet, nbet, nalp, t, *V;
  int one=1, n;

  n = A->nrows;
  Malloc(alp, msteps, double);
  Malloc(bet, msteps, double);
  Malloc(V, (msteps+1)*n, double);

  t = DDOT(&n, v, &one, v, &one);
  t = 1.0 / sqrt(t);
  DSCAL(&n, &t, v, &one);
  DCOPY(&n, v, &one, V, &one);
  double wn = 0.0; 
/*-------------------- main Lanczos loop */
  int j;
  for (j=0; j<msteps; j++) {
    // w = A*v
    matvec(A, &V[j*n], &V[(j+1)*n]);
    // w = w - bet * vold
    if (j) {
      nbet = -bet[j-1];
      DAXPY(&n, &nbet, &V[(j-1)*n], &one, &V[(j+1)*n], &one);
    }
    /* alp = w' * v */
    alp[j] = DDOT(&n, &V[(j+1)*n], &one, &V[j*n], &one);
    wn += alp[j] * alp[j];
    // w = w - alp * v
    nalp = -alp[j];
    DAXPY(&n, &nalp, &V[j*n], &one, &V[(j+1)*n], &one);
    // full reortho
    int i;
    for (i=0; i<=j; i++) {
      t = DDOT(&n, &V[(j+1)*n], &one, &V[i*n], &one);
      double mt = -t;
      DAXPY(&n, &mt, &V[i*n], &one, &V[(j+1)*n], &one);
    }
    bet[j] = DDOT(&n, &V[(j+1)*n], &one, &V[(j+1)*n], &one);
    if (bet[j]*(j+1) < orthTol*wn) {
      fprintf(stdout, "lanbounds: lucky break, j=%d, beta=%e, break\n", j, bet[j]);
      msteps = j + 1;
      break;
    }
    wn += 2.0 * bet[j];
    bet[j] = sqrt(bet[j]);
    t = 1.0 / bet[j];
    DSCAL(&n, &t, &V[(j+1)*n], &one);
  }

  double bottomBeta = bet[msteps-1];
  double *S, *ritzVal;
  Malloc(S, msteps*msteps, double);
  Malloc(ritzVal, msteps, double);
//-------------------- diagonalize tridiagonal matrix    
  SymmTridEig(ritzVal, S, msteps, alp, bet);
//-------------------- 'safe' bounds  
  double amin, amax,  x;
  amax = -INFINITY;
  amin =  INFINITY;
  for (j=0; j<msteps; j++){
    t = fabs(bottomBeta * S[(j+1)*msteps-1]);
    x = ritzVal[j]-t; 
    if (x<amin) amin = x; 
    x = ritzVal[j]+t; 
    if (x>amax) amax = x; 
  }
  *lmin = amin; 
  *lmax = amax; 
//-------------------- done 
  free(alp);
  free(bet);
  free(V);
  free(S);
  free(ritzVal);
  return 0;
}

