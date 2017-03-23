#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

/**----------------------------------------------------------------------
*
*    @param msteps   number of Lanczos steps
*    @param *v    initial vector
*
*    @param[out] *lmin, *lmax [lmin lmax] is the desired interval containing  
*    all eigenvalues of matrix pair (A,B)
*----------------------------------------------------------------------*/   
int LanBounds(int msteps, double *v, double *lmin, double *lmax) {
  double *alp, *bet, nbet, nalp, t, *V, *Z;
  int j, one=1, n;

  /* if users provided their own matvec function, matrix A will be ignored */
  if (evsldata.Amv) {
    n = evsldata.Amv->n;
  } else {
    n = evsldata.A->nrows;
  }
  msteps = min(n, msteps);
  Malloc(alp, msteps, double);
  Malloc(bet, msteps, double);
  Malloc(V, (msteps+1)*n, double);
  if (evsldata.ifGenEv) {
    /* storage for Z = B * V */
    Malloc(Z, (msteps+1)*n, double);
  } else {
    /* Z and V are the same */
    Z = V;
  }
  /* init vector */
  if (evsldata.ifGenEv) {
    /* B norm */
    matvec_B(v, Z);
    t = 1.0 / sqrt(DDOT(&n, v, &one, Z, &one));
    DSCAL(&n, &t, Z, &one);
  } else {
    /* 2-norm */
    t = 1.0 / DNRM2(&n, v, &one);
  }
  /* starting vector */
  DCOPY(&n, v, &one, V, &one);
  /* unit B-norm or 2-norm */
  DSCAL(&n, &t, V, &one);
  double wn = 0.0;
  /*-------------------- main Lanczos loop */
  for (j=0; j<msteps; j++) {
    int i;
    /* znew = A*v */
    /* vnew = A*v */
    matvec_A(&V[j*n], &Z[(j+1)*n]);
    /* znew = znew - bet * zold */
    /* vnew = vnew - bet * vold */
    if (j) {
      nbet = -bet[j-1];
      DAXPY(&n, &nbet, &Z[(j-1)*n], &one, &Z[(j+1)*n], &one);
    }
    /* alpha */
    /* alp = znew' * v */
    /* alp = vnew' * v */
    alp[j] = DDOT(&n, &Z[(j+1)*n], &one, &V[j*n], &one);
    wn += alp[j] * alp[j];
    /* znew = znew - alp * z */
    /* vnew = vnew - alp * v */
    nalp = -alp[j];
    DAXPY(&n, &nalp, &Z[j*n], &one, &Z[(j+1)*n], &one);
    /* full reortho for znew */
    for (i=0; i<=j; i++) {
      /* (znew, v) */
      /* (vnew, v) */
      t = DDOT(&n, &Z[(j+1)*n], &one, &V[i*n], &one);
      double mt = -t;
      DAXPY(&n, &mt, &Z[i*n], &one, &Z[(j+1)*n], &one);
    }
    if (evsldata.ifGenEv) {
      /* vnew = B \ znew */
      evsldata.Bsol->func(&Z[(j+1)*n], &V[(j+1)*n], evsldata.Bsol->data);
    }
    /* beta = (vnew, znew) */
    /* beta = (vnew, vnew) */
    bet[j] = DDOT(&n, &V[(j+1)*n], &one, &Z[(j+1)*n], &one);
    if (bet[j]*(j+1) < orthTol*wn) {
      fprintf(stdout, "lanbounds: lucky break, j=%d, beta=%e, break\n", j, bet[j]);
      msteps = j + 1;
      break;
    }
    wn += 2.0 * bet[j];
    bet[j] = sqrt(bet[j]);
    t = 1.0 / bet[j];
    /* vnew = vnew / bet */
    DSCAL(&n, &t, &V[(j+1)*n], &one);
    if (evsldata.ifGenEv) {
      /* znew = znew / bet */
      DSCAL(&n, &t, &Z[(j+1)*n], &one);
    }
  } /* main loop */

  double bottomBeta = bet[msteps-1];
  double *S, *ritzVal;
  Malloc(S, msteps*msteps, double);
  Malloc(ritzVal, msteps, double);
  /*-------------------- diagonalize tridiagonal matrix */
  SymmTridEig(ritzVal, S, msteps, alp, bet);
  /*-------------------- 'safe' bounds */
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
  /*-------------------- done */
  free(alp);
  free(bet);
  free(V);
  if (evsldata.ifGenEv) {
    free(Z);
  }
  free(S);
  free(ritzVal);

  return 0;
}

