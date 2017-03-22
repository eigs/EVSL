#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

/**----------------------------------------------------------------------
 *
 *    @param *A    matrix A
 *    @param msteps   number of Lanczos steps
 *    @param *v    initial vector
 *
 *    @param[out] *lmin, *lmax [lmin lmax] is the desired interval containing  
 *    all eigenvalues of A
 *    WARNING COMPLETELY NONFUNCTIONAL AT THIS POINt
 *    CURRENTLY BEING CONSTRUCTED
 *
 *----------------------------------------------------------------------*/   

int LanDos(csrMat *A, int msteps, double *v, double *lmin, double *lmax) {
    double *alp, *bet, nbet, nalp, t, *V;
    n = A->nrows;
    Malloc(alp, msteps, double);
    Malloc(bet, msteps, double);
    Malloc(V, (msteps+1)*n, double);

    double* xdos =  0;//Linear space TODO
    double* y = (double*) calloc(nvec, double); // TODO

    t = DDOT(&n, v, &one, v, &one);
    t = 1.0 / sqrt(t);
    DSCAL(&n, &t, v, &one);
    DCOPY(&n, v, &one, V, &one);
    double wn = 0.0; 
    /*-------------------- main Lanczos loop */
    int j;
    for (j=0; j<msteps; j++) {
        // w = A*v
        matvec_genev(A, &V[j*n], &V[(j+1)*n]);
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
    //theta = ritzVal

    if(false = 1) { //TODO m == 1
        double lm = ritzVal[0];
        double lM = ritzVal[msteps-1];
        double kappa = 1.25;
        int M = math(msteps, 30);
        double H = (lM - lm)/(M-1);
        double sigma = H / sqrt(8 * log(kappa));
        double sigma2 = 2 * simga^2;
        //If gaussian small than tol ignore point.
        double tol = 1e-04;
        double width = sigma * sqrt(-2 * log(tol));
    }

    for(int i = 0; i < msteps; i++) {
        double t =  ritzVal(i);
        //Todo ind = find(abs(xdos - t) < width);
        y[ind] = y[ind] + gamma2[i] * exp((-(xdos[ind] - t) * (xdos[ind] - t)) / sigma2);
    }

double* ydos = y;
double sum = 0;
for(int i = 0; i < nvec; i++) {
    sum += ydos[i];
}
for(int i = 0; i < nvec; i++) {
    ydos[i] /= (sum * (xdos[1] - xdos[0]));
}
}
