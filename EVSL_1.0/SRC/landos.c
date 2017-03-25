#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"
#include "vector.h"


/**----------------------------------------------------------------------
 *
 *    @param[in] *A    matrix A
 *    @param[in] nvec]  number of sample vectors used
 *    @param[in] msteps number of Lanczos steps
 *    @param[in] npts number of sample points used ofr the curve
 *
 *    @param[out] xdos Length-npts long vector, x-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *    @param[out] ydos Length-npts long vector, y-coordinate points for
 *    plotting the DOS. Must be preallocated.
 *
 *
 *    @warning
 *    WARNING COMPLETELY NONFUNCTIONAL AT THIS POINt
 *    CURRENTLY BEING CONSTRUCTED
 *
 *----------------------------------------------------------------------*/   

int LanDos(csrMat *A, int nvec, int msteps, int npts, double* xdos, double* ydos) {
    double *alp, *bet, nbet, nalp, t, *V;
    double *v;
    n = A->nrows;
    Malloc(alp, msteps, double);
    Malloc(bet, msteps, double);
    Malloc(V, (msteps+1)*n, double);


    //Variables that persist through iterations
    Malloc(v, n, double);
    double sigma2;
    double width;
    double* y;
    //double* ydos;
    //double* xdos;
    //Calloc(xdos, npts, double);
    //Calloc(ydos, npts, double);
    Calloc(y, npts, double);

    for(int m = 0; m <nvec; m++) {
        rand_double(n, v); // w = randn(size(A,1),1);

        t = DDOT(&n, v, &one, v, &one);
        t = 1.0 / sqrt(t); //t is now the number required to normalize vector.
        DSCAL(&n, &t, v, &one);
        DCOPY(&n, v, &one, V, &one); //v = w/norm(w);
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
        //theta = ritzVal = sorted eigenvalues

        if(m == 0) {  //On first iteration, set sigma2, width, xdos and y for future loops
            double lm = ritzVal[0]; //lm = theta(1)
            double lM = ritzVal[msteps-1]; //lM = theta(k)
            double kappa = 1.25; 
            int M = min(msteps, 30);
            double H = (lM - lm)/(M-1);
            double sigma = H / sqrt(8 * log(kappa));
            sigma2 = 2 * simga^2;
            //If gaussian small than tol ignore point.
            double tol = 1e-04;
            width = sigma * sqrt(-2 * log(tol));
            linspace(lm,LM,1,xdos);//xdos = linspace(lm,lM, npts);
            memset(y,0,npts*sizeof(y[0])); //y = zeros(size(xdos));
        }

        //Generate DOS from small gaussians centered at the ritz values
        for(int i = 0; i < msteps; i++) {
            double t =  ritzVal(i);
            //Todo ind = find(abs(xdos - t) < width);
            int* ind;
            int numind = 0;
            for(int j = 0; j < npts; j++) {  //Calculate number of elements matching pattern
                if(abs(xdos[j] - t) < width) {
                    numind++;
                }
            }
            Calloc(ind, numind, int);
            int numPlaced = 0;
            for(int j = 0; j < npts; j++) {  //Place the elements matching said pattern
                if(abs(xdos[j] - t) < width) {
                    ind[numPlaced++] = j;
                }
            }
            //ind now is = find(abs(xdos - t) < width);

            //for(int j = 0; j < numind; j++) {
            //y[ind] = y[ind] + gamma2[i] * exp((-(xdos[ind] - t) * (xdos[ind] - t)) / sigma2);

            //}
        }

        double sum = 0;
        for(int i = 0; i < npts; i++) {
            sum += ydos[i];
        }
        for(int i = 0; i < npts; i++) {
            ydos[i] /= (sum * (xdos[1] - xdos[0]));
        }
    }
    //TODO ydos = y(:);
    //TODO ydos = ydos / (sum(ydos)*(xdos(2)-xdos(1)));
}
