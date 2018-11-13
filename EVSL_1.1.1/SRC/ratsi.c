#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include "internal_header.h"

#define NBUF 2

/*
 * Subspace iteration with rational filter
 */

int RatSI(int nest, double *intv, int maxit, double tol, ratparams *rat,
	int *neig, double **eigvalout, double **eigvecout, double **resout) {

	if (tol > 1e-12) {
		printf("Tolerance is recommended to be no larger than 1e-12.\n");
	}

	// check residual every 'check' filter applications
	int check = 3;

	// misc
	int i;

	// constants
	char cT = 'T';
	char cN = 'N';
	char uplo = 'U';
	double done = 1.0;
	double dmone = -1.0;
	double dzero = 0.0;
	int one = 1;
	int nest2 = nest*nest;

	// problem properties
	int n = evsldata.n;
	const int ifGenEv = evsldata.ifGenEv;
	// interval
	double aa = intv[0];
	double bb = intv[1];

	// # converged and active eigenvectors
	int nlock = 0;
	int nact = nest;
	int nlock_new = 0;
	// # converged eigenvalues smaller than a or larger than b
	int nout = 0;
	int nout_new = 0;

	// generate random seed with current time
	srand(time(0)); 
	// genearte random initial subspace
	double *V; // subspace
	V = evsl_Malloc(n*nest, double);
	for (i = 0; i < n*nest; ++i){
		V[i] = rand() / ((double)RAND_MAX);
	}

	// memory for eiganvalue and residual
	double *eig, *res;
	eig = evsl_Malloc(nest, double);
	res = evsl_Malloc(nest, double);

	// store V(:,1:lock)'*V(:,nlock+1:nlock+nact) and the reduced matrices
	double *T;
	if (ifGenEv)
		T = evsl_Malloc(2*nest*nest, double);
	else
		T = evsl_Malloc(nest*nest, double);

	// work space and temp
	double *work, *temp, t;
	int work_size = 4*n*nest;
	if (ifGenEv)
		work_size += 2*n*nest;
	work = evsl_Malloc(work_size, double);
	temp = evsl_Malloc(2*n*nest, double);

	// subspace iteration
	int it = 0;
	int find_more = 1;
	while ( (it < maxit) && (find_more) ) {

		fprintf(stdout,"it = %d\tnl = %d\tno = %d\n",it, nlock, nout);

		int nnlock = n*nlock;
		int nnact = n*nact;
		int nnout = n*nout;

		for (i = 0; i < check; ++i) {
			// apply filter to V(:,nlock+1:nest) -> temp
			RatFiltApply(n, nact, rat, V+nnlock, temp, work);
			evsl_dcopy(&nact, temp, &one, V+nnlock, &one);

			// orth against V(:,1:nlock)
			if (nlock) {
				if (ifGenEv) {
					// block mat-vec: compute B*V(:,nlock+1:nest) -> temp
					matvec_B(V+nnlock, temp, nact);
					evsl_dgemm(&cT, &cN, &nlock, &nact, &n, &done, V, &n, temp, &n, &dzero, T, &nest);
				}
				else {
					evsl_dgemm(&cT, &cN, &nlock, &nact, &n, &done, V, &n, V+nnlock, &n, &dzero, T, &nest);
				}

				evsl_dgemm(&cN, &cN, &n, &nact, &nlock, &dmone, V, &n, T, &nest, &done, V+nnlock, &n);
			}

			// orth against itself
			// ???????????????????
			// to be replaced by a new function operates directly on the input
			// orth(V+nnlock, n, nact, temp, work);
			// evsl_dcopy(&nnact, temp, &one, V+nnlock, &one);
			orth2(V+nnlock, n, n, nact, work);

		}

		// form reduced problems
		// block mat-vec: A*V(:,nlock+1:nest) -> temp
		matvec_A(V+nnlock, temp, nact);

		// At = V'*A*V -> T(1:nact,1:nact)
		evsl_dgemm(&cT, &cN, &nact, &nact, &n, &done, V+nnlock, &n, temp, &n, &dzero, T, &nest);

		if (ifGenEv) {
			// block mat-vec: B*V -> temp 
			matvec_B(V+nnlock, temp, nact);

			// Bt = V'*B*V -> T(nest+1:nest+nact,1:nact)
			evsl_dgemm(&cT, &cN, &nact, &nact, &n, &done, V+nnlock, &n, temp, &n, &dzero, T+nest2, &nest2);

			// convert to standard problem

			// cholesky factorization of B
			int info;
			evsl_dpotrf(&uplo, &nact, T+nest2, &nest, &info);

			// convert to stardard problem
			int itype = 1;
			evsl_dsygst(&itype, &uplo, &nact, T, &nest, T+nest2, &nest, &info);
		}

		// solve reduced eig problem
		SymEigenSolver(nact, T, nest, temp, n, eig+nlock);

		// recover eigenvector
		if (ifGenEv) {
			int info;
			evsl_dtrtrs(&uplo, &cN, &cN, &nact, &nact, T+nest2, &nest, temp, &n, &info);
		}
		evsl_dgemm(&cN, &cN, &n, &nact, &nact, &done, V+nnlock, &n, temp, &n, &dzero, V+nnlock, &n);

		// compute residual
		// block mat-vec: A*V(:,nlock+1:nest) -> temp
		matvec_A(V+nnlock, temp, nact);
		for (i = 0; i < nact; ++i) {
			res[i+nlock] = sqrt(evsl_ddot(&n, temp+i*n, &one, temp+i*n, &one));
		}

		if (ifGenEv) {
			// block mat-vec: B*V(:,nlock+1:nest) -> temp+nest2
			matvec_B(V+nnlock, temp+nest2, nact);
		}
		else {
			evsl_dcopy(&nnact, V+nnlock, &one, temp+nest2, &one);
		}

		for (i = 0; i < nact; ++i) {
			t = -eig[i+nlock];

			if (ifGenEv) {
				res[i+nlock] += fabs(t)*sqrt(evsl_ddot(&n, temp+nest2+i*n, &one, temp+nest2+i*n, &one));
			}
			else {
				res[i+nlock] += fabs(t);
			}

			evsl_daxpy(&n, &t, temp+nest2+i*n, &one, temp+i*n, &one);
			res[i+nlock] = sqrt(evsl_ddot(&n, temp+i*n, &one, temp+i*n, &one))/res[i+nlock];
		}

		// collect converged eigenvalue
		nlock_new = 0;
		for (i = nlock; i < nest; ++i) {

			if (res[i] < tol) {

				if (i != nlock+nlock_new) {
					// swap i and (nlock+nlock_new)
					evsl_dswap(&n, V+n*i, &one, V+nnlock+n*nlock_new, &one);

					t = eig[i];
					eig[i] = eig[nlock+nlock_new];
					eig[nlock+nlock_new] = t;

					t = res[i];
					res[i] = res[nlock+nlock_new];
					res[nlock+nlock_new] = t;
				}

				nlock_new++;
			}
		}

		// normalize newly converged eigenvectors
		if (ifGenEv) {
			// block mat-vec: B*V(:,nlock+1:nlock+nlock_new) -> temp
			matvec_B(V+nnlock, temp, nlock_new);
			for (i = 0; i < nlock_new; ++i) {
				t = 1.0/sqrt(evsl_ddot(&n, V+nnlock+i*n, &one, temp+i*n, &one));
				evsl_dscal(&n, &t, V+nnlock+i*n, &one);
			}
		}

		// check if eigval is inside interval
		nout_new = 0;
		for (i = 0; i < nlock_new; ++i) {
			if ( ( eig[nlock+i] < aa - EVSL_DBL_EPS_MULT * DBL_EPSILON ) || 
				( eig[nlock+i] > bb + EVSL_DBL_EPS_MULT * DBL_EPSILON ) ) {
				// swap nout+nout_new and nlock+i
				evsl_dswap(&n, V+nnout+nout_new*n, &one, V+nnlock+i*n, &one);

				t = eig[nout+nout_new];
				eig[nout+nout_new] = eig[nlock+i];
				eig[nlock+i] = t;

				t = res[nout+nout_new];
				res[nout+nout_new] = res[nlock+i];
				res[nlock+i] = t;

				nout_new++;
			}
		}

		nout += nout_new;
		nlock += nlock_new;
		nact -= nlock_new;
		it++;

		// check termination
		if ( ( nout >= NBUF ) || ( nlock == nest ) ) {
			find_more = 0;
		}
	}

	// post-processing
	*neig = nlock-nout;
	int nneig = n*(*neig);

	double *eigval_out, *eigvec_out, *res_out;
	eigval_out = evsl_Malloc(*neig, double);
	eigvec_out = evsl_Malloc(nneig, double);
	res_out = evsl_Malloc(*neig, double);

	evsl_dcopy(neig, eig+nout, &one, eigval_out, &one);
	evsl_dcopy(&nneig, V+n*nout, &one, eigvec_out, &one);
	evsl_dcopy(neig, res+nout, &one, res_out, &one);

	*eigvalout = eigval_out;
	*eigvecout = eigvec_out;
	*resout = res_out;

	// free memory
	evsl_Free(V);
	evsl_Free(eig);
	evsl_Free(res);
	evsl_Free(T);
	evsl_Free(work);
	evsl_Free(temp);

	// return flag
	if ( nlock == nest ) {
		// fprintf
		return -1;
	}
	if ( it == maxit ) {
		// fprintf
		return -2;
	}

	return 0;

}