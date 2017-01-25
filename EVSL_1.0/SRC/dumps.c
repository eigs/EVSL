#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

//-------------------- miscellaneous functions for I/O 
//                     and for debugging 

void save_mtx_basic(int nrow, int ncol, int *ia, 
                    int *ja, double *a, const char *fn) {
  int i,j,nnz;
  FILE *fp = fopen(fn, "w");

  nnz = ia[nrow];
  assert(ia[0] == 0);
  fprintf(fp, "%s\n", "%%MatrixMarket matrix coordinate real general");
  fprintf(fp, "%d %d %d\n", nrow, ncol, nnz);
  for (i=0; i<nrow; i++) {
    for (j=ia[i]; j<ia[i+1]; j++) {
      fprintf(fp, "%d %d %.15e\n", i+1, ja[j]+1, a[j]);
    }
  }
  fclose(fp);
}

void savemat(csrMat *A, const char *fn) {
  fprintf(stdout, " * saving a matrix into %s\n", fn);
  save_mtx_basic(A->nrows, A->ncols, A->ia, A->ja, A->a, fn);
}

void save_vec(int n, double *x, const char fn[]) {
  fprintf(stdout, " * saving a vector into %s\n", fn);
  FILE *fp = fopen(fn, "w");
  fprintf(fp, "%s %d\n", "%", n);
  int i;
  for (i=0; i<n; i++) {
    fprintf(fp, "%.15e\n", x[i]);
  }
  fclose(fp);
}

void savedensemat(double *A, int lda, int m, int n, const char *fn) {
  fprintf(stdout, " * saving a matrix into %s\n", fn);
  FILE *fp = fopen(fn, "w");
  int i,j;
  for (i=0; i<m; i++) {
    for (j=0; j<n; j++) {
      fprintf(fp, "%.15e ", A[i+j*lda]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
