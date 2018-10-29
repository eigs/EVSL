#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "internal_header.h"

/**
 * @file dumps.c
 * @brief Miscellaneous functions used for DOS based functions
 */
//-------------------- miscellaneous functions for I/O
//                     and for debugging

/**
 * @brief Saves a matrix in MatrixMarket format
 *
 * @param[in] nrow Number of rows in matrix
 * @param[in] ncol Number of cols in matrix
 * @param[in] ia Row pointers
 * @param[in] ja Column indices
 * @param[in] a Values
 * @param[in] fn filename
 */
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

/**
 * @brief Saves a csr matrix
 * @param[in] A csr matrix to save
 * @param[in] fn filename
 */
void savemat(csrMat *A, const char *fn) {
  fprintf(stdout, " * saving a matrix into %s\n", fn);
  save_mtx_basic(A->nrows, A->ncols, A->ia, A->ja, A->a, fn);
}

/**
 * @brief Saves a vector
 * @param[in] n Vector size
 * @param[in] x Vector to save
 * @param[in] fn filename
 */
void save_vec(int n, const double *x, const char fn[]) {
  fprintf(stdout, " * saving a vector into %s\n", fn);
  FILE *fp = fopen(fn, "w");
  fprintf(fp, "%s %d\n", "%", n);
  int i;
  for (i=0; i<n; i++) {
    fprintf(fp, "%.15e\n", x[i]);
  }
  fclose(fp);
}

/**
 *
 * @brief Saves a dense matrix
 * @param[in] A Matrix to save
 * @param[in] lda leading dimension
 * @param[in] m num rows
 * @param[in] n num cols
 * @param[in] fn filename
 */
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
