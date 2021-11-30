#include <string.h>
#include "internal_header.h"
#ifdef EVSL_USING_INTEL_MKL
#include "mkl_spblas.h"
#endif

/**
 * @file spmat.c
 * @brief Sparse matrix routines
 */

/**
 * @brief convert csr to csc
 * Assume input csr is 0-based index
 * output csc 0/1 index specified by OUTINDEX      *
 * @param[in] OUTINDEX specifies if CSC should be 0/1 index
 * @param[in] nrow Number of rows
 * @param[in] ncol Number of columns
 * @param[in] job flag
 * @param[in] a Values of input matrix
 * @param[in] ia Input row pointers
 * @param[in] ja Input column indices
 * @param[out] ao Output values
 * @param[out] iao Output row pointers
 * @param[out] jao Output column indices
 */
void csrcsc(int OUTINDEX, const int nrow, const int ncol, int job,
    double *a, int *ja, int *ia,
    double *ao, int *jao, int *iao) {
  int i,k;
  for (i=0; i<ncol+1; i++) {
    iao[i] = 0;
  }
  // compute nnz of columns of A
  for (i=0; i<nrow; i++) {
    for (k=ia[i]; k<ia[i+1]; k++) {
      iao[ja[k]+1] ++;
    }
  }
  // compute pointers from lengths
  for (i=0; i<ncol; i++) {
    iao[i+1] += iao[i];
  }
  // now do the actual copying
  for (i=0; i<nrow; i++) {
    for (k=ia[i]; k<ia[i+1]; k++) {
      int j = ja[k];
      if (job) {
        ao[iao[j]] = a[k];
      }
      jao[iao[j]++] = i + OUTINDEX;
    }
  }
  /*---- reshift iao and leave */
  for (i=ncol; i>0; i--) {
    iao[i] = iao[i-1] + OUTINDEX;
  }
  iao[0] = OUTINDEX;
}

/**
 * @brief  Sort each row of a csr by increasing column
 * order
 * By double transposition
 * @param[in] A Matrix to sort
 */
void sortrow(csrMat *A) {
  /*-------------------------------------------*/
  int nrows = A->nrows;
  int ncols = A->ncols;
  int nnz = A->ia[nrows];
  // work array
  double *b;
  int *jb, *ib;
  b = evsl_Malloc(nnz, double);
  jb = evsl_Malloc(nnz, int);
  ib = evsl_Malloc(ncols+1, int);
  // double transposition
  csrcsc(0, nrows, ncols, 1, A->a, A->ja, A->ia, b, jb, ib);
  csrcsc(0, ncols, nrows, 1, b, jb, ib, A->a, A->ja, A->ia);
  // free
  evsl_Free(b);
  evsl_Free(jb);
  evsl_Free(ib);
}

/**
 * @brief  memory allocation for csr matrix
 * @param[in] nrow New number of rows
 * @param[in] ncol New number of columns
 * @param[in] nnz New number of non zeros
 * @param[in, out] csr the CSR matrix
 */
void csr_resize(int nrow, int ncol, int nnz, csrMat *csr) {
  csr->owndata = 1;
  csr->nrows = nrow;
  csr->ncols = ncol;
  csr->ia = evsl_Malloc(nrow+1, int);
  csr->ja = evsl_Malloc(nnz, int);
  csr->a = evsl_Malloc(nnz, double);
}

/**
 * @brief  memory deallocation for csr matrix
 * @param[in,out] csr Matrix to free data
 */
void free_csr(csrMat *csr) {
  /* if it does not own the data, do nothing */
  if (!csr->owndata) {
    return;
  }
  evsl_Free(csr->ia);
  evsl_Free(csr->ja);
  evsl_Free(csr->a);
}

/**
 * @brief  copy a csr matrix A into B
 * alloB: 0: will not allocate memory for B (have been alloced outside)
 *        1: will allocate memory for B (same size as A)
 *
 *  @param[in] A Source matrix
 *  @param[out] B Destination matrix
 *  @param[in] allocB Whether or not to allocate B
 */
void csr_copy(csrMat *A, csrMat *B, int allocB) {
  int nrows = A->nrows;
  int ncols = A->ncols;
  int nnz = A->ia[nrows];
  if (allocB) {
    /* allocate memory */
    csr_resize(nrows, ncols, nnz, B);
  } else {
    /* just set the sizes */
    B->nrows = nrows;
    B->ncols = ncols;
  }
  B->nnz = A->nnz;
  memcpy(B->ia, A->ia, (nrows+1)*sizeof(int));;
  memcpy(B->ja, A->ja, nnz*sizeof(int));
  memcpy(B->a,  A->a,  nnz*sizeof(double));
}

/**
 * @brief  memory deallocation for coo matrix
 * @param[in] coo Coo matrix to free
 */
void free_coo(cooMat *coo) {
  evsl_Free(coo->ir);
  evsl_Free(coo->jc);
  evsl_Free(coo->vv);
}

/**
 * @brief convert coo to csr
 * @param[in] cooidx Specify if 0 or 1 indexed
 * @param[in] coo COO matrix
 * @param[out] csr CSR matrix
 */
int cooMat_to_csrMat(int cooidx, cooMat *coo, csrMat *csr) {
  const int nnz = coo->nnz;
  csr->nnz = nnz;
  //printf("@@@@ coo2csr, nnz %d\n", nnz);
  /* allocate memory */
  csr_resize(coo->nrows, coo->ncols, nnz, csr);
  const int nrows = coo->nrows;
  /* fill (ia, ja, a) */
  int i;
  for (i=0; i<nrows+1; i++) {
    csr->ia[i] = 0;
  }
  for (i=0; i<nnz; i++) {
    int row = coo->ir[i] - cooidx;
    csr->ia[row+1] ++;
  }
  for (i=0; i<nrows; i++) {
    csr->ia[i+1] += csr->ia[i];
  }
  for (i=0; i<nnz; i++) {
    int row = coo->ir[i] - cooidx;
    int col = coo->jc[i] - cooidx;
    double val = coo->vv[i];
    int k = csr->ia[row];
    csr->a[k] = val;
    csr->ja[k] = col;
    csr->ia[row]++;
  }
  for (i=nrows; i>0; i--) {
    csr->ia[i] = csr->ia[i-1];
  }
  csr->ia[0] = 0;
  /* sort rows ? */
  sortrow(csr);
  return 0;
}

/**
 * @brief construct a csrMat copied from (ia, ja, a)
 * @param[in] nrow Number of rows
 * @param[in] ncol Number of columns
 * @param[in] ia Row pointers
 * @param[in] ja Column indices
 * @param[in] a Values
 * @param[out] A output CSR matrix
 */
int arrays_copyto_csrMat(int nrow, int ncol, int *ia, int *ja, double *a,
                         csrMat *A) {
  int nnz = ia[nrow];

  A->nnz = nnz;
  csr_resize(nrow, ncol, nnz, A);

  memcpy(A->ia, ia, (nrow+1)*sizeof(int));
  memcpy(A->ja, ja, nnz*sizeof(int));
  memcpy(A->a, a, nnz*sizeof(double));

  sortrow(A);

  return 0;
}

#if 0
/**-------------------------------------------*
 * @brief  compute the inf-norm of a csr matrix
 *-------------------------------------------*/
double dcsrinfnrm(csrMat *A){
  // computes the inf-norm of A: max abs row sum
  double ta = 0.0;
  int nrows =A->nrows, i, j;
  int *ia = A->ia;
  double *aa = A->a;
  /* for each row */
  for (i=0; i<nrows; i++) {
    /* abs row sum */
    double t = 0.0;
    for (j=ia[i]; j<ia[i+1]; j++) {
      t += fabs(aa[j]);
    }
    /* take max */
    ta = evsl_max(ta,t);
  }
  return (ta);
}
#endif

/**
 * @brief csr matrix matvec or transpose matvec, (ia, ja, a) form
 * @param[in] trans Whether or not transpose
 * @param[in] nrow Number of rows
 * @param[in] ncol Number of columns
 * @param[in] a input Values
 * @param[in] ia Row pointers
 * @param[in] ja Column indices
 * @param[in] x Input vector
 * @param[out] y Output vector
 */
void dcsrmv(char trans, int nrow, int ncol, double *a,
    int *ia, int *ja, double *x, double *y) {
  int  len, jj=nrow;
  if (trans == 'N') {
    //#pragma omp parallel for schedule(guided)
    double r;
    /*for (i=0; i<nrow; i++) {
      r = 0.0;
      for (j=ia[i]; j<ia[i+1]; j++) {
      r += a[j] * x[ja[j]];
      }
      y[i] = r;
      }
      */
    while (jj--) {
      len = *(ia+1) - *ia;
      ia++;
      r = 0.0;
      while (len--)
        r += x[*ja++]*(*a++);
      *y++ = r;
    }
  } else {
    double xi;
    int jj, len;
    jj = nrow;
    //-------------------- this is from the matvec used in FILTLAN
    //                     column oriented - gains up to 15% in time
    double *w = y;
    while (jj--)
      *w++ = 0.0;
    jj = nrow;
    //-------------------- column loop
    while (jj--){
      len = *(ia+1) - *ia;
      ia++;
      xi = *x++;
      while(len--)
        y[*ja++] += xi*(*a++);
    }
  }
}

/**
* @brief matvec for a CSR matrix, y = A*x.
* void *data points to csrMat,
* compatible form with EVSLMatvec (see struct.h)
*
* @param[in] x Input vector
* @param[out] y Output vector
* @param[in] data CSR matrix
*/
void matvec_csr(double *x, double *y, void *data) {
  csrMat *A = (csrMat *) data;
#ifdef EVSL_USING_INTEL_MKL
  char cN = 'N';
  /*
  double alp = 1.0, bet = 0.0;
  mkl_dcsrmv(&cN, &(A->nrows), &(A->ncols), &alp, "GXXCXX",
             A->a, A->ja, A->ia, A->ia+1, x, &bet, y);
  */
  mkl_cspblas_dcsrgemv(&cN, &A->nrows, A->a, A->ia, A->ja, x, y);
#else
  dcsrmv('N', A->nrows, A->ncols, A->a, A->ia, A->ja, x, y);
#endif
}


/** @brief inline function used by matadd
 * insert an element pointed by j of A (times t) to location k in C (row i)
 * */
static inline void matadd_insert(double t, csrMat *A, csrMat *C, int i, int *k,
                          int *j, int *map) {
  /* if this entry already exists in C:
   * checking if it is the first entry of this row
   * and if column indices match
   * NOTE that pointer j or k will be modified */
  if (*k > C->ia[i] && C->ja[(*k)-1] == A->ja[*j]) {
    if (map) {
      /* j maps to k-1 in C */
      map[(*j)] = (*k) - 1;
    }
    /* add to existing entry */
    C->a[(*k)-1] += t * A->a[*j];
  } else {
    if (map) {
      /* jA maps to k in C */
      map[*j] = *k;
    }
    /* create new entry */
    C->ja[*k] = A->ja[*j];
    C->a[*k] = t * A->a[*j];
    (*k)++;
  }
  (*j) ++;
}

/** @brief matrix addition C = alp * A + bet * B
 * @param[in] alp
 * @param[in] bet
 * @param[in] A
 * @param[in] B
 * @param[out] C
 * @warning the nz pattern of C will be union of those of A and B.
 * no cancellation will be considered
 * @warning A and B MUST be sorted, on output C will be sorted as well
 * @param[out] mapA (of size nnzA or null), mapB (of size nnzB or null)
 * if not null, on output mapA contains the location of each nonzero of A
 * in the CSR matrix C, i.e. mapA[i] is the position of the corresponding
 * entry in C.ja and C.a for entry in A.ja[i] and A.a[i]
 * @param[out] mapB the same as mapA
 * */
int matadd(double alp, double bet, csrMat *A, csrMat *B, csrMat *C,
           int *mapA, int *mapB) {
  int nnzA, nnzB, i, jA, jB, k;
  /* check dimension */
  if (A->nrows != B->nrows || A->ncols != B->ncols) {
    return 1;
  }
  /* nnz of A and B */
  nnzA = A->ia[A->nrows];
  nnzB = B->ia[B->nrows];
  /* alloc C [at most has nnzC = nnzA + nnzB] */
  csr_resize(A->nrows, A->ncols, nnzA+nnzB, C);
  /* nnz counter of C */
  k = 0;
  C->ia[0] = 0;
  for (i=0; i<A->nrows; i++) {
    /* open row i of A and B */
    /* merging two sorted list */
    for (jA=A->ia[i], jB=B->ia[i]; ; ) {
      if (jA < A->ia[i+1] && jB < B->ia[i+1]) {
        /* will insert the element with smaller col id */
        if (A->ja[jA] <= B->ja[jB]) {
          /* insert jA */
          matadd_insert(alp, A, C, i, &k, &jA, mapA);
        } else {
          /* instert jB */
          matadd_insert(bet, B, C, i, &k, &jB, mapB);
        }
      } else if (jA == A->ia[i+1]) {
        for (; jB < B->ia[i+1]; ) {
          /* instert jB */
          matadd_insert(bet, B, C, i, &k, &jB, mapB);
        }
        break;
      } else {
        for (; jA < A->ia[i+1]; ) {
          /* insert jA */
          matadd_insert(alp, A, C, i, &k, &jA, mapA);
        }
        break;
      }
    }
    C->ia[i+1] = k;
  }
  C->nnz = k;
  C->ja = evsl_Realloc(C->ja, k, int);
  C->a = evsl_Realloc(C->a, k, double);
  return 0;
}

/** @brief return an identity matrix of dimension n
 * @param[in] n Row/Col size
 * @param[in,out] A Matrix*/
int speye(int n, csrMat *A) {
  int i;
  csr_resize(n, n, n, A);
  for (i=0; i<n; i++) {
    A->ia[i] = A->ja[i] = i;
    A->a[i] = 1.0;
  }
  A->ia[n] = n;
  A->nnz = n;
  return 0;
}

/**
 * @brief Diagonal scaling for A such that A = D^{-1}*A*D^{-1}, i.e.,
 * A(i,j) = A(i,j) / (d(i)*d(j)), d = diag(D)
 * @param[in,out] A Coo Matrix to scale
 * @param[in] d The vector that contains d(i)
 */
void diagScalCoo(cooMat *A, double *d) {
  int i, row, col, nnz = A->nnz;
  /* diagonal scaling for A */
  for (i=0; i<nnz; i++) {
    row = A->ir[i];
    col = A->jc[i];
    A->vv[i] /= d[row] * d[col];
  }
}

/**
 * Diagonal scaling for A such that A = D^{-1}*A*D^{-1}, i.e.,
 * A(i,j) = A(i,j) / (d(i)*d(j)), d = diag(D)
 * @param[in,out] A CSR Matrix to scale
 * @param[in] d The vector that contains d(i)
 */
void diagScalCsr(csrMat *A, double *d) {
  int i, j;
  /* diagonal scaling for A */
  for (i=0; i<A->nrows; i++) {
    for (j=A->ia[i]; j<A->ia[i+1]; j++) {
      A->a[j] /= d[i] * d[A->ja[j]];
    }
  }
}

/**
 * @brief Extract the diagonal entries of csr matrix B
 *
 * @param[in]  B Matrix to extract the diagonal
 * @param[out] d preallocated vector of lengeth B.nrows
 */
void extrDiagCsr(csrMat *B, double *d) {
  int i, j, nrows = B->nrows;
  for (i=0; i<nrows; i++) {
    d[i] = 0.0;
    for (j=B->ia[i]; j<B->ia[i+1]; j++) {
      if (i == B->ja[j]) {
        d[i] = B->a[j];
      }
    }
  }
}

/**
 * @brief Extract the upper triangular part of csr matrix A
 *
 * @param[in]  A Matrix
 * @param[out] U Matrix
 */
void triuCsr(csrMat *A, csrMat *U) {
  EVSL_Int i, j;
  EVSL_Unsigned k = 0, nnzu = 0;

  /* count nnz in U */
  for (i = 0; i < A->nrows; i++) {
    for (j = A->ia[i]; j < A->ia[i+1]; j++) {
      if (i <= A->ja[j]) {
        nnzu++;
      }
    }
  }
  /* allocate memory for U */
  csr_resize(A->nrows, A->ncols, nnzu, U);
  U->ia[0] = 0;
  /* copy entries from A to U */
  for (i = 0; i < A->nrows; i++) {
    for (j = A->ia[i]; j < A->ia[i+1]; j++) {
      if (i <= A->ja[j]) {
        U->ja[k] = A->ja[j];
        U->a[k++] = A->a[j];
      }
    }
    U->ia[i+1] = k;
  }
  CHKERR(k != nnzu);
}


#ifdef EVSL_USING_CUDA_GPU
/**
* @brief matvec for a cusparse CSR matrix, y = A*x.
* void *data points to csrMat,
* compatible form with EVSLMatvec (see struct.h)
*
* @param[in] x Input vector
* @param[out] y Output vector
* @param[in] data CSR matrix
*/
void matvec_cusparse_csr(double *x, double *y, void *data) {
   csrMat *A = (csrMat *) data;
   const double alpha = 1.0;
   const double beta = 0.0;
   cusparseStatus_t cusparseStat;
#if CUSPARSE_VERSION >= 11000
   const cudaDataType        data_type  = CUDA_R_64F;
   const cusparseIndexType_t index_type = CUSPARSE_INDEX_32I;
   const cusparseIndexBase_t index_base = CUSPARSE_INDEX_BASE_ZERO;
   cusparseSpMatDescr_t      matA;

   cusparseStat =
      cusparseCreateCsr(&matA,
                        A->nrows,
                        A->ncols,
                        A->nnz,
                        A->ia,
                        A->ja,
                        A->a,
                        index_type,
                        index_type,
                        index_base,
                        data_type);

   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

   char *dBuffer = A->dBuffer;
   cusparseDnVecDescr_t vecX, vecY;

   cusparseStat = cusparseCreateDnVec(&vecX, A->ncols, x, data_type);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
   cusparseStat = cusparseCreateDnVec(&vecY, A->nrows, y, data_type);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

   if (!dBuffer)
   {
      size_t bufferSize = 0;
      cusparseStat =
         cusparseSpMV_bufferSize(evsldata.cusparseH,
                                 CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 &alpha,
                                 matA,
                                 vecX,
                                 &beta,
                                 vecY,
                                 data_type,
#if CUSPARSE_VERSION >= 11400
                                 CUSPARSE_SPMV_CSR_ALG2,
#else
                                 CUSPARSE_CSRMV_ALG2,
#endif
                                 &bufferSize);

      CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

      dBuffer = evsl_Malloc_device(bufferSize, char);

      A->dBuffer = dBuffer;
   }

   cusparseStat =
      cusparseSpMV(evsldata.cusparseH,
                   CUSPARSE_OPERATION_NON_TRANSPOSE,
                   &alpha,
                   matA,
                   vecX,
                   &beta,
                   vecY,
                   data_type,
#if CUSPARSE_VERSION >= 11400
                   CUSPARSE_SPMV_CSR_ALG2,
#else
                   CUSPARSE_CSRMV_ALG2,
#endif
                   dBuffer);

   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

   cusparseStat = cusparseDestroySpMat(matA);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
   cusparseStat = cusparseDestroyDnVec(vecX);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
   cusparseStat = cusparseDestroyDnVec(vecY);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
#else
   cusparseStat =
      cusparseDcsrmv(evsldata.cusparseH, CUSPARSE_OPERATION_NON_TRANSPOSE,
                     A->nrows, A->ncols, A->nnz,
                     &alpha, A->descr, A->a, A->ia, A->ja, x, &beta, y);
#endif

   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
}

/**
* @brief create CSR on GPU and copy from host
*/
void evsl_create_csr_gpu(csrMat *Acpu, csrMat *Agpu)
{
  int nnzA = Acpu->ia[Acpu->nrows];

  Agpu->owndata = 2;
  Agpu->nrows = Acpu->nrows;
  Agpu->ncols = Acpu->ncols;
  Agpu->nnz   = Acpu->ia[Acpu->nrows];

  cusparseStatus_t cusparseStat = cusparseCreateMatDescr(&Agpu->descr);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  cusparseSetMatIndexBase(Agpu->descr, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatType(Agpu->descr, CUSPARSE_MATRIX_TYPE_GENERAL);

  Agpu->ia = evsl_Malloc_device(Acpu->nrows + 1, int);
  Agpu->ja = evsl_Malloc_device(nnzA, int);
  Agpu->a  = evsl_Malloc_device(nnzA, double);

  cudaMemcpy(Agpu->ia, Acpu->ia, (Acpu->nrows+1)*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(Agpu->ja, Acpu->ja, nnzA*sizeof(int),            cudaMemcpyHostToDevice);
  cudaMemcpy(Agpu->a,  Acpu->a,  nnzA*sizeof(double),         cudaMemcpyHostToDevice);

  Agpu->dBuffer = NULL;
}

/**
* @brief free CSR on GPU
*/
void evsl_free_csr_gpu(csrMat *csr)
{
  if (csr->owndata) {
    evsl_Free_device(csr->ia);
    evsl_Free_device(csr->ja);
    evsl_Free_device(csr->a);
  }

  cusparseStatus_t cusparseStat = cusparseDestroyMatDescr(csr->descr);
  CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

  evsl_Free_device(csr->dBuffer);
}

/*
   if (CUSPARSE_STATUS_SUCCESS != cusparseStat) {
      if (cusparseStat == CUSPARSE_STATUS_NOT_INITIALIZED) printf("err 1\n");
      if (cusparseStat == CUSPARSE_STATUS_ALLOC_FAILED) printf("err 2\n");
      if (cusparseStat == CUSPARSE_STATUS_INVALID_VALUE) printf("err 3\n");
      if (cusparseStat == CUSPARSE_STATUS_ARCH_MISMATCH) printf("err 4\n");
      if (cusparseStat == CUSPARSE_STATUS_MAPPING_ERROR) printf("err 5\n");
      if (cusparseStat == CUSPARSE_STATUS_EXECUTION_FAILED) printf("err 6\n");
      if (cusparseStat == CUSPARSE_STATUS_INTERNAL_ERROR) printf("err 7\n");
      if (cusparseStat == CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED) printf("err 8\n");
   }
*/

#ifdef EVSL_USING_CUSPARSE_HYB
/**
* @brief matvec for a cusparse HYB matrix, y = A*x.
* void *data points to hybMat,
* compatible form with EVSLMatvec (see struct.h)
*
* @param[in] x Input vector
* @param[out] y Output vector
* @param[in] data HYB matrix
*/
void matvec_cusparse_hyb(double *x, double *y, void *data) {
   hybMat *A = (hybMat *) data;
   const double alpha = 1.0;
   const double beta = 0.0;
   cusparseStatus_t cusparseStat =
      cusparseDhybmv(evsldata.cusparseH, CUSPARSE_OPERATION_NON_TRANSPOSE,
                     &alpha, A->descr, A->hyb, x, &beta, y);

   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
}

/**
 * @brief create HYB matrix A on GPU
 *
 * @param[in] A The CSR matrix to set [on host].
 * A cusparse HYB matrix will be generated
 * */
int evsl_create_hybMat(csrMat *A, hybMat *Ahyb) {

  Ahyb->nrows = A->nrows;
  Ahyb->ncols = A->ncols;

  cusparseStatus_t cusparseStat = cusparseCreateMatDescr(&Ahyb->descr);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  cusparseSetMatIndexBase(Ahyb->descr, CUSPARSE_INDEX_BASE_ZERO);
  cusparseSetMatType(Ahyb->descr, CUSPARSE_MATRIX_TYPE_GENERAL);

  cusparseStat = cusparseCreateHybMat(&Ahyb->hyb);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  /* CSR on GPU */
  csrMat Agpu;
  evsl_create_csr_gpu(A, &Agpu);

  /* convert to hyb */
  cusparseStat = cusparseDcsr2hyb(evsldata.cusparseH, A->nrows, A->ncols,
                                  Ahyb->descr, Agpu.a, Agpu.ia, Agpu.ja,
                                  Ahyb->hyb, -1, CUSPARSE_HYB_PARTITION_AUTO);
  CHKERR(cusparseStat != CUSPARSE_STATUS_SUCCESS);

  evsl_free_csr_gpu(&Agpu);

  return 0;
}

/**
* @brief free HYB on GPU
*/
void evsl_free_hybMat(hybMat *hyb) {
   cusparseStatus_t cusparseStat = cusparseDestroyMatDescr(hyb->descr);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);

   cusparseStat = cusparseDestroyHybMat(hyb->hyb);
   CHKERR(CUSPARSE_STATUS_SUCCESS != cusparseStat);
}
#endif
#endif

