#include "def.h"
#include "blaslapack.h"
#include "struct.h"
#include "internal_proto.h"

/*-------------------------------------------------*
 * convert csr to csc
 * Assume input csr is 0-based index
 * output csc 0/1 index specified by OUTINDEX      *
 * ------------------------------------------------*/
void csrcsc(int OUTINDEX, int nrow, int ncol, int job,
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

/*-------------------------------------------*
 * Sort each row of a csr by increasing column 
 * order
 * By double transposition
 *-------------------------------------------*/
void sortrow(csrMat *A) {
  /*-------------------------------------------*/
  int nrows = A->nrows;
  int ncols = A->ncols;
  int nnz = A->ia[nrows];
  // work array
  double *b;
  int *jb, *ib;
  Malloc(b, nnz, double);
  Malloc(jb, nnz, int);
  Malloc(ib, ncols+1, int);
  // double transposition
  csrcsc(0, nrows, ncols, 1, A->a, A->ja, A->ia, b, jb, ib);
  csrcsc(0, ncols, nrows, 1, b, jb, ib, A->a, A->ja, A->ia);
  // free
  free(b);
  free(jb);
  free(ib);
}

void csr_resize(int nrow, int ncol, int nnz, csrMat *csr) {
  csr->nrows = nrow;
  csr->ncols = ncol;
  Malloc(csr->ia, nrow+1, int);
  Malloc(csr->ja, nnz, int);
  Malloc(csr->a, nnz, double);
}

void free_csr(csrMat *csr) {
  free(csr->ia);
  free(csr->ja);
  free(csr->a);
}

void free_coo(cooMat *coo) {
  free(coo->ir);
  free(coo->jc);
  free(coo->vv);
}

/*---------------------------------------------------------
 * convert coo to csr
 *---------------------------------------------------------*/
int cooMat_to_csrMat(int cooidx, cooMat *coo, csrMat *csr) {
  int nnz = coo->nnz;
  //printf("@@@@ coo2csr, nnz %d\n", nnz);
  csr_resize(coo->nrows, coo->ncols, nnz, csr);

  int i;
  for (i=0; i<coo->nrows+1; i++) {
    csr->ia[i] = 0;
  }
  for (i=0; i<nnz; i++) {
    int row = coo->ir[i] - cooidx;
    csr->ia[row+1] ++;
  }
  for (i=0; i<coo->nrows; i++) {
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
  for (i=coo->nrows; i>0; i--) {
    csr->ia[i] = csr->ia[i-1];
  }
  csr->ia[0] = 0;

  sortrow(csr);
  return 0;
}


double dcsr1nrm(csrMat *A){
  // computes the 1-norm of A 
  double ta = 0.0, t = 0.0;
  int nrows =A->nrows,  one = 1, i, k, k1, len;
  int *ia = A->ia;
  double *aa = A->a;
  k = ia[0];
  for (i=0; i<nrows;i++){
    k1 = ia[i+1];
    len = k1-k;
    t = DASUM(&len,&aa[k],&one);
    ta = max(ta,t);
    k = k1;
  }
  return (ta);
}

void dcsrgemv(char trans, int nrow, int ncol, double alp, double *a, 
    int *ia, int *ja, double *x, double bet, double *y) {
  int i,j;
  if (trans == 'N') {
    //#pragma omp parallel for schedule(guided)
    for (i=0; i<nrow; i++) {
      double r = 0.0;
      for (j=ia[i]; j<ia[i+1]; j++) {
        r += a[j] * x[ja[j]];
      }
      y[i] = bet*y[i] + alp*r;
    }
  } else {  
    for (i=0; i<ncol; i++) {
      y[i] *= bet;
    }
    for (i=0; i<nrow; i++) {
      double xi = alp * x[i];
      for (j=ia[i]; j<ia[i+1]; j++) {
        y[ja[j]] += a[j] * xi;
      }
    }
  }
}

void dcsrmv(char trans, int nrow, int ncol, double *a, 
    int *ia, int *ja, double *x, double *y) {
  int  len, jj=nrow;
  if (trans == 'N') {  
    //#pragma omp parallel for schedule(guided)
    double r ;
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
    double xi; int jj, len;
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

// y = alp*A*x + bet*y
int matvec_gen(double alp, csrMat *A, double *x, double bet, double *y) {
#ifdef EVSL_USE_MKL
  int nrows = A->nrows;
  int ncols = A->ncols;
  char cN = 'N';
  mkl_dcsrmv(&cN, &nrows, &ncols, &alp, "GXXCXX", A->a, A->ja, A->ia, A->ia+1, x, &bet, y);
#else
  dcsrgemv('N', A->nrows, A->ncols, alp, A->a, A->ia, A->ja, x, bet, y);
#endif

  return 0;
}

// y = A * x
int matvec(csrMat *A, double *x, double *y) {
#ifdef EVSL_USE_MKL
  int nrows = A->nrows;
  char cN = 'N';
  mkl_cspblas_dcsrgemv(&cN, &nrows, A->a, A->ia, A->ja, x, y);
#else
  dcsrmv('N', A->nrows, A->ncols, A->a, A->ia, A->ja, x, y);
#endif

  return 0;
}

int lapgen(int nx, int ny, int nz, cooMat *Acoo) {
  int n = nx * ny * nz;
  Acoo->nrows = n;
  Acoo->ncols = n;

  int nzmax = nz > 1 ? 7*n : 5*n;
  Malloc(Acoo->ir, nzmax, int);
  Malloc(Acoo->jc, nzmax, int);
  Malloc(Acoo->vv, nzmax, double);

  int ii, nnz=0;
  for (ii=0; ii<n; ii++) {
    double v = -1.0;
    int i,j,k,jj;
    k = ii / (nx*ny);
    i = (ii - k*nx*ny) / nx;
    j = ii - k*nx*ny - i*nx;

    if (k > 0) {
      jj = ii - nx * ny;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }
    if (k < nz-1) {
      jj = ii + nx * ny;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }

    if (i > 0) {
      jj = ii - nx;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }
    if (i < ny-1) {
      jj = ii + nx;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }

    if (j > 0) {
      jj = ii - 1;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }
    if (j < nx-1) {
      jj = ii + 1;
      Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = jj;  Acoo->vv[nnz] = v;  nnz++;
    }

    v = nz > 1 ? 6.0 : 4.0;
    Acoo->ir[nnz] = ii;  Acoo->jc[nnz] = ii;  Acoo->vv[nnz] = v;  nnz++;
  }

  Acoo->nnz = nnz;

  return 0;
}

