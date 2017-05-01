#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mmio.h"
#include "io.h"
#include "evsl.h"

#define ERR_IO  10

char *mybasename (const char *name) {
  const char *base = name;
  while (*name) {
    if (*name++ == '/') {
      base = name;
    }
  }
  return (char *) base;
}

/*-----------------------------------------------*/
int get_matrix_info( FILE *fmat, io_t *pio ){
//  char path[MAX_LINE],  MatNam[MaxNamLen], Fmt[4], ca[2], cb[2], cn_intv[2];  
  char path[MAX_LINE],  MatNam[MaxNamLen], Fmt[4], ca[100], cb[100], cn_intv[100];
  int count, n_intv;
  double a, b;
  /*-------------------- READ LINE */
  if (6 != fscanf(fmat,"%s %s %s %s %s %s\n",path,MatNam,Fmt, ca, cb, cn_intv)) {
    printf("warning: fscanf may not be successfully done\n");
  }
  /*-------------------- file pathname */
  count = strlen(path);
  memset(pio->Fname,0,MAX_LINE*sizeof(char));
  memcpy(pio->Fname,path,count*sizeof(char));
  /*-------------------- file short name */
  count = strlen(MatNam);
  memset(pio->MatNam,0,MaxNamLen*sizeof(char));
  memcpy(pio->MatNam,MatNam,count*sizeof(char));
  /*-------------------- matrix format */
  if (strcmp(Fmt,"HB")==0) 
    pio->Fmt = 1;
  else 
    if (strcmp(Fmt,"MM0")==0) 
      pio->Fmt = MM0;
    else 
      if (strcmp(Fmt,"MM1")==0) 
        pio->Fmt = MM1;
      else 
	/*-------------------- UNKNOWN_FORMAT */
	return(ERR_IO+2);
  /*-------------------- interval information */
  a = atof(ca);
  b = atof(cb);
  pio->a = a;
  pio->b = b;
  n_intv = atoi(cn_intv);
  pio->n_intv = n_intv;
  /* debug  printf(" Echo: %s %s %s \n", pio->Fname, pio->MatNam, Fmt); */
  return(0);
}

/*---------------------------------------------*
 *             READ COO Matrix Market          *
 *---------------------------------------------*/

int read_coo_MM(const char *matfile, int idxin, int idxout, cooMat *Acoo) {
  MM_typecode matcode;
  FILE *p = fopen(matfile,"r");
  int i;
  if (p == NULL) {
    printf("Unable to open mat file %s\n", matfile);
    exit(-1);
  }
  /*----------- READ MM banner */
  if (mm_read_banner(p, &matcode) != 0){
    printf("Could not process Matrix Market banner.\n");
    return 1;
  }
  if (!mm_is_valid(matcode)){
    printf("Invalid Matrix Market file.\n");
    return 1;
  }
  if ( !( (mm_is_real(matcode) || mm_is_integer(matcode)) && mm_is_coordinate(matcode) 
        && mm_is_sparse(matcode) ) ) {
    printf("Only sparse real-valued/integer coordinate \
        matrices are supported\n");
    return 1;
  }
  int nrow, ncol, nnz, nnz2, k, j;
  char line[MAX_LINE];
  /*------------- Read size */
  if (mm_read_mtx_crd_size(p, &nrow, &ncol, &nnz) !=0) {
    printf("MM read size error !\n");
    return 1;
  }
  if (nrow != ncol) {
    fprintf(stdout,"This is not a square matrix!\n");
    return 1;
  }
  /*--------------------------------------
   * symmetric case : only L part stored,
   * so nnz2 := 2*nnz - nnz of diag,
   * so nnz2 <= 2*nnz 
   *-------------------------------------*/
  if (mm_is_symmetric(matcode)){
    printf(" * * * *  matrix is symmetric * * * * \n");
    nnz2 = 2*nnz;}
  else
    nnz2 = nnz;
  /*-------- Allocate mem for COO */
  int* IR = (int *)malloc(nnz2*sizeof(int));
  int* JC = (int *)malloc(nnz2*sizeof(int));
  double* VAL = (double *)malloc(nnz2*sizeof(double));
  /*-------- read line by line */
  char *p1, *p2;
  for (k=0; k<nnz; k++) {
    if (fgets(line, MAX_LINE, p) == NULL) {return -1;};
    for( p1 = line; ' ' == *p1; p1++ );
    /*----------------- 1st entry - row index */
    for( p2 = p1; ' ' != *p2; p2++ ); 
    *p2 = '\0';
    float tmp1 = atof(p1);
    //coo.ir[k] = atoi(p1);
    IR[k] = (int) tmp1;
    /*-------------- 2nd entry - column index */
    for( p1 = p2+1; ' ' == *p1; p1++ );
    for( p2 = p1; ' ' != *p2; p2++ );
    *p2 = '\0';
    float tmp2 = atof(p1);
    JC[k] = (int) tmp2;
    //coo.jc[k]  = atoi(p1);      
    /*------------- 3rd entry - nonzero entry */
    p1 = p2+1;
    VAL[k] = atof(p1); 
  }
  /*------------------ Symmetric case */
  j = nnz;
  if (mm_is_symmetric(matcode)) {
    for (k=0; k<nnz; k++)
      if (IR[k] != JC[k]) {
        /*------------------ off-diag entry */
        IR[j] = JC[k];
        JC[j] = IR[k];
        VAL[j] = VAL[k];
        j++;
      }
    if (j != nnz2) {
      nnz2 = j;
    }
  }
  int offset = idxout - idxin;
  if (offset) {
    for (i=0; i<nnz2; i++) {
      IR[i] += offset;
      JC[i] += offset;
    }
  }
  //  printf("nrow = %d, ncol = %d, nnz = %d\n", nrow, ncol, j);
  fclose(p);
  Acoo->nrows = nrow;
  Acoo->ncols = ncol;
  Acoo->nnz = nnz2;
  Acoo->ir = IR;
  Acoo->jc = JC;
  Acoo->vv = VAL;
  return 0;
}

