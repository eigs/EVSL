#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "mmio.h"
#include "io.h"

char *mybasename (const char *name) {
  const char *base = name;
  while (*name) {
    if (*name++ == '/') {
      base = name;
    }
  }
  return (char *) base;
}

// returns number of words in str
unsigned countWords(char *str)
{
  int state = 0;
  unsigned wc = 0;  // word count
  // Scan all characters one by one
  while (*str)
  {
    // If next character is a separator, set the
    // state as OUT
     if (*str == ' ' || *str == '\n' || *str == '\t')
     {
       state = 0;
     }
     // If next character is not a word separator and
     // state is OUT, then set the state as IN and
     // increment word count
     else if (state == 0)
     {
        state = 1;
        ++wc;
     }

     // Move to next character
     ++str;
  }

  return wc;
}

/*-----------------------------------------------*/
int get_matrix_info( FILE *fmat, io_t *pio ) {
  char line[MAX_LINE];
  char path1[MAX_LINE], path2[MAX_LINE], MatNam1[MaxNamLen], MatNam2[MaxNamLen], Fmt[4], ca[100], cb[100], cn_intv[100];
  int count, n_intv;
  double a, b;
  /*-------------------- READ LINE */
  if (fgets(line, MAX_LINE, fmat) == NULL) {
    return -1;
  }
  int nfields = countWords(line);
  if (nfields != 6 && nfields != 8) {
    printf("warning: fscanf may not be successfully done\n");
    return -2;
  }
  if (nfields == 6) {
    sscanf(line,"%s %s %s %s %s %s", path1, MatNam1, Fmt, ca, cb, cn_intv);
  }
  if (nfields == 8) {
    sscanf(line,"%s %s %s %s %s %s %s %s", path1, path2, MatNam1, MatNam2, Fmt, ca, cb, cn_intv);
  }
  /*-------------------- file pathname for stiffness and mass matrices*/
  count = strlen(path1);
  memset(pio->Fname1,0,MAX_LINE*sizeof(char));
  memcpy(pio->Fname1,path1,count*sizeof(char));
  if (nfields == 8) {
    count = strlen(path2);
    memset(pio->Fname2,0,MAX_LINE*sizeof(char));
    memcpy(pio->Fname2,path2,count*sizeof(char));
  }
  /*-------------------- file short names for stiffness and mass matrices */
  count = strlen(MatNam1);
  memset(pio->MatNam1,0,MaxNamLen*sizeof(char));
  memcpy(pio->MatNam1,MatNam1,count*sizeof(char));
  if (nfields == 8) {
    count = strlen(MatNam2);
    memset(pio->MatNam2,0,MaxNamLen*sizeof(char));
    memcpy(pio->MatNam2,MatNam2,count*sizeof(char));
  }
  /*-------------------- matrix format */
  if (strcmp(Fmt,"HB")==0)
    pio->Fmt = HB;
  else if (strcmp(Fmt,"MM0")==0)
    pio->Fmt = MM0;
  else if (strcmp(Fmt,"MM1")==0)
    pio->Fmt = MM1;
  else {
     /*-------------------- UNKNOWN_FORMAT */
     return(-3);
  }
  /*-------------------- interval information */
  a = atof(ca);
  b = atof(cb);
  pio->a = a;
  pio->b = b;
  n_intv = atoi(cn_intv);
  pio->n_intv = n_intv;

  /* printf(" Echo: %s %s %d\n", pio->Fname1, pio->MatNam1, pio->Fmt); */

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
    /* printf(" * * * *  matrix is symmetric * * * * \n"); */
    nnz2 = 2*nnz;
  } else {
    nnz2 = nnz;
  }
  /*-------- Allocate mem for COO */
  int* IR = evsl_Malloc(nnz2, int);
  int* JC = evsl_Malloc(nnz2, int);
  double* VAL = evsl_Malloc(nnz2, double);
  /*-------- read line by line */
  char *p1, *p2;
  for (k=0; k<nnz; k++) {
    if (fgets(line, MAX_LINE, p) == NULL) { return -1; }
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

// parse command-line input parameters
int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv) {
  int *outint;
  double *outdouble;
  char *outchar;
  int i;
  for (i=0; i<argc; i++) {
    if (argv[i][0] != '-') {
      continue;
    }
    if (!strcmp(argname, argv[i]+1)) {
      if (type == NA) {
        return 1;
      } else {
        if (i+1 >= argc /*|| argv[i+1][0] == '-'*/) {
          return 0;
        }
        switch (type) {
          case INT:
            outint = (int *) val;
            *outint = atoi(argv[i+1]);
            return 1;
          case DOUBLE:
            outdouble = (double *) val;
            *outdouble = atof(argv[i+1]);
            return 1;
          case STR:
            outchar = (char *) val;
            sprintf(outchar, "%s", argv[i+1]);
            return 1;
          default:
            printf("unknown arg type\n");
        }
      }
    }
  }
  return 0;
}

