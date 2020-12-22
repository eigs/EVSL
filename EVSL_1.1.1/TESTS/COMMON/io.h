#ifndef IO_H
#define IO_H

#include "evsl.h"

#define MAX_MAT   100
#define MAX_LINE  1024
#define MaxNamLen 64
#define HB   1
#define MM0  2
#define MM1  3
#define UNK  4

typedef struct _io_t {
  FILE *fout;                 /* output file handle            */
  char outfile[MAX_LINE];     /* output filename               */
  char Fname1[MAX_LINE];      /* full matrix path name         */
  char MatNam1[MaxNamLen];    /* short name                    */
  char Fname2[MAX_LINE];      /* full matrix path name         */
  char MatNam2[MaxNamLen];    /* short name                    */
  char type[4];               /* type for HB matrices          */
  int Fmt;                    /* matrix format type            */
  int ndim;                   /* matrix size                   */
  int nnz;                    /* number of nonzero             */
  int n_intv;                 /* number of slices              */
  double a;                   /* [a, b]  interval of interest  */
  double b;
} io_t;

/* types of user command-line input */
typedef enum {
  INT,
  DOUBLE,
  STR,
  NA
} ARG_TYPE;

#ifdef __cplusplus
extern "C" {
#endif

int get_matrix_info( FILE *fmat, io_t *pio );
int read_coo_MM(const char *matfile, int idxin, int idxout, cooMat *Acoo);
int findarg(const char *argname, ARG_TYPE type, void *val, int argc, char **argv);

#ifdef __cplusplus
}
#endif

#endif

