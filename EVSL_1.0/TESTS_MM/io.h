#ifndef IO_H
#define IO_H

#define MAX_MAT	   100
#define MAX_LINE 1024
#define MaxNamLen 64
#define HB   1
#define MM0  2
#define MM1  3
#define UNK  4

#include "evsl.h"

typedef struct _io_t {
    FILE *fout;                 /* output file handle              */
    char outfile[MAX_LINE];     /* output filename                 */
    char Fname[MAX_LINE];       /* full matrix path name           */
    char MatNam[MaxNamLen];     /* short name                      */
    char type[4];               /* type for HB matrices            */
    int Fmt;                    /* matrix format type              */
    int ndim;                   /* matrix size                     */
    int nnz;                    /* number of nonzero               */
    int n_intv;                 /* number of slices                */
    double a;                   /* [a, b]  interval of interest    */
    double b;
} io_t;

int get_matrix_info( FILE *fmat, io_t *pio );

int read_coo_MM(const char *matfile, int idxin, int idxout, cooMat *Acoo);

#endif

