#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"

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
            break;
          case DOUBLE:
            outdouble = (double *) val;
            *outdouble = atof(argv[i+1]);
            return 1;
            break;
          case STR:
            outchar = (char *) val;
            sprintf(outchar, "%s", argv[i+1]);
            return 1;
            break;
          default:
            printf("unknown arg type\n");
        }
      }
    }
  }
  return 0;
}

