#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "evsl.h"

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
double invsqrt(const double a) { return 1 / sqrt(a); }

/*
 *-----------------------------------------------------------------------
 * Tests landos.c -- Only the  following variable are the outputs. The
 * noted  values are  the values  when the  randn_double in  landos is
 * replaced with a vector of ones.
 *-----------------------------------------------------------------------
 */
int main() {
  const double mdeg = 200;
  double *intv = (double *)malloc(sizeof(double) * 7);
  int i;

  intv[0] = -2.7395e-13;
  intv[1] = 14;
  intv[4] = 0.54793803;
  intv[5] = 2.5;
  intv[0] = 0.54793803;
  intv[1] = 2.5;
  const double tau = 1.0e-04;
  polparams pol;
  // Test v2
  pol.mu = (double *)malloc(sizeof(double) * mdeg);
  pol.cc = 5;
  lsPol(intv, mdeg, invsqrt, tau, &pol);
  free(pol.mu);
  free(intv);

  return 0;
}
