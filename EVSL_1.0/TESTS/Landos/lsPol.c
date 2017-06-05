
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "evsl.h"


#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))

double fsqrt(const double a) { 
    return sqrt(1/a);
}

double rec(const double a) {
    return 1/a;
}

/*
 *-----------------------------------------------------------------------
 * Tests landos.c -- Only the  following variable are the outputs. The
 * noted  values are  the values  when the  randn_double in  landos is
 * replaced with a vector of ones.
 *-----------------------------------------------------------------------
 */
int main() {
    double mdeg = 200;
    double* mu = (double*) malloc(sizeof(double) * mdeg);
    double c;
    double h;
    double *intv = (double *) malloc(sizeof(double) * 7);
    int *deg = (int*) malloc(sizeof(int) * 1);
    
    intv[0] = -2.7395e-13;
    intv[1] = 14;
    intv[4] = 0.54793803;
    intv[5] = 2.5;
    intv[0] = 0.54793803;
    intv[1] = 2.5;
    double tau = 1.0e-04;
    lsPol1(intv, mdeg, fsqrt, tau, mu, c, h, deg);
    int i = 0;
    for(i = 0; i < deg[0]; i++) {
	printf("mu: %f \n", mu[i]);
    }
    lsPol1(intv, mdeg, rec, tau, mu, c, h, deg);
    printf("next\n");
    for(i = 0; i < deg[0]; i++) {
	printf("mu: %f \n", mu[i]);
    }
    free(mu);
    free(deg);
    free(intv);

    return 0;
}
