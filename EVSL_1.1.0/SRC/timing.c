#include <time.h>
#include <math.h>
#ifdef USE_MKL
#include <mkl.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file timing.c
 * @brief Timer
 */

/**
 * @brief cheblan timer
 */
double evsl_timer() {
#ifdef USE_MKL
  return dsecnd();
#else
  /* POSIX C 1993 timer, requires -librt */
  struct timespec t ;
  clock_gettime (CLOCK_MONOTONIC /*CLOCK_PROCESS_CPUTIME_ID*/, &t) ;
  return ((double) (t.tv_sec) + 1e-9 * (double) (t.tv_nsec));
#endif
}

/** 
  * @brief Uses the timer to generate a seed to be used
  * for srand.. 
  */
int time_seeder() { 
  double t1,t2;
  int iseed, zero=0;
  t1   = evsl_timer();
  t1  = 1.e+09*frexp(t1, &zero);
  t1  = modf(t1, &t2);
  iseed = (int)(1.e+05*t1);
  return (iseed);
}

#ifdef __cplusplus
}
#endif
