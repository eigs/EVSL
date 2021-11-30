#include "internal_header.h"

/**
 * @file timing.c
 * @brief Timer
 */

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

#if defined(__APPLE__) && defined(__MACH__)

/* for Mac */
#include <mach/mach.h>
#include <mach/mach_time.h>

/**
 * @brief evsl timer for mac
 * @return Current time (in nanoseconds)
 */
double evsl_timer() {
  double t;
  uint64_t absNano;
  static mach_timebase_info_data_t    sTimebaseInfo;
  absNano = mach_absolute_time();
  if ( sTimebaseInfo.denom == 0 ) {
    (void) mach_timebase_info(&sTimebaseInfo);
  }
  t = (double) absNano * 1.e-09*
    (sTimebaseInfo.numer / sTimebaseInfo.denom);
  return t;
}

#elif defined(__linux__)

/* see https://stackoverflow.com/questions/40515557/compilation-error-on-clock-gettime-and-clock-monotonic/40515669 */
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif

/* for Linux */
#include <sys/time.h>
#include <time.h>
/**
 * @brief cheblan timer
 * @return Returns current time in nanoseconds
 */
double evsl_timer() {
  /* POSIX C 1993 timer, requires -librt */
  struct timespec t;
#ifdef EVSL_USING_CUDA_GPU
  cudaDeviceSynchronize();
#endif
  clock_gettime (CLOCK_MONOTONIC /*CLOCK_PROCESS_CPUTIME_ID*/, &t) ;
  return ((double) (t.tv_sec) + 1e-9 * (double) (t.tv_nsec));
  /*
  struct timeval tim;
  gettimeofday(&tim, NULL);
  double t = tim.tv_sec + tim.tv_usec/1e6;
  return(t);
  */
}

#endif

