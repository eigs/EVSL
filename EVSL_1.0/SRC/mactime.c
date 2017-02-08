#include <mach/mach.h>  
#include <mach/mach_time.h>  
#include <unistd.h>
#include <stdio.h>  
#include <math.h>

/**
 * @brief evsl timer for mac
 */
double cheblan_timer() {
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

/**
 * @brief Uses the timer to generate a seed to be used for srand.
 */
int time_seeder() {
  double t1,t2;
  int iseed, zero=0;
  t1   = cheblan_timer();
  t1  = 1.e+09*frexp(t1, &zero);
  t1  = modf(t1, &t2);
  iseed = (int)(1.e+05*t1);
  return (iseed);
}

      
