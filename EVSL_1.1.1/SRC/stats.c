#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "internal_header.h"

/**
 * @file stats.c
 * @brief Used to track various statistics (time taken by various
 * operations).
 */
/**
 * @brief Prints out stats
 * @param[in] fstats FILE to print to
 */

void StatsPrint(FILE *fstats) {
  evslStat *stats = &evslstat;
  double t_iter = stats->t_iter;
  double t_setBsv = stats->t_setBsv;
  double t_setASigBsv = stats->t_setASigBsv;
  double t_mvA = stats->t_mvA;
  double t_mvB = stats->t_mvB;
  double t_svB = stats->t_svB;
  double t_svLT = stats->t_svLT;
  double t_svASigB = stats->t_svASigB;
  double t_reorth = stats->t_reorth;
  double t_siorth = stats->t_siorth;
  double t_eig = stats->t_eig;
  double t_blas = stats->t_blas;
  double t_ritz = stats->t_ritz;
  double t_polAv = stats->t_polAv;
  double t_ratAv = stats->t_ratAv;
  size_t n_mvA = stats->n_mvA;
  size_t n_mvB = stats->n_mvB;
  size_t n_svB = stats->n_svB;
  size_t n_svLT = stats->n_svLT;
  size_t n_svASigB = stats->n_svASigB;
  size_t n_polAv = stats->n_polAv;
  size_t n_ratAv = stats->n_ratAv;
  /* memory */
  //size_t alloced       = stats->alloced;
  //size_t alloced_total = stats->alloced_total;
  //size_t alloced_max   = stats->alloced_max;
  /* time */
  fprintf(fstats, " Timing (sec):\n");
  if (t_setBsv)     { fprintf(fstats, "   Setup Solver for B       :  %f\n",  t_setBsv); }
  if (t_setASigBsv) { fprintf(fstats, "   Setup Solver for A-SIG*B :  %f\n",  t_setASigBsv); }
  if (t_iter)       { fprintf(fstats, "   Iteration time (tot)     :  %f\n",  t_iter); }

  fprintf(fstats, "   - - - - - - - - - - - - - - - - -\n");

  if (n_polAv)   { fprintf(fstats, "   Pol(A)*v                 :  %f (%8ld, avg %f)\n",  t_polAv, n_polAv, t_polAv / n_polAv); }
  if (n_ratAv)   { fprintf(fstats, "   Rat(A)*v                 :  %f (%8ld, avg %f)\n",  t_ratAv, n_ratAv, t_ratAv / n_ratAv); }
  if (n_mvA)     { fprintf(fstats, "   Matvec matrix A          :  %f (%8ld, avg %f)\n",  t_mvA, n_mvA, t_mvA / n_mvA); }
  if (n_mvB)     { fprintf(fstats, "   Matvec matrix B          :  %f (%8ld, avg %f)\n",  t_mvB, n_mvB, t_mvB / n_mvB); }
  if (n_svB)     { fprintf(fstats, "   Solve with B             :  %f (%8ld, avg %f)\n",  t_svB, n_svB, t_svB / n_svB); }
  if (n_svLT)    { fprintf(fstats, "   Solve with LT            :  %f (%8ld, avg %f)\n",  t_svLT, n_svLT, t_svLT / n_svLT); }
  if (n_svASigB) { fprintf(fstats, "   Solve with A-SIGMA*B     :  %f (%8ld, avg %f)\n",  t_svASigB, n_svASigB, t_svASigB / n_svASigB); }
  if (t_reorth)  { fprintf(fstats, "   Reorthogonalization      :  %f\n", t_reorth); }
  if (t_siorth)  { fprintf(fstats, "   SIorthogonalization      :  %f\n", t_siorth); }
  if (t_eig)     { fprintf(fstats, "   LAPACK eig               :  %f\n", t_eig); }
  if (t_blas)    { fprintf(fstats, "   Other BLAS               :  %f\n", t_blas); }
  if (t_ritz)    { fprintf(fstats, "   Compute Ritz vectors     :  %f\n", t_ritz); }
  /* if (t_sth)     { fprintf(fstats, "   Some other thing timed   :  %f\n", t_sth); } */
  /* memory */
  /*
  if (alloced_total > 1e9) {
    fprintf(fstats, " Memory (GB):\n");
    fprintf(fstats, "   Total %.2f,  Peak %.2f \n", alloced_total/1e9, alloced_max/1e9);
  } else if (alloced_total > 1e6) {
    fprintf(fstats, " Memory (MB):\n");
    fprintf(fstats, "   Total %.2f,  Peak %.2f \n", alloced_total/1e6, alloced_max/1e6);
  } else {
    fprintf(fstats, " Memory (KB):\n");
    fprintf(fstats, "   Total %.2f,  Peak %.2f \n", alloced_total/1e3, alloced_max/1e3);
  }
  if (alloced) {
    fprintf(fstats, "warning: unfreed memory %ld\n", alloced);
  }
  */
  fflush(fstats);
}

/**
 * @brief Resets stats
 * */
void StatsReset() {
  memset(&evslstat, 0, sizeof(evslStat));
}
