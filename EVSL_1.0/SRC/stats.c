#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "def.h"
#include "struct.h"
#include "internal_proto.h"

void StatsPrint(FILE *fstats) {
  evslStat *stats = &evslstat;
  double t_iter = stats->t_iter;
  double t_mvA = stats->t_mvA;
  double t_mvB = stats->t_mvB;
  double t_svB = stats->t_svB;
  double t_svASigB = stats->t_svASigB;
  size_t n_mvA = stats->n_mvA;
  size_t n_mvB = stats->n_mvB;
  size_t n_svB = stats->n_svB;
  size_t n_svASigB = stats->n_svASigB;
  /* memory */
  //size_t alloced       = stats->alloced;
  //size_t alloced_total = stats->alloced_total;
  //size_t alloced_max   = stats->alloced_max;
  /* time */
  fprintf(fstats, " Timing (sec):\n");
  fprintf(fstats, "   Iterative solver         :  %f\n",  t_iter);
  if (n_mvA)     { fprintf(fstats, "   Matvec matrix A          :  %f (%8ld, avg %f)\n",  t_mvA, n_mvA, t_mvA / n_mvA); }
  if (n_mvB)     { fprintf(fstats, "   Matvec matrix B          :  %f (%8ld, avg %f)\n",  t_mvB, n_mvB, t_mvB / n_mvB); }
  if (n_svB)     { fprintf(fstats, "   Solve with B             :  %f (%8ld, avg %f)\n",  t_svB, n_svB, t_svB / n_svB); }
  if (n_svASigB) { fprintf(fstats, "   Solve with A-SIGMA*B     :  %f (%8ld, avg %f)\n",  t_svASigB, n_svASigB, t_svASigB / n_svASigB); }
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

void StatsReset() {
  memset(&evslstat, 0, sizeof(evslStat));
}
