#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "time_it.h"

extern void mystery(uint64_t *);

#define DEFAULT_TEST           0
#define DEFAULT_N        1000000
#define DEFAULT_NTHREADS       8
#define DEFAULT_NTRIALS      100
#define DEFAULT_DENSE      FALSE

#define TRUE 1
#define FALSE 0

#define New(Type) ((Type *)(malloc(sizeof(Type))))
#define NewArray(N, Type) ((Type *)(malloc((N)*sizeof(Type))))

/***********************************************************************
 *                                                                     *
 * Stuff for timing measurements                                       *
 *                                                                     *
 ***********************************************************************/

// struct time_it_raw: used by time_it_run to record execution times
struct time_it_raw {
  int ntrials;
  struct timeval *t;
};

struct time_it_raw *time_it_create(int ntrials) {
  struct time_it_raw *tr = New(struct time_it_raw);
  tr->ntrials = ntrials;
  tr->t = NewArray(ntrials+1, struct timeval);
  return(tr);
}

void time_it_run(struct time_it_raw *tr, void (*fn)(void *), void *arg) {
  int ntrials = tr->ntrials;
  for(int j = 0; j <= ntrials; j++) {
    fn(arg);
    gettimeofday(&(tr->t[j]), NULL);
  }
}

void time_it_free(struct time_it_raw *tr) {
  free((void *)(tr->t));
  free((void *)(tr));
}

void time_it_get_stats(struct time_it_raw *tr, struct time_it_stats *s) {
  double sum_t  = 0.0;
  double sum_t2 = 0.0;
  uint n = tr->ntrials;
  for(int i = 0; i < n; i++) {
    double t =   1.0e-6*(tr->t[i+1].tv_usec - tr->t[i].tv_usec)
               + (tr->t[i+1].tv_sec - tr->t[i].tv_sec);
    sum_t  += t;
    sum_t2 += t*t;
  }
  s->mean = sum_t/n;
  s->std  = sqrt((sum_t2 - (sum_t*s->mean))/(n-1));
}
