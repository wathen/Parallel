/* timing.cu -- make some simple timing measurements.
 *   Usage: timing [test [ntrials]]
 *   Where
 *     test specifies what test to run (default 'empty')
 *       empty means time an empty kernel.
 *       sync  means run a kernel to measure the time for __syncthreads()
 *     ntrials is a positive integer (default 100).
 *   The average elapsed time and standard deviation for a kernel launch
 *     are reported.
 *   Note that sync performs 1000 __syncthreads() calls per kernel launch.
 *     This is to minimize the impact of the kernel launch time on the
 *     measurement.
 */
#include <stdio.h>
#include <math.h>
#include "time_it.h"

// DEFAULT_NTRIALS: how many time to run each kernel
//   to calculate average elapsed time and standard deviation.
#define DEFAULT_NTRIALS 100

// NSYNC: number of _syncthreads() calls per thread
//   for the SYNC_THREADS test case
#define NSYNC 1000


/************************
*  empty kernel
************************/

__global__ void empty() {}

void do_empty(void *arg) {
  empty<<<1,1>>>();
  cudaDeviceSynchronize();
}

/************************
*  syncthreds
************************/

// I tried a version with an array in shared memory
//   that each thread accesses so synchronizing would
//   serve a purpose.  It doesn't seem to have a
//   signficant effect on the runtime.

// __shared__ int sh[1024];
__global__ void sync(int m) {
  // int myId = threadIdx.x;
  // int sum = 0;
  for(int i = 0; i < m; i++) {
    __syncthreads();
    // sum += sh[(myId + i) & 1023];
  }
  // sh[myId] = sum;
}

void do_sync(void *arg) {
  sync<<<1,1024>>>(NSYNC);
  cudaDeviceSynchronize();
}

static void usage() {
  fprintf(stderr, "usage: timing [test [ntrials]]\n");
  fprintf(stderr, "  where test is either 'empty' or 'sync'\n");
  fprintf(stderr, "    'empty' indicates launch an empty kernel\n");
  fprintf(stderr,
          "    'sync'  indicates launch a kernel that performs %d %s\n",
          "1000 __syncthreads() calls\n");
  fprintf(stderr, "  and ntrials is a positive integer\n");
  exit(-1);
}

int main(int argc, char **argv) {
  int ndev;
  cudaDeviceProp prop;
  void (*what_fun)(void *);
  const char *what_name;
  int ntrials;

  // make sure we have a CUDA capable GPU on this machine
  cudaGetDeviceCount(&ndev);
  if(ndev < 1) {
    fprintf(stderr, "No CUDA device found\n");
    exit(-1);
  }
  // I'll just assume that we only have one CUDA device
  //   or that we want device 0.
  cudaGetDeviceProperties(&prop, 0);

  // determine what function to run
  if((argc <= 1) || (strcmp(argv[1], "empty") == 0)) {
    what_fun = do_empty;
    what_name = "empty kernel";
  } else if(strcmp(argv[1], "sync") == 0) {
    what_fun = do_sync;
    what_name = "synchronize threads";
  } else {
    fprintf(stderr, "timing: unknown test -- %s\n", argv[1]);
    usage();
  }

  // determine how many times to run the function
  if(argc <= 2)
    ntrials = DEFAULT_NTRIALS;
  else {
    ntrials = atoi(argv[2]);
    if(ntrials <= 0) {
      fprintf(stderr, "timing: bad value for ntrials, %d\n", ntrials);
      usage();
    }
  }

  // our data structures for time_it
  struct time_it_raw *tr = time_it_create(ntrials);
  struct time_it_stats stats;

  // call time_it_run ntrials times and record the times
  time_it_run(tr, what_fun, NULL);

  // summarize the results
  time_it_get_stats(tr, &stats);
  printf("%s: mean(T) = %10.3le, stddev(T) = %10.3le\n",
      what_name, stats.mean, stats.std);

  // clean up and return.
  time_it_free(tr);
  exit(0);
}
