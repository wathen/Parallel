/* perc.cu: percolation simulation using GPUs
 *   Usage: perc [width [height [markingProbability]]]
 *     Where
 *       width is the width of the graph.
 *       height is the height of the graph -- number of rounds to simulate.
 *       markingProbability is the probability of marking a vertex if one
 *         of its predecessors is marked.  If neither is marked, then this
 *         vertex won't get marked either.
 *     Results:
 *       Print the mean and standard deviation for the elapsed time for
 *         the simulation.  Computed based on 10 trials.
 */
#include <stdio.h>
#include <cuda.h>
#include <curand_kernel.h>
#include "time_it.h"

// CUDA_CALL and setup_kernel are from 
//   http://docs.nvidia.com/cuda/curand/device-api-overview.html#device-api-example
#define CUDA_CALL(x) do { if((x) != cudaSuccess) { \
  printf("Error at %s:%d\n",__FILE__,__LINE__); \
  return EXIT_FAILURE;}} while(0)

__global__ void setup_kernel(uint n, curandState *state)
{
  uint myId = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence 
     number, no offset */
  if(myId < n)
    curand_init(1234, myId, 0, &state[myId]);
}


// perc: a kernel that computes m rounds of percolation
__global__ void perc(uint n, uint m, float p, curandState *randState, uint *v) {
  uint myId = threadIdx.x; // my index -- assume there's only one block
  uint id1 = ((myId == 0) ? n : myId) - 1; // index of my neighbour to the left
  curandState *myRandState = &(randState[myId]); // my random number generator

  v[myId] = 1; // all initial vertices are marked (top row)
  for(int j = 1; j <= m; j++) { // for each round of percolation
    __syncthreads(); // make sure it's safe to read from v
    // should we be marked in this round?
    uint next = (v[myId] || v[id1]) && (curand_uniform(myRandState) <= p);
    __syncthreads(); // make sure it's safe to write to v.
    v[myId] = next; // update the global state so our neigbour can see it.
  }
}

// arguments to the do_perc function (below)
struct perc_arg {
  uint n,  // width of the graph
       m;  // height of the graph
  float p; // probabitity that a vertex is marked if one of its predecessors is marked
  curandState *dev_randState; // an array of n random number generators
  uint *dev_v; // write the final state here.
};

// do_perc: launch the perc kernel
//   wrapped up as a CPU function so timing_run can call us.
void do_perc(void *void_arg) {
  perc_arg *arg = (perc_arg *)(void_arg);
  perc<<<1,arg->n>>>(arg->n, arg->m, arg->p, arg->dev_randState, arg->dev_v);
  cudaDeviceSynchronize();
}

main(int argc, char **argv) {
  int ndev;
  cudaDeviceProp prop;
  uint *v;
  struct time_it_raw *tr = time_it_create(10);
  struct time_it_stats stats;
  perc_arg parg; // parameters for do_perc

  // read our command line parameters.
  parg.n = (argc > 1) ? atoi(argv[1]) : 32;  // graph width
  parg.m = (argc > 2) ? atoi(argv[2]) : 256; // graph height
  parg.p = (argc > 3) ? atof(argv[3]) : 0.5; // marking probabilty

  // make sure we have a GPU
  CUDA_CALL(cudaGetDeviceCount(&ndev));
  if(ndev < 1) {
    fprintf(stderr, "No CUDA device found\n");
    exit(-1);
  }
  CUDA_CALL(cudaGetDeviceProperties(&prop, 0));

  // Use one thread for each column in the graph.
  //   All the threads are in one block so we can use __syncthreads().
  //   If the graph is too wide for this approach, report an error and exit.
  if(parg.n > prop.maxThreadsPerBlock) {
    fprintf(stderr, "perc: array size too big, max = %d\n",
      prop.maxThreadsPerBlock);
    exit(-1);
  }

  // allocate an array for the result on the CPU
  int vsz = parg.n*sizeof(uint);
  v = (uint *)malloc(vsz);

  // allocate the result array and pseudo-random number generator states on the GPU
  CUDA_CALL(cudaMalloc((void **)(&parg.dev_v), vsz));
  CUDA_CALL(cudaMalloc((void **)(&parg.dev_randState), parg.n*sizeof(curandState)));

  // initialize the random number generators. 
  setup_kernel<<<1,parg.n>>>(parg.n, parg.dev_randState);

  // make the timing measurements.
  time_it_run(tr, do_perc, (void *)(&parg));
  // fetch the final state from the GPU
  cudaMemcpy(v, parg.dev_v, vsz, cudaMemcpyDeviceToHost);

  // count the number of marked vertices in the final generation.
  //   If the number is greater than zero, than the last row can be reached
  //   from the first row with a path where every vertex is marked.
  int sum = 0;
  for(uint i = 0; i < parg.n; i++) 
    sum += v[i];
  printf("%d vertices reachable after %d generations\n", sum, parg.m);
  time_it_get_stats(tr, &stats);

  // print the mean and standard deviation of the elapsed times for
  //   running the simulations.
  printf("perc(%u, %u, %5.3f): mean(T) = %10.3le, stddev(T) = %10.3le\n",
            parg.n, parg.m, parg.p, stats.mean, stats.std);

  // clean up
  CUDA_CALL(cudaFree(parg.dev_randState));
  CUDA_CALL(cudaFree(parg.dev_v));
  time_it_free(tr);
  free(v);
  exit(0);
}
