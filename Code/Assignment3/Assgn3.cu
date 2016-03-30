#include <stdio.h>
#include <math.h>
#include "time_it.h"
#include <iostream>
using namespace std;

#define N (512*512)
#define THREADS_PER_BLOCK 128

__global__ void logistic_cuda(unsigned int n, unsigned int m, float a, float *x, float *z) {
  unsigned int myId = blockDim.x*blockIdx.x + threadIdx.x;
  if(myId < n){
    for (int i = 1; i < m; ++i) {
      z[myId] = a*x[myId]*(1.0f - x[myId]);
      x[myId] = z[myId];
    }
  }
}

void logistic_ref(unsigned int n, unsigned int m, float a, float *x, float *z) {
  for (int j = 1; j < m; ++j) {
    for(int i = 1; i < n; ++i){
     z[i] = a*x[i]*(1.0f - x[i]);
     x[i] = z[i];
    }
  }

}
struct kernel_arg {
    float *x;
    unsigned int n;
};

__device__ void reduce_sum_dev(unsigned int n, float *x) {
  unsigned int myId = threadIdx.x;
  for(unsigned int m = n >> 1; m > 0; m = n >> 1) {
    n -= m;
    __syncthreads();
    if(myId < m)
      x[myId] += x[myId+n];
  }
}

__global__ void reduce_sum(unsigned int n, float *x) {
  reduce_sum_dev(n, x);
}

float reduce_sum_ref(unsigned int n, float *x) {
  float sum = 0.0;
  for(int i = 0; i < n; i++)
    sum += x[i];
  return(sum);
}

/******************************
*  dotprod: just like it sounds
*    Simple version: we only handle one block of threads
*******************************/

__global__ void dotprod(unsigned int n, float *x, float *z) {
  unsigned int blockBase = blockDim.x * blockIdx.x;
  unsigned int myId = blockBase + threadIdx.x;
  unsigned int m = min(blockDim.x, n - blockBase);

  if(myId < n)
    x[myId] *= x[myId];
  reduce_sum_dev(m, &(x[blockBase]));
  if((myId < n) && (threadIdx.x == 0))
    z[blockIdx.x] = x[myId];
}

float norm_ref(float *x, unsigned int n) {
  float sum = 0.0;
  for(int i = 0; i < n; i++)
    sum += x[i] * x[i];
  return(sqrt(sum));
}

/*****************************************************
*  print_vec: print the first few elements of a vector
******************************************************/

void print_vec(float *x, unsigned int n, const char *fmt, const char *who) {
  printf("%s = ", who);
  for(int i = 0; i < n; i++) {
    if(i > 0) printf(", ");
    printf(fmt, x[i]);
  }
  if(n > 10) printf(", ...");
  printf("\n");
}

void logistic(float *x, unsigned int a, unsigned int n, unsigned int m, float *z) {

  float *dev_x, *dev_z;
  int size = n*sizeof(float);
  cudaMalloc((void**)(&dev_x), size);
  cudaMalloc((void**)(&dev_z), size);
  cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);

  logistic_cuda<<<N/THREADS_PER_BLOCK , THREADS_PER_BLOCK>>>(n, m, a, dev_x, dev_z);
  cudaMemcpy(z, dev_z, sizeof(size), cudaMemcpyDeviceToHost);
}

float norm(float * x, unsigned int n) {

  float *z;
  float *dev_x, *dev_z;
  int size = n*sizeof(float);
  z = (float *) malloc(size);

  cudaMalloc((void**)(&dev_x), size);
  cudaMalloc((void**)(&dev_z), size);
  cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);

  dotprod<<<N/THREADS_PER_BLOCK , THREADS_PER_BLOCK>>>(n, dev_x, dev_z);
  reduce_sum<<<1,N/THREADS_PER_BLOCK>>>(N/THREADS_PER_BLOCK, dev_z);
  cudaMemcpy(z, dev_z, sizeof(float), cudaMemcpyDeviceToHost);
  return sqrt(z[0]);
}

void do_timing(void *void_arg) {
  struct kernel_arg *argk = (struct kernel_arg *)(void_arg);
  norm(argk->x, argk->n);
  cudaDeviceSynchronize();
}

int main(int argc, char **argv) {
  unsigned int n = N;
  unsigned int m = 10;
  float *x, *z, *z_ref;
  float a;
  cudaDeviceProp prop;
  struct kernel_arg argk;
  struct time_it_raw *tr = time_it_create(10);
  struct time_it_stats stats;

  int size = n*sizeof(float);
  x  = (float *)malloc(size);
  z  = (float *)malloc(size);
  z_ref = (float *)malloc(size);

  for(int i = 0; i < n; i++) {
    x[i] = i;
  }

  printf("The GPU is a %s\n", prop.name);
  printf("Cuda capability %d.%d.\n", prop.major, prop.minor);
  float p_norm = norm(x, n);
  z_ref[0] = norm_ref(x, n);


  printf("Parallel = %f, Sequential = %f\n\n", p_norm, z_ref[0]);
  argk.n = N;
  argk.x = x;
  time_it_run(tr, do_timing, (void *)(&argk));
  time_it_get_stats(tr, &stats);
  printf("mean(T) = %10.3e, std(T) = %10.3e\n", stats.mean, stats.std);

  a = 3.0;
  float *L;
  L  = (float*)malloc(size);
  logistic(x, a, n, m, L);
  logistic_ref(n, m, a, x, z);

  print_vec(z, min(10, N), "%5.3f", "z");
  print_vec(L, min(10, N), "%5.3f", "z");

  free(x);
  free(z_ref);
  exit(0);
}
