#include <stdio.h>
#include <math.h>
#include "time_it.h"
#include <iostream>
using namespace std;
#define N (16*16)
#define THREADS_PER_BLOCK 8

struct kernel_arg_logistic {
  float *x;
  unsigned int a;
  unsigned int n;
  unsigned int m;
  float *z;
  unsigned int blocksize;
};


__global__ void logistic_cuda(unsigned int n, unsigned int m, float a, float *x, float *z) {
  unsigned int myId = blockDim.x*blockIdx.x + threadIdx.x;
  float X, Z;
  X = x[myId];
  if(myId < n){
    for (int i = 1; i < m; ++i) {
      Z = a*X*(1.0f - X);
      X = Z;
    }
    z[myId] = Z;
  }
}

void logistic_ref(unsigned int n, unsigned int m, float a, float *x, float *z) {
  for (int j = 1; j < m; ++j) {
    for(int i = 0; i < n; ++i){
     z[i] = a*x[i]*(1.0f - x[i]);
     x[i] = z[i];
    }
  }

}

void print_vec(float *x, unsigned int n, const char *fmt, const char *who) {
  printf("%s = ", who);
  for(int i = 0; i < n; i++) {
    if(i > 0) printf(", ");
    printf(fmt, x[i]);
  }
  if(n > 10) printf(", ...");
  printf("\n");
}

void logistic(float *x, unsigned int a, unsigned int n, unsigned int m, float *z, unsigned int blocksize) {

  float *dev_x, *dev_z;
  int size = n*sizeof(float);
  cudaMalloc((void**)(&dev_x), size);
  cudaMalloc((void**)(&dev_z), size);
  cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);
  unsigned int nblks = ceil(((float)(n))/((float)(blocksize)));
  logistic_cuda<<<nblks , blocksize>>>(n, m, a, dev_x, dev_z);
  cudaMemcpy(z, dev_z, size, cudaMemcpyDeviceToHost);
}

void do_timing_logistic(void *void_arg) {
  struct kernel_arg_logistic *argk = (struct kernel_arg_logistic *)(void_arg);
  logistic(argk->x, argk->a, argk->n, argk->m, argk->z, argk->blocksize);
  cudaDeviceSynchronize();
}

int main(int argc, char *argv[] ) {
  unsigned int n = atoi(argv[1]);
  unsigned int m = atoi(argv[4]);
  unsigned int blocksize = atoi(argv[2]);
  float *x, *z, *z_ref;
  float a = atoi(argv[3]);
  cudaDeviceProp prop;
  struct kernel_arg_logistic argk;
  struct time_it_raw *tr = time_it_create(10);
  struct time_it_stats stats;

  int size = n*sizeof(float);
  x  = (float *)malloc(size);
  z  = (float *)malloc(size);
  z_ref = (float *)malloc(size);
  printf("%d\n",n );
  for(int i = 0; i < n; i++) {
    x[i] =  (float)rand() / (float)RAND_MAX;

  }

  int ndev;
  cudaGetDeviceCount(&ndev);
  cudaGetDeviceProperties(&prop, 0);
  printf("The GPU is a %s\n", prop.name);
  printf("Cuda capability %d.%d.\n", prop.major, prop.minor);


  float *L;
  L  = (float*)malloc(size);
  logistic(x, a, n, m, L, blocksize);
  logistic_ref(n, m, a, x, z);

  argk.x = x;
  argk.a = a;
  argk.n = n;
  argk.m = m;
  argk.z = z;
  argk.blocksize = blocksize;
  time_it_run(tr, do_timing_logistic, (void *)(&argk));
  time_it_get_stats(tr, &stats);
  printf("Logistic: mean(T) = %10.3e, std(T) = %10.3e\n", stats.mean, stats.std);
  printf("\n\nLogistic calculations:   # operations = %10.4i   mean time = %1.4e  time per op = %1.4e, Gflops = %5.3f\n\n", 3*n*m, stats.mean, stats.mean/(3*n*m), (3*(float)n*(float)m)/(stats.mean*pow(10,9)));

  // for(int i = 0; i < n; i++){
  //     printf("z = %5.5f,  L = %5.5f \n", z[i], L[i]);
  // }
  free(x);
  free(z_ref);
  exit(0);
}
