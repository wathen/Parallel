#include <stdio.h>
#include <math.h>
#include "time_it.h"
#include <iostream>
using namespace std;
#define N (16*16)
#define THREADS_PER_BLOCK 8

__global__ void logistic_cuda(unsigned int n, unsigned int m, unsigned int M, float a, float *x, float *z) {
  unsigned int blockBase = blockDim.x * blockIdx.x;
  unsigned int myId = blockBase + threadIdx.x;
  uint n_threadsTotal = gridDim.x * blockDim.x;
 
  float X, Z;
  if(myId < n){
    for (int j = 0; j < M; j++) {
        X = x[n_threadsTotal*j + myId];
        for (int i = 1; i < m; ++i) {
          Z = a*X*(1.0f - X);
          X = Z;
        }
        z[n_threadsTotal*j + myId] = Z;
      } 
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
struct kernel_arg_norm {
  float *x;
  unsigned int n;
  unsigned int blocksize;
  unsigned int MaxBlks;
};

struct kernel_arg_logistic {
  float *x;
  unsigned int a;
  unsigned int n;
  unsigned int m;
  float *z;
  unsigned int blocksize;
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

__global__ void dotprod(unsigned int n, unsigned int m, float *x, float *y, float *z) {
  unsigned int blockBase = blockDim.x * blockIdx.x;
  unsigned int myId = blockBase + threadIdx.x;
  uint n_threadsTotal = gridDim.x * blockDim.x;
  unsigned int M = min(blockDim.x, n - blockBase);
  float temp = 0.0;
  if(n_threadsTotal*(m-1)+myId < n)
    for (int j = 0; j < m; ++j) {
            temp = temp + x[n_threadsTotal*j + myId]*x[n_threadsTotal*j + myId];
    }
    y[myId] = temp;
  __syncthreads();
  z[blockIdx.x] = 0.0;
  for(int i = 0; i < M; i++){
    z[blockIdx.x] = z[blockIdx.x] + temp;//y[blockBase+i];
    printf("\nz = %f, y = %f", z[blockIdx.x], y[blockBase+i]);
  }

}
__global__ void dot(unsigned int n, unsigned int m, float *x, float *product, float *c) {
  unsigned int blockBase = blockDim.x * blockIdx.x;
  unsigned int myId = blockBase + threadIdx.x;
  uint n_threadsTotal = gridDim.x * blockDim.x;
  unsigned int M = min(blockDim.x, n - blockBase);
 
  product[threadIdx.x] = 0.0;
  if(n_threadsTotal*(m-1)+myId < n) { 
    for (int j = 0; j < m; ++j) {
        product[threadIdx.x] = product[threadIdx.x] + x[n_threadsTotal*j + myId]*x[n_threadsTotal*j + myId];
    }
  }
  if(myId==0) *c = 0.0;
  __syncthreads();

  if( 0 == threadIdx.x ) {
      float sum = 0.0;
      for(int j=0; j < M; j++) {
          sum += product[j];
      }
      atomicAdd(c,sum);
  }
}
float norm_ref(float *x, unsigned int n) {
  float sum = 0.0;
  for(int i = 0; i < n; i++)
    sum += x[i] * x[i];
  return(sqrt(sum));
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
  unsigned int nblks = ceil(((float)(n))/((float)(blocksize)));
  if (1024 < nblks) {
    nblks = 1024;
  }
  unsigned int M = ceil((float)(n)/(float)(nblks*blocksize));
  cudaMalloc((void**)(&dev_x), size);
  cudaMalloc((void**)(&dev_z), size);
  cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);
  logistic_cuda<<<nblks , blocksize>>>(n, m, M, a, dev_x, dev_z);
  cudaMemcpy(z, dev_z, size, cudaMemcpyDeviceToHost);
}



float norm(float * x, unsigned int n, unsigned int blocksize, unsigned int MaxBlks) {

  float *z, *zOut;
  float *dev_x, *dev_z, *dev_y;
  unsigned int nblks = ceil(((float)(n))/((float)(blocksize)));
  if (MaxBlks < nblks) {
    nblks = MaxBlks;
  }
  int size = n*sizeof(float);
  unsigned int m = ceil((float)(n)/(float)(nblks*blocksize));
  zOut = (float *) malloc(sizeof(float));
  z = (float *) malloc(size);
  cudaMalloc((void**)(&dev_x), size);
  cudaMalloc((void**)(&dev_z), size);
  cudaMalloc((void**)(&dev_y), size);
  cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);
  dot<<<nblks , blocksize>>>(n, m, dev_x, dev_y, dev_z);
  cudaMemcpy(z, dev_z, sizeof(float), cudaMemcpyDeviceToHost);
  reduce_sum<<<1,blocksize>>>(nblks, dev_z);
  cudaMemcpy(zOut, dev_z, sizeof(float), cudaMemcpyDeviceToHost);
  return sqrt(zOut[0]);
}

void do_timing_norm(void *void_arg) {
  struct kernel_arg_norm *argk = (struct kernel_arg_norm *)(void_arg);
  norm(argk->x, argk->n, argk->blocksize, argk->MaxBlks);
  cudaDeviceSynchronize();
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
  struct kernel_arg_norm argk_n;
  struct kernel_arg_logistic argk_l;
  struct time_it_raw *tr = time_it_create(10);
  struct time_it_stats stats_n;
  struct time_it_stats stats_l;

  int size = n*sizeof(float);
  x  = (float *)malloc(size);
  z  = (float *)malloc(size);
  z_ref = (float *)malloc(size);
  printf("%d\n",n );
  for(int i = 0; i < n; i++) {
    x[i] =  (float)rand() / (float)RAND_MAX;
    //x[i] = 1.0f*i;
  }

  int ndev;
  cudaGetDeviceCount(&ndev);
  cudaGetDeviceProperties(&prop, 0);
  printf("The GPU is a %s\n", prop.name);
  printf("Cuda capability %d.%d.\n", prop.major, prop.minor);
  unsigned int MaxBlks = prop.maxThreadsPerBlock;
  //unsigned int MaxBlks = 500;

  //printf("stop, %d", MaxBlks);
  float p_norm = norm(x, n, blocksize, MaxBlks);
  z_ref[0] = norm_ref(x, n);

  if (MaxBlks > blocksize) {
      MaxBlks = blocksize;
  }

  printf("\n Parallel = %5.5f,   Sequential = %5.5f",p_norm,z_ref[0]);
  printf("\n Norm error = %3.4e\n\n",1.0e-6*sqrt(n)*max(abs(p_norm-z_ref[0]), 1.0));


  argk_n.n = n;
  argk_n.x = x;
  argk_n.blocksize = blocksize;
  argk_n.MaxBlks = MaxBlks;
  time_it_run(tr, do_timing_norm, (void *)(&argk_n));
  time_it_get_stats(tr, &stats_n);


  unsigned int M = ceil(((float)(n))/((float)(MaxBlks*MaxBlks)));


  float *L;
  L  = (float*)malloc(size);
  logistic(x, a, n, m, L, blocksize);
  logistic_ref(n, m, a, x, z);

  argk_l.x = x;
  argk_l.a = a;
  argk_l.n = n;
  argk_l.m = m;
  argk_l.z = z;
  argk_l.blocksize = blocksize;
  time_it_run(tr, do_timing_logistic, (void *)(&argk_l));
  time_it_get_stats(tr, &stats_l);
  printf("\n\nLogistic calculations:   # operations = %10.4i   mean time = %1.4e  time per op = %1.4e, Gflops = %5.3f\n\n", 3*n*m, stats_l.mean, stats_l.mean/(3*n*m), (3*(float)n*(float)m)/(stats_l.mean*pow(10,9)));
  float Left_over_block = roundf((float)blocksize*( (((float)(n))/((float)(blocksize))) - floor(((float)(n))/((float)(blocksize)))));
  float nblocks = floor(((float)(n))/((float)(blocksize)));

  int NormWork = (2*(nblocks*(2*blocksize-1) + 2*Left_over_block-1))*M;
  printf("Norm calculations:       # operations = %11.4i   mean time = %1.4e  time per op = %1.4e, Gflops = %5.3f\n\n", NormWork, stats_n.mean, stats_n.mean/(NormWork),pow(stats_n.mean/(NormWork),-1)/pow(10,9));

  free(x);
  free(z_ref);
  exit(0);
}
