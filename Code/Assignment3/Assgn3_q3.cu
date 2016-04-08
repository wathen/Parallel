#include <cuda.h>
#include <curand_kernel.h>

#include <stdio.h>
#include <math.h>
#include "time_it.h"
#include <iostream>
using namespace std;
#define N (16*16)
#define THREADS_PER_BLOCK 8

typedef curandState randstate_t;


struct kernel_arg {
    float *x;
    unsigned int n;
    unsigned int blocksize;
};

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


__global__ void setup_kernel(uint n, randstate_t *state)
{
  uint myId = blockDim.x * blockIdx.x + threadIdx.x;
  /* Each thread gets same seed, a different sequence number, no offset */
  if(myId < n)
    curand_init(1234, myId, 0, &state[myId]);
}


__global__ void rand_cuda(randstate_t *state, unsigned int n, unsigned int nthreads, unsigned int m, float *x) {
    unsigned int myId = blockDim.x*blockIdx.x + threadIdx.x;
  if(myId < n) {
    for(int j = 0; j < m; j++) {  // each word generated by this thread
      x[nthreads*j + myId] = curand(state);
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

float rndm(float *x, unsigned int n, unsigned int MaxBlks, unsigned int blksize) {
    float *dev_x;
    int size = n*sizeof(float);

    randstate_t *dev_randState;
    uint sz_rnd_st = n*sizeof(dev_randState);
    HANDLE_ERROR(cudaMalloc((void **)(&dev_randState), sz_rnd_st));
    unsigned int m = ceil((float)(n)/float(MaxBlks*blksize));
    setup_kernel<<<MaxBlks, blksize>>>(n, dev_randState);

    cudaMalloc((void**)(&dev_x), size);
    rand_cuda<<<MaxBlks, blksize>>>(dev_randState, n, blksize, m, dev_x);
    cudaMemcpy(x, dev_x, size, cudaMemcpyDeviceToHost);
    return *x;
}

// void do_timing_logistic(void *void_arg) {
//   struct kernel_arg_logistic *argk = (struct kernel_arg_logistic *)(void_arg);
//   logistic(argk->x, argk->a, argk->n, argk->m, argk->z, argk->blocksize);
//   cudaDeviceSynchronize();
// }

int main(int argc, char *argv[]) {
    unsigned int n = 10000;
    unsigned int blksize = 8;
    float *x;
    int size = n*sizeof(float);
    x  = (float *)malloc(size);

    cudaDeviceProp prop;

    int ndev;
    cudaGetDeviceCount(&ndev);
    cudaGetDeviceProperties(&prop, 0);
    printf("The GPU is a %s\n", prop.name);
    printf("Cuda capability %d.%d.\n", prop.major, prop.minor);
    unsigned int MaxBlks = prop.maxThreadsPerBlock;

    print_vec(x, 10, "%f", "x");
    rndm(x, n, MaxBlks, blksize);
    print_vec(x, 10, "%f", "x");
    return 0;
}
