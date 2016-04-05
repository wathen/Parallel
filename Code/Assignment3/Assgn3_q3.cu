#include <cuda.h>
#include <curand_kernel.h>

#include <stdio.h>
#include <math.h>
#include "time_it.h"
#include <iostream>
using namespace std;
#define N (16*16)
#define THREADS_PER_BLOCK 8

struct kernel_arg {
    float *x;
    unsigned int n;
    unsigned int blocksize;
};


__global__ void rand_cuda(unsigned int n, float *x) {
    unsigned int myId = blockDim.x*blockIdx.x + threadIdx.x;
    if(myId < n){

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

void rndm(float *x, unsigned int n, unsigned int blocksize) {

    float *dev_x;
    int size = n*sizeof(float);
    cudaMalloc((void**)(&dev_x), size);
    cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);
    unsigned int nblks = ceil(((float)(n))/((float)(blocksize)));
    rand_cuda<<<nblks , blocksize>>>(n, dev_x);
    cudaMemcpy(x, dev_x, size, cudaMemcpyDeviceToHost);
}

// void do_timing_logistic(void *void_arg) {
//   struct kernel_arg_logistic *argk = (struct kernel_arg_logistic *)(void_arg);
//   logistic(argk->x, argk->a, argk->n, argk->m, argk->z, argk->blocksize);
//   cudaDeviceSynchronize();
// }

int main(int argc, char *argv[]) {
    unsigned int n = 10000;
    unsigned int blksize = 8;
    return 0;
}
