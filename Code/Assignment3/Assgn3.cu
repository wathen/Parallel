#include <stdio.h>
#include <math.h>

#define N (16*16)
#define THREADS_PER_BLOCK 16

/**************************************************************
 *  reduce_sum: compute the sum of the elements of an array
 *    Simple version: we only handle one block of threads
 ***************************************************************/

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

__global__ void dotprod(unsigned int n, float *x, float *y, float *z) {
  unsigned int blockBase = blockDim.x * blockIdx.x;
  unsigned int myId = blockBase + threadIdx.x;
  unsigned int m = min(blockDim.x, n - blockBase);

  if(myId < n)
    x[myId] *= y[myId];
  reduce_sum_dev(m, &(x[blockBase]));
  if((myId < n) && (threadIdx.x == 0))
    z[blockIdx.x] = x[myId];
}

float dotprod_ref(unsigned int n, float *x, float *y) {
  float sum = 0.0;
  for(int i = 0; i < n; i++)
    sum += x[i] * y[i];
  return(sum);
}

/********************************************************************
 *  dotprod2: just like it sounds
 *    First phase of a dotprod.  Each block computes its part of the
 *    dot-product and stores the result in the z array.  A subsequent
 *    launch of the reduce_sum kernel completes the calculation.
 ********************************************************************/

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

/*****************************************************
*  near(x, y): true if x and y are "nearly" equal
******************************************************/
int near(unsigned int n, float x, float y) {
  return(abs(x-y) < max(10.0, sqrt((float)n))*1.0e-7*max(1.0, max(abs(x), abs(y))));
}

int main(int argc, char **argv) {
  unsigned int n = N;
  unsigned int nn = n;
  float *x, *y, *z, *z_ref;
  float *dev_x, *dev_y, *dev_z;
  cudaDeviceProp prop;

  int ndev;
  cudaGetDeviceCount(&ndev);
  if(ndev < 1) {
    fprintf(stderr, "No CUDA device found\n");
    exit(-1);
  }
  cudaGetDeviceProperties(&prop, 0);

  int size = n*sizeof(float);
  x  = (float *)malloc(size);
  y  = (float *)malloc(size);
  z  = (float *)malloc(size);
  z_ref = (float *)malloc(size);

  // Use a logistic map to make some pseudo-random numbers
  // It's fast, but the distribution isn't very uniform, and
  //   the other statistical properties are lousy.  But it's
  //   fast, and that's all we need for some simple tests.
  x[0] = 0.123;
  y[0] = sqrt(0.3);
  for(int i = 1; i < n; i++) {
    x[i] = 3.8*x[i-1]*(1.0 - x[i-1]);
    y[i] = 3.9*y[i-1]*(1.0 - y[i-1]);
  }

  printf("The GPU is a %s\n", prop.name);
  printf("Cuda capability %d.%d.\n", prop.major, prop.minor);
  print_vec(x, min(10, n), "%5.3f", "x");
  print_vec(y, min(10, n), "%5.3f", "y");

  cudaMalloc((void**)(&dev_x), size);
  cudaMalloc((void**)(&dev_y), size);
  cudaMalloc((void**)(&dev_z), size);
  cudaMemcpy(dev_x, x, size, cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y, y, size, cudaMemcpyHostToDevice);


  // unsigned int blksize;
  // unsigned int nblks;
  // blksize = prop.maxThreadsPerBlock;
  // nblks = ceil(((float)(n))/((float)(blksize)));
  dotprod<<<N/THREADS_PER_BLOCK , THREADS_PER_BLOCK>>>(n, dev_x, dev_y, dev_z);
  reduce_sum<<<1,N/THREADS_PER_BLOCK>>>(N/THREADS_PER_BLOCK, dev_z);
  cudaMemcpy(z, dev_z, sizeof(float), cudaMemcpyDeviceToHost);
  z_ref[0] = dotprod_ref(n, x, y);


  for(int i = 0; i < nn; i++) { // check the result
    if(!near(n, z[i], z_ref[i])) {
      fprintf(stderr, "ERROR: i=%d: z[i] = %15.10f, z_ref[i] = %15.10f\n", z[i], z_ref[i]);
      exit(-1);
    }
  }
  print_vec(z, min(10, nn), "%5.3f", "z");
  printf("The results match!\n");

  cudaFree(dev_x);
  cudaFree(dev_y);
  cudaFree(dev_z);
  free(x);
  free(y);
  free(z);
  free(z_ref);
  exit(0);
}
