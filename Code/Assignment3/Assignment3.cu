#include <stdio.h>
#include <math.h>
#include <iostream>
#define N (512*512)
#define THREADS_PER_BLOCK 128
#define A 1
using namespace std;

__global__ void logistic(unsigned int n, float a, float *x, float *z) {
  unsigned int myId = blockDim.x*blockIdx.x + threadIdx.x;
  if(myId < n)
    z[myId] = a*x[myId]*(1 - x[myId]);
}


__device__ void reduce_sum_dev(unsigned int n, float *x) {
  unsigned int myId = threadIdx.x;
  for(unsigned int m = n >> 1; m > 0; m = n >> 1) {
    n -= m;
    __syncthreads();
    if(myId < m)
      x[myId] += x[myId+n];
  }
}

__global__ void reduce_sum(uint n, float *x) {
  reduce_sum_dev(n, x);
}

float reduce_sum_ref(uint n, float *x) {
  float sum = 0.0;
  for(int i = 0; i < n; i++)
    sum += x[i];
  return(sum);
}

/******************************
*  dotprod: just like it sounds
*    Simple version: we only handle one block of threads
*******************************/

__global__ void dotprod(uint n, float *x, float *z) {
  uint blockBase = blockDim.x * blockIdx.x;
  uint myId = blockBase + threadIdx.x;
  uint m = min(blockDim.x, n - blockBase);

  if(myId < n)
    x[myId] *= x[myId];
  reduce_sum_dev(m, &(x[blockBase]));
  if((myId < n) && (threadIdx.x == 0))
    z[blockIdx.x] = x[myId];
}

__global__ void my_dotprod(float *x, float *dot)
{
    __shared__ int product[THREADS_PER_BLOCK];

    int index = threadIdx.x + blockIdx.x * blockDim.x; //index

    product[threadIdx.x] = x[index] * x[index];


    //Make sure every thread has finished
    __syncthreads();

    if(index==0) *dot = 0; //Ask one thread to set c to zero.

    //Sum the elements serially to obtain dot product
    if( 0 == threadIdx.x ) //Every block to do c += sum
    {
        int sum = 0;
        for(int j=0; j < THREADS_PER_BLOCK; j++) {
          sum += product[j];
        }
        atomicAdd(dot,sum);
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

float norm(float * x, unsigned int n) {
    float *dot;
    float *dev_x , *dev_dot;
    int size = n * sizeof (int );
    // Allocating Cuda memory
    cudaMalloc ( (void **)& dev_x , size );
    cudaMalloc ( (void **)& dev_dot , size );
    // Allocating memory
    dot = (float *) malloc (sizeof (float ) );
    // Copying values
    cudaMemcpy (dev_x , x, size, cudaMemcpyHostToDevice );
    dotprod<<< N/THREADS_PER_BLOCK , THREADS_PER_BLOCK >>>(N, dev_x, dev_dot);
    reduce_sum<<<1,N/THREADS_PER_BLOCK>>>(N/THREADS_PER_BLOCK, dev_dot);
    cudaMemcpy ( dot, dev_dot , sizeof (float ) , cudaMemcpyDeviceToHost );
    return *dot;
}

float seq_norm(float *x, unsigned int n) {
    float y = 0.0;
    for(int i = 0; i<n; i++)
        y += x[i] * x[i];
    return (y);
}

int main(void) {
    float *x;
    float *p_norm, *s_norm;
    int size = N * sizeof(float);
    p_norm = (float *) malloc (size);
    s_norm = (float *) malloc (size);


    x = (float *) malloc ( size );

    for (int i = 0; i < N; ++i){
        x[i] = i;
    }

    print_vec(x, N, "%5.3f", "y");
    p_norm[0] = norm(x, N);
    s_norm[0] = seq_norm(x, N);

    printf("\nParallel = %f, Sequential = %f \n", p_norm, s_norm);
}

// int main( void ) {
//   float *x, *y, *dot;
//   float * dev_x , * dev_y , * dev_dot;
//   float a;
//   int size = N * sizeof (int );

//   // allocate device copies of a, b, c
//   cudaMalloc ( (void **)& dev_x , size );
//   cudaMalloc ( (void **)& dev_y , size );
//   cudaMalloc ( (void **)& dev_dot , sizeof (int ) );

//   x = (float *) malloc ( size );
//   y = (float *) malloc ( size );
//   dot = (float *) malloc (sizeof (float ) );

//   cout<<N;

//   for (int i = 0; i < N; ++i){
//     x[i] = i;
//   }
//   a = 1.0;
//   cudaMemcpy (dev_x , x, size, cudaMemcpyHostToDevice );
//   cudaMemcpy (dev_y , y, size, cudaMemcpyHostToDevice );
//   // cudaMemcpy (dev_a , a, sizeof (int ), cudaMemcpyHostToDevice );

//   // launch dot() kernel
//   dotprod<<< N/THREADS_PER_BLOCK , THREADS_PER_BLOCK >>>(dev_x, dev_dot);
//   logistic<<< N/THREADS_PER_BLOCK , THREADS_PER_BLOCK >>>(N, a, dev_x, dev_y);

//   // copy device result back to host copy of c
//   cudaMemcpy ( dot, dev_dot , sizeof (int ) , cudaMemcpyDeviceToHost );

//   // fprintf(" Ans = %15.10d\n", c);
//   // cout<< *c;

//   printf("\ndot(a,b) = %f\n", *dot); //Output result
//   print_vec(x, N, "%5.3f", "y");


//   free( x ); free( y ); free( dot );
//   cudaFree (dev_x);
//   cudaFree (dev_y);
//   cudaFree (dev_dot);
//   return 0;
// }


