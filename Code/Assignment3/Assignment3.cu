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

__global__ void dotprod(float *x, float *dot)
{
    __shared__ int product[THREADS_PER_BLOCK];

    int index = threadIdx.x + blockIdx.x * blockDim.x; //index

    product[threadIdx.x] = x[index] * x[index];

    if(index==0) *dot = 0; //Ask one thread to set c to zero.

    //Make sure every thread has finished
    __syncthreads();

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

float norm(float * x, int n) {
  float *dot;
  float *dev_x , *dev_dot;
  int size = n * sizeof (int );


  // Allocating Cuda memory
  cudaMalloc ( (void **)& dev_x , size );
  cudaMalloc ( (void **)& dev_dot , sizeof (int ) );
  // Allocating memory
  dot = (float *) malloc (sizeof (float ) );
  // Copying values
  cudaMemcpy (dev_x , x, size, cudaMemcpyHostToDevice );
  dotprod<<< N/THREADS_PER_BLOCK , THREADS_PER_BLOCK >>>(dev_x, dev_dot);
  cudaMemcpy ( dot, dev_dot , sizeof (int ) , cudaMemcpyDeviceToHost );

  return sqrt(*dot);
}


int main(void) {
  float *x;
  float Norm_out;
  int size = N * sizeof (int );

  x = (float *) malloc ( size );

  for (int i = 0; i < N; ++i){
    x[i] = i;
  }

  Norm_out = norm(x, N);
  printf("\ndot(a,b) = %f\n", Norm_out);
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


