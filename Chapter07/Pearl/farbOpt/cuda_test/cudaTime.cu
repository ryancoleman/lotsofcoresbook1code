#define __declspec(x)
#define restrict
// Rob Farber
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>
#include <cuda.h>

#ifndef REDUCT_TYPE
#define REDUCT_TYPE double
#endif

__device__
inline float G(float x) { return( x ) ;} 
#define G_ESTIMATE 0 

// helper macros to index into the example array
#define IN(i,nExamples,j)  (i*nExamples+j)
#define OUT(i,nExamples,j)  ((i+N_INPUT)*nExamples+j)
inline double getTime() { return(omp_get_wtime());}

__device__ inline void atomicAdd (double *address, double value)
 {
   unsigned long long oldval, newval, readback; 

   oldval = __double_as_longlong(*address);
   newval = __double_as_longlong(__longlong_as_double(oldval) + value);
   while ((readback=atomicCAS((unsigned long long *)address, oldval, newval)) != oldval)
     {
      oldval = readback;
      newval = __double_as_longlong(__longlong_as_double(oldval) + value);
     }
 }
#ifndef N_CONCURRENT_BLOCKS
#define N_CONCURRENT_BLOCKS (13*16)
#endif

__device__
#include "fcn.h"

template <class T, class T1, unsigned int WARP_SIZE>
__global__ void d_objFunc(T* d_param, T *d_example, int nExamples, T1 *out)
{
  __shared__ T1 ssum[WARP_SIZE];
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if(tid==0) *out=0.f;

  if(threadIdx.x < WARP_SIZE)
    ssum[threadIdx.x] = 0.;

  register T1 partial=0.f;
  while(tid < nExamples) {
    T d= myFunc(tid, d_param, d_example, nExamples, NULL);
    partial += d*d;
    //partial += tid;
    tid += blockDim.x * gridDim.x;
  }
  // sum all the partials on each multiprocessor into shared memory
  ssum[threadIdx.x & (WARP_SIZE-1)] += partial;
  __syncthreads();

  tid = blockIdx.x*blockDim.x + threadIdx.x;
  volatile T1 *smem = ssum;
  // sum all threads in each multiprocessor into a single value
  // Note: this is easy to modify for 32-bit shuffle ops with CUDA 5.0
  if(threadIdx.x < 16) smem[threadIdx.x] += smem[threadIdx.x + 16];
  if(threadIdx.x < 8) smem[threadIdx.x] += smem[threadIdx.x + 8];
  if(threadIdx.x < 4) smem[threadIdx.x] += smem[threadIdx.x + 4];
  if(threadIdx.x < 2) smem[threadIdx.x] += smem[threadIdx.x + 2];
  if(threadIdx.x < 1) smem[threadIdx.x] += smem[threadIdx.x + 1];

  // each thread puts its local sum into shared memory
  if(threadIdx.x == 0) atomicAdd(out, smem[0]);
}

REDUCT_TYPE *d_out=NULL;
float objFunc(float* param, float* d_param, float *d_example, int nExamples)
{
  REDUCT_TYPE sum=0.f;
  //int nBlocks = nExamples/nThreads + ((nExamples % nThreads>0)?1:0);
  //fprintf(stderr,"nExamples %d nBlocks %d, nThreads %d\n",nExamples,nBlocks,nThreads);
  cudaError_t err;
  err=cudaMemcpy(d_param, param, sizeof(float)*N_PARAM, cudaMemcpyDefault);
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA error (cudaMemcpy param): %s\n", cudaGetErrorString(err));
    exit(-1);
  }

  d_objFunc<float,REDUCT_TYPE,32><<<N_CONCURRENT_BLOCKS, 32>>>(d_param, d_example, nExamples,d_out);
  err=cudaGetLastError();
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  cudaMemcpy(&sum, d_out, sizeof(REDUCT_TYPE), cudaMemcpyDeviceToHost);
  err=cudaGetLastError();
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA memcpy(sum): %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  return sum;
}

int main(int argc, char* argv[])
{
  if(argc < 3) {
    fprintf(stderr,"Use: device nTests nExamples\n");
    return -1;
  }
  int device=atoi(argv[1]);
  int nTests=atoi(argv[2]);
  int nExamples=atoi(argv[3]);

  cudaSetDevice(device);
  float *d_example=NULL, *d_param=NULL;
  cudaMalloc((void**) &d_out, sizeof(REDUCT_TYPE));
  cudaMalloc((void**) &d_example, sizeof(float)*nExamples*EXAMPLE_SIZE);
  cudaMalloc((void**) &d_param, sizeof(float)*N_PARAM);

  // ********* Allocation  *********************
  // aligned allocation of the data
  float *example = (float*) memalign(64,nExamples*EXAMPLE_SIZE*sizeof(float));
  if(!example) {
    fprintf(stderr,"Not enough memory for examples!\n");
    exit(1);
  }
  // aligned allocation of the on-device parameters
  float *param =  (float*)memalign(64,N_PARAM*sizeof(float));
  if(!param) {
    fprintf(stderr,"Not enough memory for the parameters!\n");
    exit(1);
  }

  // ********* Fill data  *********************
  // generate random examples
  // fill with random data
  for(int i=0; i < nExamples*EXAMPLE_SIZE; i++)
    example[i]= (rand()/RAND_MAX);
  
  // fill with random data
  for(int i=0; i < N_PARAM; i++) param[i] = 0.1*(rand()/(double)RAND_MAX);

  cudaMemcpy((void*)d_example, (void*)example, sizeof(float)*nExamples*EXAMPLE_SIZE, cudaMemcpyHostToDevice);
  cudaError_t err =cudaGetLastError();
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA error (cudaMemcpy example): %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  
  // Perform timing runs
  float result[nTests]; // Sanity check for consistency
  double perCall[nTests];
  double startTime=getTime();
    
  for(int i=0; i < nTests; i++) {
    double callStart=getTime();
    result[i] = objFunc(param, d_param, d_example, nExamples);
    perCall[i] = getTime() - callStart;
  }
  double totalTime=getTime()-startTime;
  
  printf("Runtime %g, aveObjRuntime %g\n",totalTime, totalTime/nTests);
  printf("ave GF/s %g\n",((double)nExamples*FLOP_ESTIMATE)/(totalTime/nTests)/1.e9);
  printf("result[0] %g\n", result[0]);
  long sum=0;
  for(int i=0; i < nExamples; i++) sum += i;
  printf("sum %g\n", (double) sum);
 
  return 0;
}
