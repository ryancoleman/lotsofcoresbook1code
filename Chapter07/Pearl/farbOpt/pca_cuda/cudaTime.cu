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

__device__
inline float G(float x) { return( x ) ;} 
#define G_ESTIMATE 0 

// helper macros to index into the example array
#define IN(i,nExamples,j)  (i*nExamples+j)
#define OUT(i,nExamples,j)  ((i+N_INPUT)*nExamples+j)
inline double getTime() { return(omp_get_wtime());}

#define N_CONCURRENT_BLOCKS (13*16)
//#define N_CONCURRENT_BLOCKS (65535)
__device__
#include "fcn.h"

__global__ void d_objFunc(float* d_param, float *d_example, int nExamples, float* out)
{
  int tid = blockIdx.x*blockDim.x + threadIdx.x;
  if(tid==0) *out=0.f;

  register float partial=0.f;
  while(tid < nExamples) {
    float d= myFunc(tid, d_param, d_example, nExamples, NULL);
    partial += d*d;
    tid += blockDim.x * gridDim.x;
  }
  // each thread puts its local sum into shared memory
  tid = blockIdx.x*blockDim.x + threadIdx.x;
  atomicAdd(out, partial);
}

float *d_out=NULL;
float objFunc(float* param, float* d_param, float *d_example, int nExamples)
{
  float sum=0.f;
  //int nBlocks = nExamples/nThreads + ((nExamples % nThreads>0)?1:0);
  //fprintf(stderr,"nExamples %d nBlocks %d, nThreads %d\n",nExamples,nBlocks,nThreads);
  cudaError_t err;
  err=cudaMemcpy(d_param, param, sizeof(float)*N_PARAM, cudaMemcpyDefault);
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA error (cudaMemcpy param): %s\n", cudaGetErrorString(err));
    exit(-1);
  }

  d_objFunc<<<N_CONCURRENT_BLOCKS, 32>>>(d_param, d_example, nExamples,d_out);
  err=cudaGetLastError();
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  cudaMemcpy(&sum, d_out, sizeof(float), cudaMemcpyDeviceToHost);
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
  cudaMalloc((void**) &d_out, sizeof(float));
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
 
  return 0;
}
