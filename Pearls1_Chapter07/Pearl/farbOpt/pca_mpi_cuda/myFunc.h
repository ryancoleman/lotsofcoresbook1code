#ifndef USE_CUDA
#define __device__
#endif
#define restrict
#define __declspec(x)
// Rob Farber
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <malloc.h>
#include <math.h>
#include <omp.h>

#define MIC_DEV 0

#define	ALLOC alloc_if(1) free_if(0)
#define	FREE alloc_if(0) free_if(1)
#define	REUSE alloc_if(0) free_if(0)

// Use a struct to pass and get data from the objective function
typedef struct userData { 
  // Data information
  int nExamples;
  __declspec(align(64)) float * restrict example; 
  __declspec(align(64)) float * restrict param;
#ifdef USE_CUDA
  float *d_example;
  float *d_param;
  double *d_out;
#endif

  // Timing information
  int isWarmup;
  double timeObjFunc;
  int countObjFunc;
  double timeDataLoad;
  double minTime, maxTime;
} userData_t;

// function to measure wall clock time
inline double getTime() { return(omp_get_wtime());}

#pragma offload_attribute (push, target (mic))

// helper macros to index into the example array
#define IN(i,nExamples,j)  (i*nExamples+j)
#define OUT(i,nExamples,j)  ((i+N_INPUT)*nExamples+j)

// Define the Sigmoid
#ifdef USE_LINEAR
char *desc="generated_PCA_func LINEAR()";
__device__
inline float G(float x) { return( x ) ;} 
#define G_ESTIMATE 0 
#elif USE_TANH
char *desc="generated_func tanh()";
__device__
inline float G(float x) { return( tanhf(x) ) ;} 
#define G_ESTIMATE 7 // estimate 7 flops for G
#elif LOGISTIC
char *desc="generated func logistic()";
__device__
inline float G(float x) { return( 1.f/(1.f+expf(-x)) ) ;} 
#define G_ESTIMATE 7 // estimate flops for G
#else // Use Elliott function
char *desc="generated func Eliott activation: x/(1+fabsf(x))";
__device__
inline float G(float x) { return( x/(1.f+fabsf(x)) ) ;} 
#define G_ESTIMATE 3 // estimate flops for G
#endif

// This file defines the function to be evaluated
__device__
#include "fcn.h"

#ifdef USE_CUDA
#define N_CONCURRENT_BLOCKS (15*16)
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
  if(threadIdx.x < 16) smem[threadIdx.x] += smem[threadIdx.x + 16];
  if(threadIdx.x < 8) smem[threadIdx.x] += smem[threadIdx.x + 8];
  if(threadIdx.x < 4) smem[threadIdx.x] += smem[threadIdx.x + 4];
  if(threadIdx.x < 2) smem[threadIdx.x] += smem[threadIdx.x + 2];
  if(threadIdx.x < 1) smem[threadIdx.x] += smem[threadIdx.x + 1];

  // each thread puts its local sum into shared memory
  if(threadIdx.x == 0) atomicAdd(out, smem[0]);
}
#endif

// The offload objective function
double _objFunc(unsigned int n,  const double * restrict x,
		double * restrict grad, void * restrict my_func_data)
{
  double err;
  userData_t *uData = (userData_t *) my_func_data;

  // convert from double to float for speed
  for(int i=0; i < N_PARAM; i++) uData->param[i]=x[i];
  
  int nExamples = uData->nExamples;
  // compiler workaround
  __declspec(align(64)) float * restrict example = uData->example;
  __declspec(align(64)) float * restrict param = uData->param; 
#pragma acc data copyin(param[0:N_PARAM-1]) pcopyin(example[0:nExamples*EXAMPLE_SIZE-1])
#pragma offload target(mic:MIC_DEV) in(param:length(N_PARAM) REUSE) \
                                    out(err) in(example:length(0) REUSE)
  {
    err=0.; // initialize error here in case offload selected
    
#ifdef USE_CUDA
  cudaError_t ret;
  ret=cudaMemcpy(uData->d_param, param, sizeof(float)*N_PARAM, cudaMemcpyDefault);
  if( ret != cudaSuccess) {
    fprintf(stderr,"CUDA error (cudaMemcpy param): %s\n", cudaGetErrorString(ret));
    exit(-1);
  }

  d_objFunc<float,double,32><<<N_CONCURRENT_BLOCKS, 32>>>(uData->d_param, uData->d_example, nExamples,uData->d_out);
  ret=cudaGetLastError();
  if( ret != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s\n", cudaGetErrorString(ret));
    exit(-1);
  }
  double tmp;
  ret = cudaMemcpy(&tmp, uData->d_out, sizeof(double), cudaMemcpyDeviceToHost);
  if( ret != cudaSuccess) {
    fprintf(stderr,"CUDA memcpy(sum): %s\n", cudaGetErrorString(ret));
    exit(-1);
  }
  err=tmp;
#else
#pragma acc parallel loop num_gangs(15*16) vector_length(32) reduction(+:err)
#pragma omp parallel for reduction(+ : err)
    for(int i=0; i < nExamples; i++) {
      float d=myFunc(i, param, example, nExamples, NULL);
      err += d*d;
    }
#endif
  }

  //fprintf(stderr,"err %g\n", sqrt(err));
  return sqrt(err);
}
#pragma offload_attribute (pop)

// The optizimation library callable objective function that gathers timing information
double objFunc(unsigned int n,  const double * restrict x, 
	       double * restrict grad, void * restrict my_func_data)
{
  if(grad) {
    fprintf(stderr,"Gradient not implemented!\n");
    exit(1);
  }
  
  userData_t *uData = (userData_t *) my_func_data;

  double runTime=getTime();
  double err =  _objFunc(n,x,grad,my_func_data);
  runTime = getTime() - runTime;

  if(!uData->isWarmup) {
    // Note a maxTime of zero means this is the first call 
    if(uData->maxTime == 0.) {
      uData->maxTime = uData->minTime = runTime;
    }
    uData->maxTime = (uData->maxTime > runTime)?uData->maxTime:runTime;
    uData->minTime = (uData->minTime < runTime)?uData->minTime:runTime;
    
    uData->timeObjFunc += runTime;
    uData->countObjFunc++;
  }

  return( err );
}

// Called to free memory and report timing information
void fini(userData_t *uData)
{
  int nThreads=0;
  // Intel recommended way to get the number of threads in offload mode.
#pragma offload target(mic:MIC_DEV) out(nThreads)
  {
#pragma omp parallel
    {
#pragma omp single
      {
	nThreads = omp_get_num_threads();
      }
    }
  }
  // Ouput some information
  if(!uData->isWarmup) {
    printf("number OMP threads %d\n", nThreads);
    printf("DataLoadTime %g\n", uData->timeDataLoad);
    printf("AveObjTime %g, countObjFunc %d, totalObjTime %g\n",
	   uData->timeObjFunc/uData->countObjFunc, uData->countObjFunc, uData->timeObjFunc);
#ifdef FLOP_ESTIMATE
    printf("Estimated flops in myFunc %d, estimated average GFlop/s %g\n", FLOP_ESTIMATE,
	   (((double)uData->nExamples*FLOP_ESTIMATE)/(uData->timeObjFunc/uData->countObjFunc)/1.e9) );
    printf("Estimated maximum GFlop/s %g, minimum GFLop/s %g\n",
	   (((double)uData->nExamples*FLOP_ESTIMATE)/(uData->minTime)/1.e9),
	   (((double)uData->nExamples*FLOP_ESTIMATE)/(uData->maxTime)/1.e9) );
  }
#endif

  // free if using offload mode
  __declspec(align(64)) float * restrict example = uData->example;// compiler workaround
  __declspec(align(64)) float * restrict param = uData->param;// compiler workaround
#pragma offload target(mic:MIC_DEV) in(example: length(0) FREE) in(param : length(0) FREE)
  {} 

  // free on the host
  if(uData->example) free(uData->example); uData->example=NULL;
  if(uData->param) free(uData->param); uData->param=NULL;
}

void offloadData(userData_t *uData)
{
#ifdef USE_CUDA
  cudaError_t err;
  long exSize=sizeof(float)*EXAMPLE_SIZE*uData->nExamples;
  if( (err=cudaMalloc((void**) &uData->d_example, exSize)) != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  if( (err=cudaMalloc((void**) &uData->d_param, sizeof(float)*N_PARAM)) != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  if( (err=cudaMalloc((void**) &uData->d_out, sizeof(double))) != cudaSuccess) {
    fprintf(stderr,"CUDA error: %s\n", cudaGetErrorString(err));
    exit(-1);
  }
  err=cudaMemcpy(uData->d_example, uData->example, exSize, cudaMemcpyDefault);
  if( err != cudaSuccess) {
    fprintf(stderr,"CUDA error (cudaMemcpy example): %s\n", cudaGetErrorString(err));
    exit(-1);
  }
#endif
#ifdef __INTEL_OFFLOAD
  int nDevices =_Offload_number_of_devices();

  if(nDevices == 0) {
    fprintf(stderr,"No devices found!\n");
    exit -1;
  }

  // If necessary, perform offload transfer and allocation
  double startOffload=getTime();
  __declspec(align(64)) float * restrict example = uData->example; // compiler workaround
  __declspec(align(64)) float * restrict param = uData->param; // compiler workaround
  int Xsiz = uData->nExamples*EXAMPLE_SIZE; // compiler workaround
  // Note: the in for param just allocates memory on the device
#pragma offload target(mic:MIC_DEV) in(example: length(Xsiz) ALLOC) in(param : length(N_PARAM) ALLOC)
  {} 
  // set data load time if using offload mode
  uData->timeDataLoad = getTime() - startOffload;
#endif
}


// loads the binary file of the form:
//    nInput, nOutput, nExamples
//    Input [0] [0:nExamples]
//    Input [1] [0:nExamples]
//    ...
//    Output [0] [0:nExamples]
//    Output [1] [0:nExamples]
//    ...
void init(char*filename, userData_t *uData)
{
  FILE *fn=stdin;

#ifdef MPI_NUM_COPROC_PER_NODE
#ifdef USE_CUDA
  //fprintf(stderr,"Using GPU %d\n", mpiRank % MPI_NUM_COPROC_PER_NODE);
  cudaSetDevice(mpiRank % MPI_NUM_COPROC_PER_NODE);
#endif
#endif

  // check if reading from stdin
  if(strcmp("-", filename) != 0)
    fn=fopen(filename,"r");

  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }

  // read the header information
  double startTime=getTime();
  int32_t nInput, nOutput;
  int32_t nExamples;
  fread(&nInput,sizeof(int32_t), 1, fn);
  if(nInput != N_INPUT) {
    fprintf(stderr,"Number of inputs incorrect!\n");
    exit(1);
  }
  fread(&nOutput,sizeof(int32_t), 1, fn);
  if(nOutput != N_OUTPUT) {
    fprintf(stderr,"Number of outputs incorrect!\n");
    exit(1);
  }
  fread(&nExamples,sizeof(int32_t), 1, fn);
  if(nExamples <= 0) {
    fprintf(stderr,"Number of examples incorrect!\n");
    exit(1);
  }
  uData->nExamples = nExamples;
  
  // aligned allocation of the data
  uData->example=(float*) memalign(64,nExamples*EXAMPLE_SIZE*sizeof(float));
  if(!uData->example) {
    fprintf(stderr,"Not enough memory for examples!\n");
    exit(1);
  }
  // aligned allocation of the on-device parameters
  uData->param=(float*) memalign(64,N_PARAM*sizeof(float));
  if(!uData->param) {
    fprintf(stderr,"Not enough memory for the parameters!\n");
    exit(1);
  }

  // read the data
  for(int exIndex=0; exIndex < uData->nExamples; exIndex++) {
    for(int i=0; i < nInput; i++) 
      fread(&uData->example[IN(i,uData->nExamples, exIndex)],1, sizeof(float), fn);
    for(int i=0; i < nOutput; i++) 
      fread(&uData->example[OUT(i,uData->nExamples, exIndex)],1, sizeof(float), fn);
  }
  
  // offload the data
  double startOffload=getTime();
  __declspec(align(64)) float * restrict example = uData->example; // compiler workaround
  __declspec(align(64)) float * restrict param = uData->param; // compiler workaround
  int Xsiz = uData->nExamples*EXAMPLE_SIZE; // compiler workaround
  // Note: the in just allocates memory on the device
#pragma offload target(mic:MIC_DEV) in(example: length(Xsiz) ALLOC) in(param : length(N_PARAM) ALLOC)
  {} 

#ifdef USE_CUDA
  offloadData(uData);
#endif
  uData->timeDataLoad = getTime() - startTime;

  if(fn!=stdin) fclose(fn);
}

// no file io test for timings on big supercomputers
void init_noIO(int _nExamples, userData_t *uData)
{
  FILE *fn=stdin;

#ifdef MPI_NUM_COPROC_PER_NODE
#ifdef USE_CUDA
  //fprintf(stderr,"Using GPU %d\n", mpiRank % MPI_NUM_COPROC_PER_NODE);
  cudaSetDevice(mpiRank % MPI_NUM_COPROC_PER_NODE);
#endif
#endif

  // read the header information
  double startTime=getTime();
  int32_t nInput=N_INPUT, nOutput=N_OUTPUT;
  int32_t nExamples = _nExamples;

  if(nExamples <= 0) {
    fprintf(stderr,"Number of examples incorrect!\n");
    exit(1);
  }
  uData->nExamples = nExamples;
  
  // aligned allocation of the data
  uData->example=(float*) memalign(64,nExamples*EXAMPLE_SIZE*sizeof(float));
  if(!uData->example) {
    fprintf(stderr,"Not enough memory for examples!\n");
    exit(1);
  }
  // aligned allocation of the on-device parameters
  uData->param=(float*) memalign(64,N_PARAM*sizeof(float));
  if(!uData->param) {
    fprintf(stderr,"Not enough memory for the parameters!\n");
    exit(1);
  }

  // randomize the data
  for(int exIndex=0; exIndex < uData->nExamples; exIndex++) {
    for(int i=0; i < nInput; i++) 
	uData->example[IN(i,uData->nExamples, exIndex)] = ((float)random())/((float)RAND_MAX);
    for(int i=0; i < nOutput; i++) 
      uData->example[OUT(i,uData->nExamples, exIndex)] = ((float)random())/((float)RAND_MAX);
  }
  
  // offload the data
  double startOffload=getTime();
  __declspec(align(64)) float * restrict example = uData->example; // compiler workaround
  __declspec(align(64)) float * restrict param = uData->param; // compiler workaround
  int Xsiz = uData->nExamples*EXAMPLE_SIZE; // compiler workaround
  // Note: the in just allocates memory on the device
#pragma offload target(mic:MIC_DEV) in(example: length(Xsiz) ALLOC) in(param : length(N_PARAM) ALLOC)
  {} 

#ifdef USE_CUDA
  offloadData(uData);
#endif
  uData->timeDataLoad = getTime() - startTime;
}

