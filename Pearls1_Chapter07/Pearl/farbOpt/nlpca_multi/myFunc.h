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

  // information  for multi-Devices
  int nDevices;
  int startOffset;
  struct userData* dev;

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
inline float G(float x) { return( x ) ;} 
#define G_ESTIMATE 0 
#elif USE_TANH
char *desc="generated_func tanh()";
inline float G(float x) { return( tanhf(x) ) ;} 
#define G_ESTIMATE 7 // estimate 7 flops for G
#elif LOGISTIC
char *desc="generated func logistic()";
inline float G(float x) { return( 1.f/(1.f+expf(-x)) ) ;} 
#define G_ESTIMATE 7 // estimate flops for G
#else // Use Elliott function
char *desc="generated func Eliott activation: x/(1+fabsf(x))";
inline float G(float x) { return( x/(1.f+fabsf(x)) ) ;} 
#define G_ESTIMATE 3 // estimate flops for G
#endif

// This file defines the function to be evaluated
#include "fcn.h"

// The offload objective function
double _objFunc(unsigned int n,  const double * restrict x,
		double * restrict grad, void * restrict my_func_data)
{
  double err;
  userData_t *uData = (userData_t *) my_func_data;

  // convert from double to float for speed
  for(int i=0; i < N_PARAM; i++) uData->param[i]=x[i];
  
#ifdef __INTEL_OFFLOAD

  // **** Start of Offload section ******
  
  double partial[uData->nDevices];
  float * restrict param = uData->param;
  int loadFlag[uData->nDevices];

  // asynchronous transfer of parameters to all devices
  for(int device=0; device < uData->nDevices; device++) {
    loadFlag[device]=device;
#pragma offload_transfer target(mic:device) in(param:length(N_PARAM) REUSE) signal(loadFlag+device)
    {}
  }

  // perform the computations asynchronously
  for(int device=0; device < uData->nDevices; device++) {
    float * restrict deviceExample = uData->dev[device].example; // workaround for compiler bug
    int nDeviceExamples = uData->dev[device].nExamples; // workaround for compiler bug
#pragma offload target(mic:device) in(param:length(0) REUSE) in(deviceExample:length(0) REUSE) \
  out(err : into(partial[device])) signal(partial+device) wait(loadFlag+device)
    {
      err=0.; // initialize error here for offload
#pragma omp parallel for reduction(+ : err)
      for(int i=0; i < nDeviceExamples; i++) {
	float d=myFunc(i, param, deviceExample, nDeviceExamples, NULL);
	err += d*d;
      }
    }
  }
  err=0.;
  for(int device=0; device < uData->nDevices; device++) {
#pragma offload target(mic:device) wait(partial+device)
    {}
    err += partial[device];
  }

  // **** End of Offload section ******

#else
  err=0.; // initialize error here in case offload selected
  
#pragma omp parallel for reduction(+ : err)
  for(int i=0; i < uData->nExamples; i++) {
    float d=myFunc(i, uData->param, uData->example, uData->nExamples, NULL);
    err += d*d;
  }
#endif

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

#ifdef __INTEL_OFFLOAD
  // free example vector if using offload mode
  __declspec(align(64)) float * restrict example = uData->example;
  // free on the host
  if(uData->example) free(uData->example);
  uData->example=NULL;
  
  uData->param=NULL;
  if(uData->dev) {
    for(int i=0; i < uData->nDevices; i++) {
      if(uData->dev[i].example) {
	float *pt=uData->dev[i].example;
#pragma offload target(mic:i) in(pt : length(0) FREE)
	{} 
	free(uData->dev[i].example);
      }
      if(uData->dev[i].param) {
	float *pt=uData->dev[i].param;
#pragma offload target(mic:i) in(pt : length(0) FREE)
	{} 
      }
    }
  }
#endif

  // free on the host
  if(uData->example) free(uData->example); uData->example=NULL;
  if(uData->param) free(uData->param); uData->param=NULL;
}

// partitions data and asynchronously loads onto the devices
void offloadData(userData_t *uData)
{
#ifdef __INTEL_OFFLOAD
  double startTime=getTime();
  uData->nDevices = _Offload_number_of_devices();

  if(uData->nDevices == 0) {
    fprintf(stderr,"No devices found!\n");
    exit -1;
  }

  fprintf(stderr,"Number of devices %d\n",uData->nDevices);
  uData->dev = calloc(uData->nDevices, sizeof(userData_t));
  if(!uData->dev) {
    fprintf(stderr,"Out of memory!\n");
    exit -1;
  }

  // Partition examples across multiple devices
  for(int i=0; i < uData->nDevices; i++) 
    uData->dev[i].nExamples = (uData->nExamples/uData->nDevices);
  for(int i=0; i < (uData->nExamples % uData->nDevices); i++)
    uData->dev[i].nExamples++; // fill in non-multiples of nDevices
  
  // fill in remaining uData information
  int index=0;
  for(int i=0; i < uData->nDevices; i++) {
    uData->dev[i].startOffset = index;
    index += uData->dev[i].nExamples;
    uData->dev[i].param = uData->param;
    fprintf(stderr,"Device %d startOffset %d nExamples %d\n", i,uData->dev[i].startOffset, uData->dev[i].nExamples);
  }
  
  // asynchronous alloc and transfer

  // prepare to move blocks of data to each device. (This is faster than a sliced move)
  int loadFlag[uData->nDevices];
  for(int i=0; i < uData->nDevices; i++) {
    int Xsiz = (N_INPUT+N_OUTPUT)*EXAMPLE_SIZE* uData->dev[i].nExamples;
    // For convenience, use c99 to cast pointer to a multidimensional array
    float (*data)[uData->nExamples] = (float (*)[uData->nExamples])  uData->example;

    uData->dev[i].example=(float*) memalign(64,Xsiz*sizeof(float));
    if(!uData->dev[i].example) {
      fprintf(stderr,"Not enough memory for the device copies of the examples!\n");
      exit(1);
    }
    float *pt = uData->dev[i].example; 
    for(int row=0; row < (N_INPUT+N_OUTPUT); row++) {
      int nCol = uData->dev[i].startOffset + uData->dev[i].nExamples;
      for(int col=uData->dev[i].startOffset; col < nCol; col++) {
	*(pt++) = data[row][col];
      }
    }
    pt = uData->dev[i].example; // workaround for compiler bug
    float *param = uData->param; // workaround for compiler bug
#pragma offload_transfer target(mic:i) in(pt : length(Xsiz) ALLOC) in(param: length(N_PARAM) ALLOC) signal(loadFlag+i)
    {}
  }
  
  // wait for all async transfers to finish
  for(int i=0; i < uData->nDevices; i++) {
#pragma offload target(mic:i) wait(loadFlag+i)
    {} 
  }
  uData->timeDataLoad = getTime() - startTime;
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

  // check if reading from stdin
  if(strcmp("-", filename) != 0)
    fn=fopen(filename,"r");

  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }

  // read the header information
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
  
#ifdef __INTEL_OFFLOAD
  offloadData(uData);
#endif
  if(fn!=stdin) fclose(fn);
}

