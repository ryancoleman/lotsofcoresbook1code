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

#pragma offload_attribute (push, target (mic))

typedef struct userData { 
  // Data information
  int nExamples;
  __declspec(align(64)) float * restrict example;
  // Timing information
  double timeObjFunc;
  int countObjFunc;
  double timeDataLoad;
  double minTime, maxTime;
} userData_t;

// helper macros to index into the example array
#define IN(i,nExamples,j)  (i*nExamples+j)
#define OUT(i,nExamples,j)  ((i+N_INPUT)*nExamples+j)

inline double getTime() { return(omp_get_wtime());}

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

#include "fcn.h"

double objFunc(unsigned n,  const double * restrict x, double * restrict grad,
	       void * restrict my_func_data)
{
  double err;
  userData_t *uData = (userData_t *) my_func_data;

  if(grad) {
    //fprintf(stderr,"Gradient not implemented!\n");
    exit(1);
  }
  
  double runTime=getTime();
  int nExamples = uData->nExamples;
  __declspec(align(64)) float * restrict example = uData->example;

#pragma offload target(mic:MIC_DEV) in(x:length(N_PARAM)) out(err) in(example:length(0) REUSE)
  {
    err=0.; // initialize error here in case offload selected
    
    // convert from double to float for speed
    __declspec(align(64)) float P[N_PARAM];
    for(int i=0; i < N_PARAM; i++) P[i]=x[i];
    
#pragma omp parallel for reduction(+ : err)
    for(int i=0; i < nExamples; i++) {
      float d=myFunc(i, P, example, nExamples, NULL);
      err += d*d;
    }
  }

  runTime = getTime() - runTime;

  // Note a maxTime of zero means this is the first call 
  if(uData->maxTime == 0.) {
    uData->maxTime = uData->minTime = runTime;
  }
  uData->maxTime = (uData->maxTime > runTime)?uData->maxTime:runTime;
  uData->minTime = (uData->minTime < runTime)?uData->minTime:runTime;

  uData->timeObjFunc += runTime;
  uData->countObjFunc++;

  return sqrt(err);
}
#pragma offload_attribute (pop)

void fini(userData_t *uData)
{
  int nThreads=0;
  // The intel recommended way to get the number of threads in offload mode.
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
#endif

  // free example vector if using offload mode
  __declspec(align(64)) float * restrict example = uData->example;
#pragma offload target(mic:MIC_DEV) in(example: length(0) FREE)
  {} 
  // free on the host
  free(example); uData->example=NULL;
}

void init(char*filename, userData_t *uData)
{
  FILE *fn=stdin;

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

  // read the data
  for(int exIndex=0; exIndex < uData->nExamples; exIndex++) {
    for(int i=0; i < nInput; i++) 
      fread(&uData->example[IN(i,uData->nExamples, exIndex)],1, sizeof(float), fn);
    for(int i=0; i < nOutput; i++) 
      fread(&uData->example[OUT(i,uData->nExamples, exIndex)],1, sizeof(float), fn);
  }
  
  double startOffload=getTime();
  __declspec(align(64)) float * restrict example = uData->example;
  int Xsiz = uData->nExamples*EXAMPLE_SIZE; // single variable works around a compiler bug
#pragma offload target(mic:MIC_DEV) in(example: length(Xsiz) ALLOC)
  {} 
  uData->timeDataLoad = getTime() - startTime;

  if(fn!=stdin) fclose(fn);
}


