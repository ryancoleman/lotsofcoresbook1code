// Rob Farber
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <nlopt.h>

#include "myFunc.h"

void writeParam(char *filename, int nParam, double *x)
{
  FILE *fn=fopen(filename,"w");
  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }

  uint32_t n=nParam; // ensure size is uint32_t
  fwrite(&n,sizeof(uint32_t), 1, fn);
  for(int i=0; i < nParam; i++) {
    float tmp=x[i];
    fwrite(&tmp,sizeof(float), 1, fn);
  }
  fclose(fn);
}

int main(int argc, char* argv[])
{
  nlopt_opt opt;
  userData_t uData = {0};

  memset(&uData, 0, sizeof(userData_t));

  if(argc < 3) {
    fprintf(stderr,"Use: datafile paramFile\n");
    return -1;
  }
  init(argv[1],&uData);
  printf("myFunc %s\n", desc);
  printf("nExamples %d\n", uData.nExamples);
  printf("Number Parameters %d\n", N_PARAM);
  
  opt = nlopt_create(NLOPT_LN_PRAXIS, N_PARAM); // algorithm and dimensionality
  // NOTE: alternative optimization methods ...
  //opt = nlopt_create(NLOPT_LN_NEWUOA, N_PARAM);
  //opt = nlopt_create(NLOPT_LN_COBYLA, N_PARAM);
  //opt = nlopt_create(NLOPT_LN_BOBYQA, N_PARAM);
  //opt = nlopt_create(NLOPT_LN_AUGLAG, N_PARAM);

  nlopt_set_min_objective(opt, objFunc, (void*) &uData);
#if defined(MAX_RUNTIME) 
  fprintf(stderr,"Warning: performing a quick %d second test!\n", MAX_RUNTIME);
  nlopt_set_maxtime(opt, MAX_RUNTIME); // Use for running quick tests
#else
  fprintf(stderr,"MAX runtime %d seconds!\n", 120*60);
  nlopt_set_maxtime(opt, (120. * 60.)); // maximum runtime in seconds
#endif
  double minf; /* the minimum objective value, upon return */

  __declspec(align(64)) double x[N_PARAM];
  for(int i=0; i < N_PARAM; i++) x[i] = 0.1*(rand()/(double)RAND_MAX);

  int nExamples=uData.nExamples*EXAMPLE_SIZE;
  float* example=uData.example;
  float* param=uData.param;
#pragma acc data pcopyin(example[0:nExamples-1])
  {
    double startTime=getTime();
    int ret=nlopt_optimize(opt, x, &minf);
    printf("Optimization Time %g\n",getTime()-startTime);

    if (ret < 0) {
      printf("nlopt failed! ret %d\n", ret);
    } else {
      printf("found minimum %0.10g ret %d\n", minf,ret);
    }
  }
  writeParam(argv[2],N_PARAM, x);
  fini(&uData);
  nlopt_destroy(opt);

  return 0;
}
