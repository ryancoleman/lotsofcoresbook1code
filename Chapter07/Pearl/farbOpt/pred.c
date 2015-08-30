#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// tell myFunc that this is a prediction function
#define DO_PRED
#include "myFunc.h"

void readParam(char* filename, int nParam, float* param)
{
  FILE *fn=fopen(filename,"r");
  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }
  int parmInFile;
  fread(&parmInFile,sizeof(uint32_t), 1, fn);
  if(parmInFile != N_PARAM) {
    fprintf(stderr,"Number of inputs incorrect!\n");
    exit(1);
  }
  fread(param,sizeof(float), nParam, fn);
}

int main(int argc, char* argv[])
{
  userData_t uData = {0};

  if(argc < 3) {
    fprintf(stderr,"Use: paramFile dataFile\n");
    return -1;
  }

  printf("myFunc %s\n", desc);
  init(argv[2],&uData); 
  printf("nExamples %d\n", uData.nExamples);
  
  __declspec(align(64)) float x[N_PARAM];
  readParam(argv[1], N_PARAM, x);

  if(N_OUTPUT == 0) { // special case for autoencoders
    float pred[N_INPUT];
    for(int exIndex=0; exIndex < uData.nExamples; exIndex++) {
      float err = myFunc(exIndex, x, uData.example, uData.nExamples, pred);
      
      printf("input ");
      for(int j=0; j < N_INPUT; j++)
	printf("%g ", uData.example[IN(j,uData.nExamples,exIndex)]);
      
      printf("PredOutput ");
      for(int j=0; j < N_INPUT; j++) printf("%g ", pred[j]);
      printf("\n");
    }

  } else {

    float pred[N_OUTPUT+1];
    for(int exIndex=0; exIndex < uData.nExamples; exIndex++) {
      float err = myFunc(exIndex, x, uData.example, uData.nExamples, pred);
      
      printf("input ");
      for(int j=0; j < N_INPUT; j++)
	printf("%g ", uData.example[IN(j,uData.nExamples,exIndex)]);
      
      printf(" KnownOutput ");
      for(int j=0; j < N_OUTPUT; j++)
	printf("%g ", uData.example[OUT(j,uData.nExamples,exIndex)]);
      
      printf("PredOutput ");
      for(int j=0; j < N_OUTPUT; j++) printf("%g ", pred[j]);
      printf("\n");
    }

  }

  return 0;
}
