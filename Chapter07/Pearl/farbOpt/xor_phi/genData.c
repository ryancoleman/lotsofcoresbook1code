#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

float truthTable[] = {
    0.1, 0.1, 0.1,
    0.1, 0.9, 0.9,
    0.9, 0.1, 0.9,
    0.9, 0.9, 0.1};

int main(int argc, char *argv[])
{
  if(argc < 4) {
    fprintf(stderr,"Use: file nExamples/4 variance\n");
    return -1;
  }
  char *filename=argv[1];
  FILE *fn=stdout;

  if(strcmp("-", filename) != 0)
    fn=fopen(filename,"w");

  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }
  int32_t nDup = atoi(argv[2]);
  float variance = atof(argv[3]);

  // write header info
  int32_t nInput=2; fwrite(&nInput,sizeof(int32_t), 1, fn);
  int32_t nOutput=1; fwrite(&nOutput,sizeof(int32_t), 1, fn);
  int32_t nExamples=nDup*4; fwrite(&nExamples,sizeof(int32_t), 1, fn);
  
  int nTable = sizeof(truthTable)/sizeof(float);
  float out[nTable];
  for(int32_t i=0; i < nDup; i++) {
    for(int j=0; j < nTable; j++) {
      float r = variance*(random()/(double)RAND_MAX);
      out[j] = truthTable[j] +  r;
    }
    fwrite(out, sizeof(float), nTable, fn);
  }
  if(fn != stdout) fclose(fn);
  return 0;
}
