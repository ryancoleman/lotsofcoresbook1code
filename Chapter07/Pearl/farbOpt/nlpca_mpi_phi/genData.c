#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// get a uniform random number between -1 and 1 
inline float f_rand() {
  return 2*(rand()/((float)RAND_MAX)) -1.;
}
void genData(FILE *fn, int nVec, float xVar)
{
  float xMax = 1.1; float xMin = -xMax;
  float xRange = (xMax - xMin);

  // write header info
  uint32_t nInput=2; fwrite(&nInput,sizeof(int32_t), 1, fn);
  uint32_t nOutput=0; fwrite(&nOutput,sizeof(int32_t), 1, fn);
  uint32_t nExamples=nVec; fwrite(&nExamples,sizeof(int32_t), 1, fn);

  for(int i=0; i < nVec; i++) {
    float t = xRange * f_rand();
    float z1 = t +  xVar * f_rand();
#ifdef USE_LINEAR
    float z2 = t +  xVar * f_rand();
#else
    float z2 = t*t*t +  xVar * f_rand();
#endif
    fwrite(&z1, sizeof(float), 1, fn);
    fwrite(&z2, sizeof(float), 1, fn);
  }
}

int main(int argc, char *argv[])
{
  if(argc < 4) {
    fprintf(stderr,"Use: filename nVec variance seed\n");
    exit(1);
  }
  char *filename=argv[1];
  FILE *fn=stdout;

  if(strcmp("-", filename) != 0)
    fn=fopen(filename,"w");

  if(!fn) {
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(1);
  }
  int nVec = atoi(argv[2]);
  float variance = atof(argv[3]);
  srand(atoi(argv[4]));
  genData(fn, nVec, variance);

  if(fn != stdout) fclose(fn);
  return 0;
}
