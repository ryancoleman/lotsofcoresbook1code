#include <immintrin.h>
#include <omp.h>
#include "zorder2d.h"
#include "timer.h"

void pmat(int n, float *a) {
  int i,j,iad;
  return;
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      iad=zorder2d_c2i(j,i);
      printf("% 8.2f",a[iad]);
    }
    printf("\n");
  }
  printf("\n");
}

// transpose a z-block
void zmattr(float *ain, float *cout) {
#ifdef __MIC__
  __m512 a;
  __m512 at;
  __m512i pat={ 0, 2, 1, 3, 8, 10, 9, 11, 4, 6, 5, 7, 12, 14, 13, 15};
  a=_mm512_load_ps(ain);
  at=_mm512_castsi512_ps(_mm512_permutevar_epi32(pat,_mm512_castps_si512(a)));
  _mm512_store_ps(cout,at);
#else
  cout[0]=ain[0];
  cout[1]=ain[2];
  cout[2]=ain[1];
  cout[3]=ain[3];
  cout[4]=ain[8];
  cout[5]=ain[10];
  cout[6]=ain[9];
  cout[7]=ain[11];
  cout[8]=ain[4];
  cout[9]=ain[6];
  cout[10]=ain[5];
  cout[11]=ain[7];
  cout[12]=ain[12];
  cout[13]=ain[14];
  cout[14]=ain[13];
  cout[15]=ain[15];
#endif
}

int main(int argc, char **argv) {
  int i,j,k;
  int iad,uad,lad;
  int N=4;
  if(argc>1) {
    N=atoi(argv[1]);
  }
// for now, N must be a multiple of 4
  i=N/4;
  if(N>(i*4)) N=(i+1)*4;
  int N4=N/4;
  double l2=log(N)/log(2.0);
  int il2=(int)l2;
  if(il2<l2) il2++;
  int MSZ=pow(2,il2);
  printf("N=%d; N/4=%d; msize=%d\n",N,N4,MSZ);
// actual memory allocation must be power of 2 although N can be any multiple of 4
  float *a=(float *)_mm_malloc(MSZ*MSZ*sizeof(float),128);
  float *c=(float *)_mm_malloc(MSZ*MSZ*sizeof(float),128);
  int iv=0;
#pragma omp parallel for private(iad,i,j)
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      iad=zorder2d(j,i);
      a[iad]=j+i*N;
    }
  }
  pmat(N,a);
  tv ts;
  timer_reset(ts);
  timer_start(ts);
#pragma omp parallel for private(i,j,uad,lad)
  for(i=0; i<N4; i++) {
    for(j=0; j<N4; j++) {
      uad=zorder2d(j,i);
      lad=zorder2d(i,j);
      zmattr(&a[16*uad],&c[16*lad]);
    }
  }
  timer_stop(ts);
  pmat(N,c);
  printf("%d,%lf\n",N,timer_sec(ts));
  _mm_free(a);
  _mm_free(c);

}
