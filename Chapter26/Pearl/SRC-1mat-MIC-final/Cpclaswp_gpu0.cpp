#include "ooclu.h"

#define indx2(i,j,lda) ((i) + ((j)-1)*(lda))


#ifdef __cplusplus
extern "C"
#endif
void  Cpclaswp_gpu0(  int n, float *A, int ia, int ja, int *descA, 
                 int k1, int k2, int *ipiv, int incx )
{



int ix0 = k1;
int i1 = k1;
int i2 = k2;
int inc = 1;

int ip = 0;
int ix = 0;
int i = 0;


if (incx > 0) {
  ix0 = k1;
  i1 = k1;
  i2 = k2;
  inc = 1;
  }
else {
  ix0 = 1 + (1-k2)*incx;
  i1 = k2;
  i2 = k1;
  inc = -1;
}

ix = ix0;
if (inc > 0) {
 for( i=i1; i <= i2; i += inc) {
  ip = ipiv[ix-1];
  if (ip != i) {

    /*
    cublasDswap( n, &A[ indx2(ia-1+i,ja,lda)-1 ],lda, 
                    &A[ indx2((ia-1)+ip,ja,lda)-1 ], lda);
    */

    Cpcswap_gpu(n,  A, ia-1+i,ja, descA, descA[LLD_],
                   A, ia-1+ip,ja, descA, descA[LLD_] );
    };
  ix += incx;
 };
 }
else {
 for( i=i1; i2 <= i ; i += inc) {
  ip = ipiv[ix-1];
  if (ip != i) {
    /*
    cublasDswap( n, &A[ indx2(ia-1+i,ja,lda)-1 ],lda, 
                    &A[ indx2((ia-1)+ip,ja,lda)-1 ], lda);
    */

    Cpcswap_gpu( n,  A, ia-1+i,ja,descA, descA[LLD_],
                    A, ia-1+ip,ja,descA, descA[LLD_] );
    };
  ix += incx;
  };
}

 return;
}
