#ifndef OOCLU_C_H
#define OOCLU_C_H 1

#ifdef __cplusplus
extern "C" {
#endif
void Cpcgemm_hhd(char transA, char transB, int m, int n, int k,
        float *alpha,  float *A, int ia,int ja,int *descA,
                        float *B, int ib,int jb,int *descB,
        float *beta,   cuComplex *C, int ic,int jc,int *descC );


void pcgetrf_gpu( int *m, int *n, 
              cuComplex *A, int *ia, int *ja, int *descA,
                  int *ipiv, int *info );

void pcgetrf_gpu2( int *m, int *n,
              cuComplex *A, int *ia, int *ja, int *descA,
              float *Ah, int *iah_in, int *jah_in, int *desc_Ah,
              int *ipiv_Ah_, int *info );
void Cpcgecopy_d2h( int m, int n,
                 cuComplex *A, int ia, int ja, int *descA,
                 float *B, int ib, int jb, int *descB );

void Cpcgecopy_d2h_async( int m, int n,
                 cuComplex *A, int ia, int ja, int *descA,
                 float *B, int ib, int jb, int *descB );
void Cpcgecopy_h2d( int m, int n,
                 float *A, int ia, int ja, int *descA,
                 cuComplex *B, int ib, int jb, int *descB );

void Cpcgecopy_h2d_async( int m, int n,
                 float *A, int ia, int ja, int *descA,
                 cuComplex *B, int ib, int jb, int *descB );
void Cpclaprnt( int m, int n,  float *A, int ia, int ja, int *descA, 
                char *cmatnm );
void Cpcswap_gpu( int n, cuComplex *A, int ia,int ja,int *descA, int incA,
                    cuComplex *B, int ib,int jb,int *descB, int incB );

void  Cpclaswp_gpu(  char direct, char rowcol,
                 int nn, cuComplex *A, int ia, int ja, int *descA,
                 int k1, int k2, int *ipiv );

#ifdef __cplusplus
}
#endif



#endif
