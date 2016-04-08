#ifndef OOCLU_S_H
#define OOCLU_S_H 1

#ifdef __cplusplus
extern "C" {
#endif
void Cpsgemm_hhd(char transA, char transB, int m, int n, int k,
        float *alpha,  float *A, int ia,int ja,int *descA,
                        float *B, int ib,int jb,int *descB,
        float *beta,   float *C, int ic,int jc,int *descC );

void Cpssyrk_hhd(char uplo, char transA, int m, int n, int k,
        float *alpha,  float *A, int ia,int ja,int *descA,
        float *beta,   float *C, int ic,int jc,int *descC );

void psgetrf_gpu( int *m, int *n, 
              float *A, int *ia, int *ja, int *descA,
                  int *ipiv, int *info );

void psgetrf_gpu2( int *m, int *n,
              float *A, int *ia, int *ja, int *descA,
              float *Ah, int *iah_in, int *jah_in, int *desc_Ah,
              int *ipiv_Ah_, int *info );

void pspotrf_gpu2( char *uplo, int *m, int *n,
              float *A, int *ia, int *ja, int *descA,
              float *Ah, int *iah_in, int *jah_in, int *desc_Ah,
              int *info );

void Cpsgecopy_d2h( int m, int n,
                 float *A, int ia, int ja, int *descA,
                 float *B, int ib, int jb, int *descB );

void Cpsgecopy_d2h_async( int m, int n,
                 float *A, int ia, int ja, int *descA,
                 float *B, int ib, int jb, int *descB );
void Cpsgecopy_h2d( int m, int n,
                 float *A, int ia, int ja, int *descA,
                 float *B, int ib, int jb, int *descB );

void Cpsgecopy_h2d_async( int m, int n,
                 float *A, int ia, int ja, int *descA,
                 float *B, int ib, int jb, int *descB );
void Cpslaprnt( int m, int n,  float *A, int ia, int ja, int *descA, 
                char *cmatnm );
void Cpsswap_gpu( int n, float *A, int ia,int ja,int *descA, int incA,
                    float *B, int ib,int jb,int *descB, int incB );

void  Cpslaswp_gpu(  char direct, char rowcol,
                 int nn, float *A, int ia, int ja, int *descA,
                 int k1, int k2, int *ipiv );

#ifdef __cplusplus
}
#endif



#endif
