#ifndef OOCLU_Z_H
#define OOCLU_Z_H 1

#ifdef __cplusplus
extern "C" {
#endif
void Cpzgemm_hhd(char transA, char transB, int m, int n, int k,
        double *alpha,  double *A, int ia,int ja,int *descA,
                        double *B, int ib,int jb,int *descB,
        double *beta,   cuDoubleComplex *C, int ic,int jc,int *descC );


void pzgetrf_gpu( int *m, int *n, 
              cuDoubleComplex *A, int *ia, int *ja, int *descA,
                  int *ipiv, int *info );

void pzgetrf_gpu2( int *m, int *n,
              cuDoubleComplex *A, int *ia, int *ja, int *descA,
              double *Ah, int *iah_in, int *jah_in, int *desc_Ah,
              int *ipiv_Ah_, int *info );
void Cpzgecopy_d2h( int m, int n,
                 cuDoubleComplex *A, int ia, int ja, int *descA,
                 double *B, int ib, int jb, int *descB );

void Cpzgecopy_d2h_async( int m, int n,
                 cuDoubleComplex *A, int ia, int ja, int *descA,
                 double *B, int ib, int jb, int *descB );
void Cpzgecopy_h2d( int m, int n,
                 double *A, int ia, int ja, int *descA,
                 cuDoubleComplex *B, int ib, int jb, int *descB );

void Cpzgecopy_h2d_async( int m, int n,
                 double *A, int ia, int ja, int *descA,
                 cuDoubleComplex *B, int ib, int jb, int *descB );
void Cpzlaprnt( int m, int n,  double *A, int ia, int ja, int *descA, 
                char *cmatnm );
void Cpzswap_gpu( int n, cuDoubleComplex *A, int ia,int ja,int *descA, int incA,
                    cuDoubleComplex *B, int ib,int jb,int *descB, int incB );

void  Cpzlaswp_gpu(  char direct, char rowcol,
                 int nn, cuDoubleComplex *A, int ia, int ja, int *descA,
                 int k1, int k2, int *ipiv );

#ifdef __cplusplus
}
#endif



#endif
