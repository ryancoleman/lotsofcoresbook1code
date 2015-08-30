#ifndef OOCLU_D_H
#define OOCLU_D_H 1

#ifdef __cplusplus
extern "C" {
#endif
void Cpdgemm_hhd(char transA, char transB, int m, int n, int k,
        double *alpha,  double *A, int ia,int ja,int *descA,
                        double *B, int ib,int jb,int *descB,
        double *beta,   double *C, int ic,int jc,int *descC );

void Cpdsyrk_hhd(char uplo, char transA, int m, int n, int k,
        double *alpha,  double *A, int ia,int ja,int *descA,
        double *beta,   double *C, int ic,int jc,int *descC );

void pdpotrf_ooc2( char *uplo_in, int *n_in, 
             double *A, int *ia_in, int *ja_in, int *descA, 
	     int *memsize_in, int *info );

void pdgetrf_gpu( int *m, int *n, 
              double *A, int *ia, int *ja, int *descA,
                  int *ipiv, int *info );

void pdgetrf_gpu2( int *m, int *n,
              double *A, int *ia, int *ja, int *descA,
              double *Ah, int *iah_in, int *jah_in, int *desc_Ah,
              int *ipiv_Ah_, int *info );

void pdpotrf_gpu2( char *uplo, int *m, int *n,
              double *A, int *ia, int *ja, int *descA,
              double *Ah, int *iah_in, int *jah_in, int *desc_Ah,
              int *info );

void pdcopymatrix( double *A, int ia, int ja, int *descA);


void Cpdgecopy_d2h( int m, int n,
                 double *A, int ia, int ja, int *descA,
                 double *B, int ib, int jb, int *descB );

void Cpdgecopy_d2h_async( int m, int n,
                 double *A, int ia, int ja, int *descA,
                 double *B, int ib, int jb, int *descB );
void Cpdgecopy_h2d( int m, int n,
                 double *A, int ia, int ja, int *descA,
                 double *B, int ib, int jb, int *descB );

void Cpdgecopy_h2d_async( int m, int n,
                 double *A, int ia, int ja, int *descA,
                 double *B, int ib, int jb, int *descB );
void Cpdlaprnt( int m, int n,  double *A, int ia, int ja, int *descA, 
                char *cmatnm );
void Cpdswap_gpu( int n, double *A, int ia,int ja,int *descA, int incA,
                    double *B, int ib,int jb,int *descB, int incB );

void  Cpdlaswp_gpu(  char direct, char rowcol,
                 int nn, double *A, int ia, int ja, int *descA,
                 int k1, int k2, int *ipiv );

#ifdef __cplusplus
}
#endif



#endif
