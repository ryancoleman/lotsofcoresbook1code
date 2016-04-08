#ifndef SCALAPACK_S_H
#define SCALAPACK_S_H 1



#if defined(ADD_)
#define scalapack_psgemm  psgemm_
#define scalapack_pssyrk  pssyrk_
#define scalapack_psgetrf  psgetrf_
#define scalapack_psswap  psswap_
#define scalapack_psgeadd psgeadd_
#define scalapack_psgemm  psgemm_
#define scalapack_pstrsm  pstrsm_
#define scalapack_pslaswp pslaswp_
#define scalapack_sgebr2d cgebr2d_
#define scalapack_sgebs2d cgebs2d_
#define scalapack_pslaprnt pslaprnt_
#define scalapack_pselget pselget_
#define scalapack_pslapiv pslapiv_

#else
#define scalapack_psgemm  psgemm
#define scalapack_psgetrf  psgetrf
#define scalapack_pspotrf  pspotrf
#define scalapack_psswap  psswap
#define scalapack_psgeadd psgeadd
#define scalapack_psgemm  psgemm
#define scalapack_pstrsm  pstrsm
#define scalapack_pslaswp pslaswp
#define scalapack_sgebr2d cgebr2d
#define scalapack_sgebs2d cgebs2d
#define scalapack_pslaprnt pslaprnt
#define scalapack_pselget pselget
#define scalapack_pslapiv pslapiv
#endif





extern "C" {

void scalapack_psgetrf( int *m, int *n, 
          float *A, int *ia, int *ja, int *descA, int *ipiv, int *info );

void scalapack_pspotrf( char *uplo,  int *n, 
           float *A, int *ia, int *ja, int *descA, int *info );

void scalapack_psswap( int *n, 
  float *A, int *ia, int *ja, int *descA, int *incA,
  float *B, int *ib, int *jb, int *descB, int *incB );

void scalapack_pstrsm( char *side, char *uplo, char *trans, char *diag,
         int *m, int *n, float *alpha,
         float *A, int *ia, int *ja, int *descA,
         float *B, int *ib, int *jb, int *descB );

void scalapack_psgeadd( char *trans, int *m, int *n,
           float *alpha,
           float *A, int *ia, int *ja, int *descA,
           float *beta,
           float *B, int *ib, int *jb, int *descB );

void scalapack_pslaswp( char *direc, char *rowcol,
             int *n, float *A, int *ia, int *ja, int *descA,
             int *k1, int *k2, int *ipiv );

void scalapack_sgebs2d( int *icontxt, char *scope, char *top,
           int *m, int *n, float *A, int *lda );

void scalapack_sgebr2d( int *icontxt, char *scope, char *top,
           int *m, int *n, float *A, int *lda,
           int *rsrc, int *csrc );

void scalapack_psgemm( char *transA, char *transB, int *m, int *n, int *k,
       float *alpha,    float *A, int *ia, int *ja, int *descA,
                         float *B, int *ib, int *jb, int *descB,
       float *beta,     float *C, int *ic, int *jc, int *descC );

void scalapack_pslapiv( char *direc, char *rowcol, char *pivroc,
     int *m, int *n, float *A, int *ia, int *ja, int *descA,
     int *ipiv, int *ip, int *jp, int *descip, int *iwork );

void scalapack_pslaprnt(int *m, int *n, float *A, int *ia, int *ja, int *descA,
    int *irprnt, int *icprnt, char *cmatnm, int *nout, float *work,   int *dummy );


void scalapack_pselget( char *scope, char *top, float *alpha, 
    float *A, int *ia, int *ja, int *descA );

}

#endif
