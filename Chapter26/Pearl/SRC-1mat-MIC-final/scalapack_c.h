#ifndef SCALAPACK_C_H
#define SCALAPACK_C_H 1



#if defined(ADD_)
#define scalapack_pcgemm  pcgemm_
#define scalapack_pcgetrf  pcgetrf_
#define scalapack_pcswap  pcswap_
#define scalapack_pcgeadd pcgeadd_
#define scalapack_pcgemm  pcgemm_
#define scalapack_pctrsm  pctrsm_
#define scalapack_pclaswp pclaswp_
#define scalapack_cgebr2d cgebr2d_
#define scalapack_cgebs2d cgebs2d_
#define scalapack_pclaprnt pclaprnt_
#define scalapack_pcelget pcelget_
#define scalapack_pclapiv pclapiv_

#else
#define scalapack_pcgemm  pcgemm
#define scalapack_pcgetrf  pcgetrf
#define scalapack_pcswap  pcswap
#define scalapack_pcgeadd pcgeadd
#define scalapack_pcgemm  pcgemm
#define scalapack_pctrsm  pctrsm
#define scalapack_pclaswp pclaswp
#define scalapack_cgebr2d cgebr2d
#define scalapack_cgebs2d cgebs2d
#define scalapack_pclaprnt pclaprnt
#define scalapack_pcelget pcelget
#define scalapack_pclapiv pclapiv
#endif





extern "C" {

void scalapack_pcgetrf( int *m, int *n, 
          float *A, int *ia, int *ja, int *descA, int *ipiv, int *info );

void scalapack_pcswap( int *n, 
  float *A, int *ia, int *ja, int *descA, int *incA,
  float *B, int *ib, int *jb, int *descB, int *incB );

void scalapack_pctrsm( char *side, char *uplo, char *trans, char *diag,
         int *m, int *n, float *alpha,
         float *A, int *ia, int *ja, int *descA,
         float *B, int *ib, int *jb, int *descB );

void scalapack_pcgeadd( char *trans, int *m, int *n,
           float *alpha,
           float *A, int *ia, int *ja, int *descA,
           float *beta,
           float *B, int *ib, int *jb, int *descB );

void scalapack_pclaswp( char *direc, char *rowcol,
             int *n, float *A, int *ia, int *ja, int *descA,
             int *k1, int *k2, int *ipiv );

void scalapack_cgebs2d( int *icontxt, char *scope, char *top,
           int *m, int *n, float *A, int *lda );

void scalapack_cgebr2d( int *icontxt, char *scope, char *top,
           int *m, int *n, float *A, int *lda,
           int *rsrc, int *csrc );

void scalapack_pcgemm( char *transA, char *transB, int *m, int *n, int *k,
       float *alpha,    float *A, int *ia, int *ja, int *descA,
                         float *B, int *ib, int *jb, int *descB,
       float *beta,     float *C, int *ic, int *jc, int *descC );

void scalapack_pclapiv( char *direc, char *rowcol, char *pivroc,
     int *m, int *n, float *A, int *ia, int *ja, int *descA,
     int *ipiv, int *ip, int *jp, int *descip, int *iwork );

void scalapack_pclaprnt(int *m, int *n, float *A, int *ia, int *ja, int *descA,
    int *irprnt, int *icprnt, char *cmatnm, int *nout, float *work,   int *dummy );


void scalapack_pcelget( char *scope, char *top, float *alpha, 
    float *A, int *ia, int *ja, int *descA );

}

#endif
