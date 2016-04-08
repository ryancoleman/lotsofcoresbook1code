#ifndef SCALAPACK_D_H
#define SCALAPACK_D_H 1



#if defined(ADD_)
#define scalapack_pdgetrf  pdgetrf_
#define scalapack_pdpotrf  pdpotrf_
#define scalapack_pdswap  pdswap_
#define scalapack_pdgeadd pdgeadd_
#define scalapack_pdtradd pdtradd_
#define scalapack_pdgemm  pdgemm_
#define scalapack_pdsyrk  pdsyrk_
#define scalapack_pdtrsm  pdtrsm_
#define scalapack_pdlaswp pdlaswp_
#define scalapack_dgebr2d dgebr2d_
#define scalapack_dgebs2d dgebs2d_
#define scalapack_pdlaprnt pdlaprnt_
#define scalapack_pdelget pdelget_
#define scalapack_pdlapiv pdlapiv_

#else
#define scalapack_pdgetrf  pdgetrf
#define scalapack_pdpotrf  pdpotrf
#define scalapack_pdswap  pdswap
#define scalapack_pdgeadd pdgeadd
#define scalapack_pdtradd pdtradd
#define scalapack_pdgemm  pdgemm
#define scalapack_pdsyrk  pdsyrk
#define scalapack_pdtrsm  pdtrsm
#define scalapack_pdlaswp pdlaswp
#define scalapack_dgebr2d dgebr2d
#define scalapack_dgebs2d dgebs2d
#define scalapack_pdlaprnt pdlaprnt
#define scalapack_pdelget pdelget
#define scalapack_pdlapiv pdlapiv
#endif





extern "C" {

void scalapack_pdgetrf( int *m, int *n, 
          double *A, int *ia, int *ja, int *descA, int *ipiv, int *info );

void scalapack_pdpotrf( char *uplo,  int *n, 
           double *A, int *ia, int *ja, int *descA, int *info );

void scalapack_pdswap( int *n, 
  double *A, int *ia, int *ja, int *descA, int *incA,
  double *B, int *ib, int *jb, int *descB, int *incB );

void scalapack_pdtrsm( char *side, char *uplo, char *trans, char *diag,
         int *m, int *n, double *alpha,
         double *A, int *ia, int *ja, int *descA,
         double *B, int *ib, int *jb, int *descB );

void scalapack_pdsyrk( char *uplo, char *trans, 
         int *n, int *k, double *alpha,
         double *A, int *ia, int *ja, int *descA,
         double *beta,
         double *B, int *ib, int *jb, int *descB );

void scalapack_pdgeadd( char *trans, int *m, int *n,
           double *alpha,
           double *A, int *ia, int *ja, int *descA,
           double *beta,
           double *B, int *ib, int *jb, int *descB );

void scalapack_pdtradd( char *uplo, char *trans, int *m, int *n,
           double *alpha,
           double *A, int *ia, int *ja, int *descA,
           double *beta,
           double *B, int *ib, int *jb, int *descB );

void scalapack_pdlaswp( char *direc, char *rowcol,
             int *n, double *A, int *ia, int *ja, int *descA,
             int *k1, int *k2, int *ipiv );

void scalapack_dgebs2d( int *icontxt, char *scope, char *top,
           int *m, int *n, double *A, int *lda );

void scalapack_dgebr2d( int *icontxt, char *scope, char *top,
           int *m, int *n, double *A, int *lda,
           int *rsrc, int *csrc );

void scalapack_pdgemm( char *transA, char *transB, int *m, int *n, int *k,
       double *alpha,    double *A, int *ia, int *ja, int *descA,
                         double *B, int *ib, int *jb, int *descB,
       double *beta,     double *C, int *ic, int *jc, int *descC );


void scalapack_pdlapiv( char *direc, char *rowcol, char *pivroc,
     int *m, int *n, double *A, int *ia, int *ja, int *descA,
     int *ipiv, int *ip, int *jp, int *descip, int *iwork );

void scalapack_pdlaprnt(int *m, int *n, double *A, int *ia, int *ja, int *descA,
    int *irprnt, int *icprnt, char *cmatnm, int *nout, double *work,   int *dummy );


void scalapack_pdelget( char *scope, char *top, double *alpha, 
    double *A, int *ia, int *ja, int *descA );

}

#endif
