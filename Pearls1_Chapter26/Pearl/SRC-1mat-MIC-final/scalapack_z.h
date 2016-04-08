#ifndef SCALAPACK_Z_H
#define SCALAPACK_Z_H 1



#if defined(ADD_)
#define scalapack_pzgemm  pzgemm_
#define scalapack_pzgetrf  pzgetrf_
#define scalapack_pzswap  pzswap_
#define scalapack_pzgeadd pzgeadd_
#define scalapack_pzgemm  pzgemm_
#define scalapack_pztrsm  pztrsm_
#define scalapack_pzlaswp pzlaswp_
#define scalapack_zgebr2d zgebr2d_
#define scalapack_zgebs2d zgebs2d_
#define scalapack_pzlaprnt pzlaprnt_
#define scalapack_pzelget pzelget_
#define scalapack_pzlapiv pzlapiv_

#else
#define scalapack_pzgemm  pzgemm
#define scalapack_pzgetrf  pzgetrf
#define scalapack_pzswap  pzswap
#define scalapack_pzgeadd pzgeadd
#define scalapack_pzgemm  pzgemm
#define scalapack_pztrsm  pztrsm
#define scalapack_pzlaswp pzlaswp
#define scalapack_zgebr2d zgebr2d
#define scalapack_zgebs2d zgebs2d
#define scalapack_pzlaprnt pzlaprnt
#define scalapack_pzelget pzelget
#define scalapack_pzlapiv pzlapiv
#endif





extern "C" {

void scalapack_pzgetrf( int *m, int *n, 
          double *A, int *ia, int *ja, int *descA, int *ipiv, int *info );

void scalapack_pzswap( int *n, 
  double *A, int *ia, int *ja, int *descA, int *incA,
  double *B, int *ib, int *jb, int *descB, int *incB );

void scalapack_pztrsm( char *side, char *uplo, char *trans, char *diag,
         int *m, int *n, double *alpha,
         double *A, int *ia, int *ja, int *descA,
         double *B, int *ib, int *jb, int *descB );

void scalapack_pzgeadd( char *trans, int *m, int *n,
           double *alpha,
           double *A, int *ia, int *ja, int *descA,
           double *beta,
           double *B, int *ib, int *jb, int *descB );

void scalapack_pzlaswp( char *direc, char *rowcol,
             int *n, double *A, int *ia, int *ja, int *descA,
             int *k1, int *k2, int *ipiv );

void scalapack_zgebs2d( int *icontxt, char *scope, char *top,
           int *m, int *n, double *A, int *lda );

void scalapack_zgebr2d( int *icontxt, char *scope, char *top,
           int *m, int *n, double *A, int *lda,
           int *rsrc, int *csrc );

void scalapack_pzgemm( char *transA, char *transB, int *m, int *n, int *k,
       double *alpha,    double *A, int *ia, int *ja, int *descA,
                         double *B, int *ib, int *jb, int *descB,
       double *beta,     double *C, int *ic, int *jc, int *descC );

void scalapack_pzlapiv( char *direc, char *rowcol, char *pivroc,
     int *m, int *n, double *A, int *ia, int *ja, int *descA,
     int *ipiv, int *ip, int *jp, int *descip, int *iwork );

void scalapack_pzlaprnt(int *m, int *n, double *A, int *ia, int *ja, int *descA,
    int *irprnt, int *icprnt, char *cmatnm, int *nout, double *work,   int *dummy );


void scalapack_pzelget( char *scope, char *top, double *alpha, 
    double *A, int *ia, int *ja, int *descA );

}

#endif
