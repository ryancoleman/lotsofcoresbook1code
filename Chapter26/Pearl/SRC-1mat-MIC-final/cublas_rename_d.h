#ifndef CUBLAS_RENAME_D_H
#define CUBLAS_RENAME_D_H 1


#ifdef USE_CUBLASV2

#define CUBLAS_DGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasDgemm( cublas_get_handle(), \
       transA,transB,m,n,k,&(alpha), \
         A,lda, \
         B,ldb, &(beta),C,ldc)


#define CUBLAS_DTRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasDtrsm( cublas_get_handle(), \
      side,uplo,trans,diag,m,n, &(alpha),A,ldA,B,ldB)

#define CUBLAS_DSYRK(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC) \
   cublasDsyrk(  cublas_get_handle(), \
        uplo,trans,n,k,&(alpha),A,ldA,&(beta),C,ldC)

#else

#define CUBLAS_DGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasDgemm( \
       transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)


#define CUBLAS_DTRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasDtrsm( \
      side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB)

#define CUBLAS_DSYRK(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC) \
   cublasDsyrk(  \
        uplo,trans,n,k,alpha,A,ldA,beta,C,ldC)

#endif


#endif
