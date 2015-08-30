#ifndef CUBLAS_RENAME_Z_H
#define CUBLAS_RENAME_Z_H 1


#ifdef USE_CUBLASV2

#define CUBLAS_ZGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasZgemm( cublas_get_handle(), \
       transA,transB,m,n,k,&(alpha), \
         A,lda, \
         B,ldb, &(beta),C,ldc)


#define CUBLAS_ZTRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasZtrsm( cublas_get_handle(), \
      side,uplo,trans,diag,m,n, &(alpha),A,ldA,B,ldB)

#else

#define CUBLAS_ZGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasZgemm( \
       transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)


#define CUBLAS_ZTRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasZtrsm( \
      side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB)


#endif
#endif
