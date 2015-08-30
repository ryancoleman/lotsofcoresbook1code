#ifndef CUBLAS_RENAME_C_H
#define CUBLAS_RENAME_C_H 1

#ifdef USE_CUBLASV2

#define CUBLAS_CGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasCgemm( cublas_get_handle(), \
       transA,transB,m,n,k,&(alpha), \
         A,lda, \
         B,ldb, &(beta),C,ldc)


#define CUBLAS_CTRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasCtrsm( cublas_get_handle(), \
      side,uplo,trans,diag,m,n, &(alpha),A,ldA,B,ldB)

#else

#define CUBLAS_CGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasCgemm( \
       transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)


#define CUBLAS_CTRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasCtrsm( \
      side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB)


#endif

#endif
