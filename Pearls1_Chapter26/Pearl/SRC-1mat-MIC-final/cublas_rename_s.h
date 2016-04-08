#ifndef CUBLAS_RENAME_S_H
#define CUBLAS_RENAME_S_H 1


#ifdef USE_CUBLASV2

#define CUBLAS_SGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasSgemm( cublas_get_handle(), \
       transA,transB,m,n,k,&(alpha), \
         A,lda, \
         B,ldb, &(beta),C,ldc)

#define CUBLAS_STRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasStrsm( cublas_get_handle(), \
      side,uplo,trans,diag,m,n, &(alpha),A,ldA,B,ldB)

#define CUBLAS_SSYRK(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC) \
   cublasSsyrk(  cublas_get_handle(), \
        uplo,trans,n,k,&(alpha),A,ldA,&(beta),C,ldC)

#else

#define CUBLAS_SGEMM(transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
   cublasSgemm( \
       transA,transB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)

#define CUBLAS_STRSM(side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB) \
  cublasStrsm( \
      side,uplo,trans,diag,m,n,alpha,A,ldA,B,ldB)

#define CUBLAS_SSYRK(uplo,trans,n,k,alpha,A,ldA,beta,C,ldC) \
   cublasSsyrk(  \
        uplo,trans,n,k,alpha,A,ldA,beta,C,ldC)

#endif


#endif
