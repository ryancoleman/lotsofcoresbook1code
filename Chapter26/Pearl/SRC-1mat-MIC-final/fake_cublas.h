#ifndef FAKE_CUBLAS_H
#define FAKE_CUBLAS_H 1

#include <assert.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

#ifndef REAL_PART
#define REAL_PART 0
#endif

#ifndef IMAG_PART
#define IMAG_PART 1
#endif

typedef double complex cuDoubleComplex;

typedef unsigned int cublasStatus;

typedef unsigned int cudaError_t;

#define cudaSuccess  0

#define CUBLAS_STATUS_SUCCESS           0x00000000
#define CUBLAS_STATUS_NOT_INITIALIZED   0x00000001
#define CUBLAS_STATUS_ALLOC_FAILED      0x00000003
#define CUBLAS_STATUS_INVALID_VALUE     0x00000007
#define CUBLAS_STATUS_ARCH_MISMATCH     0x00000008
#define CUBLAS_STATUS_MAPPING_ERROR     0x0000000B
#define CUBLAS_STATUS_EXECUTION_FAILED  0x0000000D
#define CUBLAS_STATUS_INTERNAL_ERROR    0x0000000E



extern "C" {
  cublasStatus cublasAlloc( int n, int elemSize, void **ptr );

  cublasStatus cublasFree( void *ptr );

  cublasStatus cublasGetMatrix(int m, int n, int elemSize,
               void *Asrc, int lda,  void *Bdest, int ldb );

  cublasStatus cublasSetMatrix(int m, int n, int elemSize,
               void *Asrc, int lda,  void *Bdest, int ldb );

  cublasStatus cublasSetVector( int n, int elemSize,
               void *Asrc, int inc1, void *Bdest, int inc2 );


  void cublasZgemm( char transA, char transB, int m, int n, int k,
      cuDoubleComplex alpha, cuDoubleComplex *A, int ldA,
                             cuDoubleComplex *B, int ldB,
      cuDoubleComplex beta,  cuDoubleComplex *C, int ldC );

  void cublasInit();
  void cublasShutdown();

  cuDoubleComplex make_cuDoubleComplex(double x, double y );

  void cudaThreadSynchronize();

  void cublasZtrsm( char side, char uplo, char transA, char diag,
      int m, int n,
      cuDoubleComplex alpha,   cuDoubleComplex *A, int ldA,
                               cuDoubleComplex *B, int ldB );

  cudaError_t cudaMallocHost( void **ptr, size_t size );
  cudaError_t cudaFreeHost( void *ptr);
}


#endif
