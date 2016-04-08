#include "ooclu.h"


#ifdef __cplusplus
extern "C" {
#endif

  void zgemm_( char *transA, char *transB, int *m, int *n, int *k,
      double *alpha,  double *A, int *ldA,  double *B, int *ldB,
      double *beta,   double *C, int *ldC );

  void ztrsm_( char *side, char *uplo, char *transA, char *diag,
      int *m, int *n,
      double *alpha,  double *A, int *ldA,
                      double *B, int *ldB );

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" 
#endif
cublasStatus cublasAlloc( int n, int elemSize, void ** ptr )
{
  const int use_calloc = FALSE;
  if (use_calloc) {
    *ptr = (void *) calloc( n, elemSize );
  }
  else {
    size_t nbytes;
    nbytes = n;
    nbytes *= elemSize;
    *ptr = (void *) malloc( nbytes );
  };

  return(  CUBLAS_STATUS_SUCCESS );
}


#ifdef __cplusplus
extern "C" 
#endif
cublasStatus cublasFree( void *ptr )
{
  free( ptr );
  return(  CUBLAS_STATUS_SUCCESS );

}



#ifdef __cplusplus
extern "C" 
#endif

cublasStatus cublasGetMatrix( int m, int n, int elemSize,
        void *A, int ldA, void *B, int ldB )
{
  /*
   * array A is the source, B is the destination
   */

  int i,j, ipA, ipB;
  char *pA = 0;
  char *pB = 0;
  char *src = 0;
  char *dest = 0;
  size_t nbytes; 

  pA = (char *) A;
  pB = (char *) B;
  nbytes = m;
  nbytes *= elemSize;

  for (j=1; j <= n; j++) {
    src = pA;
    dest = pB;
    memcpy( (void *) dest, (void *) src, nbytes );

    pA += ldA * elemSize;
    pB += ldB * elemSize;
  };
  return( CUBLAS_STATUS_SUCCESS );
}

#ifdef __cplusplus
extern "C" 
#endif
cublasStatus cublasSetMatrix( int m, int n, int elemSize,
        void *A, int ldA, void *B, int ldB )
{
  return( cublasGetMatrix(m,n,elemSize, A,ldA, B, ldB ) );
}


#ifdef __cplusplus
extern "C" 
#endif
cublasStatus cublasSetVector( int n, int elemSize,
       void *A, int inc1,  void *B, int inc2 )
{
  char *src = 0;
  char *dest = 0;
  int i;
  size_t nbytes;

  if ((inc1 == 1) && (inc2 == 1)) {
      nbytes = elemSize;
      nbytes *= n;
      src = (char *) A;
      dest = (char *)  B;
      memcpy( (void *) dest, (void *) src, nbytes );
      }
  else {
    src =  (char *) A;
    dest =  (char *) B;
    nbytes = elemSize;
    for( i=0; i < n; i++) {
        memcpy( (void *) dest, (void *) src, nbytes );
        src += inc1 * nbytes;
        dest += inc2 * nbytes;
    };
  };
  return (CUBLAS_STATUS_SUCCESS);
}


#ifdef __cplusplus
extern "C" 
#endif
void
cublasDgemm( char transA, char transB, int m, int n, int k,
    double alpha, double *A, int ldA,
                           double *B, int ldB,
    double beta,  double *C, int ldC )
{
  double zalpha_[REAL_PART+IMAG_PART+1];
  double zbeta_[REAL_PART+IMAG_PART+1];
  double *zalpha = &(zalpha_[0]);
  double *zbeta = &(zbeta_[0]);

  zalpha[REAL_PART] = creal(alpha);
  zalpha[IMAG_PART] = cimag(alpha);
  zbeta[REAL_PART] = creal(beta);
  zbeta[IMAG_PART] = cimag(beta);

  zgemm_( &transA, &transB, &m, &n, &k,
        zalpha, (double *) A, &ldA,  (double *) B, &ldB,
        zbeta,  (double *) C, &ldC );
}

#ifdef __cplusplus
extern "C" 
#endif
void cublasDtrsm( char side, char uplo, char trans, char diag,
    int m, int n,
    double alpha,  double *A, int ldA,
                            double *B, int ldB )
{
  double zalpha_[REAL_PART+IMAG_PART+1];
  double *zalpha = &(zalpha_[0]);

  zalpha[REAL_PART] = creal(alpha);
  zalpha[IMAG_PART] = cimag(alpha);

  ztrsm_( &side, &uplo, &trans, &diag,
      &m, &n,
      zalpha, (double *) A, &ldA, (double *) B, &ldB );

}

#ifdef __cplusplus
extern "C" 
#endif
void cublasInit()
{
  return;
}

#ifdef __cplusplus
extern "C" 
#endif
void cublasShutdown()
{
  return;
}


#ifdef __cplusplus
extern "C" 
#endif
double make_double(double x, double y )
{
  return x + _Complex_I * y;
}


cudaError_t cudaMallocHost( void **ptr, size_t nbytes )
{
  *ptr = malloc( nbytes );
  assert( ptr != 0 );
  return( cudaSuccess );
}

cudaError_t cudaFreeHost( void *ptr )
{
  free( ptr );
  return( cudaSuccess );
}
