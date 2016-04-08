#ifndef OOC_OFFLOAD
#define OOC_OFFLOAD 1

#include <stdint.h>
#include <mkl.h>

#ifndef MIN
#define MIN(x,y)  (((x) < (y)) ? (x) : (y) )
#endif

#ifndef MAX
#define MAX(x,y)  (((x) > (y)) ? (x) : (y) )
#endif

#ifndef MOD
#define MOD(x,y)  ((x) % (y))
#endif

#ifndef IDX2F
#define IDX2F(i,j,lld)  (( (size_t)(i) + (size_t)((j)-1)*(size_t)(lld) ) - 1)
#endif

void offload_init(int *myrank, int *mydevice);
extern "C" void offload_init_(int *myrank, int *mydevice);
void offload_destroy();
extern "C" void offload_destroy_();

/* memory management */
intptr_t offload_Alloc(size_t size);
void offload_Free(void *p);

/* helper function */
void offload_dSetVector(int n, double *x, int incx, double *y, int incy);
void offload_dGetVector(int n, double *x, int incx, double *y, int incy);
void offload_dSetMatrix(int rows, int cols, double *A, int lda, double *B, int ldb);
void offload_dGetMatrix(int rows, int cols, double *A, int lda, double *B, int ldb);
void offload_Sync();

/* level-1 function */
void offload_dcopy(int n, const double *x, int incx, double *y, int incy);

/* level-3 function */
void offload_dgemm(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k,
                   const double *alpha, const double *a, const MKL_INT *lda, const double *b, const MKL_INT *ldb,
                   const double *beta, double *c, const MKL_INT *ldc);

void offload_dsyrk(const char *uplo, const char *trans, const MKL_INT *n, const MKL_INT *k,
                   const double *alpha, const double *a, const MKL_INT *lda, const double *beta,
                   double *c, const MKL_INT *ldc);

void offload_dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                   const MKL_INT *m, const MKL_INT *n, const double *alpha, const double *a, const MKL_INT *lda,
                   double *b, const MKL_INT *ldb);

/* matrix factorization */
void offload_dpotrf( const char* uplo, const MKL_INT* n, double* a, const MKL_INT* lda, 
                     MKL_INT* info );


#endif
