#ifndef OOCLU_H
#define OOCLU_H 1

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef USE_MAGMA
    #include "magma.h"
    #include "magma_lapack.h"
/*
 * note magma may use the old cublas v1 interface
 */
#endif


#ifdef USE_MIC
            #include "ooc_offload.h"
#else
            #include "cublas_rename.h"
    #ifdef USE_CUBLASV2
            #include <cuda_runtime.h>
            #include <cuda_runtime_api.h>
            #include "cublas_v2.h"
            #define cublasStatus cublasStatus_t
    #else
        #ifdef USE_FAKE_CUBLAS
            #include "fake_cublas.h"
            #include "cublasOps.h"
        #else
            #include <cuda_runtime_api.h>
            #include <cublas.h>
            #include "cublasOps.h"
        #endif
    #endif
#endif


#include "scalapack.h"

//#include "PBblacs.h" //use when under mvapich2 approach
#ifndef USE_MIC
    #ifdef USE_CUBLASV2

        #define CUBLAS_MALLOC(dY,isize,elemSize) { \
            size_t nbytes; \
            cudaError_t ierr; \
            nbytes = ((size_t) isize) * ( (size_t) elemSize ); \
            assert( dY == 0 ); \
            ierr = cudaMalloc( (void **) &(dY), nbytes ); \
            assert( ierr == cudaSuccess ); \
        }

        #define CUBLAS_FREE(dY) { \
            cudaError_t ierr; \
            if (dY != 0) { \
                ierr = cudaFree( (void *) dY ); \
                assert( ierr == cudaSuccess ); \
            }; \
            dY = 0; \
        }

    #else

        #define CUBLAS_MALLOC(dY,isize,elemSize) { \
            cublasStatus cu_status; \
            assert( dY == 0 ); \
            cu_status  = cublasAlloc( isize, elemSize, (void **) &(dY) ); \
            CHKERR(cu_status); \
        }

        #define CUBLAS_FREE(dY) { \
            cublasStatus cu_status; \
            if (dY != 0) { \
            cu_status  = cublasFree(dY); \
            CHKERR(cu_status); \
            }; \
            dY = 0; \
        }

    #endif
#endif

#ifdef USE_PROFILE
    #include "profinit.h"
    #define cudaThreadSynchronize() 

    #define PROFSTART(name) { char buf[] = name; profstart(buf); }
    #define PROFEND(name) { char buf[] = name; profend(buf); }
    #define PROFINIT() { profinit(); }
    #define PROFSTAT() { profstat(); }
#else
    #define PROFSTART(name) {}
    #define PROFEND(name) {}
    #define PROFINIT() {}
    #define PROFSTAT() {}
#endif

#ifndef TRUE
#define TRUE (1 == 1)
#endif

#ifndef FALSE
#define FALSE (1 == 0)
#endif


#ifndef IDX1F
#define IDX1F(i)  ((i)-1)
#endif

#ifndef IDX2F
#define IDX2F(i,j,lld)  (( (size_t)(i) + (size_t)((j)-1)*(size_t)(lld) ) - 1)
#endif

#ifndef MIN
#define MIN(x,y)  (((x) < (y)) ? (x) : (y) )
#endif

#ifndef MAX
#define MAX(x,y)  (((x) > (y)) ? (x) : (y) )
#endif

#ifndef MOD
#define MOD(x,y)  ((x) % (y))
#endif

#ifndef USE_MIC
    #ifndef PRINT_ERR
        #define PRINT_ERR(cu_status) { \
                switch(cu_status)  \
                { \
                case CUBLAS_STATUS_NOT_INITIALIZED: \
                  { printf("CUBLAS_STATUS_NOT_INITIALIZED\n");  break; } \
                case CUBLAS_STATUS_ALLOC_FAILED:                       \
                  { printf("CUBLAS_STATUS_ALLOC_FAILED\n"); break; }   \
                case CUBLAS_STATUS_INVALID_VALUE: \
                  { printf("CUBLAS_STATUS_INVALID_VALUE\n");  break; } \
                case CUBLAS_STATUS_MAPPING_ERROR: \
                  { printf("CUBLAS_STATUS_MAPPING_ERROR\n");  break; } \
                case CUBLAS_STATUS_EXECUTION_FAILED: \
                  { printf("CUBLAS_STATUS_EXECUTION_FAILED\n");  break; } \
                case CUBLAS_STATUS_INTERNAL_ERROR: \
                  { printf("CUBLAS_STATUS_INTERNAL_ERROR\n");  break; } \
                default: \
                  { printf("unknown error\n"); } \
                }; }
    #endif

    #ifndef CHKERR
        #define CHKERR(cu_status)  { \
            if (cu_status != CUBLAS_STATUS_SUCCESS) {  \
                PRINT_ERR(cu_status); \
                }; \
            assert( cu_status == CUBLAS_STATUS_SUCCESS );  \
        }

    #endif
#endif






#ifdef __cplusplus
extern "C" {
#endif

int Cindxl2g( int indxloc, int nb, int iproc, int isrcproc, int nprocs );

int Cindxg2p(int indxglob, int nb, int iproc, int isrcproc, int nprocs );


void setup_desc( int m, int n, int ia, int ja, int *descA,  
                  int *isize, int *descB );

void Cdescinit( int *desc, int m, int n, int mb, int nb,
                int irsrc, int icsrc, int ictxt, int lld, int *info);

void Cdescset( int *desc, int m, int n, int mb, int nb,
                int irsrc, int icsrc, int ictxt, int lld );

int Cnumroc( int n, int nb, int iproc, int isrcproc, int nprocs );

int Cnumroc2( int ia, int n, int nb, int iproc, int isrcproc, int nprocs );

void Cinfog1l( int gindx, int nb, int nprocs, int myroc, int isrcproc,
               int *lindx, int *rocsrc );

void local_extent( int m, int n, int ia, int ja, int *descA,
              int *msize, int *nsize,
              int *lrA1, int *lcA1,  int *lrA2, int *lcA2 );


void Cinfog2l( int grindx, int gcindx, int *desc, int nprow, int npcol,
               int myrow, int mycol, 
               int *lrindx, int *lcindx, int *rsrc, int *csrc );


void Cpilaprnt( int m, int n,  int *A, int ia, int ja, int *descA, char *cmatnm );

int Ciafirst(int ia,int mb,int myprow,int rsrc,int nprow);

int Cialast(int ia,int mb,int myprow,int rsrc,int nprow);


void *MallocHost( size_t nbytes );
void FreeHost( void *ptr );

#ifdef __cplusplus
}
#endif

#ifdef USE_MIC
    #include "ooclu_d.h"
#else
    #include "ooclu_z.h"
    #include "ooclu_c.h"
    #include "ooclu_d.h"
    #include "ooclu_s.h"
#endif

#endif
