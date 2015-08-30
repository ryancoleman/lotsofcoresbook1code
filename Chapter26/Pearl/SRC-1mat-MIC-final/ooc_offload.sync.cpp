#ifdef USE_MIC

#include <offload.h>
#include <mkl_blas.h>
#include <malloc.h>
#include "ooc_offload.h"
#define MAXSTREAM 60 // not implemented
#define BUFFERSIZE 16*1024*1024/sizeof(double) // 16M. NOT OPTIMIZED

#define da(i,j)  (((double*)aptr) + IDX2F(i,j,lda))
#define db(i,j)  (((double*)bptr) + IDX2F(i,j,ldb))

static int MYDEVICE;
static int MEMSIZE; // not implemented
static int *CURRENTSTREAM; // not implemented
static int STREAM[MAXSTREAM]; // not implemented
__attribute__((target(mic))) static double DBUFFER_[BUFFERSIZE];
__attribute__((target(mic))) static double *DBUFFER;

void offload_init(int *myrank, int *mydevice){
    int i;
    intptr_t ptr;
    MYDEVICE = _Offload_number_of_devices();
    if(MYDEVICE){
        MYDEVICE = *myrank % _Offload_number_of_devices();
    }
    *mydevice = MYDEVICE;
    for(i = 0; i < MAXSTREAM; i++){
        STREAM[i] = 0;
    }
//  CURRENTSTREAM = STREAM;
    DBUFFER = DBUFFER_;
    #pragma offload target(mic:MYDEVICE) nocopy(DBUFFER_:alloc_if(1) free_if(0)) out(ptr)
    {
        ptr = (intptr_t)DBUFFER_;
    }
//  printf("buffer %zd %zd\n", ptr, ptr + BUFFERSIZE - 1);

}

extern "C" void offload_init_(int *myrank, int *mydevice){
    offload_init(myrank, mydevice);
}

void offload_destroy(){
    #pragma offload_transfer target(mic:MYDEVICE) nocopy(DBUFFER:alloc_if(0) free_if(1))
}

extern "C" void offload_destroy_(){
    offload_destroy();    
}



intptr_t offload_Alloc(size_t size){
    intptr_t ptr;
    #pragma offload target(mic:MYDEVICE) out(ptr)
    {
        ptr = (intptr_t) memalign(64, size);
    }
//  printf("%zd %zd\n", ptr, ptr+size-1);
    return ptr;
}

void offload_Free(void* p){
    intptr_t ptr = (intptr_t)p;
    #pragma offload target(mic:MYDEVICE) in(ptr)
    {
        free((void*)ptr);
    }
}

void offload_dSetVector(int n, double *x, int incx, double *y, int incy){
/*
 *  copy x at host to y at device
 *  incx is the index increment of x, incy is the index increment of y
 *  n elements are copied
 *  algorithm works for negative values of incx and incy, but gives undefined behavior
 */

//  assert(n >= 0);
    // copy x to DBUFFER, offload transfer in to DBUFFER, copy to y
    int incB = 1;
    int start = 0;
    int end = start + BUFFERSIZE - 1;
        end = MIN(end, n - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    double *xstart;
    intptr_t yptr = (intptr_t)y;

    for(start = 0; start < n; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, n - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        xstart = x + start*incx;

        dcopy(&length, xstart, &incx, DBUFFER, &incB);
        #pragma offload target(mic:MYDEVICE) in(DBUFFER:length(length) alloc_if(0) free_if(0)) \
                                             in(yptr,incy,length,incB)
        {
            double *ystart = ((double*)yptr) + start*incy;
            dcopy(&length, DBUFFER, &incB, ystart, &incy);
        }
    }
}

void offload_dGetVector(int n, double *x, int incx, double *y, int incy){
/*
 *  copy x at device to y at host
 *  incx is the index increment of x, incy is the index increment of y
 *  n elements are copied
 *  algorithm works for negative values of incx and incy, but gives undefined behavior
 */
//  assert(n >= 0);
  
    // copy x to DBUFFER, offload transfer out to DBUFFER, copy to y
    int incB = 1;
    int start = 0;
    int end = start + BUFFERSIZE - 1;
        end = MIN(end, n - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    double *ystart;
    intptr_t xptr = (intptr_t)x;
    
    for(start = 0; start < n; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, n - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        ystart = y + start*incy;

        #pragma offload target(mic:MYDEVICE) out(DBUFFER:length(length) alloc_if(0) free_if(0)) \
                                             in(xptr,incx,length,incx,incB)
        {
            double *xstart = ((double*)xptr) + start*incx;
            dcopy(&length, xstart, &incx, DBUFFER, &incB);
        }
        dcopy(&length, DBUFFER, &incB, ystart, &incy);
    }
}

void offload_dSetMatrix(int rows, int cols, double *a, int lda, double *b, int ldb){
/*
 * a is at host, b is at device
 */
//  assert(rows >= 0); assert(cols >= 0); assert(ia >= 1); assert(ja >= 1); assert(lda >= 1);
//  assert(ib >= 1); assert(jb >= 1); assert(ldb >= 1);

    // copy a to DBUFFER, offload transfer in to DBUFFER, copy to y
    int inca = 1;
    int incB = 1;
    int start = 0;
    size_t end = rows;
           end = end*cols - 1;
           end = MIN(end, start + BUFFERSIZE - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    int clength;
    int ia = 1; 
    int ja = 1;
    int ib = 1;
    int jb = 1;
    int filled = 0;
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;
        
    for(start = 0; start < rows*cols; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, rows*cols - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        filled = 0;
        while(filled < length){
            if(length - filled < rows - ia + 1){
                clength = length - filled;
                dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                ia = ia + clength;
            }else{
                clength = rows - ia + 1;
                dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                ia = 1;
                ja = ja + 1;
            }
            filled = filled + clength;
        }
        #pragma offload target(mic:MYDEVICE) in(DBUFFER:length(length) alloc_if(0) free_if(0))
        {
            int incb = 1;
            filled = 0;
            while(filled < length){
                if(length - filled < rows - ib + 1){
                    clength = length - filled;
                    dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                    ib = ib + clength;
                }else{
                    clength = rows - ib + 1;
                    dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                    ib = 1;
                    jb = jb + 1;
                }
                filled = filled + clength;
            }
        }
    }
}

void offload_dGetMatrix(int rows, int cols, double *a, int lda, double *b, int ldb){
/*
 * a is at device, b is at host
 */
//  assert(rows >= 0); assert(cols >= 0); assert(ia >= 1); assert(ja >= 1); assert(lda >= 1);
//  assert(ib >= 1); assert(jb >= 1); assert(ldb >= 1);

    // copy a to DBUFFER, offload transfer out to DBUFFER, copy to y
    int incB = 1;
    int incb = 1;
    int start = 0;
    size_t end = rows;
           end = end*cols - 1;
           end = MIN(end, start + BUFFERSIZE - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    int clength;
    int ia = 1; 
    int ja = 1;
    int ib = 1;
    int jb = 1;
    int filled = 0;
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;
        
    for(start = 0; start < rows*cols; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, rows*cols - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        #pragma offload target(mic:MYDEVICE) out(DBUFFER:length(length) alloc_if(0) free_if(0))
        {
            int inca = 1;
            filled = 0;
            while(filled < length){
                if(length - filled < rows - ia + 1){
                    clength = length - filled;
                    dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                    ia = ia + clength;
                }else{
                    clength = rows - ia + 1;
                    dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                    ia = 1;
                    ja = ja + 1;
                }
                filled = filled + clength;
            }
        }
        filled = 0;
        while(filled < length){
            if(length - filled < rows - ib + 1){
                clength = length - filled;
                dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                ib = ib + clength;
            }else{
                clength = rows - ib + 1;
                dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                ib = 1;
                jb = jb + 1;
            }
            filled = filled + clength;
        }
    }
}

void offload_dcopy(int n, const double *x, int incx, double *y, int incy){
/*
 *  perform dcopy on the device. x,y pre-exist on the device
 */
    intptr_t xptr = (intptr_t)x;
    intptr_t yptr = (intptr_t)y;
    #pragma offload target(mic:MYDEVICE)
    {
        dcopy(&n, (double*)xptr, &incx, (double*)yptr, &incy);
    }
}

void offload_dgemm(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k,
                   const double *alpha, const double *a, const MKL_INT *lda, const double *b, const MKL_INT *ldb,
                   const double *beta, double *c, const MKL_INT *ldc){
/*
 * perform dgemm on the device. a,b,c pre-exist on the device
 */
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;
    intptr_t cptr = (intptr_t)c;
    #pragma offload target(mic:MYDEVICE) in(transa,transb,m,n,k:length(1)) \
                                         in(alpha,lda,ldb,beta,ldc:length(1)) 
    {
        dgemm(transa,transb,m,n,k,alpha,(double*)aptr,lda,(double*)bptr,ldb,beta,(double*)cptr,ldc); 
    }
}

void offload_dsyrk(const char *uplo, const char *trans, const MKL_INT *n, const MKL_INT *k,
                   const double *alpha, const double *a, const MKL_INT *lda, const double *beta,
                   double *c, const MKL_INT *ldc){
/*
 * perform dsyrk on the device. a,c pre-exist on the device
 */
    intptr_t aptr = (intptr_t)a;
    intptr_t cptr = (intptr_t)c;
    #pragma offload target(mic:MYDEVICE) in(uplo,trans,n,k:length(1)) \
                                         in(alpha,lda,beta,ldc:length(1)) in(aptr,cptr)
    {
        dsyrk(uplo,trans,n,k,alpha,(double*)aptr,lda,beta,(double*)cptr,ldc);
    }
}

void offload_dtrsm(const char *side, const char *uplo, const char *transa, const char *diag,
                   const MKL_INT *m, const MKL_INT *n, const double *alpha, const double *a, const MKL_INT *lda,
                   double *b, const MKL_INT *ldb){
/*
 * perform dtrsm on the device. a,b pre-exist on the device
 */
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;
    #pragma offload target(mic:MYDEVICE) in(side,uplo,transa,diag,m,n,alpha,lda,ldb:length(1)) 
    {
        dtrsm(side,uplo,transa,diag,m,n,alpha,(double*)aptr,lda,(double*)bptr,ldb);
    }
}

void offload_dpotrf( const char* uplo, const MKL_INT* n, double* a, const MKL_INT* lda, 
                     MKL_INT* info ){
/*
 * perform dpotrf on the device. a pre-exists on the device
 */
    intptr_t aptr = (intptr_t)a;
    #pragma offload target(mic:MYDEVICE) in(uplo,n,lda:length(1)) out(info:length(1))
    {
        dpotrf(uplo,n,(double*)aptr,lda,info);
    }
}

#endif
