/*
Copyright (c) The University of Tennessee.  All rights reserved.


$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution.

- Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or any other intellectual property rights of third parties.  The copyright holders disclaim any liability to any recipient for claims brought against recipient by any third party for infringement of that parties intellectual property rights.
*/
#ifdef USE_MIC
#include <mpi.h>
#include "profinit.h"
#include <offload.h>
#include <mkl_blas.h>
#include <malloc.h>
#include "ooc_offload.h"
#define MAXSTREAM 60 // not implemented
#define BUFFERSIZE 16*1024*1024/sizeof(double) // 16M. NOT OPTIMIZED

#define da(i,j)  (((double*)aptr) + IDX2F(i,j,lda))
#define db(i,j)  (((double*)bptr) + IDX2F(i,j,ldb))

static int MYRANK;
static int MYDEVICE;
static int MEMSIZE; // not implemented
static int *CURRENTSTREAM; // not implemented
static int STREAM[MAXSTREAM]; // not implemented
__attribute__((target(mic))) static double DBUFFER_[BUFFERSIZE];
__attribute__((target(mic))) static double *DBUFFER;

void offload_init(int *myrank, int *mydevice){
    int i;
    intptr_t ptr;
    MYRANK = *myrank;
    MYDEVICE = _Offload_number_of_devices();
    if(MYDEVICE){
        MYDEVICE = MYRANK % _Offload_number_of_devices();
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
    offload_touch((void*)ptr, size);
    return ptr;
}

void offload_touch(void* p, size_t size){
    intptr_t ptr = (intptr_t) p;
    #pragma offload target(mic:MYDEVICE)
    {
        char* C = (char*) ptr;
        double* D;
        size_t i, iend;
        double B[8];
        iend = size % 8;
        for(i = 0; i < iend; i++){
            B[i] = C[i];
        }
        D = (double*) (C+i);
        iend = (size / 8) % 8;
        for(i = 0 ; i < iend; i++){
            B[i] = D[i];
        }
        iend = size / 8;
        for( ; i < iend; i = i + 8){
            B[0] = D[i];
            B[1] = D[i+1];
            B[2] = D[i+2];
            B[3] = D[i+3];
            B[4] = D[i+4];
            B[5] = D[i+5];
            B[6] = D[i+6];
            B[7] = D[i+7];
        }
    }
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
    size_t start = 0;
    size_t end = start + BUFFERSIZE - 1;
           end = MIN(end, n - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    double *xstart;
    intptr_t yptr = (intptr_t)y;

//  printf("offload_dSetVector start\n");
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
//  printf("offload_dSetVector end\n");
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
    size_t start = 0;
    size_t end = start + BUFFERSIZE - 1;
           end = MIN(end, n - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    double *ystart;
    intptr_t xptr = (intptr_t)x;
    
//  printf("offload_dGetVector start\n");
    for(start = 0; start < n; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, n - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        ystart = y + start*incy;

        #pragma offload target(mic:MYDEVICE) out(DBUFFER:length(length) alloc_if(0) free_if(0)) in(xptr,incx,length,incy,incB)
//        #pragma offload target(mic:MYDEVICE) out(DBUFFER:length(length) alloc_if(0) free_if(0)) in(xptr,incx,length,incx,incB)  //original version, has two "incx", but it compiles well on Beacon
        {
            double *xstart = ((double*)xptr) + start*incx;
            dcopy(&length, xstart, &incx, DBUFFER, &incB);
        }
        dcopy(&length, DBUFFER, &incB, ystart, &incy);
    }
//  printf("offload_dGetVector end\n");
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
    size_t start = 0;
    size_t end = rows;
           end = end*cols - 1;
           end = MIN(end, start + BUFFERSIZE - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    int clength;
    int ia = 1; 
    int ja = 1;
    int ib = 1;
    int jb = 1;
    int ra, rb;
    int filled = 0;
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;

//  printf("offload_dSetMatrix start\n");
//  printf("rows %d cols %d a %zd lda %d b %zd ldb %d\n", rows, cols, aptr, lda, bptr, ldb);
//  printf("offload_dSetMatrix skip\n"); return;
    for(start = 0; start < rows*cols; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, rows*cols - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        if(lda == rows){
            dcopy(&length, da(ia,ja), &inca, DBUFFER, &incB);
            ra = length % rows;
            ia = ia + ra;
            if(ia > rows){
                ia = ia - rows;
                ja = ja + (length - ra)/rows + 1;
            }else{
                ja = ja + (length - ra)/rows;
            }
        }else{
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
        #pragma offload target(mic:MYDEVICE) in(DBUFFER:length(length) alloc_if(0) free_if(0))
        {
            int incb = 1;
            if(ldb == rows){
                dcopy(&length, DBUFFER, &incB, db(ib,jb), &incb);
                rb = length % rows;
                ib = ib + rb;
                if(ib > rows){
                    ib = ib - rows;
                    jb = jb + (length - rb)/rows + 1;
                }else{
                    jb = jb + (length - rb)/rows;
                }
            }else{
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
//  printf("offload_dSetMatrix end\n");
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
    size_t start = 0;
    size_t end = rows;
           end = end*cols - 1;
           end = MIN(end, start + BUFFERSIZE - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    int clength;
    int ia = 1; 
    int ja = 1;
    int ib = 1;
    int jb = 1;
    int ra, rb;
    int filled = 0;
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;

//  printf("offload_dGetMatrix start\n");
    for(start = 0; start < rows*cols; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, rows*cols - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        #pragma offload target(mic:MYDEVICE) out(DBUFFER:length(length) alloc_if(0) free_if(0))
        {
            int inca = 1;
            if(lda == rows){
                dcopy(&length, da(ia,ja), &inca, DBUFFER, &incB);
                ra = length % rows;
                ia = ia + ra;
                if(ia > rows){
                    ia = ia - rows;
                    ja = ja + (length - ra)/rows + 1;
                }else{
                    ja = ja + (length - ra)/rows;
                }
            }else{
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
        }
        if(ldb == rows){
            dcopy(&length, DBUFFER, &incB, db(ib,jb), &incb);
            rb = length % rows;
            ib = ib + rb;
            if(ib > rows){
                ib = ib - rows;
                jb = jb + (length - rb)/rows + 1;
            }else{
                jb = jb + (length - rb)/rows;
            }
        }else{
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
//  printf("offload_dGetMatrix end\n");
}

void offload_dtrSetMatrix(char uplo, int rows, int cols, double *a, int lda, double *b, int ldb){
/*
 * a is at host, b is at device, for trapezoidal matrices
 */
//  assert(rows >= 0); assert(cols >= 0); assert(ia >= 1); assert(ja >= 1); assert(lda >= 1);
//  assert(ib >= 1); assert(jb >= 1); assert(ldb >= 1);

    // copy a to DBUFFER, offload transfer in to DBUFFER, copy to y
    int is_lower = (uplo == 'L')||(uplo == 'l');
    int inca, incb;
    int incB = 1;
    size_t start = 0;
    size_t end = rows;
           end = end*cols - 1;
           end = MIN(end, start + BUFFERSIZE - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    int clength;
    size_t total;
    int ia = 1; 
    int ja = 1;
    int ib = 1;
    int jb = 1;
    int ra, rb;
    int filled = 0;
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;
        
    if(is_lower){
        cols = MIN(rows, cols);
        inca = 1; incb = 1;
        total = 2*rows - cols;
        total = total*cols + cols;
        total = total/2; 
    }else{
        rows = MIN(rows, cols);
        inca = lda; incb = ldb;
        total = 2*cols - rows;
        total = total*rows + rows;
        total = total/2;
    }
    for(start = 0; start < total; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, total - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        filled = 0;
        if(is_lower){
            while(filled < length){
                if(length - filled < rows - ia + 1){
                    clength = length - filled;
                    dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                    ia = ia + clength;
                }else{
                    clength = rows - ia + 1;
                    dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                    ja = ja + 1;
                    ia = ja;
                }
                filled = filled + clength;
            }
        }else{
            while(filled < length){
                if(length - filled < cols - ja + 1){
                    clength = length - filled;
                    dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                    ja = ja + clength;
                }else{
                    clength = cols - ja + 1;
                    dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                    ia = ia + 1;
                    ja = ia;
                }
                filled = filled + clength;
            }
        }
        #pragma offload target(mic:MYDEVICE) in(DBUFFER:length(length) alloc_if(0) free_if(0))
        {
            filled = 0;
            if(is_lower){
                while(filled < length){
                    if(length - filled < rows - ib + 1){
                        clength = length - filled;
                        dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                        ib = ib + clength;
                    }else{
                        clength = rows - ib + 1;
                        dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                        jb = jb + 1;
                        ib = jb;
                    }
                    filled = filled + clength;
                }
                while(filled < length){
                    if(length - filled < cols - ib + 1){
                        clength = length - filled;
                        dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                        jb = jb + clength;
                    }else{
                        clength = rows - ib + 1;
                        dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                        ib = ib + 1;
                        jb = ib;
                    }
                    filled = filled + clength;
                }
            }
        }
    }
}

void offload_dtrGetMatrix(char uplo, int rows, int cols, double *a, int lda, double *b, int ldb){
    // for trapezoidol matrices
    // copy a to DBUFFER, offload transfer out to DBUFFER, copy to y
    
    int is_lower = (uplo == 'L')||(uplo == 'l');
    int incB = 1;
    int inca, incb;
    size_t start = 0;
    size_t end = rows;
           end = end*cols - 1;
           end = MIN(end, start + BUFFERSIZE - 1);
    int length = MIN(end - start + 1, BUFFERSIZE);
    int clength;
    size_t total;
    int ia = 1; 
    int ja = 1;
    int ib = 1;
    int jb = 1;
    int filled = 0;
    intptr_t aptr = (intptr_t)a;
    intptr_t bptr = (intptr_t)b;

    if(is_lower){
        cols = MIN(rows, cols);
        inca = 1; incb = 1;
        total = 2*rows - cols;
        total = total*cols + cols;
        total = total/2; 
    }else{
        rows = MIN(rows, cols);
        inca = lda; incb = ldb;
        total = 2*cols - rows;
        total = total*rows + rows;
        total = total/2;
    }
    for(start = 0; start < total; start = end + 1){ 
        end = start + BUFFERSIZE - 1;
        end = MIN(end, total - 1);
    
        length = MIN(end - start + 1, BUFFERSIZE);
        #pragma offload target(mic:MYDEVICE) out(DBUFFER:length(length) alloc_if(0) free_if(0))
        {
            filled = 0;
            if(is_lower){
                while(filled < length){
                    if(length - filled < rows - ia + 1){
                        clength = length - filled;
                        dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                        ia = ia + clength;
                    }else{
                        clength = rows - ia + 1;
                        dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                        ja = ja + 1;
                        ia = ja;
                    }
                    filled = filled + clength;
                }
            }else{
                while(filled < length){
                    if(length - filled < cols - ja + 1){
                        clength = length - filled;
                        dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                        ja = ja + clength;
                    }else{
                        clength = cols - ja + 1;
                        dcopy(&clength, da(ia,ja), &inca, DBUFFER+filled, &incB);
                        ia = ia + 1;
                        ja = ia;
                    }
                    filled = filled + clength;
                }
            }
        }

        filled = 0;
        if(is_lower){
            while(filled < length){
                if(length - filled < rows - ib + 1){
                    clength = length - filled;
                    dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                    ib = ib + clength;
                }else{
                    clength = rows - ib + 1;
                    dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                    jb = jb + 1;
                    ib = jb;
                }
                filled = filled + clength;
            }
            while(filled < length){
                if(length - filled < cols - ib + 1){
                    clength = length - filled;
                    dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                    jb = jb + clength;
                }else{
                    clength = rows - ib + 1;
                    dcopy(&clength, DBUFFER+filled, &incB, db(ib,jb), &incb);
                    ib = ib + 1;
                    jb = ib;
                }
                filled = filled + clength;
            }
        }
    }
}

void offload_dcopy(int n, const double *x, int incx, double *y, int incy){
/*
 *  perform dcopy on the device. x,y pre-exist on the device
 */
    intptr_t xptr = (intptr_t)x;
    intptr_t yptr = (intptr_t)y;
//  printf("offload_dcopy start\n");
    #pragma offload target(mic:MYDEVICE)
    {
        dcopy(&n, (double*)xptr, &incx, (double*)yptr, &incy);
    }
//  printf("offload_dcopy end\n");
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
//  printf("offload_dgemm start\n");
    #pragma offload target(mic:MYDEVICE) in(transa,transb,m,n,k:length(1)) \
                                         in(alpha,lda,ldb,beta,ldc:length(1)) in(aptr,bptr,cptr)
    {
        dgemm(transa,transb,m,n,k,alpha,(double*)aptr,lda,(double*)bptr,ldb,beta,(double*)cptr,ldc); 
    }
//  printf("offload_dgemm end\n");
}

void offload_dsyrk(const char *uplo, const char *trans, const MKL_INT *n, const MKL_INT *k,
                   const double *alpha, const double *a, const MKL_INT *lda, const double *beta,
                   double *c, const MKL_INT *ldc){
/*
 * perform dsyrk on the device. a,c pre-exist on the device
 */
    intptr_t aptr = (intptr_t)a;
    intptr_t cptr = (intptr_t)c;
//  printf("offload_dsyrk start\n");
    #pragma offload target(mic:MYDEVICE) in(uplo,trans,n,k:length(1)) \
                                         in(alpha,lda,beta,ldc:length(1)) in(aptr,cptr)
    {
        dsyrk(uplo,trans,n,k,alpha,(double*)aptr,lda,beta,(double*)cptr,ldc);
    }
//  printf("offload_dsyrk end\n");
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
