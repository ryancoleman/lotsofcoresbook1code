#ifndef CUBLAS_RENAME_H
#define CUBLAS_RENAME_H 1

extern void cublas_sync_stream();
#ifdef USE_CUBLASV2
#include "cublas_v2.h"

extern cublasHandle_t cublas_get_handle();
extern cudaStream_t cublas_get_stream();
#endif


#include "cublas_rename_z.h"
#include "cublas_rename_c.h"
#include "cublas_rename_d.h"
#include "cublas_rename_s.h"





#endif
