#ifndef USE_MIC

#include "ooclu.h"
#include <cuda_runtime_api.h>

#include <assert.h>
#include <mpi.h>

#ifdef USE_CUBLASV2
#include "cublas_v2.h"
static cublasHandle_t default_handle;
static  cudaStream_t default_stream;
cublasHandle_t  cublas_get_handle() { return default_handle; };
cudaStream_t cublas_get_stream() { return default_stream; };
void cublas_sync_stream() {
	cudaError_t istatus = cudaStreamSynchronize( cublas_get_stream() );
	assert( istatus == cudaSuccess );
}
#else
void cublas_sync_stream() {
	/*
	 * do nothing
	 */
}
#endif



#ifdef __cplusplus
extern "C"
#endif
void cublasinit_(int *IAM, int *NPROCS)
{

      /*
       * select the GPU device
       */
		int number_of_gpus = 0;
		cudaError_t istat;

		istat = cudaGetDeviceCount( &number_of_gpus );
		assert( istat == cudaSuccess );
		assert( number_of_gpus >= 1 );


		int my_device = (*IAM % number_of_gpus);
		istat = cudaSetDevice(my_device);
		assert( istat == cudaSuccess );

                // Get the name of the proccessor
                char processor_name[MPI_MAX_PROCESSOR_NAME]; 
                int name_len;
                MPI_Get_processor_name(processor_name, &name_len); 

		printf("I am from processor %s MPI task%d of %d and my device is %d\n",processor_name,*IAM,*NPROCS,my_device);fflush(stdout);

#ifdef USE_CUBLASV2
	{
		cublasStatus_t cuStatus;
		cudaError_t ierr;



		/*
		   the GPU is set, now call cublasCreate()
		 */
		cuStatus = cublasCreate( &default_handle);
		assert( cuStatus == CUBLAS_STATUS_SUCCESS );

		ierr = cudaStreamCreate( &default_stream );
		assert(ierr == cudaSuccess);

		cuStatus = cublasSetStream(default_handle, default_stream);
		assert( cuStatus == CUBLAS_STATUS_SUCCESS );
	}
#else
	cublasInit();
#endif

}


#ifdef __cplusplus
extern "C"
#endif
void cublasshutdown_()
{

#ifdef USE_CUBLASV2
	{
		cudaError_t ierr = cudaStreamDestroy( default_stream );
		assert( ierr == cudaSuccess );

		cublasStatus_t cuStatus = cublasDestroy( default_handle );
		assert( cuStatus == CUBLAS_STATUS_SUCCESS );
	}
#else
	cublasShutdown();
#endif

}
#endif
