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
#include "ooclu.h"
void *MallocHost( size_t nbytes )
{

  void *p = 0;

#ifdef USE_CUDA_MALLOC_HOST
  cudaError_t cuda_status;

  cuda_status = cudaMallocHost( (void **) &p, nbytes );
  assert( cuda_status == cudaSuccess );
#else
  p = malloc( nbytes );
#endif


  return( p );
}


void FreeHost( void *p )
{
#ifdef USE_CUDA_MALLOC_HOST
  cudaError_t cuda_status;

  cuda_status = cudaFreeHost( p );
  assert( cuda_status == cudaSuccess );
#else
  free(p);
#endif

}
