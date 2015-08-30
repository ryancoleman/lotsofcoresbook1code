//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, FLORIAN WENDE, ZIB
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the ZIB nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FLORIAN WENDE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
//////////////////////////////////////////////////////////////////////////////////


#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <omp.h>
#include <sched.h>

using namespace std;

#if !defined GROUP_SIZE
#define GROUP_SIZE ( 4 )
#endif
//#define BATCH_OUTPUT

int main( int argc, char **argv ) {

  MPI::Init( argc, argv );

  uint32_t
    groupId = MPI::COMM_WORLD.Get_rank(),
    numGroups = MPI::COMM_WORLD.Get_size(),
    *groupCoreId = (uint32_t *) new uint32_t[GROUP_SIZE*256];
    
#pragma omp target device( 0 ) \
  map( from : groupCoreId[0:GROUP_SIZE*256] )
  {

#pragma omp parallel num_threads( GROUP_SIZE )
    {

      uint32_t
	ompId = omp_get_thread_num();
      cpu_set_t
	cpuMask;
      
      sched_getaffinity( 0, sizeof( cpu_set_t ), &cpuMask );
      
      for( uint32_t i=0; i<256; i++ ) {
	
	if( CPU_ISSET( i, &cpuMask ) )
	  groupCoreId[ompId*256+i] = i;
	else
	  groupCoreId[ompId*256+i] = 0xFFFFFFFF;
	
      }
      
    }
    
  }

  if( groupId == 0 ) {

    for( uint32_t g=0; g<numGroups; g++ ) {

      for( uint32_t i=0; i<GROUP_SIZE; i++ ) {

	// print mask
#if !defined BATCH_OUTPUT
	cout << "thread(hostId=" << g << ",phiId=" << i << "):";
#endif
	for( uint32_t j=0; j<256; j++ ) {
	  if( groupCoreId[i*256+j] != 0xFFFFFFFF ) {
#if defined BATCH_OUTPUT
	    cout << "# group\tthread\tcoreid " << endl << g << "\t" << i << "\t" << j << endl;
#else
	    cout << " " << j;
#endif	
	  }
	}

#if !defined BATCH_OUTPUT
	cout << endl;
#endif
	
      }

#if defined BATCH_OUTPUT
      cout << endl;
      cout << endl;
#endif

      if( g == ( numGroups-1 ) )
	break;
      
      // collect masks from other ranks
      MPI::COMM_WORLD.Recv( groupCoreId, GROUP_SIZE*256, MPI::UNSIGNED, g+1, 1 );

    }
    
  } else {

    MPI::COMM_WORLD.Send( groupCoreId, GROUP_SIZE*256, MPI::UNSIGNED, 0, 1 );

  }
  
  delete [] groupCoreId;
  groupCoreId = NULL;

  MPI::Finalize();

  exit( 0 );

}
