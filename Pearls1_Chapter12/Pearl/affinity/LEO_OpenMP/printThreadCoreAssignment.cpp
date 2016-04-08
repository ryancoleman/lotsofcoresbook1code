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

#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <omp.h>
#include <sched.h>

using namespace std;

#if !defined GROUPS
#define GROUPS ( 4 )
#endif

#if !defined GROUP_SIZE
#define GROUP_SIZE ( 4 )
#endif

//#define BATCH_OUTPUT

int main() {

  uint32_t
    coreId[GROUPS*GROUP_SIZE*256];

#pragma omp parallel num_threads( GROUPS )
  {
    
    uint32_t
      groupId = omp_get_thread_num();
    uint32_t
      *groupCoreId = (uint32_t *) new uint32_t[GROUP_SIZE*256];
    
    //#pragma omp critical // fix!
    {
    
#pragma offload target( mic : 0 )
      {

#pragma omp parallel num_threads( GROUP_SIZE )
	{
	  ; // just create threads
	}

      }

    }
    
#pragma offload target( mic : 0 ) \
  out( groupCoreId : length( GROUP_SIZE*256 ) )
    {

#pragma omp parallel num_threads( GROUP_SIZE )
      {

	uint32_t
	  ompId = omp_get_thread_num();
	cpu_set_t
	  cpuMask;
	
	sched_getaffinity( 0, sizeof( cpu_set_t ), &cpuMask );

	for( uint32_t i=0; i<256; i++ ) {

	  if( CPU_ISSET( i, &cpuMask ) ) {
	    groupCoreId[ompId*256+i] = i;
	  } else {
	    groupCoreId[ompId*256+i] = 0xFFFFFFFF;
	  }

	}

      }

    }
    
    for( uint32_t i=0; i<( GROUP_SIZE*256 ); i++ )
      coreId[groupId*( GROUP_SIZE*256 )+i] = groupCoreId[i];
    
    delete [] groupCoreId;
    groupCoreId = NULL;

  }

  for( uint32_t i=0; i<( GROUPS*GROUP_SIZE ); i++ ) {

#if !defined BATCH_OUTPUT
    cout << "thread(hostId=" << i/GROUP_SIZE << ",phiId=" << i%GROUP_SIZE << "):";
#endif
    for( uint32_t j=0; j<256; j++ ) {
      if( coreId[( i/GROUP_SIZE )*( GROUP_SIZE*256 )+( i%GROUP_SIZE )*256+j] != 0xFFFFFFFF ) {
#if defined BATCH_OUTPUT
	cout << "# group\tthread\tcoreid " << endl << i/GROUP_SIZE << "\t" << i << "\t" << j << endl;
#else
	cout << " " << j;
#endif	
      }
    }

#if defined BATCH_OUTPUT
    if( ( ( i+1 )%GROUP_SIZE ) == 0 ) {

      cout << endl;
      cout << endl;

    }
#else
    cout << endl;
#endif

  }

  exit( 0 );

}
