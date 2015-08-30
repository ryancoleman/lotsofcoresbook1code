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
#include <iomanip>
#include <cstdlib>
#include <cstdint>
#include <omp.h>
#include <sched.h>
#include <time.h>

using namespace std;

#if !defined LENGTH
#define LENGTH ( 1 )
#endif

//#define BATCH_OUTPUT
#define MIN_MINIMIZATION_RUNS ( 20 )
#define MAX_MINIMIZATION_RUNS ( 200 )
#define ACCEPTANCE ( 0.75 )
#define ITERATIONS ( 5000 )

#define MIN( X, Y ) ( ( X ) < ( Y ) ? ( X ) : ( Y ) )

double get_time_stamp();

// main function
int main( int argc, char **argv ) {
  
  MPI::Init( argc, argv );
  
  uint32_t
    mpiId = (uint32_t)MPI::COMM_WORLD.Get_rank(),
    numRanks = (uint32_t)MPI::COMM_WORLD.Get_size();

  // data package
  uint8_t
    *array = (uint8_t *) new uint8_t[LENGTH];

  for( uint32_t i=0; i<LENGTH; i++ )
    array[i] = (uint8_t)i;
  
  // create buffer on Xeon Phi and apply thread pinning
#pragma offload target( mic : 0 ) \
  in( array : length( LENGTH ) alloc_if( true ) free_if( false ) )
  {
    cpu_set_t
      cpuMask;
    
    CPU_ZERO( &cpuMask );
    CPU_SET( 1+4*mpiId, &cpuMask );
    sched_setaffinity( 0, sizeof( cpu_set_t ), &cpuMask );
  }

  uint32_t
    currentIteration = 0;
  uint64_t
    temp_1 = 0;
  double
    temp_2,
    temp_3,
    temp_4 = 1.0E100,
    temp_5,
    acceptanceRatio = 0.0;
  
  // iterate offload copy so that...
  //   - ...at least MIN_MINIMIZATION_RUNS runs have been performed, or
  //   - ...at most MAX_MINIMIZATION_RUNS runs have been performed, or
  //   - ...more than 75% of the data transfers fall into the same time frame
  while( ( currentIteration < MAX_MINIMIZATION_RUNS ) &&
	 ( currentIteration < MIN_MINIMIZATION_RUNS || acceptanceRatio < ACCEPTANCE ) ) {

    currentIteration++;
    
    double
      start[ITERATIONS],
      stop[ITERATIONS];
  
    // perform ITERATIONS offload copies to the Phi and take time stamps  
    MPI::COMM_WORLD.Barrier();
    
    for( uint32_t i=0; i<ITERATIONS; i++ ) { 
      
      start[i] = get_time_stamp();
      
#pragma offload_transfer target( mic : 0 ) \
  in( array : length( LENGTH ) alloc_if( false ) free_if( false ) )
      
      stop[i] = get_time_stamp();
      
    }
    
    double
      minStartTime,
      maxStopTime,
      startTime, // contains maximum start time
      stopTime;  // contains minimum end time
    
    // total execution time across all threads: maxStopTime-minStartTime
    MPI::COMM_WORLD.Reduce( &start[0], &minStartTime, 1, MPI::DOUBLE, MPI::MIN, 0 );
    MPI::COMM_WORLD.Reduce( &stop[ITERATIONS-1], &maxStopTime, 1, MPI::DOUBLE, MPI::MAX, 0 );
    
    // time frame in which all threads transfer data packages: stopTime-startTime
    MPI::COMM_WORLD.Reduce( &start[0], &startTime, 1, MPI::DOUBLE, MPI::MAX, 0 );
    MPI::COMM_WORLD.Reduce( &stop[ITERATIONS-1], &stopTime, 1, MPI::DOUBLE, MPI::MIN, 0 );
    
    // determine number of overlapping data transfers
    uint64_t
      overlappingTransfers = 0;
    
    for( uint32_t i=0; i<ITERATIONS; i++ )
      if( start[i] >= startTime && stop[i] <= stopTime )
	overlappingTransfers++;
    
    if( mpiId == 0 ) {
      
      for( uint32_t r=1; r<numRanks; r++ ) {
	
	MPI::COMM_WORLD.Recv( start, ITERATIONS, MPI::DOUBLE, r, 10 );
	MPI::COMM_WORLD.Recv( stop, ITERATIONS, MPI::DOUBLE, r, 20 );
	
	for( uint32_t i=0; i<ITERATIONS; i++ )
	  if( start[i] >= startTime && stop[i] <= stopTime )
	    overlappingTransfers++;
	
      }

      // assume 100% concurrency and distribute all messages evenly across
      // all threads, and estimate package transfertime
      if( overlappingTransfers > 0 &&
	  temp_4 > ( ( stopTime-startTime )*numRanks*1.0E6/overlappingTransfers ) ) {

	temp_1 = overlappingTransfers;
	temp_2 = ( stopTime-startTime )*1.0E3;
	temp_3 = ( maxStopTime-minStartTime )*1.0E6/ITERATIONS;
	temp_4 = ( stopTime-startTime )*numRanks*1.0E6/overlappingTransfers;
	temp_5 = overlappingTransfers*LENGTH*1.0E-9/( stopTime-startTime );

      }
      
      acceptanceRatio = (double)( temp_1 )/( numRanks*ITERATIONS );
      
    } else {
      
      MPI::COMM_WORLD.Send( start, ITERATIONS, MPI::DOUBLE, 0, 10 );
      MPI::COMM_WORLD.Send( stop, ITERATIONS, MPI::DOUBLE, 0, 20 );
      
    }
    
    MPI::COMM_WORLD.Bcast( &acceptanceRatio, 1, MPI::DOUBLE, 0 );

    // for larges packages a single thread already saturates the PCIe bandwidth 
    // -> just one run necessary
    if( LENGTH > 128*1024 )
      break;

  }

  // free buffers on Xeon Phi
#pragma offload target( mic : 0 ) \
  in( array : length( 0 ) alloc_if( false ) free_if( true ) )
  {
    ;
  }
  
  // release host memory
  delete [] array;

  // print results
  if( mpiId == 0 ) {

#if defined BATCH_OUTPUT
    cout << LENGTH << " " << numRanks << " " << ::setprecision( 6 ) << MIN( temp_3, temp_4 ) << " " << temp_5 << endl;
#else
    cout << "# overlapping transfers: " << temp_1 << endl;
    cout << "# overlapping transfer time: " << temp_2 << " ms" << endl;
    cout << "# naive transfertime per package: " << temp_3 << " us" << endl;
    cout << "# packagesize[bytes], threadcount, transfertime per package[us], bandwidth[GiB/s]" << endl;
    cout << LENGTH << "\t" << numRanks << "\t" << ::setprecision( 4 ) << MIN( temp_3, temp_4 ) << "\t\t" << temp_5 << endl;
#endif

  }
  
  MPI::Finalize();
	    
  return 0;
  
}

// timer
double get_time_stamp() {

  timespec t;
  clock_gettime( CLOCK_REALTIME, &t ); // system-wide clock
  return (double)( t.tv_sec+t.tv_nsec*1.0E-9 );

}
