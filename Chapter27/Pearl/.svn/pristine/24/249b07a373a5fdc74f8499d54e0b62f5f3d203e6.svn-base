/*
 * Copyright (c) 2007-2014, A. N. Yzelman,   Utrecht University 2007-2011;
 *                                                    KU Leuven 2011-2014.
 *                          R. H. Bisseling, Utrecht University 2007-2014.
 * 
 * This file is part of the Sparse Library.
 * 
 * This library was developed under supervision of Prof. dr. Rob H. Bisseling at
 * Utrecht University, from 2007 until 2011. From 2011-2014, development continued 
 * at KU Leuven, where Prof. dr. Dirk Roose contributed significantly to the ideas 
 * behind the newer parts of the library code.
 * 
 *     The Sparse Library is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by the
 *     Free Software Foundation, either version 3 of the License, or (at your
 *     option) any later version.
 * 
 *     The Sparse Library is distributed in the hope that it will be useful, but
 *     WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 *     or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 *     for more details.
 * 
 *     You should have received a copy of the GNU General Public License along
 *     with the Sparse Library. If not, see <http://www.gnu.org/licenses/>.
 */


/*
 * File created by:
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2014.
 */


#define FIXED_SEED 13925

/**Define OUTPUT_Z to have the output vector printed to stdout */
//#define OUTPUT_Z

#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>
#include <unistd.h>
#ifndef _NO_LIBNUMA
 #include <numa.h>
#endif

#include "Matrix.hpp"
#include "oICRS.hpp"
#include "FBICRS.hpp"
#include "RDBHilbert.hpp"

#ifndef RDBH_NO_COLLECT
#include "TS.hpp"
#endif

//The below is good for getting total user time,
//but will thus combine clocks for multiple threads
//(i.e., not useful for parallel SpMV timing)
//#define CLOCK_ID CLOCK_PROCESS_CPUTIME_ID

//This is a real-time clock, but may be less
//accurate with respect to CLOCK_PROCESS_CPUTIME_ID
#define CLOCK_ID CLOCK_MONOTONIC

//WARNING: a different clock you might be tempted to
//         use is CLOCK_THREAD_CPUTIME_ID.
//That clock would report time used by the master
//thread. This should be fine for sequential SpMVs,
//but for parallel (threaded) SpMVs that clock will
//FAIL since it does NOT count waits (the numbers
//will thus be waaay too low).
//
//In summary: never use CLOCK_THREAD_CPUTIME_ID here.

/**
 * When the CSB library is locally available, define WITH_CSB
 * flag for a row-distributed CSB-backed parallel scheme (14).
 */
//#define WITH_CSB

#ifdef WITH_CSB
	#include "RDCSB.hpp"
#endif

double checksum( double* z, unsigned long int m ) {
	if( m==0 ) return 0.0;
	double sum = z[ 0 ];
	for( unsigned long int i=1; i<m; i++ )
		sum += z[ i ];
	return sum / static_cast< double >( m );
}

Matrix< double >* selectMatrix( const int p, const int q, const std::string file ) {
	double zero = 0.0;
	size_t pos = file.find_last_of( '.' );
	std::string ext = file.substr( pos + 1, file.length() );
	std::cout << "Matrix file extension: \"" << ext << "\"." << std::endl;
	ULI m, n;
	std::vector< Triplet< double > > input;
	if( ext.compare( "trp" ) == 0 ) {
		std::cout << "Expecting binary triplet format, reading in..." << std::endl;
		input = Triplet< double >::load( file, m, n );
	} else if( ext.compare( "crs" ) == 0 || ext.compare( "csr" ) == 0 ) {
		std::cout << "Expecting a text-based CRS format, reading in..." << std::endl;
		input = Triplet< double >::loadCRS( file, m, n );
	}
	if( input.size() > 0 ) { //if the file was parsed into a vector of triplets, then follow the below
		if( p == 1 && q == 1 ) {
			return new RDBHilbert< double, FBICRS< double, LI, ICRS< double, uint16_t >, 15 > >( input, m, n, zero );
		} else if( p == 1 && q == 8 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 1, 8, int16_t > > >( input, m, n, zero );
		} else if( p == 2 && q == 4 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 2, 4, int16_t > > >( input, m, n, zero );
		} else if( p == 4 && q == 2 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 4, 2, int16_t > > >( input, m, n, zero );
		} else if( p == 8 && q == 1 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 8, 1, int16_t > > >( input, m, n, zero );
		} else {
			std::cerr << "Error: undefined case for p (" << p << ") and q (" << q << ")!" << std::endl;
			exit( EXIT_FAILURE );
		}
	} else /*if( ext.compare( "mtx" ) == 0 )*/ {
		std::cout << "Matrix-market format expected, reading in..." << std::endl;
		if( p == 1 && q == 1 ) {
			return new RDBHilbert< double, FBICRS< double, LI, ICRS< double, uint16_t >, 15 > >( file, zero );
		} else if( p == 1 && q == 8 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 1, 8, int16_t > > >( file, zero );
		} else if( p == 2 && q == 4 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 2, 4, int16_t > > >( file, zero );
		} else if( p == 4 && q == 2 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 4, 2, int16_t > > >( file, zero );
		} else if( p == 8 && q == 1 ) {
			return new RDBHilbert< double, FBICRS< double, LI, oICRS< double, 8, 1, int16_t > > >( file, zero );
		} else {
			std::cerr << "Error: undefined case for p (" << p << ") and q (" << q << ")!" << std::endl;
			exit( EXIT_FAILURE );
		}
	}
}

int main( int argc, char** argv ) {
	struct timespec start, stop; 
	double time;

#ifndef NDEBUG
	std::cout << "-->WARNING: COMPILED *WITH* ASSERTIONS!<--" << std::endl;
#endif
	
	if( argc<=3 ) {
		std::cout << "Usage: " << argv[0] << " <mtx> <p> <q>" << std::endl << std::endl;
		std::cout << "calculates Ax=y and reports average time taken as well as the mean of y." << std::endl;
		std::cout << "with\t\t <mtx> filename of the matrix A in matrix-market or binary triplet format." << std::endl;
		std::cout << "    \t\t <p> the row-wise block length." << std::endl;
		std::cout << "    \t\t <q> the column-wise block length." << std::endl;
		std::cout << "(This program is written specifically for the mic; pq should always equal 8!)" << std::endl;
		std::cout << "Note: binary triplet format is machine-dependent. ";
		std::cout << "Take care when using the same binary files on different machine architectures." << std::endl;
		return EXIT_FAILURE;
	}

	std::string file = std::string( argv[1] );
	int p = atoi( argv[2] );
	int q = atoi( argv[3] );
	//pin master thread to a single core
	cpu_set_t mask;
	CPU_ZERO( &mask );
	CPU_SET ( 0, &mask );
	if( pthread_setaffinity_np( pthread_self(), sizeof( mask ), &mask ) != 0 ) {
		std::cerr << "Error setting main thread affinity!" << std::endl;
		exit( EXIT_FAILURE );
	}

	std::cout << argv[0] << " called with matrix input file " << file << "." << std::endl;

	clock_gettime( CLOCK_ID, &start );
	Matrix< double >* matrix = selectMatrix( p, q, file );
	clock_gettime( CLOCK_ID, &stop );
	time  = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
	if( matrix == NULL ) {
		std::cerr << "Error during sparse scheme loading, exiting." << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Matrix dimensions: " << matrix->m() << " times " << matrix->n() << "." << std::endl;
	std::cout << "Datastructure loading time: " << time << " ms." << std::endl << std::endl;

	srand( FIXED_SEED );
	double *x = (double*) _mm_malloc( matrix->n() * sizeof( double ), 64 );

	//initialise input vector
	for( unsigned long int j = 0; j < matrix->n(); j++ ) {
		x[ j ] = rand()/((double)RAND_MAX+1);
	}

	//do verification
	double* z = matrix->mv( x );
#ifndef RDBH_NO_COLLECT
	//build verification matrix
	TS< double > * verificationMatrix = new TS< double >( file, 0.0 );
	//do verified SpMV
	double* c = verificationMatrix->mv( x );
	//check for errors
	double checkMSE = 0;
	unsigned long int max_e_index = 0;
	double max_e = fabs( z[0] - c[0] );
	for( unsigned long int j=0; j<matrix->m(); j++ ) {
		double curdiff = fabs( z[j] - c[j] );
		if( curdiff > max_e ) {
			max_e = curdiff;
			max_e_index = j;
		}
		curdiff  *= curdiff;
		curdiff  /= (double)(matrix->m());
		checkMSE += curdiff;
	}
	//report
	std::cout << "Verification step. MSE = " << checkMSE << ", max abs error = " << max_e << " while comparing y[ " << max_e_index << " ] = " <<  z[max_e_index] << " and c[ " << max_e_index << " ] = " <<  c[max_e_index] << std::endl;
	//clean up
	_mm_free( c );
	delete verificationMatrix;
#else
	std::cout << "Verification step skipped. To enable, disable the RDBH_NO_COLLECT flag (see the RDBHilbert.hpp file)." << std::endl;
#endif

	//prep time measurement
	time = 0;
	//warm up
	matrix->zax( x, z );
	//now with timer
	matrix->zax( x, z, 1, CLOCK_ID, &time );
	std::cout << "Single SpMV time = " << time << " ms." << std::endl;
	//derive number of SpMVs to get to a one-minute runtime, approx.
	const size_t N = static_cast< size_t >( 60000.0 / time ) + 1;
	//split in two
	const size_t rep1 = N >= 10 ? 10 : 3;
	const size_t rep2 = N / rep1 + 1;

	//Run ~N instances
	clock_gettime( CLOCK_ID, &start );
	double times[ rep1 ];
	for( unsigned long int run = 0; run < rep1; run++ ) {
		//sleep a little
		sleep( 1 );
		//reset time
		time = 0;
		//"prefetch"
		matrix->zax( x, z );
		//`real' benchmark
		matrix->zax( x, z, rep2, CLOCK_ID, &time );
		//record time
		times[ run ] = time / static_cast< double >(rep2);
	}
	clock_gettime( CLOCK_ID, &stop);
	time = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
	std::cout << "In-place results, averaged over " << N << " ~= " << rep1 << "*" << rep2 << " runs, ran in " << time << "ms.:" << std::endl;

	//calculate time statistics
	double meantime, mintime, vartime;
	meantime = vartime = 0.0;
	mintime = times[ 0 ];
	for( unsigned long int run = 0; run < rep1; run++ ) {
		if( times[ run ] < mintime ) mintime = times[ run ];
		meantime += times[ run ] / static_cast< double >( rep1 );
	}
	for( unsigned long int run = 0; run < rep1; run++ ) {
		vartime += ( times[ run ] - meantime ) * ( times[ run ] - meantime ) / static_cast< double >( rep1 - 1 );
	}
	vartime = sqrt( vartime );

	//output statistics
	std::cout << "Time  = " << meantime << " (average), \t" <<  mintime << " (fastest), \t" << vartime << " (stddev) ms. " << std::endl;
	const double avgspeed = static_cast< double >( 2*matrix->nzs() ) / meantime / 1000000.0;
	const double minspeed = static_cast< double >( 2*matrix->nzs() ) / mintime / 1000000.0;
	const double varspeed = fabs( avgspeed - static_cast< double >( 2*matrix->nzs() ) / (meantime - vartime) / 1000000.0 );
	std::cout << "Speed = " << avgspeed << " (average), \t"
					<< minspeed << " (fastest), \t"
					<< varspeed << " (variance) Gflop/s." << std::endl;
	const size_t memuse1 = matrix->bytesUsed() + sizeof( double ) * 2 * matrix->nzs();
	const double avgmem1 = static_cast< double >( 1000*memuse1 ) / meantime / 1073741824.0;
	const double minmem1 = static_cast< double >( 1000*memuse1 ) / mintime / 1073741824.0;
	const double varmem1 = fabs( avgmem1 - static_cast< double >( 1000*memuse1 ) / (meantime-vartime) / 1073741824.0 );
	std::cout << "        " << avgmem1 << " (average), \t"
					<< minmem1 << " (fastest), \t"
					<< varmem1 << " (variance) Gbyte/s (upper bound)." << std::endl;
	const size_t memuse2 = matrix->bytesUsed() + sizeof( double ) * ( matrix->m() + matrix->n() );
	const double avgmem2 = static_cast< double >( 1000*memuse2 ) / meantime / 1073741824.0;
	const double minmem2 = static_cast< double >( 1000*memuse2 ) / mintime / 1073741824.0;
	const double varmem2 = fabs( avgmem2 - static_cast< double >( 1000*memuse2 ) / (meantime-vartime) / 1073741824.0 );
	std::cout << "        " << avgmem2 << " (average), \t"
					<< minmem2 << " (fastest), \t"
					<< varmem2 << " (variance) Gbyte/s (lower bound)." << std::endl;

	//clean up
	_mm_free( x );
	_mm_free( z );
	delete matrix;

	return EXIT_SUCCESS;
}

