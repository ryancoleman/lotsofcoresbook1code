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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2009.
 */


#define FIXED_SEED 13925

/**Define OUTPUT_Z to have the output vector printed to stdout */
//#define OUTPUT_Z

#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <unistd.h>
#ifndef _NO_LIBNUMA
#include <numa.h>
#endif

#include "SparseMatrix.hpp"
#include "TS.hpp"
#include "CRS.hpp"
#include "ICRS.hpp"
#include "oICRS.hpp"
#include "ZZ_CRS.hpp"
#include "ZZ_ICRS.hpp"
#include "SVM.hpp"
#include "HTS.hpp"
#include "BICRS.hpp"
#include "CBICRS.hpp"
#include "CCSWrapper.hpp"
#include "Hilbert.hpp"
#include "BlockHilbert.hpp"
#include "BetaHilbert.hpp"
#include "RDBHilbert.hpp"
#include "BisectionHilbert.hpp"
#include "RDScheme.hpp"
#include "McCRS.hpp"
#include "CompressedHilbert.hpp"

#ifdef WITH_MKL
 #include "MKLCRS.hpp"
#endif
#ifdef WITH_CUDA
 #include "CuHyb.hpp"
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
 * When using one of the parallel SpMV schemes, it may be useful
 * to allocate x using numa_interleaved_alloc, instead of using
 * the default local allocation policy. Define INTERLEAVE_X to
 * enable this whenever appropriate.
 */
#define INTERLEAVE_X

//no libnuma, no interleaving
#ifdef _NO_LIBNUMA
 #undef INTERLEAVE_X
#endif

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

Matrix< double >* selectMatrix( const int scheme, const int ccs, const std::string file ) {
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
		if( ccs ) {
			switch( scheme ) {
				case 0: { return new CCSWrapper< double, TS< double >, ULI >     ( input, m, n, zero ); }
				case 1: { return new CCSWrapper< double, CRS< double >, ULI >    ( input, m, n, zero ); }
				case 2: { return new CCSWrapper< double, ICRS< double >, ULI >   ( input, m, n, zero ); }
				case 3: { return new CCSWrapper< double, ZZ_CRS< double >, LI >  ( input, m, n, zero ); }
				case 4: { return new CCSWrapper< double, ZZ_ICRS< double >, LI > ( input, m, n, zero ); }
				case 5: { return new CCSWrapper< double, SVM< double >, ULI >    ( input, m, n, zero ); }
				case 6: { return new CCSWrapper< double, HTS< double >, ULI >    ( input, m, n, zero ); }
				case 7: { return new CCSWrapper< double, BICRS< double >, ULI >  ( input, m, n, zero ); }
				case 8: { return new CCSWrapper< double, Hilbert< double >, ULI >( input, m, n, zero ); }
				case 9: { return new CCSWrapper< double, BlockHilbert< double >, LI >( input, m, n, zero ); }
				case 10:{ return new CCSWrapper< double, BisectionHilbert< double >, LI >( input, m, n, zero ); }
				case 11:{ return CBICRS_factory< double >::getCBICCS( input, m, n, zero ); }
				case 12:{ return new CCSWrapper< double, BetaHilbert< double >, ULI >( input, m, n, zero ); }
				case 13:{ return new CCSWrapper< double, RDBHilbert< double >, ULI >( input, m, n, zero ); }
#ifdef WITH_CSB
				case 14:{ return new CCSWrapper< double, RDCSB< double >, ULI >( input, m, n, zero ); }
#endif
				case 15:{ return new CCSWrapper< double, RDScheme< double, Hilbert< double > >, ULI >( input, m, n, zero ); }
				case 16:{ return new CCSWrapper< double, McCRS< double >, ULI >( input, m, n, zero ); }
				case 17:{ return new CCSWrapper< double, RDScheme< double, CompressedHilbert< double > >, ULI >( input, m, n, zero ); }
#ifdef WITH_MKL
				case 18:{ return new CCSWrapper< double, MKLCRS< double >, ULI >( input, m, n, zero ); }
#endif
				case 19:{ return new CCSWrapper< double, oICRS< double >, ULI >( input, m, n, zero ); }
#ifdef WITH_CUDA
				case 20:{ return new CCSWrapper< double, CuHyb, ULI >( input, m, n, zero ); }
#endif
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		} else {
			switch( scheme ) {
				case 0: { return new TS< double >     ( input, m, n, zero ); }
				case 1: { return new CRS< double >    ( input, m, n, zero ); }
				case 2: { return new ICRS< double >   ( input, m, n, zero ); }
				case 3: { return new ZZ_CRS< double > ( input, m, n, zero ); }
				case 4: { return new ZZ_ICRS< double >( input, m, n, zero ); }
				case 5: { return new SVM< double >    ( input, m, n, zero ); }
				case 6: { return new HTS< double >    ( input, m, n, zero ); }
				case 7: { return new BICRS< double >  ( input, m, n, zero ); }
				case 8: { return new Hilbert< double >( input, m, n, zero ); }
				case 9: { return new BlockHilbert< double >( input, m, n, zero ); }
				case 10:{ return new BisectionHilbert< double >( input, m, n, zero ); }
				case 11:{ return CBICRS_factory< double >::getCBICRS( input, m, n, zero ); }
				case 12:{ return new BetaHilbert< double >( input, m, n, zero ); }
				case 13:{ return new RDBHilbert< double >( input, m, n, zero ); }
#ifdef WITH_CSB
				case 14:{ return new RDCSB< double >( input, m, n, zero ); }
#endif
				case 15:{ return new RDScheme< double, Hilbert< double > >( input, m, n, zero ); }
				case 16:{ return new McCRS< double >( input, m, n, zero ); }
				case 17:{ return new RDScheme< double, CompressedHilbert< double > >( input, m, n, zero ); }
#ifdef WITH_MKL
				case 18:{ return new MKLCRS< double >( input, m, n, zero ); }
#endif
				case 19:{ return new oICRS< double >( input, m, n, zero ); }
#ifdef WITH_CUDA
				case 20:{ return new CuHyb( input, m, n, zero ); }
#endif
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		}
	} else /*if( ext.compare( "mtx" ) == 0 )*/ {
		std::cout << "Matrix-market format expected, reading in..." << std::endl;
		if( ccs ) {
			switch( scheme ) {
				case 0: { return new CCSWrapper< double, TS< double >, ULI >     ( file, zero ); }
				case 1: { return new CCSWrapper< double, CRS< double >, ULI >    ( file, zero ); }
				case 2: { return new CCSWrapper< double, ICRS< double >, ULI >   ( file, zero ); }
				case 3: { return new CCSWrapper< double, ZZ_CRS< double >, LI > ( file, zero ); }
				case 4: { return new CCSWrapper< double, ZZ_ICRS< double >, LI >( file, zero ); }
				case 5: { return new CCSWrapper< double, SVM< double >, ULI >    ( file, zero ); }
				case 6: { return new CCSWrapper< double, HTS< double >, ULI >    ( file, zero ); }
				case 7: { return new CCSWrapper< double, BICRS< double >, ULI >  ( file, zero ); }
				case 8: { return new CCSWrapper< double, Hilbert< double >, ULI >( file, zero ); }
				case 9: { return new CCSWrapper< double, BlockHilbert< double >, LI >( file, zero ); }
				case 10:{ return new CCSWrapper< double, BisectionHilbert< double >, LI >( file, zero ); }
				case 11:{ return CBICRS_factory< double >::getCBICCS( file, zero ); }
				case 12:{ return new CCSWrapper< double, BetaHilbert< double >, ULI >( file, zero ); } 
				case 13:{ return new CCSWrapper< double, RDBHilbert< double >, ULI >( file, zero ); } 
#ifdef WITH_CSB
				case 14:{ return new CCSWrapper< double, RDCSB< double >, ULI >( file, zero ); }
#endif
				case 15:{ return new CCSWrapper< double, RDScheme< double, Hilbert< double > >, ULI >( file, zero ); }
				case 16:{ return new CCSWrapper< double, McCRS< double >, ULI >( file, zero ); }
				case 17:{ return new CCSWrapper< double, RDScheme< double, CompressedHilbert< double > >, ULI >( file, zero ); }
#ifdef WITH_MKL
				case 18:{ return new CCSWrapper< double, MKLCRS< double >, ULI >( file, zero ); }
#endif
				case 19:{ return new CCSWrapper< double, oICRS< double >, ULI >( file, zero ); }
#ifdef WITH_CUDA
				case 20:{ return new CCSWrapper< double, CuHyb, ULI >( file, zero ); }
#endif
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
		} else {
			switch( scheme ) {
				case 0: { return new TS< double >     ( file, zero ); }
				case 1: { return new CRS< double >    ( file, zero ); }
				case 2: { return new ICRS< double >   ( file, zero ); }
				case 3: { return new ZZ_CRS< double > ( file, zero ); }
				case 4: { return new ZZ_ICRS< double >( file, zero ); }
				case 5: { return new SVM< double >    ( file, zero ); }
				case 6: { return new HTS< double >    ( file, zero ); }
				case 7: { return new BICRS< double >  ( file, zero ); }
				case 8: { return new Hilbert< double >( file, zero ); }
				case 9: { return new BlockHilbert< double >( file, zero ); }
				case 10:{ return new BisectionHilbert< double >( file, zero ); }
				case 11:{ return CBICRS_factory< double >::getCBICRS( file, zero ); }
				case 12:{ return new BetaHilbert< double >( file, zero ); }
				case 13:{ return new RDBHilbert< double >( file, zero ); }
#ifdef WITH_CSB
				case 14:{ return new RDCSB< double >( file, zero ); }
#endif
				case 15:{ return new RDScheme< double, Hilbert< double > >( file, zero ); }
				case 16:{ return new McCRS< double >( file, zero ); }
				case 17:{ return new RDScheme< double, CompressedHilbert< double > >( file, zero ); }
#ifdef WITH_MKL
				case 18:{ return new MKLCRS< double >( file, zero ); }
#endif
				case 19:{ return new oICRS< double >( file, zero ); }
#ifdef WITH_CUDA
				case 20:{ return new CuHyb( file, zero ); }
#endif
				default: {
					std::cerr << "Invalid scheme ID, matrix not loaded into sparse matrix structure!" << std::endl;
					return NULL;
				}
			}
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
		std::cout << "Usage: " << argv[0] << " <mtx> <scheme> <x> <REP1> <REP2>" << std::endl << std::endl;
		std::cout << "calculates Ax=y and reports average time taken as well as the mean of y." << std::endl;
		std::cout << "with\t\t <mtx> filename of the matrix A in matrix-market or binary triplet format." << std::endl;
		std::cout << "    \t\t <scheme> number of a sparse scheme to use, see below." << std::endl;
		std::cout << "    \t\t <x> 0 for taking x to be the 1-vector, 1 for taking x to be random (fixed seed)." << std::endl;
		std::cout << "    \t\t <REP1> (optional, default is 1) number of repititions of the entire experiment." << std::endl;
		std::cout << "    \t\t <REP2> (optional, default is 1) number of repititions of the in-place SpMV multiplication, per experiment." << std::endl;
		std::cout << std::endl << "Possible schemes:" << std::endl;
		std::cout << " 0: TS (triplet scheme)" << std::endl;
		std::cout << " 1: CRS (also known as CSR)" << std::endl;
		std::cout << " 2: ICRS (Incremental CRS)" << std::endl;
		std::cout << " 3: ZZ-CRS (Zig-zag CRS)" << std::endl;
		std::cout << " 4: ZZ-ICRS (Zig-zag ICRS)" << std::endl;
		std::cout << " 5: SVM (Sparse vector matrix)" << std::endl;
		std::cout << " 6: HTS (Hilbert-ordered triplet scheme)" << std::endl;
		std::cout << " 7: BICRS (Bi-directional Incremental CRS)" << std::endl;
		std::cout << " 8: Hilbert (Hilbert-ordered triplets backed by BICRS)" << std::endl;
		std::cout << " 9: Block Hilbert (Sparse matrix blocking, backed by Hilbert and HBICRS)" << std::endl;
		std::cout << "10: Bisection Hilbert (Sparse matrix blocking by bisection, backed by Hilbert and HBICRS)" << std::endl;
		std::cout << "11: CBICRS (Compressed Bi-directional Incremental CRS)" << std::endl;
		std::cout << "12: Beta Hilbert (known as Block CO-H+ in the paper by Yzelman & Roose, 2012: parallel compressed blocked Hilbert with BICRS)" << std::endl;
		std::cout << "13: Row-distributed Beta Hilbert (known as Row-distributed block CO-H in the paper by Yzelman & Roose, 2012: same as 12, but simpler distribution)" << std::endl;
#ifdef WITH_CSB
		std::cout << "14: Row-distributed CSB (Uses CSB sequentially within the row-distributed scheme of 13)" << std::endl;
#endif
		std::cout << "15: Row-distributed Hilbert (Parallel row-distributed Hilbert scheme, see also 8)" << std::endl;
		std::cout << "16: Row-distributed parallel CRS (using OpenMP, known as OpenMP CRS in the paper by Yzelman & Roose, 2012)" << std::endl;
		std::cout << "17: Row-distributed SpMV using compressed Hilbert indices." << std::endl;
#ifdef WITH_MKL
		std::cout << "18: Intel MKL SpMV based on the CRS data structure." << std::endl;
#endif
		std::cout << "19: Optimised ICRS." << std::endl;
#ifdef WITH_CUDA
		std::cout << "20: CUDA CuSparse HYB format." << std::endl;
#endif
		std::cout << std::endl << "The in-place Ax=y calculation is preceded by a quasi pre-fetch." << std::endl;
		std::cout << "Add a minus sign before the scheme number to enable use of the CCS wrapper (making each CRS-based structure CCS-based instead)" << std::endl;
		std::cout << "Note: binary triplet format is machine-dependent. ";
		std::cout << "Take care when using the same binary files on different machine architectures." << std::endl;
		return EXIT_FAILURE;
	}

	std::string file = std::string( argv[1] );
	int scheme = atoi( argv[2] );
	int ccs    = scheme < 0 ? 1 : 0;
	if( ccs ) scheme = -scheme;
	int x_mode = atoi( argv[3] );
	unsigned long int rep1 = 1;
	unsigned long int rep2 = 1;
	if( argc >= 5 )
		rep1 = static_cast< unsigned long int >( atoi( argv[4] ) );
	if( argc >= 6 )
		rep2 = static_cast< unsigned long int >( atoi( argv[5] ) );

	if( scheme != 16 && scheme != -16 && //pin master thread to a single core
		scheme != 18 && scheme != -18 ) { //but not when OpenMP is used (otherwise serialised computations)
		cpu_set_t mask;
		CPU_ZERO( &mask );
		CPU_SET ( 0, &mask );
		if( pthread_setaffinity_np( pthread_self(), sizeof( mask ), &mask ) != 0 ) {
			std::cerr << "Error setting main thread affinity!" << std::endl;
			exit( 1 );
		}
	} else {
		omp_set_num_threads( MachineInfo::getInstance().cores() );
	}

#ifdef WITH_MKL
	if( scheme == 18 ) {
		mkl_set_num_threads( MachineInfo::getInstance().cores() );
	}
#endif
	std::cout << argv[0] << " called with matrix input file " << file << ", scheme number ";
	std::cout << scheme << " and x being " << (x_mode?"random":"the 1-vector") << "." << std::endl;
	std::cout << "Number of repititions of in-place zax is " << rep2 << std::endl;
	std::cout << "Number of repititions of the " << rep2 << " in-place zax(es) is " << rep1 << std::endl;

	Matrix< double >* checkm = new TS< double >( file );
	clock_gettime( CLOCK_ID, &start);
	Matrix< double >* matrix = selectMatrix( scheme, ccs, file );
	clock_gettime( CLOCK_ID, &stop);
	time  = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
	if( matrix == NULL ) {
		std::cerr << "Error during sparse scheme loading, exiting." << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Matrix dimensions: " << matrix->m() << " times " << matrix->n() << "." << std::endl;
	std::cout << "Datastructure loading time: " << time << " ms." << std::endl << std::endl;

	srand( FIXED_SEED );
	double* x = NULL;
#ifdef INTERLEAVE_X
	if( scheme == 13 || scheme == 14 || scheme == 15 || scheme == 16 || scheme == 17 || scheme == 18 )
		x = (double*) numa_alloc_interleaved( matrix->n() * sizeof( double ) );
	else
#endif
		x = (double*) _mm_malloc( matrix->n() * sizeof( double ), 64 );

	//initialise input vector
	for( unsigned long int j=0; j<matrix->n(); j++ ) {
		x[ j ] = x_mode?(rand()/(double)RAND_MAX):1.0;
	}

	//do one trial run, also for verification
	double* c = checkm->mv( x );
	clock_gettime( CLOCK_ID, &start );
	double* z = matrix->mv( x );
	clock_gettime( CLOCK_ID, &stop);
	time = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
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
#ifdef OUTPUT_Z
	for( unsigned long int j=0; j<matrix->m(); j++ ) {
		std::cout << z[ j ] << std::endl;
	}
#endif
	std::cout << "out-of-place z=Ax: mean= " << checksum( z, matrix->m() ) << ", ";
	std::cout << "MSE = " << checkMSE << ", ";
	std::cout << "max abs error = " << max_e << " while comparing y[ " << max_e_index << " ] = " <<  z[max_e_index] << " and c[ " << max_e_index << " ] = " <<  c[max_e_index] << ", ";
	std::cout << "time= " << time << " ms." << std::endl;
#ifdef RDBH_NO_COLLECT
	if( scheme == 13 ) {
		std::cout << "WARNING: MSE and max abs error are not correct for the Row-distributed Beta Hilbert scheme; please see the RDBHilbert.hpp file, and look for the RDBH_NO_COLLECT flag." << std::endl;
	}
#else
	if( scheme == 13 ) {
		std::cout << "WARNING: timings are pessimistic for the Row-distributed Beta Hilbert scheme; each spmv a (syncing) collect is executed to write local data to the global output vector as required by this library. To get the correct timings, turn this collect off via the RDBH_NO_COLLECT flag in the RDBHilbert.hpp file. Note that this causes the verification process to fail, since all data is kept in private local output subvectors." << std::endl;
	}
#endif
	double *times = new double[ rep1 ];

	//Run rep*rep instances
	for( unsigned long int run = 0; run < rep1; run++ ) {
		sleep( 1 );
		time = 0.0;
		//"prefetch"
		matrix->zax( x, z );
		matrix->zax( x, z, rep2, CLOCK_ID, &time );
		time /= static_cast<double>( rep2 );
		times[ run ] = time;
	}

	//calculate statistics
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

	std::cout << "In-place:" << std::endl;
	std::cout << "Mean  = " << checksum( z, matrix->m() ) << std::endl;
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

	delete [] times;
#ifdef INTERLEAVE_X
	if( scheme == 13 || scheme == 14 || scheme == 15 || scheme == 16 || scheme == 17 || scheme == 18 ) {
		numa_free( x, matrix->n() * sizeof( double ) );
	} else
#endif
		_mm_free( x );

	if( scheme == 12 || scheme == 13 || scheme == 14 || scheme == 15 || scheme == 16 || scheme == 17 || scheme == 18 ) {
#ifdef _NO_LIBNUMA
		_mm_free( z );
#else
		numa_free( z, matrix->m() * sizeof( double ) );
#endif
	} else {
		_mm_free( z );
	}
	_mm_free( c );
	delete matrix;
	delete checkm;

	return EXIT_SUCCESS;
}

