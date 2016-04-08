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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2010.
 * 
 * Reads in Matrix-Market format plus meta-data after partitioning (by e.g.,
 * Mondriaan), and loads the doubly SBD format into a hierarchical BICRS scheme
 * for enhanced SpMV. Will use ICRS and ICCS formats within the hierarchical
 * scheme.
 *
 * Warning: this code uses assertions. Define the NDEBUG flag for production
 * use.
 */


#include<iostream>
#include<fstream>
#include<limits>
#include<string>
#include<vector>
#include<map>

#include<assert.h>

#include "Triplet.hpp"
#include "FileToVT.hpp"
#include "HBICRS.hpp"
#include "MinCCS.hpp"
#include "MinCRS.hpp"
#include "BlockCRS.hpp"
#include "U2.hpp"
#include "Duck.hpp"
#include "SepLast.hpp"
#include "util.hpp"

//#define _DEBUG
//#define NDEBUG
//#define _INFO
//#define _SBD2TRP_COUT
//#define OUTPUT_Z
//#define x_mode 0 //input one-vector (useful for checksumming)
#define x_mode 1 //random input vector

#define SUCCESS EXIT_SUCCESS
#define FAILURE EXIT_FAILURE

unsigned long int count( std::vector< std::vector< Triplet< double > > > &hierarchy ) {
	unsigned long int nnz = 0;
	for( unsigned long int i=0; i<hierarchy.size(); i++ )
		nnz += hierarchy[ i ].size();
	return nnz;
}

int main( int nargs, char** args ) {
	//timer structs
	struct timespec start, stop;
	double time;

	if( nargs != 5 && nargs != 7 ) {
		std::cout << "Usage: " << args[0] << " <matrix file> <row file> <column file> <repeats> <block order> <incremental>" << std::endl;
		std::cout << "where \t<matrix file> is a matrix in 2D SBD form stored in Matrix-Market format" << std::endl;
		std::cout << " \t<row file> determines row boundaries of the 2D SBD blocks" << std::endl;
		std::cout << " \t<column file> determines column boundaries of the 2D SBD blocks" << std::endl;
		std::cout << " \t<repeats> will apply repeat * repeat SpMVs" << std::endl;
		std::cout << " \t<block order> an integer representing a specific block order (see below)" << std::endl;
		std::cout << " \t<incremental> 1 for ICRS/ICCS, 0 for CRS/CCS" << std::endl;
		std::cout << std::endl << "OR" << std::endl << std::endl;
		std::cout << "Usage: " << args[0] << " <matrix file> <repeats> <block order> <incremental>" << std::endl;
		std::cout << "where \t: <matrix file> is a matrix in 2D SBD form in Extended Matrix-Market format" << std::endl;
		std::cout << " \t<repeats> how many times to repeat the SpMV" << std::endl;
		std::cout << " \t<block order> an integer representing a specific block order (see below)" << std::endl;
		std::cout << " \t<incremental> 1 for ICRS/ICCS, 0 for CRS/CCS" << std::endl;
		std::cout << std::endl;
		std::cout << "Possible block orders:" << std::endl;
		std::cout << "0: CRS" << std::endl;
		std::cout << "1: U2 block ordering, cache-inefficiencies in the row direction" << std::endl;
		std::cout << "2: Duck block ordering, cache-inefficiencies in the column direction" << std::endl;
		std::cout << "3: Separators last ordering, cache-inefficiencies in the row direction" << std::endl;
		std::cout << "4: Zig-zag CCS block ordering" << std::endl;
		std::cout << "5: Zig-zag CRS block ordering" << std::endl;
		std::cout << std::endl;
		std::cout << "Use, for example, the Mondriaan package to obtain 2D SBD forms for arbitrary sparse matrices; see:" << std::endl;
		std::cout << "http://www.math.uu.nl/people/bisselin/Mondriaan" << std::endl;
		std::cout << "For more information on SBD form and cache-oblivious SpMV, see the relevant paper(s) at:" << std::endl;
		std::cout << "http://www.math.uu.nl/people/yzelman/publications" << std::endl;
		return SUCCESS;
	}

	unsigned long int repeat = 0;
	unsigned long int incr   = 2; //set to invalid
	if( nargs == 5 ) {
		std::cout << "Extended Matrix-Market format support not yet implemented; sorry." << std::endl;
		return FAILURE;
	} else {
		repeat = atoi( args[ 4 ] );
	}

	char blockorder = 99;

	//set block mode
	if( args[nargs-2][0]=='0' ) blockorder = 0;
	else if( args[nargs-2][0]=='1' ) blockorder = 1;
	else if( args[nargs-2][0]=='2' ) blockorder = 2;
	else if( args[nargs-2][0]=='3' ) blockorder = 3;
	else if( args[nargs-2][0]=='4' ) blockorder = 4;
	else if( args[nargs-2][0]=='5' ) blockorder = 5;
	else {
		std::cerr << "Undefined block ordering given (" << args[nargs-2][0] << ")! Exiting." << std::endl;
		return FAILURE;
	}
	incr = atoi( args[ nargs-1 ] );

	if( blockorder == 0 ) {
#ifdef _INFO
		std::cout << "Info: CRS block order will be applied" << std::endl;
#endif
	}
	else if( blockorder == 1 ) {
#ifdef _INFO
		std::cout << "Info: U2 block order will be applied" << std::endl;
#endif
	}
	else if( blockorder == 2 ) {
#ifdef _INFO
		std::cout << "Info: Duck block order will be applied" << std::endl;
#endif
	} else if( blockorder == 3 ) {
#ifdef _INFO
		std::cout << "Info: Separators are ordered last, with most misses in the row-direction" << std::endl;
#endif
	} else if( blockorder == 4 ) {
#ifdef _INFO
		std::cout << "Info: Minimal CCS block ordering will be applied" << std::endl;
#endif
	} else if( blockorder == 5 ) {
#ifdef _INFO
		std::cout << "Info: Minimal CRS block ordering will be applied" << std::endl;
#endif
	}

	std::vector< Triplet< double > > naive;
	std::vector< unsigned long int > rowBounds;
	std::vector< unsigned long int > colBounds;
	std::vector< unsigned long int > rowHierarchy;
	std::vector< unsigned long int > colHierarchy;
	std::vector< std::vector< unsigned long int >* > rowFile;
	std::vector< std::vector< unsigned long int >* > colFile;
	rowFile.push_back( &rowBounds );
	rowFile.push_back( &rowHierarchy );
	colFile.push_back( &colBounds ); 
	colFile.push_back( &colHierarchy );
	unsigned long int nnz;
	ULI m, n;

	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start );

	naive = FileToVT::parse( std::string( args[1] ), m, n, nnz );
	if( parseULIVector( std::string( args[2] ), rowFile ) == FAILURE || 
		parseULIVector( std::string( args[3] ), colFile ) == FAILURE ) {
		std::cerr << "I/O Failure" << std::endl;
		return FAILURE;
	}

#ifdef _INFO
	std::cout << std::string(args[1]) << " read in; m=" << m << ", n=" << n << ", nnz=" << nnz << std::endl;
#endif

	//sanity checking using assertions, transforming to 0-based; O(m+n+nnz) operation
	std::vector< unsigned long int >::iterator it = rowBounds.begin();
	for( ; it!=rowBounds.end(); ++it ) {
		--(*it); //0-based
		assert( *it <= m );
	}
	it = colBounds.begin();
	for( ; it!=colBounds.end(); ++it ) {
		--(*it); //0-based
		assert( *it <= n );
	}

#ifdef _INFO
	std::cout << "Inducing block order..." << std::endl;
#endif

	BlockOrderer< double > *order;
	switch( blockorder ) {
	case 0:
		order = new BlockCRS< double >();
		break;
	case 1:
		order = new U2< double >();
		break;
	case 2:
		order = new Duck< double >();
		break;
	case 3:
		order = new SepLast< double >();
		break;
	case 4:
		order = new MinCCS< double >();
		break;
	case 5:
		order = new MinCRS< double >();
		break;
	default:
		std::cerr << "Invalid or unimplemented blockorder " << static_cast< int >( blockorder ) << std::endl;
		exit( EXIT_FAILURE );
	}

	std::vector< signed char > hierarchy_datatype;
	std::vector< std::vector< Triplet< double > > > hierarchy = 
		order->induce( naive, rowHierarchy, colHierarchy, rowBounds, colBounds, &hierarchy_datatype );
	assert( hierarchy_datatype.size() == hierarchy.size() );


	//happy asserting
#ifndef NDEBUG
	assert( count( hierarchy ) == nnz );
#endif

	//clear naive datastructure
	naive.clear();

#ifdef _INFO
	std::cout << "Number of blocks found:    " << hierarchy.size() << std::endl;
	std::cout << "Total number of nonzeroes: " << count( naive ) << std::endl;
	std::cout << "Loading into HBICRS..." << std::endl;
#endif

	//feed data to constructor
	if( !incr ) {
		//change Incremental CRS/CCS to normal CRS/CCS
		for( unsigned long int i=0; i<hierarchy_datatype.size(); ++i ) {
			if( hierarchy_datatype[i] == -2 ) hierarchy_datatype[i] = -1;
			else if( hierarchy_datatype[i] == 2 ) hierarchy_datatype[i] = 1;
		}
	}
	HBICRS< double > matrix( hierarchy, &(hierarchy_datatype[0]), m, n, 0.0 );
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop );
	time  = (stop.tv_sec-start.tv_sec)*1000;
	time += (stop.tv_nsec-start.tv_nsec)/1000000.0;

	//clear hierarchy structure
	hierarchy.clear();
	//clear block orderer
	delete order;

	std::cout << "Buildtime = " << time << std::endl;

#ifdef _INFO
	std::cout << "Going for " << repeat << " SpMVs..." << std::endl;
#endif

	//construct  x
	double *x = new double[ n ];
	for( unsigned long int i=0; i<n; i++ ) x[ i ] = x_mode?rand():1.0;

	double *times = new double[ repeat ];

	//Initialise
	double *z = matrix.mv( x );
	for( unsigned long int run = 0; run < repeat; run++ ) {
		//"prefetch"
		matrix.zax( x, z );
		//do it
		clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start );
		for( unsigned long int i=0; i<repeat; i++ )
			matrix.zax( x, z );
		clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop );
		time  = (stop.tv_sec-start.tv_sec)*1000;
		time += (stop.tv_nsec-start.tv_nsec)/1000000.0;
		time /= static_cast< double >( repeat );
		times[ run ] = time;
	}
	
	//calculate statistics, checksum
	double meantime, mintime, vartime;
	meantime = vartime = 0.0;
	double mean = 0.0;
	for( unsigned long int i=0; i<m; i++ ) {
#ifdef OUTPUT_Z
		std::cout << z[i] << std::endl;
#endif
		mean += z[ i ] / static_cast< double >( m );
	}
	mintime = times[ 0 ];
	for( unsigned long int run = 0; run < repeat; run++ ) {
		if( times[ run ] < mintime ) mintime = times[ run ];
		meantime += times[ run ] / static_cast< double >( repeat );
	}
	for( unsigned long int run = 0; run < repeat; run++ ) {
		vartime += ( times[ run ] - meantime ) * ( times[ run ] - meantime ) / static_cast< double >( repeat - 1 );
	}

	//free up memory
	delete [] times;
	delete [] x;
	delete [] z;

	std::cout << "Mean = " << mean << std::endl;
	std::cout << "Time = " << meantime << " (average), \t" <<  mintime << " (fastest), \t" << vartime << " (variance) ms. " << std::endl;
#ifdef _INFO
	std::cout << "Done!" << std::endl;
#endif
}

