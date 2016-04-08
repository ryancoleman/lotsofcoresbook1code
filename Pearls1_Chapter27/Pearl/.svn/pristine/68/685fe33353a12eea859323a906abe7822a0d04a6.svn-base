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
 * Transforms Matrix-Market format plus meta-data after partitioning (by e.g.,
 * Mondriaan), to a binary triplet format optimised for BICRS use. 
 * 
 * Warning: this code uses assertions. Define the NDEBUG flag for production
 * use.
*/

//debug statements
//#define DEBUG

//more debug statements
//#define _DEBUG

//don't do any sanity checking whatsoever
//#define NDEBUG

//print info to standard out
#define _INFO

//print nonzero coordinates to standard out
//#define _SBD2TRP_COUT

#include "util.hpp"
#include "MinCCS.hpp"
#include "BlockCRS.hpp"
#include "U2.hpp"
#include "Duck.hpp"
#include "SepLast.hpp"
#include "MinCRS.hpp"
#include "FileToVT.hpp"

#define SUCCESS EXIT_SUCCESS
#define FAILURE EXIT_FAILURE

int main( int nargs, char** args ) {

	if( nargs != 4 && nargs != 6 ) {
		std::cout << "Usage: " << args[0] << " <matrix file> <row file> <column file> <output file> <block order>" << std::endl;
		std::cout << "where \t<matrix file> is a matrix in 2D SBD form stored in Matrix-Market format" << std::endl;
		std::cout << " \t<row file> determines row boundaries of the 2D SBD blocks" << std::endl;
		std::cout << " \t<column file> determines column boundaries of the 2D SBD blocks" << std::endl;
		std::cout << " \t<output file> output file name (binary triplet format)" << std::endl;
		std::cout << " \t<block order> an integer representing a specific block order (see below)" << std::endl;
		std::cout << std::endl << "OR" << std::endl << std::endl;
		std::cout << "Usage: " << args[0] << " <matrix file> <output file> <block order>" << std::endl;
		std::cout << "where \t: <matrix file> is a matrix in 2D SBD form in Extended Matrix-Market format" << std::endl;
		std::cout << " \t<output file> output file name (binary triplet format)" << std::endl;
		std::cout << " \t<block order> an integer representing a specific block order (see below)" << std::endl;
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

	if( nargs == 4 ) {
		std::cout << "Extended Matrix-Market format support not yet implemented; sorry." << std::endl;
		return FAILURE;
	}

	char blockorder = 99;

	//set block mode
	if( args[nargs-1][0]=='0' ) blockorder = 0;
	else if( args[nargs-1][0]=='1' ) blockorder = 1;
	else if( args[nargs-1][0]=='2' ) blockorder = 2;
	else if( args[nargs-1][0]=='3' ) blockorder = 3;
	else if( args[nargs-1][0]=='4' ) blockorder = 4;
	else if( args[nargs-1][0]=='5' ) blockorder = 5;
	else {
		std::cerr << "Undefined block ordering given (" << args[nargs-1][0] << ")! Exiting." << std::endl;
		return FAILURE;
	}

	//print info
	if( blockorder == 0 ) {
#ifdef _INFO
		std::cout << "Info: CRS block order will be applied" << std::endl;
#endif
	}
	else if( blockorder == 1 ) {
#ifdef _INFO
		std::cout << "Info: U2 block ordering will be applied" << std::endl;
#endif
	}
	else if( blockorder == 2 ) {
#ifdef _INFO
		std::cout << "Info: Duck block ordering will be applied" << std::endl;
#endif
	} else if( blockorder == 3 ) {
#ifdef _INFO
		std::cout << "Info: Separators are ordered last, with most misses in the row-direction" << std::endl;
#endif
	} else if( blockorder == 4 ) {
#ifdef _INFO
		std::cout << "Info: Minimal CCS ordering will be applied" << std::endl;
#endif
	} else if( blockorder == 5 ) {
#ifdef _INFO
		std::cout << "Info: Zig-zag CRS block ordering will be applied" << std::endl;
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

	//new method: using BlockOrderer classes
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

	std::vector< std::vector< Triplet< double > > > ordered = order->induce( naive, rowHierarchy, colHierarchy, rowBounds, colBounds, NULL );
	//clear old vector
	naive.clear();

	//flatten results
	std::vector< std::vector< Triplet< double > > >::iterator ordered_it = ordered.begin();
	for( ; ordered_it!=ordered.end(); ++ordered_it ) naive.insert( naive.end(), ordered_it->begin(), ordered_it->end() );

	//delete orderer
	delete order;

#ifdef _INFO
	std::cout << "Writing results..." << std::endl;
#endif
	//write away results
	Triplet< double >::save( std::string( args[4] ), &(naive[0]), m, n, naive.size() );

#ifdef _SBD2TRP_COUT
	//write to cout as well
	for( unsigned long int i=0; i<naive.size(); i++ )
		std::cout << naive[i].i() << "\t" << naive[i].j() << std::endl;
#endif

#ifdef _INFO
	std::cout << "Done!" << std::endl;
#endif
}

