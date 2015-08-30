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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2011.
 *
 * Multicore fully 2D distributed sparse matrix--vector multiplication.
 * Can perform 2D cache-oblivious multiplications using HBICRS.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>
#include <sstream>
#include <cmath>
#include <pthread.h>

#ifndef _NO_LIBNUMA
 #include <numa.h>
#endif

#include "Matrix.hpp"
#include "Triplet.hpp"
#include "Upscaler.hpp"

/** How  many times to repeat SpMVs for timing */
#define REPEAT_EXP 1000


/** USER SHOULD DEFINE ***ONE*** OF THE FOLLOWING FOUR FLAGS */

/** If defined, this program performs a 2D-SBD cache-oblivious SpMV on each core, backed by CBICRS */
#define _MCDMV_BICRS

/** If defined, this program performs a 2D-SBD cache-oblivious SpMV on each core */
//#define _MCDMV_SBD

/** If defined, this program performs a BetaHilbert parallel multiply on each local part (hybrid) */
//#define _MCDMV_BH

/** If defined, this program performs a row-distributed BetaHilbert parallel multiply on each local part (hybrid) */
//#define _MCDMV_1DBH

/** ---------------------------------------------------------------- */


/** USER SHOULD DEFINE ***ONE*** OF THE FOLLOWING THREE FLAGS */

/** If defined, this program processes inter-part communication in-order of storage in the SBD tree */
//#define _MCDMV_REMOTE_INORDER

/** If defined, inter-part communication is executed in row-major order */
#define _MCDMV_REMOTE_ROW_MAJOR

/** If defined, this program performs a BetaHilbert parallel multiply for inter-part areas of the SpMV -- 
 *  use this or MCDMV_REMOTE_1DBH (next flag) when usingg _MCDMV_BH or _MCDMV_1DBH, for scalability */
//#define _MCDMV_REMOTE_BH

/** If defined, this program performs a row-distributed BetaHilbert parallel multiply on inter-part areas --
 *  use this or MCDMV_REMOTE_BH (previous flag) when using _MCDMV_BH or _MCDMV_1DBH, for scalability */
//#define _MCDMV_REMOTE_1DBH

/** ---------------------------------------------------------------- */

//apply compression on index data used during fan-in/fan-out
#define _MCDMV_COMPRESSED_COMMUNICATION

//sort separator ranges according to vector distribution
//#define _MCDMV_VECTOR_SORT

//take x=1 instead of random
//#define ONE_X

//define for spin locks instead of mutex-based barriers
#define SPINLOCK

//enable to print out local output vector
//#define PRINTY

//include files based on the above defs
#if defined _MCDMV_SBD || defined _MCDMV_BICRS
#include "MinCRS.hpp"
#endif
#if defined _MCDMV_SBD || defined _MCDMV_BICRS
#include "CBICRS.hpp"
#endif
#ifdef _MCDMV_SBD
#include "HBICRS.hpp"
#endif
#if defined _MCDMV_BH || defined _MCDMV_REMOTE_BH
#include "BetaHilbert.hpp"
#endif
#if defined _MCDMV_1DBH || defined _MCDMV_REMOTE_1DBH
#include "RDBHilbert.hpp"
#endif
#if defined _MCDMV_REMOTE_INORDER || defined _MCDMV_REMOTE_ROW_MAJOR
#include "CBICRS.hpp"
#endif
#if defined _MCDMV_REMOTE_ROW_MAJOR
#include "ICRS.hpp"
#endif

//set convenience flags
#if (!defined _MCDMV_SBD && !defined _MCDMV_BICRS) || (!defined _MCDMV_REMOTE_INORDER && !defined _MCDMV_REMOTE_ROW_MAJOR)
#define _MCDMV_HYBRID
#endif

#ifdef _MCDMV_COMPRESSED_COMMUNICATION
class fanQuadlet {
public:
	unsigned long int remoteP;
	unsigned long int localStart;
	unsigned long int remoteStart;
	unsigned long int length;
	fanQuadlet( const unsigned long int _P, const unsigned long int _lS, const unsigned long int _rS, const unsigned long int _l ):
		remoteP( _P ), localStart( _lS ), remoteStart( _rS ), length( _l ) {}
};
typedef std::vector< fanQuadlet > fanInContainerType;
typedef std::vector< fanQuadlet > fanOutContainerType;
#else
class fanPair {
public:
	unsigned long int P;
	unsigned long int i;
	fanPair() { P = ULONG_MAX; i = ULONG_MAX; }
	fanPair( unsigned long int _P, unsigned long int _i ): P( _P ), i( _i ) {}
};
typedef std::map< unsigned long int, fanPair > fanInContainerType;
typedef std::multimap< unsigned long int, fanPair > fanOutContainerType;
#endif

class MVData {
public:
	//temporary storage used for construction
	std::vector< Triplet< double > > local;
	//distributed temporary structures
	std::vector< std::vector< Triplet< double > > > main_matrices;
	std::vector< std::vector< Triplet< double > > > remote_separators;
	//remember lower and higher bounds for local matrix parts
	std::vector< unsigned long int > lo_i;
	std::vector< unsigned long int > lo_j;
	std::vector< unsigned long int > hi_i;
	std::vector< unsigned long int > hi_j;
	//fan-in/out data
	std::vector< fanInContainerType > fanIn;
	std::vector< fanOutContainerType > fanOut;
	//the various local input / output vectors
	double **localX, **localY;
	//EMM file storage
	unsigned int *Pstart;
	unsigned long int *rowboundaries,		*colboundaries;
	unsigned long int *upscaled_rowboundaries,	*upscaled_colboundaries;
	unsigned long int *rowhierarchy,		*colhierarchy;
	unsigned long int *input_vector_lengths,	*output_vector_lengths;
	unsigned long int *localM,			*localN;
	unsigned long int *row_perm,			*col_perm;
	unsigned long int *output_vector_distribution,  *input_vector_distribution;
	std::vector< std::vector< unsigned long int > > local_rowboundaries;
	std::vector< std::vector< unsigned long int > > local_colboundaries;
	std::vector< std::vector< unsigned long int > > local_rowhierarchy;
	std::vector< std::vector< unsigned long int > > local_colhierarchy;
	unsigned long int **row2global, **col2global;
	unsigned long int **row2index, **col2index;
	unsigned long int **row2proc, **col2proc;
	unsigned long int **inv_row_perm;
	unsigned long int **inv_col_perm;
	unsigned long int m, n, nz; //inhereted now
	//unsigned long int nz;
	unsigned long int input_length, output_length;
	unsigned int P, Pref;

	MVData() {
		localX = localY = NULL;
		row_perm = col_perm = rowboundaries = colboundaries = rowhierarchy = colhierarchy = localM = localN =
		input_vector_lengths = output_vector_lengths = input_vector_distribution = output_vector_distribution = NULL;
		row2global = col2global = this->row2proc = this->col2proc = this->row2index = this->col2index = NULL;
		Pstart = NULL;
		inv_row_perm = new unsigned long int*[3];
		inv_col_perm = new unsigned long int*[3];
		for( int i=0; i<3; i++ )
			inv_row_perm[i] = inv_col_perm[i] = NULL;
	}

	~MVData() {
		for( unsigned long int s=0; s<P; s++ ) {
			if( localX != NULL     && localX[ s ]      != NULL ) delete [] localX    [ s ];
			if( localY != NULL     && localY[ s ]      != NULL ) delete [] localY    [ s ];
			if( row2global != NULL && row2global [ s ] != NULL ) delete [] row2global[ s ];
			if( col2global != NULL && col2global [ s ] != NULL ) delete [] col2global[ s ];
			if( row2index != NULL  && row2index [ s ]  != NULL ) delete [] row2index [ s ];
			if( col2index != NULL  && col2index [ s ]  != NULL ) delete [] col2index [ s ];
			if( row2proc != NULL   && row2proc [ s ]   != NULL ) delete [] row2proc  [ s ];
			if( col2proc != NULL   && col2proc [ s ]   != NULL ) delete [] col2proc  [ s ];
		}
		if( localX != NULL ) delete [] localX;
		if( localY != NULL ) delete [] localY;
		if( Pstart != NULL ) delete [] Pstart;
		if( row2global != NULL ) delete [] row2global;
		if( col2global != NULL ) delete [] col2global;
		if( row2proc != NULL ) delete [] row2proc;
		if( col2proc != NULL ) delete [] col2proc;
		if( row2index != NULL ) delete [] row2index;
		if( col2index != NULL ) delete [] col2index;
		if( output_vector_distribution != NULL ) delete [] output_vector_distribution;
		if( input_vector_distribution != NULL ) delete [] input_vector_distribution;
		if( rowboundaries != NULL ) delete [] rowboundaries;
		if( colboundaries != NULL ) delete [] colboundaries;
		if( upscaled_rowboundaries != NULL ) delete [] upscaled_rowboundaries;
		if( upscaled_colboundaries != NULL ) delete [] upscaled_colboundaries;
		if( rowhierarchy != NULL ) delete [] rowhierarchy;
		if( colhierarchy != NULL ) delete [] colhierarchy;
		if( localM != NULL ) delete [] localM;
		if( localN != NULL ) delete [] localN;
		if( output_vector_lengths != NULL ) delete [] output_vector_lengths;
		if( input_vector_lengths != NULL ) delete [] input_vector_lengths;
		if( row_perm != NULL ) delete [] row_perm;
		if( col_perm != NULL ) delete [] col_perm;
		for( int i=0; i<3; i++ ) {
			if( inv_row_perm[i] != NULL ) delete [] inv_row_perm[i];
			if( inv_col_perm[i] != NULL ) delete [] inv_col_perm[i];
		}
		delete [] inv_row_perm;
		delete [] inv_col_perm;
	}
};

MVData matrix;
double *checkX = NULL, *checkY = NULL, *checkY2 = NULL;

void printproc( const unsigned int proc ) { std::cout << "(" << proc << ")"; }

void parseVC( std::fstream &file, std::string array_name, unsigned long int **array, unsigned long int *lengths, char *chars, const unsigned long int P ) {
	std::string line;
	unsigned long int uli;
	std::cout << "Will read " << array_name << std::endl;
	if( file.peek() != '%' ) {
		std::cout << "Not at header line!" << std::endl;
		file >> line;
		std::cout << line << std::endl;
		file >> line;
		std::cout << line << std::endl;
		exit( 1 );
	}
	file.getline( chars, 500 ); //skip header
	std::cout << "Reading " << chars << std::endl;
	file >> uli; //skip number of vectors (=ret.P)
	assert( uli == P );
	unsigned long int ai;
	for( unsigned int k=0; k<P; k++ ) { //read the P vectors
		file >> ai; //get vector size
		lengths[ k ] = ai;
		array[ k ] = new unsigned long int[ ai ];
		for( unsigned int i=0; i<ai; i++ ) {
			file >> (array[ k ])[ i ];
			(array[ k ])[i]--; //corect base
		}
	}
}

unsigned long int *p_translate = NULL;

//used for gathering global timing info
double *global_time = NULL;

/** @param fn File to parse
    @param s  Local processor ID */
void parse( const std::string &fn, MVData &ret ) {
	unsigned long int temp;

	//open file
	std::cout << "Opening file " << fn << std::endl;
	std::fstream file( fn.c_str(), std::ios::in );
	if( !file.is_open() ) {
		std::cerr << "Error while opening file" << std::endl;
		exit( 1 );
	}
	
	//skip comments
	std::string line;
	char *chars = new char[500];
	file.getline( chars, 500 );
	const bool pattern_matrix = ( strstr( chars, "pattern" ) != NULL );
	if( pattern_matrix )
		std::cout << "Input is a pattern matrix" << std::endl;
	else
		std::cout << "Input is a general matrix" << std::endl;
	while( file.peek() == '%' ) {
		file.getline( chars, 500 );
		std::cout << "Skipping " << chars << std::endl;
	}

	//read parameters
	unsigned long int uli, nnz;
	file >> ret.m;
	file >> ret.n;
	file >> nnz; //full nnz
	ret.nz = nnz;
	file >> ret.Pref;
	ret.P = ret.Pref; //true P (<=Pref) will be set by main

	//read Pstart
	ret.Pstart = new unsigned int[ ret.P+1 ];
	for( unsigned int k=0; k<=ret.P; k++ ) {
		file >> ret.Pstart[ k ];
//		std::cout << "PStart: " << ret.Pstart[ k ] << std::endl;
	}

	//read nonzeroes
	std::cout << "Will read " << nnz << " nonzeroes " << std::endl;
	double av = 1.0; unsigned long int ai, aj;
	for( unsigned int k=0; k<ret.P; k++ ) {
		for( unsigned int i=ret.Pstart[k]; i<ret.Pstart[k+1]; i++ ) {
			file >> ai; file >> aj;
			if( !pattern_matrix )
				file >> av;
//			std::cout << "Read " << ai << ", " << aj << ", " << av << std::endl;
			ret.local.push_back( Triplet< double >( ai-1, aj-1, av ) ); //also correct for base
			ret.local[ ret.local.size() - 1 ].meta = k;
		}
	}
	file.getline( chars, 500 ); //skip EOL

	//skip remainder
	while( file.peek() != '%' )
		file.getline( chars, 500 );

	//read PAQ header
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	file.getline( chars, 500 ); //skip size headers

	//skip current matrix PAQ
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read Row boundaries header
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;

	//read row boundaries
	file >> temp;
	if( temp != 2 * ret.Pref ) {
		std::cerr << "Row-wise boundary/hierarchy is incomplete! Aborting..." << std::endl;
		exit( 1 );
	}
	ret.rowboundaries = new unsigned long int[ temp ];
	for( unsigned int i=0; i<temp; i++ ) {
		file >> (ret.rowboundaries[i]);
		ret.rowboundaries[i]--; //correct base
	}
	file.getline( chars, 500 ); //go past EOL
	std::cout << "Read " << temp << " row boundaries" << std::endl;
	std::cout << "Pref is " << ret.Pref << std::endl;

	//read Row hierarchy header
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;

	//read row hierarchy
	file >> temp;
	assert( 2 * ret.Pref - 1 == temp );
	ret.rowhierarchy = new unsigned long int[ temp ];
	for( unsigned int i=0; i<temp; i++ ) {
		file >> (ret.rowhierarchy[ i ]);
	}
	file.getline( chars, 500 ); //go past EOL
	std::cout << "Read " << temp << " row hierarchy entries" << std::endl;

	//read Col boundaries header
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;

	//skip Col boundaries
	file >> uli; //half should equal ret.Pref
	if( uli/2 != ret.Pref ) {
		std::cerr << "Error: Refinements in row direction does not equal that in column direction!" << std::endl;
		std::cerr << "(The column-wise boundary/hierarchy is incomplete)" << std::endl;
		exit( 1 );
	}
	ret.colboundaries = new unsigned long int[ uli ];
	for( unsigned int j=0; j<uli; j++ ) {
		file >> (ret.colboundaries[j]);
		ret.colboundaries[j]--; //correct base
	}
	file.getline( chars, 500 ); //go past EOL
	std::cout << "Read " << uli << " column boundaries" << std::endl;

	//read Col hierarchy header
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;

	//read column hierarchy
	file >> temp;
	assert( ret.Pref-1 == temp/2 );
	ret.colhierarchy = new unsigned long int[ temp ];
	for( unsigned int i=0; i<temp; i++ ) {
		file >> (ret.colhierarchy[ i ]);
	}
	file.getline( chars, 500 ); //go past EOL
	std::cout << "Read " << temp << " column hierarchy entries" << std::endl;

	//read Local-A header
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> ret.m;
	file >> ret.n;
	file >> nnz;
	file >> ret.P;

	//read Pstart
	for( unsigned int k=0; k<=ret.P; k++ ) {
		file >> ret.Pstart[ k ];
	}

	//read nonzeroes
	std::cout << "Will read " << nnz << " nonzeroes " << std::endl;
	for( unsigned int k=0; k<ret.P; k++ ) {
		for( unsigned int i=ret.Pstart[k]; i<ret.Pstart[k+1]; i++ ) {
			file >> ai; file >> aj;
			if( !pattern_matrix )
				file >> av;
			//NOTE: cannot use local-A because upscaling causes trouble for p>2, and even more trouble for p>4.
			//	so ignoring local view
		}
	}
	file.getline( chars, 500 ); //go past EOL

	//read row2global
	ret.row2global = new unsigned long int*[ ret.P ];
	ret.output_vector_lengths = new unsigned long[ ret.P ];
	ret.localM = new unsigned long[ ret.P ];
	parseVC( file, "row2global", ret.row2global, ret.localM, chars, ret.P );
	file.getline( chars, 500 ); //go past EOL

	std::cout << "Row dimension of local matrices: ";
	for( unsigned long int i=0; i<ret.P; i++ )
		std::cout << ret.localM[ i ] << " ";
	std::cout << std::endl;

	//read col2global
	ret.col2global = new unsigned long int*[ ret.P ];
	ret.input_vector_lengths = new unsigned long[ ret.P ];
	ret.localN = new unsigned long[ ret.P ];
	parseVC( file, "col2global", ret.col2global, ret.localN, chars, ret.P );
	file.getline( chars, 500 ); //go past EOL

	std::cout << "Column dimension of local matrices: ";
	for( unsigned long int i=0; i<ret.P; i++ )
		std::cout << ret.localN[ i ] << " ";
	std::cout << std::endl;

	//read Row-permutation and skip global-A if there
	file.getline( chars, 500 );
	if( strstr( chars, "Global-A" ) != NULL ) { //skip global-A if there
		std::cout << "Skipping " << chars << std::endl;
		while( file.peek() != '%' )
			file.getline( chars, 500 );
		file.getline( chars, 500 );
	}
	std::cout << "Reading " << chars << std::endl;
	file >> uli;
	ret.row_perm = new unsigned long int[ uli ];
	assert( uli == ret.m );
	for( unsigned int i=0; i<ret.m; i++ ) {
		file >> ret.row_perm[ i ];
		ret.row_perm[ i ]--;
	}
	file.getline( chars, 500 ); //skip EOL
	std::cout << "Last row permutation index read: " << ret.row_perm[ ret.m-1 ] << std::endl;

	//invert row permutation
	unsigned long int *dummy = ret.row_perm;
	ret.row_perm = new unsigned long int[ ret.m ];
	for( unsigned int k=0; k<ret.m; k++ ) ret.row_perm[ dummy[ k ] ] = k;
	delete [] dummy;

	//read Column-permutation
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> uli;
	assert( uli == ret.n );
	ret.col_perm = new unsigned long int[ uli ];
	for( unsigned int i=0; i<ret.n; i++ ) {
		file >> ret.col_perm[ i ];
		ret.col_perm[ i ]--;
	}
	file.getline( chars, 500 ); //skip EOL
	std::cout << "Last column permutation index read: " << ret.col_perm[ ret.n-1 ] << std::endl;

	//invert column permutation
	dummy = ret.col_perm;
	ret.col_perm = new unsigned long int[ ret.n ];
	for( unsigned int k=0; k<ret.n; k++ ) ret.col_perm[ dummy[ k ] ] = k;
	delete [] dummy;

	//read Input-vector distribution
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> uli; //should be ret.n
	ret.input_vector_distribution = new long unsigned int[ uli ];
	file >> uli; //should be ret.P
	for( unsigned int i=0; i<ret.n; i++ ) {
		file >> uli;
		file >> ret.input_vector_distribution[ uli-1 ];
		ret.input_vector_distribution[ uli-1 ]--; //correct base
	}
	file.getline( chars, 500 ); //skip EOL

	//read Output-vector distribution
	if( file.peek() != '%' ) {
		std::cout << "Not at header line!" << std::endl;
		file.getline( chars, 500 );
		std::cout << chars << std::endl;
		file.getline( chars, 500 );
		std::cout << chars << std::endl;
		exit( 1 );
	}
	file.getline( chars, 500 );
	std::cout << "Reading output_vector_distribution " << chars << std::endl;
	file >> uli; //should be ret.m
	ret.output_vector_distribution = new long unsigned int[ uli ];
	file >> uli; //should be ret.P
	for( unsigned int i=0; i<ret.m; i++ ) {
		file >> uli;
		file >> ret.output_vector_distribution[ uli-1 ];
		ret.output_vector_distribution[ uli-1 ]--; //correct base
	}
	file.getline( chars, 500 ); //skip EOL

	//read outputvectorlengths
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> uli; //should be ret.P
	assert( uli == ret.P );
	for( unsigned int i=0; i<ret.P; i++ ) {
		file >> uli;
		ret.output_vector_lengths[ i ] = uli;
	}
	file.getline( chars, 500 ); //skip EOL

	//read row2proc
	ret.row2proc = new unsigned long int*[ ret.P ];
	parseVC( file, "LocalRow2Processor", ret.row2proc, ret.localM, chars, ret.P );
	file.getline( chars, 500 ); //go past EOL

	//read row2index
	ret.row2index = new unsigned long int*[ ret.P ];
	parseVC( file, "LocalRow2Index", ret.row2index, ret.localN, chars, ret.P );
	file.getline( chars, 500 ); //go past EOL

	//read outputvectorlengths
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> uli; //should be ret.P
	assert( uli == ret.P );
	for( unsigned int i=0; i<ret.P; i++ ) {
		file >> uli;
		ret.input_vector_lengths[ i ] = uli;
	}
	file.getline( chars, 500 ); //skip EOL

	//read col2proc
	ret.col2proc = new unsigned long int*[ ret.P ];
	parseVC( file, "LocalCol2Processor", ret.col2proc, ret.localN, chars, ret.P );
	file.getline( chars, 500 ); //go past EOL

	//read col2index
	ret.col2index = new unsigned long int*[ ret.P ];
	parseVC( file, "LocalCol2Index", ret.col2index, ret.localN, chars, ret.P );
	file.getline( chars, 500 ); //go past EOL

	//done
	file.getline( chars, 500 ); //trigger EOF
	if( !file.eof() ) {
		std::cout << "Warning: input file not at end!" << std::endl;
		std::cout << "Next lines are: " << chars << std::endl;
		file >> line;
		std::cout << line;
	} else {
		std::cout << "End of file; parse OK." << std::endl;
	}
	delete [] chars;
}

std::string fn;

#ifdef SPINLOCK
unsigned char *condition   = NULL;
#endif
pthread_t      *threads    = NULL;
pthread_mutex_t sync_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  sync_cond  = PTHREAD_COND_INITIALIZER;
unsigned int sync_val = 0;
unsigned int P = 0;
#ifdef _MCDMV_HYBRID
unsigned int R = 0;
std::vector< std::vector< unsigned short int > > sub_p_translate;
#endif

void bsp_sync( const unsigned int bsp_id ) {
#ifdef SPINLOCK
	const unsigned char sync_number = ++(condition[ bsp_id ]);
	bool go = false;
	while( !go  ) {
		go = true;
		for( unsigned int s = 0; s < P; ++s ) {
			if( condition[ s ] != sync_number && condition[ s ] != (unsigned char)(sync_number + 1) ) {
				go = false;
				break;
			}
		}
	}
#else
        pthread_mutex_lock( &sync_mutex );
        sync_val++;
        if( sync_val < P ) pthread_cond_wait( &sync_cond, &sync_mutex );
        else {
                sync_val = 0;
                pthread_cond_broadcast( &sync_cond );
        }
        pthread_mutex_unlock( &sync_mutex );
#endif
}

/** Sorts a permutation range according to a given distribution (so that it become contiguous) */
/*void applySubPermutation( const unsigned long int lo, const unsigned long int hi, const unsigned long int P,
				unsigned long int * const vector_distribution,
				unsigned long int * const perm ) {
	//find row sub-permutation
	std::vector< unsigned long int > subperm, Pcount, Pind;
	//init counting sort
	subperm.resize( hi - lo, ULONG_MAX );
	Pcount .resize( P, 0ul );
	for( unsigned long int i = lo; i < hi; ++i )
		Pcount[ vector_distribution[ i ] ]++;
	Pind.push_back( 0ul );
	//counting sort phase 1
	for( unsigned long int s = 1; s < P; ++s )
		Pcount[ s ] += Pcount[ s - 1 ];
	for( unsigned long int s = 0; s < P; ++s )
		Pind.push_back( Pcount[ s ] );
	//counting sort phase 2
	for( unsigned long int i = lo; i < hi; ++i )
		subperm[ i-lo ] = Pind[ vector_distribution[ i ] ]++;
	//check permutation
	for( unsigned long int i = lo; i < hi; ++i )
		assert( subperm[ i-lo ] != ULONG_MAX );
	for( unsigned long int s = 0; s < P; ++s )
		assert( Pind[ s ] == Pcount[ s ] );
	//apply permutation on perm
	std::vector< unsigned long int > inverse, cache;
	inverse.resize( hi-lo, ULONG_MAX );
	for( unsigned long int i = lo; i < hi; ++i )
		inverse[ subperm[ i - lo ] ] = i;
	for( unsigned long int i = lo; i < hi; ++i )
		assert( inverse[ i - lo ] != ULONG_MAX );
	for( unsigned long int i = lo; i < hi; ++i )
		cache.push_back( perm[ inverse[ i - lo ] ] );
	for( unsigned long int i = lo; i < hi; ++i ) 
		perm[ i ] = cache[ i - lo ];
	//apply permutation on distribution
	cache.clear();
	for( unsigned long int i = lo; i < hi; ++i )
		cache.push_back( vector_distribution[ inverse[ i - lo ] ] );
	for( unsigned long int i = lo; i < hi; ++i )
		vector_distribution[ i ] = cache[ i - lo ];
}*/

//The below is an alternative for the above, and is applied after upscaling. This is parallelisable.
//the goal is to permute LOCAL rows and columns so that in global view everything is sorted by
//input or output vector distribution.
void postPermuteKernel( const unsigned long int P, const unsigned long int s,
				const unsigned int long lo,
				const unsigned long int size,
				const unsigned long int * const distribution,
				unsigned long int * const toGlobal,
				std::vector< std::map< unsigned long int, unsigned long int > > &global2local,
				std::vector< Triplet< double > > &rmatrix,
				std::vector< Triplet< double > > &lmatrix,
				const bool rowwise ) {
	//input check
	assert( toGlobal != NULL );
	assert( distribution != NULL );
	assert( global2local.size() == P );
	//project local matrices indices to global indices
	for( unsigned long int k = 0; k < rmatrix.size(); ++k )
		if( rowwise ) {
			const unsigned long int local = rmatrix[ k ].i();
			rmatrix[ k ].setRowPosition   ( toGlobal[ local ] );
		} else {
			const unsigned long int local = rmatrix[ k ].j();
			rmatrix[ k ].setColumnPosition( toGlobal[ local ] );
		}
	for( unsigned long int k = 0; k < lmatrix.size(); ++k )
		if( rowwise ) {
			const unsigned long int local = lmatrix[ k ].i();
			lmatrix[ k ].setRowPosition   ( toGlobal[ local ] );
		} else {
			const unsigned long int local = lmatrix[ k ].j();
			lmatrix[ k ].setColumnPosition( toGlobal[ local ] );
		}
	//do stable sort on range 0 to size[ s ]
	std::vector< unsigned long int > permute, count, id;
	count.resize( P, 0ul );
	for( unsigned long int k = 0; k < size; ++k )
		count[ distribution[ toGlobal[ lo + k ] ] ]++;
	id.push_back( 0 );
	for( unsigned long int k = 0; k < P; ++k )
		id.push_back( count[ k ] );
	//make cumulative
	for( unsigned long int k = 1; k <=P; ++k )
		id[ k ] += id[ k - 1];
	for( unsigned long int k = 0; k < size; ++k )
#ifdef NDEBUG
		permute.push_back( lo + id[ distribution[ toGlobal[ lo + k ] ] ]++ );
#else
		permute.push_back( id[ distribution[ toGlobal[ lo + k ] ] ]++ ); 
	//check permutation (entire range should be covered)
	std::vector< bool > check;
	check.resize( size, false );
	for( unsigned long int k = 0; k < size; ++k ) {
		assert( check[ permute[ k ] ] == false );
		check[ permute[ k ] ] = true;
	}
	for( unsigned long int k = 0; k < size; ++k )
		assert( check[ k ] );
	//let permute correspond to the real local index
	for( unsigned long int k = 0; k < size; ++k )
		permute[ k ] += lo;
#endif
	//do permute
	std::vector< unsigned long int > replace2Global;
	for( unsigned long int k = 0; k < size; ++k ) {
		const unsigned long int replacement = permute[ k ];
		const unsigned long int global = toGlobal[ replacement ];
		replace2Global.push_back( global );
		//update globalToLocal maps
		const std::map< unsigned long int, unsigned long int >::iterator it = 
			global2local[ s ].find( global );
		if( it != global2local[s].end() ) {
			assert( it->second == permute[ k ] );
			it->second = lo + k; //update with new location
			assert( global2local[ s ].find( global )->second == lo + k );
		}
	}
	//copy back into actual arrays
	for( unsigned long int k = 0; k < size; ++k ) {
		toGlobal[ lo + k ] = replace2Global[ k ];
		assert( global2local[ s ][ toGlobal[ lo + k ] ] == lo + k );
	}
	//transform global matrix entries back to (permuted) local ones
	for( unsigned long int k = 0; k < rmatrix.size(); ++k )
		if( rowwise ) {
			const unsigned long int global = rmatrix[ k ].i();
			rmatrix[ k ].setRowPosition   ( global2local[ s ][ global] );
		} else {
			const unsigned long int global = rmatrix[ k ].j();
			rmatrix[ k ].setColumnPosition( global2local[ s ][ global ] );
		}
	for( unsigned long int k = 0; k < lmatrix.size(); ++k )
		if( rowwise ) {
			const unsigned long int global = lmatrix[ k ].i();
			lmatrix[ k ].setRowPosition   ( global2local[ s ][ global ] );
		} else {
			const unsigned long int global = lmatrix[ k ].j();
			lmatrix[ k ].setColumnPosition( global2local[ s ][ global ] );
		}
}


void prepare( MVData &matrix ) {
	//upscale boundaries
	std::cout << "Upscaling from " << matrix.Pref << " to " << matrix.P << std::endl;

	//deduce and print p-translate
	unsigned long int *proc2proc = new unsigned long int[ matrix.Pref ];
	for( unsigned long int i=0; i<matrix.Pref; i++ ) proc2proc[ i ] = i * matrix.P / matrix.Pref; //scale block-by-block
	//NOTE: the above transation does assume the original blocks are ordered (top-left to bottom-right of the matrix)!
	std::cout << "Triplet p-translate: ";
	for( unsigned int i=0; i<matrix.Pref; i++ )
		std::cout << proc2proc[ i ] << " ";
	std::cout << std::endl;

	//optimising vector assignments through extra permutations
	/* Applying permutations at this stage causes errors in Upscaler,
	 * so turned off for now. Will work on post-upscaled data instead.
	std::cout << "Refining separator permutations" << std::endl;
	for( unsigned long int c = 0; c < 2*matrix.Pref - 1; ++c ) {
		if( c % 2 == 0 ) continue; //disregard pure blocks
		applySubPermutation( matrix.rowboundaries[ c ], matrix.rowboundaries[ c+1 ], matrix.Pref, matrix.output_vector_distribution, matrix.row_perm );
		applySubPermutation( matrix.colboundaries[ c ], matrix.colboundaries[ c+1 ], matrix.Pref,  matrix.input_vector_distribution, matrix.col_perm );
	} */

	//translate to upscaled processor assignments
	std::cout << "Translating vector assignments" << std::endl;
	for( unsigned long int i=0; i<matrix.m; i++ )
		matrix.output_vector_distribution[ i ] = proc2proc[ matrix.output_vector_distribution[ i ] ];
	for( unsigned long int j=0; j<matrix.n; j++ )
		matrix.input_vector_distribution[ j ] = proc2proc[ matrix.input_vector_distribution[ j ] ];

	std::cout << "Initialising upscaler\n" << std::endl;
	std::vector< unsigned long int > stdrb, stdcb, stdrh, stdch;
	for( unsigned long int i=0; i<2*matrix.Pref - 1; i++ ) {
		stdrb.push_back( matrix.rowboundaries[ i ] );
		stdcb.push_back( matrix.colboundaries[ i ] );
		stdrh.push_back( matrix.rowhierarchy[ i ] );
		stdch.push_back( matrix.colhierarchy[ i ] );
		//simultaneously allocate local objects
		if( i < matrix.P ) {
			matrix.main_matrices.push_back( std::vector< Triplet< double > >() );
			matrix.remote_separators.push_back( std::vector< Triplet< double > >() );
			matrix.local_rowhierarchy.push_back( std::vector< unsigned long int >() );
			matrix.local_colhierarchy.push_back( std::vector< unsigned long int >() );
			matrix.local_rowboundaries.push_back( std::vector< unsigned long int >() );
			matrix.local_colboundaries.push_back( std::vector< unsigned long int >() );
		}
	}
	stdrb.push_back( matrix.rowboundaries[ 2*matrix.Pref - 1 ] );
	stdcb.push_back( matrix.colboundaries[ 2*matrix.Pref - 1 ] );
	
	//delete unused file input
	for( unsigned long int s=0; s<matrix.Pref; ++s ) {
		delete [] matrix.row2global[ s ];
		delete [] matrix.col2global[ s ];
		delete [] matrix.row2index [ s ];
		delete [] matrix.col2index [ s ];
		delete [] matrix.row2proc  [ s ];
		delete [] matrix.col2proc  [ s ];
		matrix.row2global[ s ] = matrix.col2global[ s ] = matrix.row2index[ s ] = matrix.col2index[ s ] = matrix.row2proc[ s ] = matrix.col2proc[ s ] = NULL;
	}
	
	//initialise inverse maps
	std::vector< std::map< unsigned long int, unsigned long int > > rowGlobalToLocal, colGlobalToLocal;
	for( unsigned long int s=0; s<matrix.P; ++ s ) {
		rowGlobalToLocal.push_back( std::map< unsigned long int, unsigned long int >() ) ;
		colGlobalToLocal.push_back( std::map< unsigned long int, unsigned long int >() );
	}

	//do upscaling, processor by processor
	Upscaler< double > upscaler( matrix.local, matrix.P, stdrh, stdch, stdrb, stdcb, matrix.row_perm, matrix.col_perm, proc2proc );
	for( unsigned long int s=0; s<matrix.P; ++s ) {
		std::vector< unsigned long int > tempRow2Glob, tempCol2Glob;
		upscaler.getSubTree( s, matrix.main_matrices[ s ], matrix.remote_separators[ s ],
			matrix.local_rowhierarchy[ s ], matrix.local_rowboundaries[ s ],
			matrix.local_colboundaries[ s ], tempRow2Glob, tempCol2Glob,
			rowGlobalToLocal[ s ], colGlobalToLocal[ s ], matrix.P, matrix.Pref );

		//hierarchy should be 1-based for use with BlockOrderer classes
		for( unsigned long int i=0; i < matrix.local_rowhierarchy[ s ].size(); ++i )
			if( matrix.local_rowhierarchy[ s ][ i ] > 0 )
				matrix.local_rowhierarchy[ s ][ i ]++;

		//copy row hierarchy to column hierarchy as well (assume symmetric hierarchies)
		matrix.local_colhierarchy[ s ] = matrix.local_rowhierarchy[ s ];

		//copy row/col to global data to matrix
		matrix.row2global[ s ] = new unsigned long int[ tempRow2Glob.size() ];
		for( unsigned long int i=0; i<tempRow2Glob.size(); i++ )
			matrix.row2global[ s ][ i ] = tempRow2Glob[ i ];
		matrix.col2global[ s ] = new unsigned long int[ tempCol2Glob.size() ];
		for( unsigned long int j=0; j<tempCol2Glob.size(); j++ )
			matrix.col2global[ s ][ j ] = tempCol2Glob[ j ];

		//remember local sizes
		matrix.localM[ s ] = tempRow2Glob.size();
		matrix.localN[ s ] = tempCol2Glob.size();

		//hierarchy outs:
		std::cout << "\nResults for processor " << s << ": " << std::endl;
		std::cout << "Row boundaries:    [ ";
		for( unsigned long int i=0; i<matrix.local_rowboundaries[ s ].size(); i++ )
			std::cout << matrix.local_rowboundaries[ s ][ i ] << " ";
		std::cout << "]" << std::endl;
		std::cout << "Column boundaries: [ ";
		for( unsigned long int i=0; i<matrix.local_colboundaries[ s ].size(); i++ )
			std::cout << matrix.local_colboundaries[ s ][ i ] << " ";
		std::cout << "]" << std::endl;
		std::cout << "Hierarchy:         [ ";
		for( unsigned long int i=0; i<matrix.local_rowhierarchy[ s ].size(); i++ )
			std::cout << matrix.local_rowhierarchy[ s ][ i ] << " ";
		std::cout << "]\n" << std::endl;
	}

	//if we do permute here, the data to modify are:
	//row2global, col2global,
	//rowGlobalToLocal, colGlobalToLocal,
	//input_vector_distribution, output_vector_distribution
	//
	//global/local indicators are:
	//local_rowboundaries, local_colboundaries
	//
	//data to use are:
	//localM, localN
	//
	//now do permute, if requested:
#ifdef _MCDMV_VECTOR_SORT
	for( unsigned long int s = 0; s < P; ++s ) {
		postPermuteKernel( P, s, 0, matrix.local_rowboundaries[ s ][ 0 ],
			matrix.output_vector_distribution, matrix.row2global[ s ],
			rowGlobalToLocal, matrix.remote_separators[ s ],
			matrix.main_matrices[ s ], true );
		postPermuteKernel( P, s, matrix.local_rowboundaries[ s ].back(),
			( matrix.localM[ s ] - matrix.local_rowboundaries[ s ].back() ),
			matrix.output_vector_distribution, matrix.row2global[ s ],
			rowGlobalToLocal, matrix.remote_separators[ s ],
			matrix.main_matrices[ s ], true );
		postPermuteKernel( P, s, 0, matrix.local_colboundaries[ s ][ 0 ],
			matrix.input_vector_distribution, matrix.col2global[ s ],
			colGlobalToLocal, matrix.remote_separators[ s ],
			matrix.main_matrices[ s ], false );
		postPermuteKernel( P, s, matrix.local_colboundaries[ s ].back(),
			( matrix.localN[ s ] - matrix.local_colboundaries[ s ].back() ),
			matrix.input_vector_distribution, matrix.col2global[ s ],
			colGlobalToLocal, matrix.remote_separators[ s ],
			matrix.main_matrices[ s ], false );
	}
#endif

	//allocate fanIn/fanOut objects for all processors
	for( unsigned long int s=0; s<matrix.P; s++ ) {
		matrix.fanIn .push_back( fanInContainerType() );
		matrix.fanOut.push_back( fanOutContainerType() );
	}

	//now rebuild row/col2index, row/col2proc
	for( unsigned long int s=0; s<matrix.P; s++ ) {
		//build fanIn for processor s
		for( unsigned long int j=0; j<matrix.localN[s]; ++j ) {
			//do not process local area;
			//that is, jump over local area
			if( j == matrix.local_colboundaries[ s ][ 0 ] ) {
				j = matrix.local_colboundaries[ s ].back();
				//be sure we are still in range
				if( j == matrix.localN[ s ] )
					break;
			}
			//get global column index from local one
			const unsigned long int globalJ = matrix.col2global[ s ][ j ];
			//get the owner of that global column index
			const unsigned long int remoteP = matrix.input_vector_distribution[ globalJ ];
			//get the corresponding local index
			const unsigned long int remoteI = colGlobalToLocal[ remoteP ][ globalJ ];
			assert( globalJ < matrix.n );
			assert( remoteP < matrix.P );
			assert( remoteI < matrix.localN[ remoteP ] );
			//if we are the remote processor, skip
			if( remoteP == s ) assert( j == remoteI );
			if( remoteP == s ) continue;
			//in range and non-local processor; fan-in required
#ifdef _MCDMV_COMPRESSED_COMMUNICATION
			//This loop differs from the uncompressed strategy (below). That is because
			//access must be consecutive in both the source and destination vector (otherwise a quadlet 
			//cannot encode the required communication).
			
			//find sequence length
			unsigned long int inc = 1;
			while( j + inc < matrix.localN[ s ] &&
				matrix.input_vector_distribution[ matrix.col2global[ s ][ j + inc ] ] == remoteP &&
				colGlobalToLocal[ remoteP ][ matrix.col2global[ s ][ j + inc ] ] == remoteI + inc ) ++inc;
			//found range, store communication
			matrix.fanIn[ s ].push_back( fanQuadlet( remoteP, j, remoteI, inc ) );
			j += inc - 1; //one less since the for-loop will increment one
#else
			matrix.fanIn[ s ][ j ] = fanPair( remoteP, remoteI );
#endif
		}
		//build fanOut for remote processors
		for( unsigned long int i=0; i<matrix.localM[s]; ++i ) {
			//do not process local area
			if( i == matrix.local_rowboundaries[ s ][ 0 ] ) {
				i = matrix.local_rowboundaries[ s ].back();
				if( i == matrix.localM[ s ] )
					break;
			}
			const unsigned long int globalI = matrix.row2global[ s ][ i ];
			const unsigned long int remoteP = matrix.output_vector_distribution[ globalI ];
			const unsigned long int remoteI = rowGlobalToLocal[ remoteP ][ globalI ];
			assert( globalI < matrix.m );
			assert( remoteP < matrix.P );
			assert( remoteI < matrix.localM[ remoteP ] );
			//if we are the remote processor, skip
			if( remoteP == s ) continue;
#ifdef _MCDMV_COMPRESSED_COMMUNICATION
			//This loop differs from the uncompressed strategy (below). That is because
			//access must be consecutive in both the source and destination vector (otherwise a quadlet 
			//cannot encode the required communication).
			
			//this is the start of a remote section
			//find length, make sure destination also remains contiguous
			unsigned long int inc = 1;
			while( i + inc < matrix.localM[s] &&
				matrix.output_vector_distribution[ matrix.row2global[ s ][ i + inc ] ] == remoteP &&
				rowGlobalToLocal[ remoteP ][ matrix.row2global[ s ][ i + inc ] ] == remoteI + inc ) ++inc;
			matrix.fanOut[ remoteP ].push_back( fanQuadlet( s, remoteI, i, inc ) );
			i += inc - 1; //one less since the for-loop will increment one
#else
			if( remoteP == s ) continue;
			//remote processor should get the value we have
			fanPair pair( s, i );
			matrix.fanOut[ remoteP ].insert( std::pair< unsigned long int, fanPair >( remoteI, pair ) );
#endif
		}
	}

#ifdef _MCDMV_COMPRESSED_COMMUNICATION
	//assertion time!
	std::vector< unsigned long int > numFanOut;
	unsigned long int totalCompressed, totalUncompressed; //for compression statistics
	totalCompressed = totalUncompressed = 0ul;
	numFanOut.resize( P, 0ul );
	for( unsigned long int s = 0; s < matrix.P; ++s ) {
		unsigned long int numFanIn = 0;
		for( unsigned long int i = 0; i < matrix.localN[ s ]; ++i ) {
			if( matrix.input_vector_distribution[ matrix.col2global[ s ][ i ] ] == s ) continue;
			numFanIn++;
		}
		totalUncompressed += numFanIn;
		totalCompressed   += matrix.fanIn[ s ].size();
		for( unsigned long int i = 0; i < matrix.fanIn[ s ].size(); ++i ) {
			numFanIn -= matrix.fanIn[ s ][ i ].length;
			//each request should be in locally in range
			assert( matrix.fanIn[ s ][ i ].remoteStart +
				matrix.fanIn[ s ][ i ].length - 1 < 
				matrix.localN[ matrix.fanIn[ s ][ i ].remoteP ] );
			//and also remotely
			assert( matrix.fanIn[ s ][ i ].localStart +
				matrix.fanIn[ s ][ i ].length - 1 < 
				matrix.localN[ s ] );
			for( unsigned long int k = 0; k < matrix.fanIn[ s ][ i ].length; ++k ) {
				//check if local2global is the same as remote2global
				assert( matrix.col2global[ s ][ matrix.fanIn[ s ][ i ].localStart + k ] ==
					matrix.col2global[ matrix.fanIn[ s ][ i ].remoteP ][ matrix.fanIn[ s ][ i ].remoteStart + k ] );
				//check if remoteP is really the owner of the global indices
				assert( matrix.input_vector_distribution[ matrix.col2global[ s ][ matrix.fanIn[ s ][ i ].localStart + k ] ] == matrix.fanIn[ s ][ i ].remoteP );
			}
		}
		assert( numFanIn == 0 );
		for( unsigned long int i = 0; i < matrix.localM[ s ]; ++i ) {
			if( matrix.output_vector_distribution[ matrix.row2global[ s ][ i ] ] == s ) continue;
			numFanOut[ s ]++;
		}
		totalUncompressed += numFanOut[ s ];
		totalCompressed   += matrix.fanOut[ s ].size();
		for( unsigned long int i = 0; i < matrix.fanOut[ s ].size(); ++i ) {
			numFanOut[ matrix.fanOut[ s ][ i ].remoteP ] -= matrix.fanOut[ s ][ i ].length;
			assert( matrix.fanOut[ s ][ i ].remoteStart +
				matrix.fanOut[ s ][ i ].length - 1 <
				matrix.localM[ matrix.fanOut[ s ][ i ].remoteP ] );
			assert( matrix.fanOut[ s ][ i ].localStart +
				matrix.fanOut[ s ][ i ].length - 1 <
				matrix.localM[ s ] );
			for( unsigned long int k = 0; k < matrix.fanOut[ s ][ i ].length; ++k ) {
				assert( matrix.row2global[ s ][ matrix.fanOut[ s ][ i ].localStart + k ] ==
					matrix.row2global[ matrix.fanOut[ s ][ i ].remoteP ][ matrix.fanOut[ s ][ i ].remoteStart + k ] );
				assert( matrix.output_vector_distribution[ matrix.row2global[ s ][ matrix.fanOut[ s ][ i ].localStart + k ] ] == s );
			}
		}
	}
	for( unsigned long int s = 0; s < matrix.P; ++s )
		assert( numFanOut[ s ] == 0 );
	std::cout << "Communication compression is " << (totalUncompressed / ((double)totalCompressed) * 100.0 ) << "%" << std::endl;
#endif

	//allocate arrays needed during parallel SpMV
	matrix.localX = new double*[ matrix.P ];
	matrix.localY = new double*[ matrix.P ];

#ifndef NDEBUG
	//happy asserting
	std::cout << "Sanity checks...";
	unsigned long int nonzeroes = 0;
	for( unsigned long int s=0; s<matrix.P; s++ ) {
		nonzeroes += matrix.main_matrices[ s ].size();
		nonzeroes += matrix.remote_separators[ s ].size();
		for( unsigned long int k=0; k<matrix.main_matrices[ s ].size(); k++ ) {
			const unsigned long int local_i = matrix.main_matrices[ s ][ k ].i();
			const unsigned long int local_j = matrix.main_matrices[ s ][ k ].j();
			assert( matrix.main_matrices[ s ][ k ].meta == s );
			assert( local_i < matrix.localM[ s ] );
			assert( local_j < matrix.localN[ s ] );
			assert( matrix.row2global[ s ][ local_i ] < matrix.m );
			assert( matrix.col2global[ s ][ local_j ] < matrix.n );
		}
		for( unsigned long int k=0; k<matrix.remote_separators[ s ].size(); k++ ) {
			const unsigned long int local_i = matrix.remote_separators[ s ][ k ].i();
			const unsigned long int local_j = matrix.remote_separators[ s ][ k ].j();
			assert( matrix.remote_separators[ s ][ k ].meta == s );
			assert( local_i < matrix.localM[ s ] );
			assert( local_j < matrix.localN[ s ] );
			assert( matrix.row2global[ s ][ local_i ] < matrix.m );
			assert( matrix.col2global[ s ][ local_j ] < matrix.n );
		}
	}
	assert( nonzeroes == matrix.local.size() );

	//check if all global rows are mapped
	std::vector< unsigned long int > check;
	for( unsigned long int i = 0; i < matrix.m; ++i )
		check.push_back( 1 ); //default a row is OK
	for( unsigned long int k = 0; k < matrix.local.size(); ++k )
		check[ matrix.row_perm[ matrix.local[ k ].i() ] ] = 0; //nonempty rows temporarily flagged not-OK
	for( unsigned long int s = 0; s< matrix.P; s++ )
		for( unsigned long int i = 0; i < matrix.localM[ s ]; i++ ) {
			const unsigned long int global = matrix.row2global[ s ][ i ];
			if( matrix.output_vector_distribution[ global ] == s ) {
				assert( check[ global ] == 0 );
				check[ global ] = 2; //flag OK if a parallel process is mapping this
			}
		}
	for( unsigned long int i = 0; i < matrix.m; ++i )
		assert( check[ i ] > 0 ); //everything should be OK (1 if empty, 2 if nonempty and mapped)

	//check if all global columns are mapped, see above
	check.clear();
	for( unsigned long int i = 0; i < matrix.n; ++i )
		check.push_back( 1 );
	for( unsigned long int k = 0; k < matrix.local.size(); ++k )
		check[ matrix.col_perm[ matrix.local[ k ].j() ] ] = 0;
	for( unsigned long int s = 0; s< matrix.P; s++ )
		for( unsigned long int i = 0; i < matrix.localN[ s ]; i++ ) {
			const unsigned long int global = matrix.col2global[ s ][ i ];
			if( matrix.input_vector_distribution[ global ] == s ) {
				assert( check[ global ] == 0 );
				check[ global ] = 2;
			}
		}
	for( unsigned long int i = 0; i < matrix.n; ++i ) assert( check[ i ] > 0 );

	//check if all local rows are used
	for( unsigned long int s = 0; s< matrix.P; s++ ) {
		check.clear();
		for( unsigned long int i = 0; i < matrix.localM[ s ]; ++i )
			check.push_back( 0 );
		for( unsigned long int k = 0; k < matrix.main_matrices[ s ].size(); ++k )
			check[ matrix.main_matrices[ s ][ k ].i() ] = 1;
		for( unsigned long int k = 0; k < matrix.remote_separators[ s ].size(); ++k )
			check[ matrix.remote_separators[ s ][ k ].i() ] = 1;
		for( unsigned long int i = 0; i < matrix.localM[ s ]; ++i )
			assert( check[ i ] );
	}

	//check if all local columns are used
	for( unsigned long int s = 0; s< matrix.P; s++ ) {
		check.clear();
		for( unsigned long int i = 0; i < matrix.localN[ s ]; ++i )
			check.push_back( 0 );
		for( unsigned long int k = 0; k < matrix.main_matrices[ s ].size(); ++k )
			check[ matrix.main_matrices[ s ][ k ].j() ] = 1;
		for( unsigned long int k = 0; k < matrix.remote_separators[ s ].size(); ++k )
			check[ matrix.remote_separators[ s ][ k ].j() ] = 1;
		for( unsigned long int i = 0; i < matrix.localN[ s ]; ++i )
			assert( check[ i ] );
	}	

	std::cout << " done." << std::endl;

#endif

	delete [] proc2proc;
	return;
}

void spmv( const unsigned int s, const MVData &data,
		Matrix< double > *A, Matrix< double > *S,
		double *x, double *y ) {

	assert( s < data.P );
	assert( data.local_colboundaries[ s ][ 0 ] < data.localN[ s ] );

	fanInContainerType::const_iterator it = data.fanIn[ s ].begin();
	for( ; it != data.fanIn [ s ].end(); ++it ) {
#ifdef _MCDMV_COMPRESSED_COMMUNICATION
		for( unsigned long int i=0; i<it->length; ++i )
			x[ it->localStart + i ] = data.localX[ it->remoteP ][ it->remoteStart + i ];
#else
	//fan-in (could be done in another thread)
		const unsigned long int remoteP = it->second.P;
		const unsigned long int remoteI = it->second.i;
		assert( it->first < data.localN[ s ] );
		assert( remoteP < data.P );
		assert( remoteI < data.localN[ remoteP ] );
		x[ it->first ] = data.localX[ remoteP ][ remoteI ];
#endif
	}

	//local multiply
	A->zax( x, y );
	S->zax( x, y );
	bsp_sync( s );
	//fan-out
	fanOutContainerType::const_iterator m_it = data.fanOut[ s ].begin();
	for( ; m_it != data.fanOut[ s ].end(); ++m_it ) {
	//each element of the entire buffered region can require multiple gets
	//hence the use of a multimap (a sorted vector storing triplets is possible too,
	//however memory usage should be a bit better with a multimap, at the cost of an
	//increased build time though)
#ifdef _MCDMV_COMPRESSED_COMMUNICATION
		for( unsigned long int i=0; i<m_it->length; ++i ) {
			assert( m_it->localStart + i < data.localM[ s ] );
			assert( m_it->remoteP < data.P );
			assert( m_it->remoteStart + i < data.localM[ m_it->remoteP ] );
			y[ m_it->localStart + i ] += data.localY[ m_it->remoteP ][ m_it->remoteStart + i ];
		}
#else
		const unsigned long int remoteP = m_it->second.P;
		const unsigned long int remoteI = m_it->second.i;
		assert( m_it->first < data.localM[ s ] );
		assert( remoteP < data.P );
		assert( remoteI < data.localM[ remoteP ] );
		//execute get
		y[ m_it->first ] += data.localY[ remoteP ][ remoteI ];
#endif
	}
	//no sync, since if this is done, the local data has arrived and local computation
	//can immediately continue.
}

void *parallel( void *threadid ) {
	struct timespec start, stop;
	time_t t0, t1;
	double proc_time;
	unsigned int s = (unsigned int) ((long)threadid);

#ifndef _NO_LIBNUMA
	//set kernel to local thread allocation if it wasn't already the case
	numa_set_localalloc();
#endif

	if( s == 0 ) std::cout << "Starting up " << (REPEAT_EXP) << " syncs for benchmarking" << std::endl;
	clock_gettime( CLOCK_MONOTONIC, &start );
	for( unsigned int i=0; i<REPEAT_EXP; i++ ) bsp_sync( s );
	clock_gettime( CLOCK_MONOTONIC, &stop );
	proc_time = (stop.tv_sec-start.tv_sec) * 1000;
	proc_time+= (stop.tv_nsec-start.tv_nsec)/1000000.0;
	if( s == 0 ) { printproc(s); std::cout << ": local time estimate for l = " << (proc_time/((double)REPEAT_EXP)) << " ms." << std::endl; }
	bsp_sync( s );
	clock_gettime( CLOCK_MONOTONIC, &start );
	printproc(s); std::cout << ": Building main matrix of size " << matrix.main_matrices[ s ].size() << std::endl;
	Matrix< double > *local = NULL;
	//now build block orderer
#if defined _MCDMV_BICRS || defined _MCDMV_SBD
	std::vector< signed char > hierarchy_datatype;
	MinCRS< double > order;
	std::vector< std::vector< Triplet< double > > > hierarchy = order.induce( matrix.main_matrices[s],
		matrix.local_rowhierarchy[s],  matrix.local_colhierarchy[s],
		matrix.local_rowboundaries[s], matrix.local_colboundaries[s],
		&hierarchy_datatype );
#endif
#ifdef _MCDMV_BICRS
	matrix.main_matrices[ s ].clear();
	for( unsigned long int i=0; i < hierarchy.size(); ++i )
		for( unsigned long int k=0; k<hierarchy[i].size(); ++k )
			matrix.main_matrices[ s ].push_back( hierarchy[i][k] );
	local = CBICRS_factory< double >::getCBICRS( matrix.main_matrices[ s ], matrix.localM[ s ], matrix.localN[ s ] );
#elif defined _MCDMV_SBD
	//now build block orderer
	std::vector< signed char > hierarchy_datatype;
	MinCRS< double > order;
	std::vector< std::vector< Triplet< double > > > hierarchy = order.induce( matrix.main_matrices[s],
		matrix.local_rowhierarchy[s],  matrix.local_colhierarchy[s],
		matrix.local_rowboundaries[s], matrix.local_colboundaries[s],
		&hierarchy_datatype );
#ifndef NDEBUG
	unsigned long int newnnz = 0;
	for( unsigned long int i=0; i<hierarchy.size(); i++ )
		newnnz += hierarchy[ i ].size();
	assert( newnnz == matrix.main_matrices[ s ].size() );
	assert( hierarchy_datatype.size() == hierarchy.size() );
#endif
	local = new HBICRS< double >( hierarchy, &(hierarchy_datatype[0]), matrix.localM[ s ], matrix.localN[ s ] );
#elif defined _MCDMV_BH
	local = new BetaHilbert< double >( matrix.main_matrices[ s ], matrix.localM[ s ], matrix.localN[ s ], 0, &(sub_p_translate[ s ]) );
#elif defined _MCDMV_1DBH
	local = new RDBHilbert< double >( matrix.main_matrices[ s ], matrix.localM[ s ], matrix.localN[ s ], 0, &(sub_p_translate[ s ]) );
#else
	std::cerr << "No underlying scheme defined (See McDMV.cpp, _MCDMV_*)" << std::endl;
	exit( 1 );
#endif

	printproc(s); std::cout << ": Building remote matrix of size " << matrix.remote_separators[ s ].size() << std::endl;
	Matrix< double > *remote = NULL;
#ifdef _MCDMV_REMOTE_INORDER
	remote = CBICRS_factory< double >::getCBICRS( matrix.remote_separators[ s ], matrix.localM[ s ], matrix.localN[ s ] );
#elif defined _MCDMV_REMOTE_ROW_MAJOR
	qsort( &(matrix.remote_separators[ s ][ 0 ]), matrix.remote_separators[ s ].size(), sizeof( Triplet< double > ), &(ICRS< double >::compareTriplets) );
	remote = CBICRS_factory< double >::getCBICRS( matrix.remote_separators[ s ], matrix.localM[ s ], matrix.localN[ s ] );
#elif defined _MCDMV_REMOTE_BH
	remote = new BetaHilbert< double >( matrix.remote_separators[ s ], matrix.localM[ s ], matrix.localN[ s ], 0, &(sub_p_translate[ s ]) );
#elif defined _MCDMV_REMOTE_1DBH
	remote = new RDBHilbert< double >( matrix.remote_separators[ s ], matrix.localM[ s ], matrix.localN[ s ], 0, &(sub_p_translate[ s ]) );
#else
	std::cerr << "No underlying scheme defined (See McDMV.cpp, _MCDMV_REMOTE_*)" << std::endl;
	exit( 1 );
#endif

	std::cout << "Thread " << s << ": ratio local to remote number of nonzeroes is " << (((double)matrix.main_matrices[ s ].size())/((double)matrix.remote_separators[ s ].size())) << std::endl;

	matrix.main_matrices[s].clear();
	matrix.remote_separators[s].clear();
	//local input/output vectors
	const unsigned long int rowBuffer_lo = matrix.local_rowboundaries[ s ][ 0 ];
	const unsigned long int rowBuffer_hi = matrix.local_rowboundaries[ s ].back();
	const unsigned long int colBuffer_lo = matrix.local_colboundaries[ s ][ 0 ];
	const unsigned long int colBuffer_hi = matrix.local_colboundaries[ s ].back();
	std::cout << "Creating local input/output vectors, including buffers" << std::endl;
	std::cout << "(local input  vector length including buffer: " << matrix.localN[s] << ", buffer size: ";
	std::cout << (colBuffer_lo + matrix.localN[s] - colBuffer_hi ) << ";" << std::endl;
	std::cout << " local output vector length including buffer: " << matrix.localM[s] << ", buffer size: ";
	std::cout << (rowBuffer_lo + matrix.localM[s] - rowBuffer_hi ) << ".)" << std::endl;
	double *x = matrix.localX[ s ] = new double[ matrix.localN[s] ];
	double *y = matrix.localY[ s ] = new double[ matrix.localM[s] ];
	
	std::cout << "Intialising local input/output vectors" << std::endl;
	for( unsigned long int i=0; i<matrix.localM[s]; i++ ) y[ i ] = 0.0;
	for( unsigned long int j=0; j<matrix.localN[s]; j++ ) {
		//if not owner, skip this element
		if( j < colBuffer_lo || j >= colBuffer_hi )
			if( matrix.input_vector_distribution[ matrix.col2global[ s ][ j ] ] != s )
				continue;
		//otherwise initialise
		x[ j ] = checkX[ matrix.col2global[ s ][ j ] ];
		//(note that valgrind will detect unitialised accesses, so use that to determine
		// if something here goes awry)
	}
	//here we go
	clock_gettime( CLOCK_MONOTONIC, &stop );
	proc_time = stop.tv_sec-start.tv_sec;
	proc_time+= (stop.tv_nsec-start.tv_nsec)/1000000000.0;
	printproc(s); std::cout << ": time taken for build = " << proc_time << " seconds." << std::endl;
	bsp_sync( s );
	if( s == 0 ) {
		std::cout << "Will do 1 parallel SpMV in the next superstep for verification...";
		fflush( stdout );
	}
	bsp_sync( s );
	spmv( s, matrix, local, remote, x, y );
	if( s == 0 ) {
		std::cout << "done. Verification results:" << std::endl;
	}

#ifdef PRINTY
	bsp_sync( s );
	for( unsigned long int r =0; r < matrix.P; ++r ) {
		for( unsigned long int i = 0; s == r && i < matrix.localM[ s ]; i++ ) 
			std::cout << "--- " << s << " " << i << " " << y[ i ] << std::endl;
		bsp_sync( s );
	}
#else
	bsp_sync( s );
#endif

	if( s == 0 ) {
#ifndef NDEBUG
		//check whether rows are touched at most once
		std::vector< unsigned short int > check;
		for( unsigned long int i = 0; i < matrix.m; ++i )
			check.push_back( 0 );
#endif
		double checkMSE   = 0.0;
		double checksumTS = 0;
		double checksum2D = 0;
		unsigned long int max_e_index = 0;
		unsigned long int max_e_proc  = 0;
		double max_e = 0.0;
		for( unsigned long int s = 0; s < matrix.P; s++ ) {
			for( unsigned long int i = 0; i < matrix.localM[ s ]; i++ ) {
				if( matrix.output_vector_distribution[ matrix.row2global[ s ][ i ] ] == s ) {
#ifndef NDEBUG
					assert( check[ matrix.row2global[ s ][ i ] ] < 2 );
					check[ matrix.row2global[ s ][ i ] ]++;
#endif
					assert( matrix.localY[ s ][ i ] );
					assert( checkY[ matrix.row2global[ s ][ i ] ] );
					double curdiff = fabs( matrix.localY[ s ][ i ] - checkY[ matrix.row2global[ s ][ i ] ] );
#ifdef _DEBUG
					std::cout << matrix.localY[ s ][ i ] << " \t " << checkY[ matrix.row2global[ s ][ i ] ];
					if( curdiff > 1e-15 )
						std::cout << "*** " << matrix.output_vector_distribution[ matrix.row2global[ s ][ i ] ];
#endif
					if( curdiff > max_e ) {
						max_e = curdiff;
						max_e_index = i;
						max_e_proc  = s;
					}
#ifdef _DEBUG
					std::cout << std::endl;
#endif
					curdiff *= curdiff;
					curdiff /= (double)matrix.m;
					checkMSE += curdiff;
					checksum2D += matrix.localY[ s ][ i ] / (double)matrix.m;
					checksumTS += checkY[ matrix.row2global[ s ][ i ] ] / (double)matrix.m;
				}
			}
		}
#ifndef NDEBUG
		for( unsigned long int i = 0; i < matrix.m; ++i )
			assert( check[ i ] < 2 );
#endif
		std::cout << "MSE = " << checkMSE << ", max abs error = " << max_e << " at index " << max_e_index << " while comparing " << matrix.localY[ max_e_proc][ max_e_index ];
		std::cout << " to " << checkY[ matrix.row2global[ max_e_proc][ max_e_index ] ] << std::endl;
		std::cout << "Parallel mean = " << checksum2D << ", sequential mean = " << checksumTS << std::endl;
		std::cout << "Will do 1 parallel SpMV in the next superstep to get the caches `hot'...";
		fflush( stdout );
	}
	bsp_sync( s );
	spmv( s, matrix, local, remote, x, y );
	if( s == 0 ) {
		std::cout << "done" << std::endl;
		std::cout << "Main experiment, " << REPEAT_EXP << " SpMV multiplications, will start in the next superstep...";
		fflush( stdout );
	}
	bsp_sync( s );
	clock_gettime( CLOCK_MONOTONIC, &start );
	t0 = time( NULL );
	for( int i=0; i<REPEAT_EXP; i++ )
		spmv( s, matrix, local, remote, x, y );
	t1 = time( NULL );
	clock_gettime( CLOCK_MONOTONIC, &stop );
	if( s==0 ) {
		std::cout << "done" << std::endl;
		proc_time = (stop.tv_sec-start.tv_sec) * 1000;
		proc_time+= (stop.tv_nsec-start.tv_nsec)/1000000.0;
		printproc(s); std::cout << ": local wall-clock time (GNU)  taken for SpMV = " << (proc_time/((double)REPEAT_EXP)) << " ms." << std::endl;
		printproc(s); std::cout << ": local wall-clock time (ANSI) taken for SpMV = " << ((double)(t1-t0)/((double)REPEAT_EXP)*1000.0) << " ms." << std::endl;
	}

	global_time[ s ]  = (stop.tv_sec-start.tv_sec) * 1000;
	global_time[ s ] += (stop.tv_nsec-start.tv_nsec)/1000000.0;
	global_time[ s ] /= (double)REPEAT_EXP;

	bsp_sync( s );

	//delete local stuff
	delete remote;
	delete local;

	pthread_exit( NULL );
}

int main( int argc, char** argv ) {

#ifndef _MCDMV_HYBRID
	if( argc < 3 ) {
		std::cout << "Usage: " << argv[0] << " <.emm file> <P> <ID_1> <ID_2> ... <ID_P>" << std::endl;
		return 0;
	}
	P = atoi( argv[ 2 ] );
#ifdef SPINLOCK
	condition = new unsigned char[ P ];
	for( size_t s = 0; s < P; ++s )
		condition[ s ] = 0;
#endif
	if( argc != 3+(int)P ) {
		std::cout << "Usage: " << argv[0] << " <.emm file> <P> <ID_1> <ID_2> ... <ID_P>" << std::endl;
		return 0;
	}
#else
	if( argc < 4 ) {
		std::cout << "Usage: " << argv[ 0 ] << " <.emm file> <P> <R> <ID_1> <ID_2> ... <ID_P> <ID_11> <ID_12> ... <ID_1R> <ID_21 > ... <ID_PR>" << std::endl;
		std::cout << std::endl << "(Hybrid scheme: executes P*R threads on manually given pinnings)" << std::endl;
		return 0;
	}
	P = atoi( argv[ 2 ] );
	R = atoi( argv[ 3 ] );
	std::cout << "Will run hybrid on a " << P << " x " << R << " processor grid." << std::endl;
	if( static_cast< unsigned long int >( argc ) != 4 + P + P*R ) {
		std::cout << "Usage: " << argv[ 0 ] << " <.emm file> <P> <R> <ID_1> <ID_2> ... <ID_P> <ID_11> <ID_12> ... <ID_1R> <ID_21 > ... <ID_PR>" << std::endl;
		std::cout << std::endl << "(Hybrid scheme: executes P*R threads on manually given pinnings)" << std::endl;
		return 0;
	}
#endif

	global_time = (double*) malloc( P * sizeof( double ) );

	fn = std::string( argv[1] );

	p_translate = new unsigned long int[ P ];
	for( unsigned int i=0; i<P; i++ ) {
#ifdef _MCDMV_HYBRID
		sub_p_translate.push_back( std::vector< unsigned short int >() );
		p_translate[ i ] = (unsigned int)atoi( argv[ 4 + i ] );
		std::cout << "Master process " << i << " pins to " << p_translate[ i ] << std::endl;
		for( unsigned int j=0; j<R; j++ ) {
			sub_p_translate[ i ].push_back( (unsigned short int)atoi( argv[ 4 + P + i*R + j ] ) );
			std::cout << "Slave proceses " << j << " of master process " << i << " pins to " << sub_p_translate[ i ][ j ] << std::endl;
		}
#else
		p_translate[i] = (unsigned int)atoi( argv[ 3+i ] );
#endif
	}

	std::cout << "Parsing input file" << std::endl;
	parse( fn, matrix );

	std::cout << "Preparing verification vectors" << std::endl;
	checkX  = new double[ matrix.n ];
	checkY  = new double[ matrix.m ];
	checkY2 = new double[ matrix.m ];
	for( unsigned long int i = 0; i < matrix.m; ++i ) checkY[ i ] = checkY2[ i ] = 0.0;
#ifdef ONE_X
	for( unsigned long int j = 0; j < matrix.n; ++j ) checkX[ j ] = 1.0;
#else
	for( unsigned long int j = 0; j < matrix.n; ++j ) checkX[ j ] = rand()/((double)RAND_MAX);
#endif

	std::cout << "Preparing matrix" << std::endl;
	matrix.P = P; //set true P
	prepare( matrix );

	std::cout << "Doing one sequential SpMV for verification purposes, ";
	std::cout << "sequential input/output vector sizes correspond to " << ((matrix.m+matrix.n)*4/1024.0/1024.0) << " MBs of data" << std::endl;
	for( unsigned long int k = 0; k < matrix.local.size(); ++k )
		checkY[ matrix.row_perm[ matrix.local[ k ].i() ] ] += checkX[ matrix.col_perm[ matrix.local[ k ].j() ] ] * matrix.local[ k ].value;
	//done, clear variables
	matrix.local.clear();
	//check
	for( unsigned long int i = 0; i < matrix.m; ++i ) checkY2[ i ] = 0.0;
	for( unsigned long int s = 0; s < matrix.P; s++ ) {
		for( unsigned long int k = 0; k < matrix.main_matrices[ s ].size(); ++k ) {
			checkY2[ matrix.row2global[ s ][ matrix.main_matrices[ s ][ k ].i() ] ] +=
				checkX[ matrix.col2global[ s ][ matrix.main_matrices[ s ][ k ].j() ] ] *
				matrix.main_matrices[ s ][ k ].value;
		}
		for( unsigned long int k = 0; k < matrix.remote_separators[ s ].size(); ++k ) {
			checkY2[ matrix.row2global[ s ][ matrix.remote_separators[ s ][ k ].i() ] ] +=
				checkX[ matrix.col2global[ s ][ matrix.remote_separators[ s ][ k ].j() ] ] * matrix.remote_separators[ s ][ k ].value;
		}
	}
	double max_e = 0, MSE = 0;
	unsigned long int max_e_i = 0;
	for( unsigned long int i = 0; i < matrix.m; ++i ) {
		double curdiff = fabs( checkY2[ i ] - checkY[ i ] );
		if( curdiff > max_e ) {
			max_e = curdiff;
			max_e_i = i;
		}
		curdiff *= curdiff;
		curdiff /= (double)matrix.m;
		MSE += curdiff;
	}
	std::cout << "MSE between serialised SpMV and sequential SpMV: " << MSE << " with the max abs of " << max_e << " while comparing " << checkY2[ max_e_i ] << " and " << checkY[ max_e_i ] << std::endl;

	std::cout << "Spawning threads" << std::endl;
	threads = (pthread_t*)malloc( ( P ) * sizeof( pthread_t ) );
	for( unsigned long int i=0; i<P; i++ ) {
		//set fixed affinity for threads
		cpu_set_t mask;
		CPU_ZERO( &mask );
		CPU_SET ( p_translate[ i ], &mask );
#ifdef _MCDMV_HYBRID
		for( unsigned long int s=0; s<R; s++ ) {
			CPU_SET( sub_p_translate[ i ][ s ], &mask );
		}
#endif
		//prepare attributes
		pthread_attr_t attr;
		pthread_attr_init( &attr );
		//set fixed affinity in attribute, so that it starts binded immediately
		pthread_attr_setaffinity_np( &attr, sizeof( cpu_set_t ), &mask );
		//spawn thread
		pthread_create( &threads[i], &attr, parallel, (void *)i );
		//free attr
		pthread_attr_destroy( &attr );
	}

	//wait for exit
	for( unsigned int i=0; i<P; i++ )
		pthread_join( threads[i], NULL );

	//get statistics
	double maxavg = global_time[ 0 ];
	for( unsigned int s=1; s<P; ++s ) {
		std::cout << "Process " << s << " reports " << global_time[ s ] << " milliseconds taken." << std::endl;
		if( maxavg < global_time[ s ] )
			maxavg = global_time[ s ];
	}
	std::cout << "Maximum average time taken is " << maxavg << " ms." << std::endl;

	std::cout << "Cleanup" << std::endl;
	delete [] p_translate;
	delete [] global_time;
	delete [] checkX;
	delete [] checkY;
	delete [] checkY2;
#ifdef SPINLOCK
	delete []  condition;
#endif

	free(threads);
	return 0;
}

