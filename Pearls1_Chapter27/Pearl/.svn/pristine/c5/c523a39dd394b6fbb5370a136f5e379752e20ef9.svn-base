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
 */


#include <cstdlib>
#include <iostream>
#include <pthread.h>
#include <fstream>
#include <vector>
#include <time.h>
#include <string.h>

#include "Triplet.hpp"
#include "CRS.hpp"
#include "ICRS.hpp"
#include "Hilbert.hpp"
#include "CCSWrapper.hpp"

//user defines:
#define REPEAT_EXP 10000
#define USE_CBICRS
#ifdef USE_CBICRS
 #define MAINMATRIXTYPE CBICRS
#else
 #define MAINMATRIXTYPE BICRS //if you are changing this, also change the corresponding include below!
#endif

//dependent includes
#ifdef USE_CBICRS
 #include "CBICRS.hpp"
#else
 #include "BICRS.hpp" //if mainmatrixtype is changed, also change this include accordingly
#endif

class MVData {
public:
	std::vector< Triplet< double > > local;
	std::vector< std::vector< Triplet< double > > > main_matrices;
	std::vector< std::vector< std::vector< Triplet< double > > > > separated_matrices;
	unsigned int *Pstart;
	unsigned long int *rowboundaries,  *colboundaries;
	unsigned long int *row2global, *col2global;
	unsigned long int *row2proc,   *col2proc;
	unsigned long int *row2index,  *col2index; //inhereted now
	unsigned long int *row_perm, *col_perm;
	unsigned long int **inv_row_perm;
	unsigned long int **inv_col_perm;
	unsigned long int m, n, nz; //inhereted now
	unsigned long int input_length, output_length;
	unsigned int P, Pref;

	MVData() {
		row_perm = col_perm = rowboundaries = colboundaries =
		this->row2global = this->col2global = 
		this->row2proc = this->col2proc = 
		this->row2index = this->col2index = NULL;
		Pstart = NULL;
		inv_row_perm = new unsigned long int*[3];
		inv_col_perm = new unsigned long int*[3];
		for( int i=0; i<3; i++ )
			inv_row_perm[i] = inv_col_perm[i] = NULL;
	}

	~MVData() {
		if( Pstart != NULL ) delete [] Pstart;
		if( rowboundaries != NULL ) delete [] rowboundaries;
		if( colboundaries != NULL ) delete [] colboundaries;
		if( row_perm != NULL ) delete [] row_perm;
		if( col_perm != NULL ) delete [] col_perm;
		if( row2global != NULL ) delete [] row2global;
		if( row2proc != NULL ) delete [] row2proc;
		if( row2index != NULL ) delete [] row2index;
		if( col2global != NULL ) delete [] col2global;
		if( col2proc != NULL ) delete [] col2proc;
		if( col2index != NULL ) delete [] col2index;
		
		for( int i=0; i<3; i++ ) {
			if( inv_row_perm[i] != NULL ) delete [] inv_row_perm[i];
			if( inv_col_perm[i] != NULL ) delete [] inv_col_perm[i];
		}
		delete [] inv_row_perm;
		delete [] inv_col_perm;
	}
};

MVData matrix;

void printproc( const unsigned int proc ) { std::cout << "(" << proc << ")"; }

void parseVC( std::fstream &file, std::string array_name, unsigned long int **array, char *chars,
		const unsigned int s, const unsigned long int P, const char skip=1 ) {
	std::string line;
	unsigned long int uli;
	std::cout << "Will read " << array_name << std::endl;
	if( skip )
		file.getline( chars, 500 ); //go past EOL
	if( file.peek() != '%' ) {
		std::cout << "Not at header line!" << std::endl;
		file >> line;
		std::cout << line << std::endl;
		file >> line;
		std::cout << line << std::endl;
		exit( 1 );
	}
	file.getline( chars, 500 ); //skip header
	std::cout << "Skipping " << chars << std::endl;
	file >> uli; //skip number of vectors (=ret.P)
	std::cout << "Skipping " << uli << std::endl;
	unsigned long int ai, aj;
	for( unsigned int k=0; k<P; k++ ) { //read the P vectors
		file >> ai; //get vector size
		std::cout << k << "th vector is of size " << ai << std::endl;
		if( k==s ) { //this is our vector
			*array = new unsigned long int[ ai ];
			for( unsigned int i=0; i<ai; i++ ) {
				file >> (*array)[ i ];
				(*array)[i]--; //corect base
			}
		} else	
			for( unsigned int i=0; i<ai; i++ )
				file >> aj; //read(skip) element
	}
}

unsigned long int *p_translate = NULL;

/** @param fn File to parse
    @param s  Local processor ID */
void parse( const std::string &fn, MVData &ret ) {
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
	file >> ret.P;

	//read Pstart
	ret.Pstart = new unsigned int[ ret.P+1 ];
	for( unsigned int k=0; k<=ret.P; k++ ) {
		file >> ret.Pstart[ k ];
	}

	//read nonzeroes
	std::cout << "Will read " << nnz << " nonzeroes " << std::endl;
	double av = 1.0; unsigned long int ai, aj;
	for( unsigned int k=0; k<ret.P; k++ ) {
		for( unsigned int i=ret.Pstart[k]; i<ret.Pstart[k+1]; i++ ) {
			file >> ai; file >> aj;
			if( !pattern_matrix )
				file >> av;
			ret.local.push_back( Triplet< double >( ai-1, aj-1, av ) ); //also correct for base
			ret.local[ ret.local.size() - 1 ].meta = k;
		}
	}
	std::cout << "Read a " << ret.m << " X " << ret.n << " matrix of " << ret.local.size() << " nonzeroes (==" << nnz << ")." << std::endl;
	std::cout << "First entry is: (" << ret.local[ 0 ].i() << "," << ret.local[ 0 ].j() << ")=" << ret.local[ 0 ].value << " at " << ret.local[ 0 ].meta<< std::endl;
	std::cout << "Last entry is: (" << ret.local[ ret.local.size()-1 ].i() << "," << ret.local[ ret.local.size()-1 ].j() << ")=" << ret.local[ ret.local.size()-1 ].value << " at " << ret.local[ ret.local.size()-1 ].meta << std::endl;
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
	file >> ret.Pref;
	ret.rowboundaries = new unsigned long int[ ret.Pref ];
	for( unsigned int i=0; i<ret.Pref; i++ ) {
		file >> (ret.rowboundaries[i]);
		ret.rowboundaries[i]--; //correct base
	}
	std::cout << "Read " << ret.Pref << " row boundaries" << std::endl;
	ret.Pref /= 2; //correct P value (is double the number of boundaries)
	std::cout << "Pref is " << ret.Pref << std::endl;

	//read Row hierarchy header
	file.getline( chars, 500 ); //go past EOL
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;

	//skip Row hierarchy
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read Col boundaries header
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;

	//skip Col boundaries
	file >> uli; //half should equal ret.Pref
	if( uli/2 != ret.Pref ) {
		std::cerr << "Error: Refinements in row direction does not equal that in column direction!" << std::endl;
		exit( 1 );
	}
	ret.colboundaries = new unsigned long int[ uli ];
	for( unsigned int j=0; j<uli; j++ ) {
		file >> (ret.colboundaries[j]);
		ret.colboundaries[j]--; //correct base
	}
	std::cout << "Read " << uli << " column boundaries" << std::endl;

	//read Col hierarchy header
	file.getline( chars, 500 ); //go past EOL
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;

	//skip Col hierarchy
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read Local-A header
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;

	//skip current matrix
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//skip row2global
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//skip col2global
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read/skip Row-permutation
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
	for( unsigned int i=0; i<ret.m; i++ ) {
		file >> ret.row_perm[ i ];
		ret.row_perm[ i ]--;
	}
	std::cout << "Last row permutation index read: " << ret.row_perm[ ret.m-1 ] << std::endl;
	//invert
	unsigned long int *dummy = ret.row_perm;
	ret.row_perm = new unsigned long int[ ret.m ];
	for( unsigned int k=0; k<ret.m; k++ ) ret.row_perm[ dummy[ k ] ] = k;
	delete [] dummy;
	//read/skip Column-permutation
	file.getline( chars, 500 ); //skip EOL
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> uli;
	ret.col_perm = new unsigned long int[ uli ];
	for( unsigned int i=0; i<ret.n; i++ ) {
		file >> ret.col_perm[ i ];
		ret.col_perm[ i ]--;
	}
	std::cout << "Last column permutation index read: " << ret.col_perm[ ret.n-1 ] << std::endl;
	//invert
	dummy = ret.col_perm;
	ret.col_perm = new unsigned long int[ ret.n ];
	for( unsigned int k=0; k<ret.n; k++ ) ret.col_perm[ dummy[ k ] ] = k;
	delete [] dummy;
	//read/skip Input-vector
	file.getline( chars, 500 ); //skip EOL
	file.getline( chars, 500 );
	std::cout << "Reading " << chars << std::endl;
	file >> uli; //should be ret.n
	ret.col2proc = new long unsigned int[ uli ];
	std::cout << uli << "==" << ret.m << "?" << std::endl;
	file >> uli; //should be ret.P
	for( unsigned int i=0; i<ret.n; i++ ) {
		file >> uli;
		file >> ret.col2proc[ uli-1 ];
		ret.col2proc[ uli-1 ]--; //correct base
	}
	//read Output-vector
	file.getline( chars, 500 ); //skip EOL
	if( file.peek() != '%' ) {
		std::cout << "Not at header line!" << std::endl;
		file.getline( chars, 500 );
		std::cout << chars << std::endl;
		file.getline( chars, 500 );
		std::cout << chars << std::endl;
		exit( 1 );
	}
	file.getline( chars, 500 );
	std::cout << "Reading row2proc " << chars << std::endl;
	file >> uli; //should be ret.m
	ret.row2proc = new long unsigned int[ uli ];
	file >> uli; //should be ret.P
	for( unsigned int i=0; i<ret.m; i++ ) {
		file >> uli;
		file >> ret.row2proc[ uli-1 ];
		ret.row2proc[ uli-1 ]--; //correct base
	}

	//done
	delete [] chars;
	if( !file.eof() ) {
		std::cout << "Warning: input file not at end!" << std::endl;
		file >> line;
		std::cout << "Next line is: " << line << std::endl;
	} else {
		std::cout << "End of file; parse OK." << std::endl;
	}
}

std::string fn;

pthread_t      *threads    = NULL;
pthread_mutex_t sync_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  sync_cond  = PTHREAD_COND_INITIALIZER;
unsigned int sync_val = 0;
unsigned int P = 0;

void bsp_sync() {
        pthread_mutex_lock( &sync_mutex );
        sync_val++;
        if( sync_val < P ) pthread_cond_wait( &sync_cond, &sync_mutex );
        else {
                sync_val = 0;
                pthread_cond_broadcast( &sync_cond );
        }
        pthread_mutex_unlock( &sync_mutex );
}

void prepare( MVData &matrix ) {
	//upscale boundaries
	std::cout << "Upscaling from " << matrix.Pref << " to " << matrix.P << std::endl;
	const unsigned int boundary_stepsize = matrix.Pref / matrix.P - 1;
	std::cout << "Upscaled from: "; for( unsigned int i=0; i<2*matrix.Pref; i++ ) std::cout << matrix.rowboundaries[ i ] << " "; std::cout << std::endl;
	std::cout << "Upscaled from: "; for( unsigned int i=0; i<2*matrix.Pref; i++ ) std::cout << matrix.colboundaries[ i ] << " "; std::cout << std::endl;
	for( unsigned int i=1, k=1; k<2*matrix.P-1; i+=2, k+=2 ) {
		for( unsigned int j=0; j<boundary_stepsize; j++ ) //skip this many per P
			i+=2;
		matrix.rowboundaries[ k ] = matrix.rowboundaries[ i ];
		matrix.rowboundaries[k+1] = matrix.rowboundaries[i+1];
		matrix.colboundaries[ k ] = matrix.colboundaries[ i ];
		matrix.colboundaries[k+1] = matrix.colboundaries[i+1];
	}
	matrix.rowboundaries[ 2*matrix.P-1 ] = matrix.rowboundaries[ 2*matrix.Pref-1 ]; //copy final index
	matrix.colboundaries[ 2*matrix.P-1 ] = matrix.colboundaries[ 2*matrix.Pref-1 ]; //copy final index
	std::cout << "Upscaled to: "; for( unsigned int i=0; i<2*matrix.P; i++ ) std::cout << matrix.rowboundaries[ i ] << " "; std::cout << std::endl;
	std::cout << "Upscaled to: "; for( unsigned int i=0; i<2*matrix.P; i++ ) std::cout << matrix.colboundaries[ i ] << " "; std::cout << std::endl;

	//scale large number of processors back to the smaller number of processors
	unsigned int *proc2proc = new unsigned int[ matrix.Pref ];
	for( unsigned int i=0; i<matrix.Pref; i++ ) proc2proc[ i ] = matrix.P; //invalid
	std::cout << "Looping over " << matrix.nz << " elements" << std::endl;
	for( unsigned int k=0; k<matrix.nz; k++ ) {
		const Triplet< double > cur = matrix.local[ k ];
		for( unsigned int i=0; i<2*matrix.P-1; i+= 2 ) { //check if nonzero is in pure block
			if( (matrix.row_perm[ cur.i() ] >= matrix.rowboundaries[i] && matrix.row_perm[ cur.i() ] < matrix.rowboundaries[i+1]) &&
			    (matrix.col_perm[ cur.j() ] >= matrix.colboundaries[i] && matrix.col_perm[ cur.j() ] < matrix.colboundaries[i+1]) ) {
				if( (unsigned int)cur.meta >= matrix.Pref ) {
					std::cerr << cur.meta << " is larger than the number of processors!" << std::endl;
					exit( 1 );
				}
				if( proc2proc[ cur.meta ] == matrix.P ) proc2proc[ cur.meta ] = i/2;
				else if( proc2proc[ cur.meta ] != i/2 ) {
					std::cerr << "There is a processor in more than one pure block: " << cur.meta << " is in " << proc2proc[ cur.meta ] << " and " << (i/2) << std::endl;
					std::cerr << "Caused by nonzero at " << cur.i() << ", " << cur.j() << " which fits in pure category " << (i/2) << std::endl;
					std::cerr << "This nonzero resides at processor " << cur.meta << " and should be upscaled to processor " << (i/2) << " while previous nonzeroes from processor ";
					std::cerr << cur.meta << " were upscaled to processor " << proc2proc[ cur.meta ] << std::endl;
					exit( 1 );
				}
				break;
			}
		}
	}
	std::cout << "Triplet p-translate: "; for( unsigned int i=0; i<matrix.Pref; i++ ) std::cout << proc2proc[ i ] << " "; std::cout << std::endl;
	std::cout << "Translating vector assignments" << std::endl;
	for( unsigned int i=0; i<matrix.m; i++ ) matrix.row2proc[ i ] = proc2proc[ matrix.row2proc[ i ] ]; 
	for( unsigned int i=0; i<matrix.n; i++ ) matrix.col2proc[ i ] = proc2proc[ matrix.col2proc[ i ] ]; 
	std::cout << "Permuting vector assignments" << std::endl;
	unsigned long int *dummy = matrix.row2proc;
	matrix.row2proc = new unsigned long int[ matrix.m ];
	for( unsigned long int i=0; i<matrix.m; i++ )
		matrix.row2proc[ matrix.row_perm[ i ] ] = dummy[ i ];
	delete [] dummy;
	dummy = matrix.col2proc;
	matrix.col2proc = new unsigned long int[ matrix.n ];
	for( unsigned long int j=0; j<matrix.n; j++ )
		matrix.col2proc[ matrix.col_perm[ j ] ] = dummy[ j ];
	delete [] dummy;

	std::cout << "Optimising vector permutation" << std::endl;
	for( unsigned int i=1; i<2*matrix.P-1; i+= 2 ) { //loop over separator bounds
		unsigned int *Pcount = new unsigned int[ P ]; //perform counting sort
		for( unsigned int k=0; k<matrix.P; k++ ) Pcount[ k ] = 0;
		//count, rows first
		for( unsigned int k=matrix.rowboundaries[i]; k<matrix.rowboundaries[i+1]; k++ ) //loop over indices in separator
			Pcount[ matrix.row2proc[ k ] ]++;
		for( unsigned int k=1; k<matrix.P; k++ ) //cumulative
			Pcount[ k ] += Pcount[k-1];
		unsigned int *Pplace = new unsigned int[ P ];
		Pplace[ 0 ] = 0;
		for( unsigned int k=1; k<matrix.P; k++ ) //determine start locations
			Pplace[ k ] = Pcount[ k-1 ];
		//sort 
		delete [] Pcount;
		Pcount = new unsigned int[ matrix.rowboundaries[ i+1 ] - matrix.rowboundaries[ i ] ];
		for( unsigned int k=matrix.rowboundaries[i]; k<matrix.rowboundaries[i+1]; k++ )
			Pcount[ k-matrix.rowboundaries[i] ] = Pplace[ matrix.row2proc[ k ] ]++;
		//Pcount contains new subpermutation, apply
		dummy = new unsigned long int[ matrix.rowboundaries[ i+1 ] - matrix.rowboundaries[ i ] ]; //first row2proc
		for( unsigned int k=matrix.rowboundaries[i]; k<matrix.rowboundaries[i+1]; k++ )
			dummy[ Pcount[ k-matrix.rowboundaries[i] ] ] = matrix.row2proc[ k ];
		for( unsigned int k=matrix.rowboundaries[i]; k<matrix.rowboundaries[i+1]; k++ )
			matrix.row2proc[ k ] = dummy[ k - matrix.rowboundaries[i] ];
		delete [] dummy;
		//now row_perm
		for( unsigned int k=0; k<matrix.m; k++ )
			if( matrix.row_perm[ k ] >= matrix.rowboundaries[i] && matrix.row_perm[ k ] < matrix.rowboundaries[i+1] )
				matrix.row_perm[ k ] = Pcount[ matrix.row_perm[ k ] - matrix.rowboundaries[ i ] ] + matrix.rowboundaries[ i ];
		delete [] Pcount;
		//count, now cols
		Pcount = new unsigned int[ P ]; //perform counting sort
		for( unsigned int k=0; k<matrix.P; k++ ) Pcount[ k ] = 0;
		for( unsigned int k=matrix.colboundaries[i]; k<matrix.colboundaries[i+1]; k++ ) //loop over indices in separator
			Pcount[ matrix.col2proc[ k ] ]++;
		for( unsigned int k=1; k<matrix.P; k++ ) //cumulative
			Pcount[ k ] += Pcount[k-1];
		Pplace[ 0 ] = 0;
		for( unsigned int k=1; k<matrix.P; k++ ) //determine start locations
			Pplace[ k ] = Pcount[ k-1 ];
		//sort 
		delete [] Pcount;
		Pcount = new unsigned int[ matrix.colboundaries[ i+1 ] - matrix.colboundaries[ i ] ];
		for( unsigned int k=matrix.colboundaries[i]; k<matrix.colboundaries[i+1]; k++ )
			Pcount[ k-matrix.colboundaries[i] ] = Pplace[ matrix.col2proc[ k ] ]++;
		//Pcount contains new subpermutation, apply
		dummy = new unsigned long int[ matrix.colboundaries[ i+1 ] - matrix.colboundaries[ i ] ]; //first col2proc
		for( unsigned int k=matrix.colboundaries[i]; k<matrix.colboundaries[i+1]; k++ )
			dummy[ Pcount[ k-matrix.colboundaries[i] ] ] = matrix.col2proc[ k ];
		for( unsigned int k=matrix.colboundaries[i]; k<matrix.colboundaries[i+1]; k++ )
			matrix.col2proc[ k ] = dummy[ k - matrix.colboundaries[i] ];
		delete [] dummy;
		//now col_perm
		for( unsigned int k=0; k<matrix.n; k++ )
			if( matrix.col_perm[ k ] >= matrix.colboundaries[i] && matrix.col_perm[ k ] < matrix.colboundaries[i+1] )
				matrix.col_perm[ k ] = Pcount[ matrix.col_perm[ k ] - matrix.colboundaries[ i ] ] + matrix.colboundaries[ i ];
		delete [] Pcount;
		delete [] Pplace;

	}
	
	std::cout << "Translating matrix" << std::endl;
	for( unsigned int i=0; i<matrix.P; i++ ) {
		matrix.main_matrices.push_back( std::vector< Triplet< double > >() );
		matrix.separated_matrices.push_back( std::vector< std::vector< Triplet< double > > >() );
		for( unsigned int j=0; j<matrix.P; j++ )
			matrix.separated_matrices[i].push_back( std::vector< Triplet< double > >() );
	}
	unsigned long int new_nnz = 0;
	//loop over all entries
	for( unsigned int k=0; k<matrix.nz; k++ ) {
		const Triplet< double > cur = matrix.local[ k ];
		//get permuted entry
		Triplet< double > newt = Triplet< double >( matrix.row_perm[ cur.i() ], matrix.col_perm[ cur.j() ], cur.value );
		newt.meta = proc2proc[ cur.meta ];
		char not_pure = 1;
		//check if permuted entry is in pure category
		for( unsigned int i=0; i<2*matrix.P-1; i+= 2 ) { //loop over row-wise pure blocks
			if( newt.i() >= matrix.rowboundaries[i] && newt.i() < matrix.rowboundaries[i+1] ) {
				//it is, add to main matrix
				matrix.main_matrices[newt.meta].push_back( newt );
				new_nnz++;
				not_pure = 0;
			}
		}
		//it is a separator entry
		if( not_pure ) {
			//check if it is owned by the same processor as the pure block is
			if( matrix.row2proc[ newt.i() ] == newt.meta ) {
				matrix.main_matrices[newt.meta].push_back( newt );
				new_nnz++;
			} else { //if it isn't, add to approprite separated matrix
				matrix.separated_matrices[ newt.meta ][ matrix.row2proc[ newt.i() ] ].push_back( newt );
				new_nnz++;
			}
		}
	}
	assert( new_nnz == matrix.nz );
	for( unsigned int i=0; i<matrix.P; i++ ) {
		assert( matrix.separated_matrices[i][i].size() == 0 );
	}

	//done!
	matrix.local.clear();
	delete [] proc2proc;
	return;
}

double *x = NULL;
double *y = NULL;

void spmv( const unsigned int s,
		Matrix< double > *A,
		CCSWrapper< double, ICRS< double >, ULI > **Seps,
		const unsigned int numSeps ) {
	//superstep 1, local mv
	A->zax( x, y );
	bsp_sync();
	//supersteps 2,3,...,P
	for( unsigned int i=0; i<numSeps; i++ ) {
		Seps[i]->zax( x, y );
		bsp_sync();
	}
	return;
}

void *parallel( void *threadid ) {
	struct timespec start, stop;
	time_t t0, t1;
	double proc_time;
	unsigned int s = (unsigned int) ((long)threadid);

	if( s == 0 ) std::cout << "Starting up " << (REPEAT_EXP) << " syncs for benchmarking" << std::endl;
	clock_gettime( CLOCK_MONOTONIC, &start );
	for( unsigned int i=0; i<REPEAT_EXP; i++ ) bsp_sync();
	clock_gettime( CLOCK_MONOTONIC, &stop );
	proc_time = (stop.tv_sec-start.tv_sec) * 1000;
	proc_time+= (stop.tv_nsec-start.tv_nsec)/1000000.0;
	if( s == 0 ) { printproc(s); std::cout << ": local time estimate for l = " << (proc_time/((double)REPEAT_EXP)) << " ms." << std::endl; }
	bsp_sync();
	clock_gettime( CLOCK_MONOTONIC, &start );
	printproc(s); std::cout << ": Doing initialisation" << std::endl;
	printproc(s); std::cout << ": Building data structures" << std::endl;

#ifdef USE_CBICRS
	Matrix< double > *local = CBICRS_factory< double >::getCBICRS( matrix.main_matrices[s], matrix.m, matrix.n, 0 );
#else
	Matrix< double > *local = new MAINMATRIXTYPE< double >( matrix.main_matrices[s], matrix.m, matrix.n, 0 );
#endif
	matrix.main_matrices[s].clear();
	typedef CCSWrapper< double, ICRS< double >, ULI > SEP_DATA_T;
	typedef SEP_DATA_T* P_SEP_DATA_T;
	P_SEP_DATA_T *Seps = new P_SEP_DATA_T[matrix.P-1];
	for( unsigned int i=1; i<matrix.P; i++ ) {
		Seps[i-1] = new SEP_DATA_T( matrix.separated_matrices[s][(s+i)%matrix.P], matrix.m, matrix.n, 0 );
	}
	//make available
	clock_gettime( CLOCK_MONOTONIC, &stop );
	proc_time = stop.tv_sec-start.tv_sec;
	proc_time+= (stop.tv_nsec-start.tv_nsec)/1000000000.0;
	printproc(s); std::cout << ": time taken for build = " << proc_time << " seconds." << std::endl;
	if( s == 0 ) std::cout << "Main experiment, " << REPEAT_EXP << " SpMV multiplications, will start in the next superstep" << std::endl;
	bsp_sync();
	clock_gettime( CLOCK_MONOTONIC, &start );
	t0 = time( NULL );
	for( int i=0; i<REPEAT_EXP; i++ )
		spmv( s, local, &(Seps[0]), matrix.P-1 );
	t1 = time( NULL );
	clock_gettime( CLOCK_MONOTONIC, &stop );
	if( s==0 ) {
		proc_time = (stop.tv_sec-start.tv_sec) * 1000;
		proc_time+= (stop.tv_nsec-start.tv_nsec)/1000000.0;
		printproc(s); std::cout << ": local wall-clock time (GNU)  taken for SpMV = " << (proc_time/((double)REPEAT_EXP)) << "ms." << std::endl;
		printproc(s); std::cout << ": local wall-clock time (ANSI) taken for SpMV = " << ((double)(t1-t0)/((double)REPEAT_EXP)*1000.0) << "ms." << std::endl;
	}
	for( unsigned int i=0; i<matrix.P-1; i++ )
		delete Seps[i];
	delete [] Seps;
	delete local;
	pthread_exit( NULL );
}

int main( int argc, char** argv ) {

	if( argc < 3 ) {
		std::cout << "Usage: " << argv[0] << " <.emm file> <P> <ID_1> <ID_2> ... <ID_P>" << std::endl;
		return 0;
	}

	P = atoi( argv[2] );
	fn = std::string( argv[1] );

	if( argc != 3+(int)P ) {
		std::cout << "Usage: " << argv[0] << " <.emm file> <P> <ID_1> <ID_2> ... <ID_P>" << std::endl;
		return 1;
	}

	p_translate = new unsigned long int[ P ];
	for( unsigned int i=0; i<P; i++ )
		p_translate[i] = (unsigned int)atoi( argv[ 3+i ] );

	std::cout << "Parsing input file" << std::endl;
	parse( fn, matrix );
	std::cout << "Preparing matrix" << std::endl;
	matrix.P = P;
	prepare( matrix );
	std::cout << "Initialising input/output vectors, corresponding to " << ((matrix.m+matrix.n)*4/1024.0/1024.0) << " MB of data" << std::endl;
	x = new double[ matrix.n ];
	y = new double[ matrix.m ];
	for( unsigned int i=0; i<matrix.m; i++ ) y[i] = 0.0;
	for( unsigned int j=0; j<matrix.n; j++ ) x[j] = 1.0;

	std::cout << "Spawning threads" << std::endl;
	threads = (pthread_t*)malloc( ( P ) * sizeof( pthread_t ) );
	for( unsigned long int i=0; i<P; i++ ) {
		//set fixed affinity for threads
		cpu_set_t mask;
		CPU_ZERO( &mask );
		CPU_SET ( p_translate[ i ], &mask );
		pthread_attr_t attr;
		pthread_attr_init( &attr );
		//set fixed affinity in attribute, so that it starts binded immediately
		pthread_attr_setaffinity_np( &attr, sizeof( cpu_set_t ), &mask );
		pthread_create( &threads[i], &attr, parallel, (void *)i );
		//free attr
		pthread_attr_destroy( &attr );
	}

	//wait for exit
	for( unsigned int i=0; i<P; i++ )
		pthread_join( threads[i], NULL );

	std::cout << "Cleanup" << std::endl;
	delete [] p_translate;
	delete [] x;
	delete [] y;
	free(threads);
	return 0;
}

