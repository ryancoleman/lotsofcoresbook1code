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

#include "Triplet.hpp"
#include "ICRS.hpp"
#include "CCSWrapper.hpp"

#define REPEAT_EXP 10

class DMVData { //distributed MV data
public:
	std::vector< Triplet< double > > data;
	unsigned long int *row2global, *col2global, *row2proc, *col2proc, *row2index, *col2index;
	unsigned long int m, n;
	double *x_buffer;
	double *y_buffer;

	void allocate() {
		y_buffer = new double[ m ];
		x_buffer = new double[ n ];
	}

	void prepare() {
		for( unsigned int i=0; i<m; i++ ) y_buffer[ i ] = 0.0;
	}

	DMVData() {
		row2global = col2global = row2proc = col2proc = row2index = col2index = NULL;
		x_buffer = y_buffer = NULL;
		m = n = 0;
	}

	~DMVData() {
		if( row2global != NULL ) delete [] row2global;
		if( col2global != NULL ) delete [] col2global;
		if( row2proc != NULL ) delete [] row2proc;
		if( col2proc != NULL ) delete [] col2proc;
		if( row2index != NULL ) delete [] row2index;
		if( col2index != NULL ) delete [] col2index;
		if( x_buffer != NULL ) delete [] x_buffer;
		if( y_buffer != NULL ) delete [] y_buffer;
	}
};

class MVData {
public:
	std::vector< Triplet< double > > local;
	DMVData sepA1; //fanout with neighbour only
	DMVData sepA2; //fanin  with neighbour only
	DMVData sepA3; //fanout and fanin with neighbour only
	DMVData sepB1; //fanout with all processors
	DMVData sepB2; //fanin  with all processors
	DMVData sepB3; //fanout and fanin with all processors
	unsigned long int *rowboundaries,  *colboundaries;
	unsigned long int *row2global, *col2global;
	unsigned long int *row2proc,   *col2proc;
	unsigned long int *row2index,  *col2index; //inhereted now
	unsigned long int **inv_row_perm;
	unsigned long int **inv_col_perm;
	unsigned long int m,n,nz; //inhereted now
	//unsigned long int nz;
	unsigned long int input_length, output_length;
	unsigned int P, Pref;
	unsigned long int r_localcount, r_subset1, r_subset2;
	unsigned long int c_localcount, c_subset1, c_subset2;

	MVData() {
		rowboundaries = colboundaries = this->row2global = this->col2global = this->row2proc = this->col2proc = this->row2index = this->col2index = NULL;
		inv_row_perm = new unsigned long int*[3];
		inv_col_perm = new unsigned long int*[3];
		for( int i=0; i<3; i++ )
			inv_row_perm[i] = inv_col_perm[i] = NULL;
	}

	~MVData() {
		if( row2global != NULL ) delete [] row2global;
		if( col2global != NULL ) delete [] col2global;
		if( row2proc != NULL ) delete [] row2proc;
		if( col2proc != NULL ) delete [] col2proc;
		if( row2index != NULL ) delete [] row2index;
		if( col2index != NULL ) delete [] col2index;
		if( rowboundaries != NULL ) delete [] rowboundaries;
		if( colboundaries != NULL ) delete [] colboundaries;
		if( sepA1.row2global != NULL ) {
			delete [] sepA3.row2global;
			delete [] sepA3.col2global;
			delete [] sepA3.row2proc;
			delete [] sepA3.col2proc;
			delete [] sepA3.row2index;
			delete [] sepA3.col2index;
			delete [] sepB3.row2global;
			delete [] sepB3.col2global;
			delete [] sepB3.row2proc;
			delete [] sepB3.col2proc;
			delete [] sepB3.row2index;
			delete [] sepB3.col2index;
		}
		for( int i=0; i<3; i++ ) {
			if( inv_row_perm[i] != NULL ) delete [] inv_row_perm[i];
			if( inv_col_perm[i] != NULL ) delete [] inv_col_perm[i];
		}
		delete [] inv_row_perm;
		delete [] inv_col_perm;
	}
};

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
void parse( const std::string &fn, const unsigned int s, MVData &ret ) {
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
	while( file.peek() == '%' ) {
		file.getline( chars, 500 );
		std::cout << "Skipping " << chars << std::endl;
	}

	//read parameters
	unsigned long int uli, nnz;
	file >> ret.m;
	file >> ret.n;
	file >> nnz; //full nnz
	file >> ret.P;

	//skip current matrix A
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read PAQ header
	file.getline( chars, 500 );
	std::cout << "File reader at " << chars << std::endl;

	//skip current matrix PAQ
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read Row boundaries header
	file.getline( chars, 500 );
	std::cout << "File reader at " << chars << std::endl;

	//read row boundaries
	file >> ret.Pref;
	ret.rowboundaries = new unsigned long int[ ret.Pref ];
	for( unsigned int i=0; i<ret.Pref; i++ ) {
		file >> (ret.rowboundaries[i]);
		ret.rowboundaries[i]--; //correct base
	}
	std::cout << "Read " << ret.Pref << " row boundaries" << std::endl;
	ret.Pref /= 2; //correct P value (is double the number of boundaries)

	//read Row hierarchy header
	file.getline( chars, 500 ); //go past EOL
	file.getline( chars, 500 );
	std::cout << "File reader at " << chars << std::endl;

	//skip Row hierarchy
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read Col boundaries header
	file.getline( chars, 500 );
	std::cout << "File reader at " << chars << std::endl;

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
	std::cout << "File reader at " << chars << std::endl;

	//skip Col hierarchy
	while( file.peek() != '%' ) {
		file.getline( chars, 500 );
	}

	//read Local-A header
	file.getline( chars, 500 );
	std::cout << "File reader at " << chars << std::endl;

	//read parameter line
	file.getline( chars, 500 ); //should be the same as read from main distributed A

	//determine which nnz's to read
	unsigned long int from, to;
	from = to = 0;
	for( unsigned int k=0; k<=ret.P; k++ ) {
		file >> uli;
		if( k==s ) {
			from = uli;
			file >> to;
			k++;
		}
	}
	
	//read nonzeroes
	printproc(s);std::cout << "Will read nonzeroes from " << from << " till " << to << std::endl;
	from--; to--; //correct base
	double av; unsigned long int ai, aj;
	ret.m = ret.n = 0;
	for( unsigned int i=0; i<nnz; i++ ) {
		file >> ai; file >> aj; file >> av;
		if( i >= from && i < to ) {
			if( ai > ret.m ) ret.m = ai;
			if( aj > ret.n ) ret.n = aj;
			ret.local.push_back( Triplet< double >( ai-1, aj-1, av ) ); //also correct for base
		}
	}
	printproc(s);std::cout << "Read a " << ret.m << " X " << ret.n << " matrix of " << ret.local.size() << " nonzeroes (==" << (to-from) << ")." << std::endl;
	printproc(s);std::cout << "First entry is: (" << ret.local[ 0 ].i() << "," << ret.local[ 0 ].j() << ")=" << ret.local[ 0 ].value << std::endl;
	printproc(s);std::cout << "Last entry is: (" << ret.local[ ret.local.size()-1 ].i() << "," << ret.local[ ret.local.size()-1 ].j() << ")=" << ret.local[ ret.local.size()-1 ].value << std::endl;

	//read row2global
	parseVC( file, "row2global", &(ret.row2global), chars, s, ret.P );
	//read col2global
	parseVC( file, "col2global", &(ret.col2global), chars, s, ret.P );
	//skip EOL
	file.getline( chars, 500 );
	//read global-A header
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	//skip global-A
	while( file.peek() != '%' )
		file.getline( chars, 500 );
	//read/skip Row-permutation
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;	
	while( file.peek() != '%' )
		file.getline( chars, 500 );
	//read/skip Column-permutation
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	while( file.peek() != '%' )
		file.getline( chars, 500 );
	//read/skip Input-vector
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	while( file.peek() != '%' )
		file.getline( chars, 500 );
	//read/skip Output-vector
	file.getline( chars, 500 );
	std::cout << "Skipping " << chars << std::endl;
	while( file.peek() != '%' )
		file.getline( chars, 500 );
	//read local output vector length
	file.getline( chars, 500 ); //skip header
	file >> uli; //skip P
	for( unsigned int k=0; k<ret.P; k++ ) {
		file >> uli;
		if (k==s) ret.output_length = uli;
	}
	std::cout << "Read local number of rows: " << ret.m << std::endl;
	//read row2proc
	parseVC( file, "row2proc", &(ret.row2proc), chars, s, ret.P );
	//read row2index
	parseVC( file, "row2index", &(ret.row2index), chars, s, ret.P );
	//read local n length
	file.getline( chars, 500 ); //skip EOL 
	file.getline( chars, 500 ); //skip header
	file >> uli; //skip P
	for( unsigned int k=0; k<ret.P; k++ ) {
		file >> uli;
		if (k==s) ret.input_length = uli;
	}
	std::cout << "Read local number of columns: " << ret.n << std::endl;
	//read col2proc
	parseVC( file, "col2proc", &(ret.col2proc), chars, s, ret.P );
	//read col2index
	parseVC( file, "col2index", &(ret.col2index), chars, s, ret.P );
	//skip EOL
	file.getline( chars, 500 );
	file.getline( chars, 500 );

	//done
	delete [] chars;
	if( !file.eof() ) {
		printproc( s );
		std::cout << "Warning: input file not at end!" << std::endl;
		file >> line;
		printproc( s );
		std::cout << "Next line is: " << line << std::endl;
	} else {
		printproc( s );
		std::cout << "End of file; parse OK." << std::endl;
	}
}

std::string fn;

pthread_t      *threads    = NULL;
pthread_mutex_t sync_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  sync_cond  = PTHREAD_COND_INITIALIZER;
unsigned int sync_val = 0;
unsigned int P = 0;
pthread_mutex_t sync_mutex1= PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  sync_cond1 = PTHREAD_COND_INITIALIZER;
unsigned int sync_val1 = 0;
pthread_mutex_t sync_mutex2= PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  sync_cond2 = PTHREAD_COND_INITIALIZER;
unsigned int sync_val2 = 0;

void bsp_sync1() {
        pthread_mutex_lock( &sync_mutex1 );
        sync_val1++;
        if( sync_val1 < 2 ) pthread_cond_wait( &sync_cond1, &sync_mutex1 );
        else {
                sync_val1 = 0;
                pthread_cond_broadcast( &sync_cond1 );
        }
        pthread_mutex_unlock( &sync_mutex1 );
}

void bsp_sync2() {
        pthread_mutex_lock( &sync_mutex2 );
        sync_val2++;
        if( sync_val2 < 2 ) pthread_cond_wait( &sync_cond2, &sync_mutex2 );
        else {
                sync_val2 = 0;
                pthread_cond_broadcast( &sync_cond2 );
        }
        pthread_mutex_unlock( &sync_mutex2 );
}

void bsp_sync() {
        pthread_mutex_lock( &sync_mutex );
        sync_val++;
//      printf( "Sync counter after inc: %d\n", sync_val );
        if( sync_val < P ) pthread_cond_wait( &sync_cond, &sync_mutex );
        else {
                sync_val = 0;
                pthread_cond_broadcast( &sync_cond );
        }
//      printf( "Sync counter after leave: %d\n", sync_val );
        pthread_mutex_unlock( &sync_mutex );
}

unsigned int *inv_p_translate = NULL;

void check( DMVData &matrix ) {
	for( unsigned long int i=0; i<matrix.data.size(); i++ ) {
		if( matrix.data[i].i() > matrix.m )
			std::cerr << "Element outside matrix!" << std::endl;
		else if( matrix.data[i].j() > matrix.n )
			std::cerr << "Element outside matrix!" << std::endl;
		else continue;
		exit( 1 );
	}
}

void localise( DMVData &matrix, unsigned long int *row2proc, unsigned long int *col2proc,
		unsigned long int *row2index, unsigned long int *col2index ) {
	std::vector< Triplet< double > > replacement;
	unsigned int r_min = matrix.data.size() > 0 ? matrix.data[0].i() : 0;
	unsigned int c_min = matrix.data.size() > 0 ? matrix.data[0].j() : 0;
	for( unsigned int i=1; i<matrix.data.size(); i++ ) {
		if( r_min > matrix.data[i].i() ) r_min = matrix.data[i].i();
		if( c_min > matrix.data[i].j() ) c_min = matrix.data[i].j();
	}
	for( unsigned int i=0; i<matrix.data.size(); i++ ) {
		if( matrix.data[i].i() >= matrix.m || matrix.data[i].j() >= matrix.n ) {
			std::cerr << "Error: matrix to be localised is invalid! (" << matrix.data[i].i() << "," << matrix.data[i].j() << ") while matrix is " << matrix.m << " X " << matrix.n << std::endl;
			exit( 1 );
		}
		replacement.push_back( Triplet< double >( matrix.data[i].i() - r_min, matrix.data[i].j() - c_min, matrix.data[i].value ) );
	}
	matrix.row2proc = new unsigned long int [ matrix.m-r_min ];
	matrix.row2index = new unsigned long int [ matrix.m-r_min ];
	for( unsigned int i=r_min; i<matrix.m; i++ ) {
		matrix.row2proc[i-r_min]  = row2proc[ i ];
		matrix.row2index[i-r_min] = row2index[ i ];
	}
	matrix.col2proc = new unsigned long int [ matrix.n-c_min ];
	matrix.col2index = new unsigned long int [ matrix.n-c_min ];
	for( unsigned int i=c_min; i<matrix.n; i++ ) {
		matrix.col2proc[i-c_min]  = col2proc[ i ];
		matrix.col2index[i-c_min] = col2index[ i ];
	}
	matrix.data = replacement;
	matrix.m   -= r_min;
	matrix.n   -= c_min;
	check( matrix );
}

unsigned long int prepare( MVData &matrix, const unsigned int s ) {
	//determine direct neighbour
	//take modulo 2 neighbour
	const unsigned int real_id   = inv_p_translate[ s ];
	const unsigned int neighbour = p_translate[ (real_id%2==1) ? 2*(real_id/2) : 2*(real_id/2) + 1 ];
	/*
	//determine automagically (bad idea wrt to affinity on NUMA)
	unsigned int *neighbours = new unsigned int[ matrix.P ];
	for( unsigned int i=0; i<matrix.P; i++ )
		neighbours[i] = 0;
	for( unsigned int i=0; i<matrix.m; i++ )
		neighbours[ matrix.row2proc[i] ]++;
	unsigned int max = (s==0) ? 1 : 0;
	for( unsigned int i=1; i<matrix.P; i++ )
		if( i!=s && neighbours[ i ] > neighbours[ max ] ) max = i;
	const unsigned int neighbour = max;
	delete [] neighbours;*/

	//upscale boundaries
	const unsigned int boundary_stepsize = matrix.Pref / matrix.P - 1;
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
	printproc(s);std::cout << "Locally upscaled to: "; for( unsigned int i=0; i<2*matrix.P; i++ ) std::cout << matrix.rowboundaries[ i ] << " "; std::cout << std::endl;
	printproc(s);std::cout << "Locally upscaled to: "; for( unsigned int i=0; i<2*matrix.P; i++ ) std::cout << matrix.colboundaries[ i ] << " "; std::cout << std::endl;

	//create row categories
	unsigned char *row_categories = new unsigned char[ 2*matrix.P-1 ];
	unsigned int  *row2boundary   = new unsigned int[ matrix.m ];
	for( unsigned int i=0; i<2*matrix.P-1; i++ ) row_categories[i] = 0; //unused category
	for( unsigned int i=0; i<matrix.m; i++ ) {
		//get row2boundary
		const unsigned long int global_i = matrix.row2global[i];
		long int k=-1;
		for( unsigned int j=0; j<2*matrix.P-1; j++ )
			if( matrix.rowboundaries[j]<=global_i && global_i<=matrix.rowboundaries[j+1] )
				k=j;
		if( k==-1 ) {
			std::cerr << "Error: local row " << i << " is outside the matrix!" << std::endl;
			exit( 1 );
		}
		if( k < 0 || k >= 2*(signed long int)matrix.P-1 ) {
			std::cerr << "Error: category index out of range!" << std::endl;
			exit( 1 );
		}
		row2boundary[ i ] = k;

		//infer boundary category
		if( (k%2)==0 ) { //pure block
			if( matrix.row2proc[ i ] != s ) {
				std::cerr << "Error: local row in non-local location!" << std::endl;
				exit( 1 );
			} else if( row_categories[ k ] > 1 ) {
				std::cerr << "Error: pure block has non-pure category!" << std::endl;
				exit( 1 );
			} else
				row_categories[ k ] = 1; //in my pure set 
		} else { //separator
			if( row_categories[ k ] == 0 ) row_categories[ k ] = 2; //local separator by default
			if( matrix.row2proc[ i ] != s && matrix.row2proc[ i ] != neighbour )
				row_categories[ k ] = 3; //this is a full separator
		}
	}
	//check number of pure categories
	unsigned int pure_categories = 0;
	for( unsigned int i=0; i<2*matrix.P-1; i+=2 )
		if( row_categories[ i ] == 1 ) pure_categories++;
	if( pure_categories > 1 ) {
		std::cerr << "Error: there are more than one pure row categories!" << std::endl;
		exit( 1 );
	}
	
	//create column categories
	unsigned char *col_categories = new unsigned char[ 2*matrix.P-1 ];
	unsigned int  *col2boundary   = new unsigned int[ matrix.n ];
	for( unsigned int i=0; i<2*matrix.P-1; i++ ) col_categories[i] = 0; //unused category
	for( unsigned int i=0; i<matrix.n; i++ ) {
		//get col2boundary
		const unsigned long int global_i = matrix.col2global[i];
		long int k=-1;
		for( unsigned int j=0; j<2*matrix.P-1; j++ )
			if( matrix.colboundaries[j]<=global_i && global_i<=matrix.colboundaries[j+1] )
				k=j;
		if( k==-1 ) {
			std::cerr << "Error: local column " << i << " is outside the matrix!" << std::endl;
			exit( 1 );
		}
		col2boundary[ i ] = k;
		//infer boundary category
		if( k%2==0 ) { //pure block
			if( matrix.col2proc[ i ] != s ) {
				std::cerr << "Error: local column in non-local location!" << std::endl;
				exit( 1 );
			} else if ( col_categories[ k ] > 1 ) {
				std::cerr << "Error: pure block has non-pure category!" << std::endl;
				exit( 1 );
			} else
				col_categories[ k ] = 1; //in my pure set 
		} else { //separator
			if( col_categories[ k ] == 0 ) col_categories[ k ] = 2; //local separator by default
			if( matrix.col2proc[ i ] != s && matrix.col2proc[ i ] != neighbour )
				col_categories[ k ] = 3; //this is a full separator
		}
	}

	//check number of pure categories
	pure_categories = 0;
	for( unsigned int i=0; i<2*matrix.P-1; i+=2 )
		if( col_categories[ i ] == 1 ) pure_categories++;
	if( pure_categories > 1 ) {
		std::cerr << "Error: there are more than one pure col categories!" << std::endl;
		exit( 1 );
	}
	
	//perform permutations, splits
	std::vector< Triplet< double > > replacement;
	std::vector< Triplet< double > >::iterator it = matrix.local.begin();
	//save old size, for applying permutations later on
	unsigned long int old_m = matrix.m;
	unsigned long int old_n = matrix.n;
	matrix.m = matrix.n = 0; //reset local matrix dimensions
	for( ; it != matrix.local.end(); ++it ) {
		if( it->i() >= old_m || it->j() >= old_n ) {
			std::cerr << "Invalid nonzero location (" << it->i() << "," << it->j() << ")!" << std::endl;
			exit( 1 );
		}
		//apply permutation FIXME perm
		unsigned long int new_i = it->i();
		unsigned long int new_j = it->j();

		//determine nonzero category
		const unsigned char row_category = row_categories[ row2boundary[ new_i ] ];
		const unsigned char col_category = col_categories[ col2boundary[ new_j ] ];
		if( ( row_category == 0 && col_category == 0 ) || 
		    ( row_category >  3 && col_category >  3 ) ) { //error
			std::cerr << "Unknown categories encountered!" << std::endl;
			exit( 1 );
		}
		//move to correct submatrix
		Triplet< double > newentry = Triplet< double >( new_i, new_j, it->value );
		if( row_category == 1 && col_category == 1 ) { //purely local
			replacement.push_back( newentry );
			if( new_i > matrix.m ) matrix.m = new_i;
			if( new_j > matrix.n ) matrix.n = new_j;
		} else if( row_category == 1 ) {
			if( col_category == 2 ) { //pure fanin with neighbour
				matrix.sepA1.data.push_back( newentry );
				if( new_i > matrix.sepA1.m ) matrix.sepA1.m = new_i;
				if( new_j > matrix.sepA1.n ) matrix.sepA1.n = new_j;
			} else if( col_category == 3 ) { //pure fanin with all processors
				matrix.sepB1.data.push_back( newentry );
				if( new_i > matrix.sepB1.m ) matrix.sepB1.m = new_i;
				if( new_j > matrix.sepB1.n ) matrix.sepB1.n = new_j;
			}
		} else if( col_category == 1 ) {
			if( row_category == 2 ) {
				matrix.sepA2.data.push_back( newentry ); //pure fanout with neighbour
				if( new_i > matrix.sepA2.m ) matrix.sepA2.m = new_i;
				if( new_j > matrix.sepA2.n ) matrix.sepA2.n = new_j;
			} else if( row_category == 3 ) {
				matrix.sepB2.data.push_back( newentry ); //pure fanout with all processors
				if( new_i > matrix.sepB2.m ) matrix.sepB2.m = new_i;
				if( new_j > matrix.sepB2.n ) matrix.sepB2.n = new_j;
			}
		} else if( row_category == 2 && col_category == 2 ) {
			matrix.sepA3.data.push_back( newentry ); //mixed fanin/fanout with neighbour
			if( new_i > matrix.sepA3.m ) matrix.sepA3.m = new_i;
			if( new_j > matrix.sepA3.n ) matrix.sepA3.n = new_j;
		} else {
			matrix.sepB3.data.push_back( newentry ); //mixed fanin/fanout with all processors
			if( new_i > matrix.sepB3.m ) matrix.sepB3.m = new_i;
			if( new_j > matrix.sepB3.n ) matrix.sepB3.n = new_j;
		}
	}
	//replace local
	matrix.local = replacement;
	//correct counts
	matrix.m++; matrix.n++;
	matrix.sepA1.m++; matrix.sepA2.m++; matrix.sepA3.m++;
	matrix.sepA1.n++; matrix.sepA2.n++; matrix.sepA3.n++;
	matrix.sepB1.m++; matrix.sepB2.m++; matrix.sepB3.m++;
	matrix.sepB1.n++; matrix.sepB2.n++; matrix.sepB3.n++;

	//localise boundary matrices
	printproc(s);std::cout << ": Localising A1" << std::endl;
	localise( matrix.sepA1, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	printproc(s);std::cout << ": Localising A2" << std::endl;
	localise( matrix.sepA2, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	printproc(s);std::cout << ": Localising A3" << std::endl;
	localise( matrix.sepA3, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	printproc(s);std::cout << ": Localising B1" << std::endl;
	localise( matrix.sepB1, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	printproc(s);std::cout << ": Localising B2" << std::endl;
	localise( matrix.sepB2, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	printproc(s);std::cout << ": Localising B3" << std::endl;
	localise( matrix.sepB3, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	//localise main matrix
	printproc(s);std::cout << ": Localising A" << std::endl;
	printproc(s);std::cout << ": Allocating inter-processor separator buffer" << std::endl;
	matrix.sepA1.allocate();
	matrix.sepA2.allocate();
	matrix.sepA3.allocate();
	matrix.sepB1.allocate();
	matrix.sepB2.allocate();
	matrix.sepB3.allocate();
	DMVData temp;
	temp.data = matrix.local; //copy my data
	temp.m    = matrix.m;
	temp.n    = matrix.n;
	localise( temp, matrix.row2proc, matrix.col2proc, matrix.row2index, matrix.col2index );
	//copy back
	matrix.local = temp.data;
	matrix.m     = temp.m;
	matrix.n     = temp.n;
	delete [] matrix.row2proc;
	delete [] matrix.col2proc;
	delete [] matrix.row2index;
	delete [] matrix.col2index;
	matrix.row2proc = temp.row2proc;
	matrix.col2proc = temp.col2proc;
	matrix.row2index = temp.row2index;
	matrix.col2index = temp.col2index;
	temp.row2proc = temp.col2proc = temp.row2index = temp.col2index = NULL;

	//done, clean up
	delete [] row2boundary;
	delete [] col2boundary;
	delete [] row_categories;
	delete [] col_categories;
	return neighbour;
}

MVData **shared_pool = NULL;
double **shared_x = NULL;
double **shared_y = NULL;

void spmv( const unsigned int s, const unsigned int t,
		ICRS< double > &A, ICRS< double > &S1,
		CCSWrapper< double, ICRS< double >, ULI > &S2,
		ICRS< double > &S3, ICRS< double > &S4,
		CCSWrapper< double, ICRS< double >, ULI > &S5,
		ICRS< double > &S6,
		DMVData &D1, DMVData &D2,
		DMVData &D3, DMVData &D4,
		DMVData &D5, DMVData &D6 ) {
	//superstep 1, first set shared vars
	double *x = shared_x[s];
	double *y = shared_y[s];
	//fanin with direct get
	for( unsigned long int j=0; j<D1.n; j++ )
		D1.x_buffer[ j ] = shared_x[ D1.col2proc[j] ][ D1.col2index[ j ] ];
	for( unsigned long int j=0; j<D2.n; j++ )
		D2.x_buffer[ j ] = shared_x[ D2.col2proc[j] ][ D2.col2index[ j ] ];
	for( unsigned long int j=0; j<D3.n; j++ )
		D3.x_buffer[ j ] = shared_x[ D3.col2proc[j] ][ D3.col2index[ j ] ];
	//local mv
//	printproc(s);std::cout << "MV on " << A.m() << " by " << A.n() << " matrix; ";
//	printproc(s);std::cout << "Corresponding vector sizes: " << data->output_length << ", " << data->input_length << std::endl;
	A.zax( x, y );
	//fan-in mv's
	D1.prepare();
	D2.prepare();
	D3.prepare();
	S1.zax( D1.x_buffer, D1.y_buffer );
	S2.zax( D2.x_buffer, D2.y_buffer );
	S3.zax( D3.x_buffer, D3.y_buffer );
	//local translation
	for( unsigned long int j=0; j<D1.m; j++ )
		y[ D1.row2index[ j ] ] += D1.y_buffer[ j ];
	for( unsigned long int j=0; j<D2.m; j++ )
		y[ D2.row2index[ j ] ] += D2.y_buffer[ j ];
	for( unsigned long int j=0; j<D3.m; j++ )
		y[ D3.row2index[ j ] ] += D3.y_buffer[ j ];
	//end superstep (ensure local MVs are done)
	if( s==p_translate[0] || s==p_translate[1] )
		bsp_sync1();
	else
		bsp_sync2();
	//fanout communication starts now with neighbour only
	for( unsigned long int j=0; j<shared_pool[t]->sepA2.m; j++ ) {
		if( shared_pool[t]->sepA2.row2proc[j] == s )
			y[ shared_pool[t]->sepA2.row2index[j] ] += shared_pool[t]->sepA2.y_buffer[j];
	}
	for( unsigned long int j=0; j<shared_pool[t]->sepA3.m; j++ ) {
		if( shared_pool[t]->sepA3.row2proc[j] == s )
			y[ shared_pool[t]->sepA3.row2index[j] ] += shared_pool[t]->sepA3.y_buffer[j];
	}
	//direct get from all
	for( unsigned long int j=0; j<D4.n; j++ )
		D4.x_buffer[ j ] = shared_x[ D4.col2proc[j] ][ D4.col2index[ j ] ];
	for( unsigned long int j=0; j<D5.n; j++ )
		D5.x_buffer[ j ] = shared_x[ D5.col2proc[j] ][ D5.col2index[ j ] ];
	for( unsigned long int j=0; j<D6.n; j++ )
		D6.x_buffer[ j ] = shared_x[ D6.col2proc[j] ][ D6.col2index[ j ] ];
	D4.prepare();
	D5.prepare();
	D6.prepare();
	S4.zax( D4.x_buffer, D4.y_buffer );
	S5.zax( D5.x_buffer, D5.y_buffer );
	S6.zax( D6.x_buffer, D6.y_buffer );
	//local copy
	for( unsigned long int j=0; j<D4.m; j++ )
		if( D4.row2proc[j]==s )
			y[ D4.row2index[j] ] += D4.y_buffer[ j ];
	for( unsigned long int j=0; j<D5.m; j++ )
		if( D5.row2proc[j]==s )
			y[ D5.row2index[j] ] += D5.y_buffer[ j ];
	for( unsigned long int j=0; j<D6.m; j++ )
		if( D6.row2proc[j]==s )
			y[ D6.row2index[j] ] += D6.y_buffer[ j ];
	bsp_sync(); //sync with all
	for( unsigned int k=0; k<P; k++ ) {
		if( k==s ) continue;
		for( unsigned long int j=0; j<shared_pool[k]->sepB1.m; j++ ) {
			if( shared_pool[k]->sepB1.row2proc[j] == s )
				y[ shared_pool[k]->sepB1.row2index[j] ] += shared_pool[k]->sepB1.y_buffer[j];
		}
		for( unsigned long int j=0; j<shared_pool[k]->sepB2.m; j++ ) {
			if( shared_pool[k]->sepB2.row2proc[j] == s )
				y[ shared_pool[k]->sepB2.row2index[j] ] += shared_pool[k]->sepB2.y_buffer[j];
		}
		for( unsigned long int j=0; j<shared_pool[k]->sepB3.m; j++ ) {
			if( shared_pool[k]->sepB3.row2proc[j] == s )
				y[ shared_pool[k]->sepB3.row2index[j] ] += shared_pool[k]->sepB3.y_buffer[j];
		}
	}
}

void *parallel( void *threadid ) {
	MVData matrix;
	struct timespec start, stop;
	double time;
	unsigned int s = (unsigned int) ((long)threadid);

	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start );
	for( unsigned int i=0; i<REPEAT_EXP*REPEAT_EXP; i++ ) bsp_sync();
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop );
	time = (stop.tv_sec-start.tv_sec) * 1000;
	time+= (stop.tv_nsec-start.tv_nsec)/1000000.0;
	if( s==0 ) { printproc(s); std::cout << ": local time estimate for l = " << (time/((double)REPEAT_EXP*REPEAT_EXP)) << "ms." << std::endl; }
	bsp_sync();
	printproc(s); std::cout << ": Executing parse" << std::endl;
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start );
	parse( fn, s, matrix );
	printproc(s); std::cout << ": Doing vector-based permutations" << std::endl;
	const unsigned long int neighbour = prepare( matrix, s );
	printproc(s); std::cout << ": Building data structures" << std::endl;
	printproc(s); std::cout << "Creating matrix with m=" << matrix.m << ", n=" << matrix.n << ", nz=" << matrix.local.size() << std::endl;
	ICRS< double > local( matrix.local, matrix.m, matrix.n, 0 );
	printproc(s); std::cout << "Creating matrix with m=" << matrix.sepA1.m << ", n=" << matrix.sepA1.n << ", nz=" << matrix.sepA1.data.size() << std::endl;
	ICRS< double > sepA1( matrix.sepA1.data, matrix.sepA1.m, matrix.sepA1.n, 0 );
	printproc(s); std::cout << "Creating matrix with m=" << matrix.sepA2.m << ", n=" << matrix.sepA2.n << ", nz=" << matrix.sepA2.data.size() << std::endl;
	CCSWrapper< double, ICRS< double >, ULI > sepA2( matrix.sepA2.data, matrix.sepA2.m, matrix.sepA2.n, 0 );
	printproc(s); std::cout << "Creating matrix with m=" << matrix.sepA3.m << ", n=" << matrix.sepA3.n << ", nz=" << matrix.sepA3.data.size() << std::endl;
	ICRS< double > sepA3( matrix.sepA3.data, matrix.sepA3.m, matrix.sepA3.n, 0 );
	printproc(s); std::cout << "Creating matrix with m=" << matrix.sepB1.m << ", n=" << matrix.sepB1.n << ", nz=" << matrix.sepB1.data.size() << std::endl;
	ICRS< double > sepB1( matrix.sepB1.data, matrix.sepB1.m, matrix.sepB1.n, 0 );
	printproc(s); std::cout << "Creating matrix with m=" << matrix.sepB2.m << ", n=" << matrix.sepB2.n << ", nz=" << matrix.sepB2.data.size() << std::endl;
	CCSWrapper< double, ICRS< double >, ULI > sepB2( matrix.sepB2.data, matrix.sepB2.m, matrix.sepB2.n, 0 );
	printproc(s); std::cout << "Creating matrix with m=" << matrix.sepB3.m << ", n=" << matrix.sepB3.n << ", nz=" << matrix.sepB3.data.size() << std::endl;
	ICRS< double > sepB3( matrix.sepB3.data, matrix.sepB3.m, matrix.sepB3.n, 0 );
	//make available
	shared_pool[ s ] = &matrix;
	std::cout << "Creating local vectors" << std::endl;
	shared_x[ s ] = new double[ matrix.input_length  ];
	shared_y[ s ] = new double[ matrix.output_length ];
	for( unsigned int i=0; i<matrix.input_length;  i++ ) shared_x[s][i] = 1.0;
	for( unsigned int i=0; i<matrix.output_length; i++ ) shared_y[s][i] = 0.0;
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop );
	time = stop.tv_sec-start.tv_sec;
	time+= (stop.tv_nsec-start.tv_nsec)/1000000000.0;
	printproc(s); std::cout << ": time taken for build = " << time << " seconds." << std::endl;
	bsp_sync();
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &start );
	for( int i=0; i<REPEAT_EXP; i++ )
		spmv( s, neighbour, local, sepA1, sepA2, sepA3, sepB1, sepB2, sepB3,
			matrix.sepA1, matrix.sepA2, matrix.sepA3, matrix.sepB1, matrix.sepB2, matrix.sepB3 );
	bsp_sync();
	clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &stop );
	if( s==0 ) {
		time = (stop.tv_sec-start.tv_sec) * 1000;
		time+= (stop.tv_nsec-start.tv_nsec)/1000000.0;
		printproc(s); std::cout << ": local time taken for SpMV = " << (time/((double)REPEAT_EXP)) << "ms." << std::endl;
	}
//	printproc(s); std::cout << ": Exiting parallel code" << std::endl;
	delete [] shared_x[ s ];
	delete [] shared_y[ s ];
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

	shared_pool = new MVData*[ P ];
	shared_x = new double*[ P ];
	shared_y = new double*[ P ];
	p_translate = new unsigned long int[ P ];
	for( unsigned int i=0; i<P; i++ )
		p_translate[i] = (unsigned int)atoi( argv[ 3+i ] );
	inv_p_translate = new unsigned int[ P ];
	for( unsigned int i=0; i<P; i++ )
		inv_p_translate[ p_translate[ i ] ] = i;

	threads = (pthread_t*)malloc( ( P ) * sizeof( pthread_t ) );
	for( unsigned int i=0; i<P; i++ )
		pthread_create( &threads[p_translate[i]], NULL, parallel, (void *)p_translate[i] );
	//wait for exit
	for( unsigned int i=0; i<P; i++ )
		pthread_join( threads[i], NULL );

	delete [] shared_pool;
	delete [] shared_x;
	delete [] shared_y;
	delete [] p_translate;
	delete [] inv_p_translate;
	free(threads);
	return 0;
}

