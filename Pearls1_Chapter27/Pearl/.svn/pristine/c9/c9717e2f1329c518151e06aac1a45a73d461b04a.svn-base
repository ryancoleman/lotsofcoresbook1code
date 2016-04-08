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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2007.
 */


#include "Triplet.hpp"
#include "SparseMatrix.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>
#include <cmath>

//#define _DEBUG

#ifndef _H_ICRS
#define _H_ICRS

#ifdef _DEBUG
#include<iostream>
#endif

/**
 *  The *incremental* compressed row storage sparse matrix data structure.
 */
template< typename T, typename _i_value=ULI >
class ICRS: public SparseMatrix< T, ULI > {

  private:

  protected:

	/** Array containing the actual nnz non-zeros. */
	T* ds;

	/** Start position, row */
	ULI r_start;

	/** Start position, column */
	ULI c_start;

	/** Array containing the column jumps. */
	_i_value* c_ind;

	/** Array containing the row jumps. */
	_i_value* r_ind;

	/** Remembers the number of bytes allocated. */
	size_t bytes;

  public:

	/** Fill-in field for interoperability with oICRS. Fill-in is always 0 for regular ICRS. */
	static const size_t fillIn = 0;

	/** Comparison function used for sorting input data. */
	static int compareTriplets( const void * left, const void * right ) {
                const Triplet< T > one = *( (Triplet< T > *)left );
                const Triplet< T > two = *( (Triplet< T > *)right );
                if( one.i() < two.i() )
                        return -1;
                if ( one.i() > two.i() )
                        return 1;
                if ( one.j() < two.j() )
                        return -1;
                if ( one.j() > two.j() )
                        return 1;
                return 0;
        }

	/** Base constructor. */
	ICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	ICRS( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid ICRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows of the matrix to be stored.
	 *  @param number_of_cols The number of columns of the matrix to be stored.
	 *  @param zero Which element is considered to be the zero element.
	 */
	ICRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ):
		SparseMatrix< T, ULI >( number_of_nonzeros, number_of_cols, number_of_rows, zero ) {
		ds = new T[ this->nnz ];
		c_ind = new _i_value[ this->nnz ];
		r_ind = new _i_value[ this->nnz ];
		bytes = sizeof( ULI ) * 2 + sizeof( _i_value ) * 2 * this->nnz + sizeof( T ) * this->nnz;
	}

	/**
	 *  Copy constructor.
	 *  @param toCopy Reference to the ICRS datastructure to copy.
	 */
	ICRS( ICRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->noc = toCopy.noc;
		this->nor = toCopy.nor;
		this->r_start = toCopy.r_start;
		this->c_start = toCopy.c_start;
		ds = new T[ this->nnz ];
		c_ind = new _i_value[ this->nnz - 1 ];
		r_ind = new _i_value[ this->nnz - 1 ];
		for( ULI i=0; i<this->nnz; i = i + 1 ) {
			ds[ i ] = toCopy.ds[ i ];
			c_ind[ i ]= toCopy.c_ind[ i ];
			r_ind[ i ] = toCopy.r_ind[ i ];
		}
		bytes = sizeof( ULI ) * 2 + sizeof( _i_value ) * 2 * this->nnz + sizeof( T ) * this->nnz;
	}

	/**
	 *  Constructor which transforms a collection of input triplets to CRS format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indices. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *
	 *  @param input The input collection of triplets (i,j,val).
	 *  @param m The number of rows of the input matrix.
	 *  @param n The number of columns of the input matrix.
	 *  @param zero Which element is considered zero.
	 */
	ICRS( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero = 0 ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		this->zero_element = zero;

		this->nor = m;
		this->noc = n;
	
		if( log2(this->noc) > sizeof( _i_value ) * 8 )
			std::cerr << "Warning: the matrix with " << this->noc << " columns cannot be represented within " << (sizeof( _i_value )*8) << "-bit index values this ICRS instance uses!" << std::endl;

		if( log2(this->nor) > sizeof( _i_value ) * 8 )
			std::cerr << "Warning: the matrix with " << this->nor << " rows cannot be represented within " << (sizeof( _i_value )*8) << "-bit index values this ICRS instance uses!" << std::endl;

		if( m==0 || n==0 || input.size() == 0 ) { //empty matrix
			this->nor = this->noc = this->nnz = 0;
			ds = NULL;
			r_ind = NULL;
			c_ind = NULL;
			return;
		}

		typename std::vector< Triplet< T > >::iterator in_it;

		//WARNING: noc*nor typically overflows on 32 bits!
		//         immediately start recording differences
		//         instead of trying to directly calculate
		//         the index as i*noc+j.

		//Complexity compiler-package dependent. Probably O(nnz log(nnz)) average, O(nnz^2) worst-case.
		//for standard C++ sort:
		//sort( input.begin(), input.end(), compareTriplets );
		//for C quicksort:
		qsort( &input[ 0 ], input.size(), sizeof( Triplet< T > ), &compareTriplets );

		//filtering out zeros is skipped for now.
		
		//Count the number of row jumps
		std::vector< ULI > r_ind_temp;
		typename std::vector< Triplet< T > >::iterator it = input.begin();
		ULI prev_row = (*it).i();
		r_ind_temp.push_back( prev_row );
		it++;

		//O(nnz log(nnz)); push_back on a vector uses amortised log(n) array growing algorithm on inserts
		for( ; it!=input.end(); it++ ) {
			assert( (*it).i() < this->nor );
			assert( (*it).j() < this->noc );
			assert( (*it).i() >= prev_row );
			if( (*it).i() > prev_row ) {
				assert( (*it).i() - prev_row < m );
				r_ind_temp.push_back( (*it).i() - prev_row );
				prev_row = (*it).i();
			}
		}
		this->nnz = input.size();

		//allocate arrays
		const unsigned long int allocsize =  this->nnz;
		ds    = new T[ allocsize ];
		c_ind = new _i_value[ allocsize ];
		r_ind = new _i_value[ r_ind_temp.size() ];

		//record #bytes used
		bytes = sizeof( ULI ) * 2 + sizeof( _i_value ) * ( allocsize + r_ind_temp.size() ) + sizeof( T ) * allocsize;

		//set last entry
		c_ind[ allocsize - 1 ] = this->noc;
		//r_ind does not have to be set; altough the last element is read, it is actually never used.

		//copy row-jump vector
		//O(m) worst case
		r_start = r_ind_temp[ 0 ];
		for( ULI i=1; i<r_ind_temp.size(); i++ ) {
			r_ind[i-1] = r_ind_temp[ i ];
			if( static_cast< ULI >( r_ind[i-1] ) != r_ind_temp[ i ] ) {
				std::cerr << "Row increment too large to store in this ICRS instance!" << std::endl;
				exit( 1 );
			}
		}

		//make ICRS
		prev_row = r_start;
		ULI prev_col = c_start = input[ 0 ].j(); //now r_- and c_-start have been set
		//O(nnz)
		unsigned long int check_jumps = 0;
		unsigned long int check_row   = r_ind_temp[0];
		assert( r_ind_temp[ 0 ] == r_start );
		ds[ 0 ] = input[ 0 ].value;
		for( ULI i=1; i<this->nnz; ++i ) {
			const Triplet< T > cur = input[ i ];
			const ULI currow = cur.i();
			const ULI curcol = cur.j();
			if( currow == prev_row ) {
				c_ind[i-1] = curcol - prev_col;
				if( static_cast< ULI >( c_ind[i-1] ) != curcol - prev_col ) {
					std::cerr << "Column increment too large to store in this ICRS instance!" << std::endl;
					exit( 1 );
				}
				assert( currow == check_row );
			} else {
				assert( currow > prev_row );
				c_ind[i-1] = this->noc + ( curcol - prev_col );
				if( static_cast< ULI >( c_ind[i-1] ) - this->noc + prev_col != curcol ) {
					std::cerr << "Overflowed column increment too large to store in this ICRS instance!" << std::endl;
					exit( 1 );
				}
				check_row += r_ind[ check_jumps++ ];
				assert( currow == check_row );
				prev_row = currow;
			}
			ds[ i ] = cur.value;
			prev_col = curcol;

#ifdef _DEBUG
			std::cout << currow << "," << curcol << "(" << cur.value << ") maps to " << c_ind[ i ] << std::endl;
#endif
		}

		//append with zeroes
		for( unsigned long int i=this->nnz; i<allocsize; ++i ) {
			c_ind[ i - 1 ] = 0;
			ds[ i ] = 0;
		}

		//assert row jumps is equal to r_ind_temp's size
		assert( check_jumps == r_ind_temp.size()-1 );
		
		//clear temporary r_ind vector
		r_ind_temp.clear();
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = this->r_start;
		col = this->c_start;
	}

	/** Gets starting position (first nonzero location) */
	void getStartingPos( ULI &row_start, ULI &column_start ) {
		row_start = this->r_start;
		column_start = this->c_start;
	}

	/**
	 *  Sets starting position of matrix multiplication.
	 *  (Useful for example when the input/output vectors will be
	 *   shifted before passed on to this class' zax method.)
	 *
	 *  @param row_start    New row start location
	 *  @param column_start New column start location
	 */
	void setStartingPos( const ULI row_start, const ULI column_start ) {
		assert( row_start <= this->r_start );
		assert( column_start <= this->c_start );
		this->r_start = row_start;
		this->c_start = column_start;
	}

	/**
	 *  In-place z=xA function. Adapted from the master thesis of Joris Koster, Utrecht University.
	 *
	 *  @param pDataX Pointer to array x to multiply by the current matrix (Ax).
	 *  @param pDataZ Pointer to result array. Must be pre-allocated and its elements set to zero for correct results.
	 */
        virtual void zxa( const T*__restrict__ pDataX, T*__restrict__ pDataZ ) {
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) return;
		         T *__restrict__ pDataA    = ds;
		const    T *__restrict__ pDataAend = ds + this->nnz;
		//const    T *__restrict__ const pDataXend = pDataX + this->nor;
		const T * const pDataZend = pDataZ + this->noc;
		//const    T *pDataZend = z + nor; //unused

		_i_value *__restrict__ pIncRow   = r_ind;
		_i_value *__restrict__ pIncCol   = c_ind;

		//go to first position
		pDataZ += this->c_start;
		pDataX += this->r_start;
		while( pDataA < pDataAend ) {
			while( pDataZ < pDataZend ) {
				*pDataZ += *pDataA * *pDataX;
				pDataA++;
				pDataZ  += *pIncCol;
				pIncCol++;
			}
			pDataZ -= this->noc;
			pDataX += *pIncRow++;
		}
        }

	/**
	 *  In-place z=Ax function. Adapted from the master thesis of Joris Koster, Utrecht University.
	 *
	 *  @param pDataX Pointer to array x to multiply by the current matrix (Ax).
	 *  @param pDataZ Pointer to result array. Must be pre-allocated and its elements set to zero for correct results.
	 */
        virtual void zax( const T*__restrict__ pDataX, T*__restrict__ pDataZ ) {
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) return;
		//go to first position
		      T *__restrict__ pDataA    = ds;
		const T * const       pDataAend = ds     + this->nnz;
		const T * const       pDataXend = pDataX + this->noc;
#ifndef NDEBUG
		const T * const pDataXst  = pDataX;
		const T * const pDataZst  = pDataZ;
		const T * const pDataZend = pDataZ + this->nor;
#endif

		_i_value *__restrict__ pIncRow   = r_ind;
		_i_value *__restrict__ pIncCol   = c_ind;

		//go to first column
		assert( r_start < this->nor );
		assert( c_start < this->noc );
		assert( this->nnz > 0 );
		pDataX += c_start;
		pDataZ += r_start;
		while( pDataA < pDataAend ) {
			while( pDataX < pDataXend ) {
				assert( pDataA < pDataAend );
				assert( pDataX >= pDataXst );
				assert( pDataX < pDataXend );
				assert( pDataZ >= pDataZst );
				assert( pDataZ < pDataZend );
				assert( pDataX + this->noc >= pDataXend ); //otherwise pDataX is before the start of x!
				assert( pDataA + this->nnz >= pDataAend );
				assert( pDataZ + this->nor >= pDataZend );

				*pDataZ += *pDataA++ * *pDataX;
				 pDataX += *pIncCol++;
			}
			pDataX -= this->noc;
			//jump to correct row
			pDataZ += *pIncRow++;
		}
        }

	/** Base deconstructor. */
	~ICRS() {
		if( ds    != NULL ) delete [] ds;
		if( c_ind != NULL ) delete [] c_ind;
		if( r_ind != NULL ) delete [] r_ind;
	}

	virtual size_t bytesUsed() {
		return bytes;
	}

};

#endif

