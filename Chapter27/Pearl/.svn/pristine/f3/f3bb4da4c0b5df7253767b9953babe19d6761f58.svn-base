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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2008.
 */


#include "alignment.hpp"
#include "Triplet.hpp"
#include "SparseMatrix.hpp"

#include <map>
#include <set>
#include <cmath>
#include <vector>
#include <assert.h>
#include <algorithm>

#if defined __INTEL_COMPILER && defined __KNC__
 #include <stdint.h>
 #include <immintrin.h>
 #include <smmintrin.h>
 #include <zmmintrin.h>
#endif

//#define _DEBUG

#ifndef _H_oICRS
#define _H_oICRS

#ifdef _DEBUG
#include<iostream>
#endif

/**
 *  The *incremental* compressed row storage sparse matrix data structure.
 */
template< typename T, size_t block_length_row=1, size_t block_length_column=8, typename _i_value=LI >
class oICRS: public SparseMatrix< T, ULI > {

  private:

  protected:

	/** Combined blocking size. */
	static const size_t block_length = block_length_row * block_length_column;

	/** Array containing the actual nnz non-zeros. */
	T* ds;

	/** Array containing the column jumps. */
	_i_value* c_ind;

	/** Array containing the row jumps. */
	_i_value* r_ind;

	/** Remembers the number of bytes allocated. */
	size_t bytes;

	/** The number of nonzeroes allocated (may differ from the actual number of nonzeroes). */
	size_t allocsize;

	/** The number of row jumps plus one; i.e., the length of r_ind. */
	size_t jumpsize;

	/**
	 * Utility function: allocate memory areas according @allocsize and @jumpsize.
	 */
	void allocate() {
#ifdef __INTEL_COMPILER
		ds    = (T*) _mm_malloc( this->allocsize * sizeof( T ), _SL_ALIGN_DOUBLE );
		c_ind = (_i_value*) _mm_malloc( this->allocsize * sizeof( _i_value ), _SL_ALIGN_INT16 );
		r_ind = (_i_value*) _mm_malloc( (this->jumpsize + 1) * sizeof( _i_value ), _SL_ALIGN_INT16 );
#else
		ds    = new T[ this->allocsize ];
		c_ind = new _i_value[ this->allocsize ];
		r_ind = new _i_value[ this->jumpsize + 1 ];
#endif
		bytes = sizeof( _i_value ) * (this->jumpsize + this->allocsize + 1) + sizeof( T ) * this->allocsize;
	}

	/**
	 * Utility function: free memory areas according @allocsize and @jumpsize.
	 */
	void deallocate() {
#ifdef __INTEL_COMPILER
		_mm_free( ds );
		_mm_free( c_ind );
		_mm_free( r_ind );
#else
		if( ds    != NULL ) delete [] ds;
		if( c_ind != NULL ) delete [] c_ind;
		if( r_ind != NULL ) delete [] r_ind;
#endif
	}

	/**
	 * Helper function for oICRS construction.
	 * 
	 *  @param blockingSize Pads to this boundary (sensible values are block_length_row or block_length)
	 */
	static void addPadding( std::vector< _i_value > &row_increments, const std::map< size_t, size_t > &row2local, const size_t blockingSize, const size_t prev_row, const size_t m ) {
		if( row_increments.size() % blockingSize != 0 ) {
			//derive the current column index, to prepare for reverse iteration
			size_t base = row2local.rbegin()->first;
			//build a set of row indices in the full current block_length rows
			std::set< size_t > current;
			//get current index in row increments vector
			size_t index = row_increments.size() - 1;
			//walk through current block, in reverse
			while( index % block_length != block_length - 1 ) {
				//DBG
				//std::cout << "Increment at index " << index << " is " << row_increments[ index ] << " and current base equals " << base << std::endl;
				//add this row index
				current.insert( base );
				//derive previous index base
				base -= row_increments[ index ];
				//decrement index
				--index;
			}
			//derive padding index; closest unused row index not in current set
			size_t padding_index = *(current.begin());
			if( padding_index > 0 ) {
				--padding_index;
			}
			const size_t firstguess = padding_index;
			//DBG
			//std::cout << "First guess for padding index: " << firstguess << std::endl;
			//keep incrementing padding index until it is no longer in current row set
			while( current.find( padding_index ) != current.end() ) {
				++padding_index;
			}
			//check for overflow
			if( padding_index >= m ) {
				//reset index
				padding_index = firstguess;
				//search in the other direction
				while( current.find( padding_index ) != current.end() && padding_index < m ) {
					--padding_index;
				}
				//sanity check
				if( padding_index >= m ) {
					std::cerr << "Matrix too small to be vectorised! Using this matrix in computations may crash due to illegal reads (by one data element after the output vector end). Overallocating the output vector will ensure proper, error-free execution." << std::endl;
				}
			}
			//DBG
			//std::cout << "Final guess for padding index: " << padding_index << ", blocking on boundaries of " << blockingSize << std::endl;
			//add the padding
			row_increments.push_back( padding_index - prev_row );
			//DBG
			//std::cout << "Pushed back: " << row_increments[ row_increments.size() - 1 ] << std::endl;
			//pad the remainder
			while( row_increments.size() % blockingSize != 0 ) {
				//stay on that row
				row_increments.push_back( 0 );
				//DBG
				//std::cout << "Pushed back: " << row_increments[ row_increments.size() - 1 ] << std::endl;
				//sanity check
				assert( row2local.find( padding_index ) == row2local.end() );
			}
		}
	}

	/**
	 *  Helper method for oBICRS constructor.
	 *
	 *  Takes an input array of nnz triplets. This function looks at all triplets from index k on,
	 *  takes the first block_length_row different rows that are encountered,
	 *  and adds the nonzeroes on those row with a maximum of eight per row.
	 *
	 *  If a new row index is encountered but block_length_rows were already added, this function
	 *  stops iterating.
	 *
	 *  If a new column index on a given row is encountered, but that row already contains
	 *  block_length_col elements, then this function also stops iterating.
	 *
	 *  The above two constraints ensure the nonzeroes are always processed in the intended order.
	 *
	 *  @param row_indices an array of block_length_row sets that will contain the nonzero indices that will be added to each block.
	 *  @param row_increments Keep track of row increments.
	 *  @param prev_row The global index of the last previously added row.
	 *  @param input the input triplet array.
	 *  @param k From which k on to start iterating.
	 *  @param m The number of matrix rows.
	 *
	 *  @return The k at which iteration stopped (exclusive).
	 */
	static size_t prepareBlockIteration(
		std::vector< size_t > * const row_indices,
		std::vector< _i_value > &row_increments,
		_i_value &prev_row,
		const std::vector< Triplet< double > > &input,
		const size_t k,
		const size_t m
	) {
		std::map< size_t, size_t > row2local; //used to map global row indices to the block row indices
		size_t k_walk = k + 1;                //the nonzero iterator
		row2local[ input[ k ].i() ] = 0;      //initially, only the row of the first nonzero of this block is in the map.
		row_indices[ 0 ].push_back( k );      //initially, only the first nonzero index is in this block.
		//walk through input as long as we can add nonzeroes to the current block.
		while( k_walk < input.size() ) {
			//current row
			const size_t currow = input[ k_walk ].i();
			//check if the new row is already mapped
			std::map< size_t, size_t >::const_iterator it = row2local.find( currow );
			if( it == row2local.end() ) {
				//new row will be at this block-local location
				const size_t location = row2local.size();
				//if there is no more room for a new row, exit the loop
				if( location == block_length_row ) {
					break;
				}
				//record next nonzero index on this row 
				row_indices[ location ].push_back( k_walk );
				//add to map
				row2local[ currow ] = location;
			} else {
				//add this nonzero to local row set
				row_indices[ it->second ].push_back( k_walk );
			}
			//go to next nonzero
			++k_walk;
		}
		//iterate through all local rows.
		std::map< size_t, size_t >::const_iterator cit = row2local.begin();
		//sanity check; there should be at least one row
		assert( cit != row2local.end() );
		//add increments for each local row
		for( ; cit != row2local.end(); ++cit ) {
			//DBG
			//std::cout << "Current row is " << cit->first << " previous row was " << prev_row << std::endl;
			//put row_incremenet with relative index
			row_increments.push_back( cit->first - prev_row );
			//DBG
			//std::cout << "Pushed back " << row_increments[ row_increments.size() - 1 ] << std::endl;
			//update prev_row
			prev_row = cit->first;
		}

		//add padding increments
		addPadding( row_increments, row2local, block_length_row, prev_row, m );

		//check if we are at the end
		if( k_walk == input.size() ) {
			//add final padding to row_increments
			addPadding( row_increments, row2local, block_length, prev_row, m );
		}

		//DBG
		/*std::cout << "blk_indices = {" << std::endl;
		for( size_t blockrow = 0; blockrow < block_length_row; ++blockrow ) {
			std::cout << "                 ";
			std::vector< size_t >::const_iterator it = row_indices[ blockrow ].begin();
			for( ; it != row_indices[ blockrow ].end(); ++it ) {
				std::cout << *it << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "              }" << std::endl;*/
		//return current nonzero that did not get included into this block
		return k_walk;
	}

	void postProcessRowIncrements( std::vector< _i_value > &r_ind ) {
		//check if we vectorise in row direction
		if( block_length_row == 1 ) {
			//no, so no post-process
			return;
		}
		if( block_length_row == 2 ) {
			//assume regular, non-blocked in z-direction, is faster
			return;
		}

		//sanity check
		assert( r_ind.size() % block_length == 0 );

		//make each block of block_length increments cumulative
		size_t base  = r_ind[ 0 ];
		size_t index = r_ind[ 0 ];
		for( size_t i = 0; i < r_ind.size(); ) {
			//make suitable for vectorisation on output vector
			index += r_ind[ i ];
			r_ind[ i ] = index - base;
			base = index;
			//DBG
			//std::cout << "At row base " << base << ", relative increments are:" << std::endl;
			for( size_t j = 1; j < block_length; ++j ) {
				index += r_ind[ i + j ];
				//DBG
				//std::cout << "\t" << r_ind[ i + j ] << "-->" << (index - base) << " (index presumed to be " << index << ")" << std::endl;
				r_ind[ i + j ] = index - base;
			}
#if defined __INTEL_COMPILER && defined __KNC__
			//permute block for Xeon Phi friendliness, in case of block_length_row > 1
			if( block_length_row > 1 && block_length_row < block_length ) {
				//commit to cyclic distribution with p = block_length_row
				_i_value temp[ block_length ];
				//copy to temporary array
				for( size_t j = 0; j < block_length; ++j ) {
					temp[ j ] = r_ind[ i + j ];
				}
				//commit
				for( size_t s = 0; s < block_length_row; ++s ) {
					for( size_t t = s; t < block_length; t += block_length_row ) {
						r_ind[ i++ ] = temp[ t ];
					}
				}
			} else {
				i += block_length;
			}
#else
			i += block_length;
#endif
		}
	}

  public:

	/** Fill-in (explicitly added nonzeroes to enable vectorisation). */
	size_t fillIn;

	/** Comparison function used for sorting input data. */
	static int compareTriplets( const void * left, const void * right ) {
                const Triplet< T > one      = *( (Triplet< T > *)left );
                const Triplet< T > two      = *( (Triplet< T > *)right );
                if( one.i() < two.i() )
                        return -1;
                if( one.i() > two.i() )
                        return 1;
                if( one.j() < two.j() )
                        return -1;
                if( one.j() > two.j() )
                        return 1;
                return 0;
        }

	/** Comparison function used for sorting vectorised blocks; in-row column ordering */
	static int pColumnSort( const void * left, const void * right ) {
		const Triplet< T >* one = *( (Triplet< T >**)left  );
		const Triplet< T >* two = *( (Triplet< T >**)right );
		if( one->j() < two->j() )
			return -1;
		if( one->j() > two->j() )
			return  1;
		return 0;
	}

	/** Comparison function used for sorting vectorised blocks; row ordering */
	static int pRowSort( const void * left, const void * right ) {
		const Triplet< T >* one = *( (Triplet< T >**)left  ); //takes the first triplet of the row
		const Triplet< T >* two = *( (Triplet< T >**)right );
		if( one->j() < two->j() )
			return -1;
		if( one->j() > two->j() )
			return  1;
		return 0;
	}

	/** Base constructor. */
	oICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	oICRS( std::string file, T zero = 0 ) {
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
	oICRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ):
		SparseMatrix< T, ULI >( number_of_nonzeros, number_of_cols, number_of_rows, zero ) {
		this->allocsize = this->jumpsize = this->nnz;
		allocate();
	}

	/**
	 *  Copy constructor.
	 *  @param toCopy Reference to the ICRS datastructure to copy.
	 */
	oICRS( oICRS< T >& toCopy ) {
		//rather straightforward implementation
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->noc = toCopy.noc;
		this->nor = toCopy.nor;
		this->r_start = toCopy.r_start;
		this->c_start = toCopy.c_start;
		this->allocsize = toCopy.allocsize;
		this->jumpsize  = toCopy.jumpsize;
		allocate();
		for( size_t i = 0; i < this->allocsize; ++i ) {
			ds[ i ] = toCopy.ds[ i ];
			c_ind[ i ]= toCopy.c_ind[ i ];
		}
		for( size_t i = 0; i < this->jumpsize; ++i ) {
			r_ind[ i ] = toCopy.r_ind[ i ];
		}
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
	oICRS( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero = 0 ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		//set superclass fields
		this->zero_element = zero;
		this->nor = m;
		this->noc = n;
	
		//sanity checks
		if( log2(this->noc) > sizeof( _i_value ) * 8 ) {
			std::cerr << "Warning: the matrix with " << this->noc << " columns cannot be represented within " << (sizeof( _i_value )*8) << "-bit index values this oICRS instance uses! Attempting to continue anyway..." << std::endl;
		}
		if( log2(this->nor) > sizeof( _i_value ) * 8 ) {
			std::cerr << "Warning: the matrix with " << this->nor << " rows cannot be represented within " << (sizeof( _i_value )*8) << "-bit index values this oICRS instance uses! Attempting to continue anyway..." << std::endl;
		}
		if( m==0 || n==0 || input.size() == 0 ) { //empty matrix
			this->nor = this->noc = this->nnz = 0;
			ds    = NULL;
			r_ind = NULL;
			c_ind = NULL;
			return;
		}

		//WARNING: noc*nor typically overflows on 32 bits!
		//         immediately start recording differences
		//         instead of trying to directly calculate
		//         the index as i*noc+j.

		//TODO this is what makes this an ICRS implementation; removing the sort will result in BICRS.
		//Complexity dependens on qsort implementation. Probably O(nnz log(nnz)) average, O(nnz^2) worst-case.
		//for standard C++ sort:
		//sort( input.begin(), input.end(), compareTriplets );
		//for C quicksort:
		qsort( &input[ 0 ], input.size(), sizeof( Triplet< T > ), &compareTriplets );

		//filtering out zeros is skipped for now (TODO).
		
		//column indices per row in each block
		std::vector< size_t > indices[ block_length_row ];

		//vector aliases of the internal flat arrays.
		//The flat arrays will be allocated taking alignment into
		//account; this is why no C++ vectors are used internally.
		std::vector< T > ds;
		std::vector< _i_value > c_ind;
		std::vector< _i_value > r_ind;
		//std::vector< _i_value > r_ind2;

		//nonzero iterator
		size_t k = 0;

		//first block row tracker
		bool first_row = false; //set to false, since the very first row of blocks should not trigger an overflow

		//row tracker
		_i_value prev_row = 0;

		//column tracker
		_i_value prev_col = 0;

		//indicates when there are no more nonzeroes to add
		bool depleted;
		
		//start getting block rows
		do {
			//get next block indices
			k = prepareBlockIteration( indices, r_ind, prev_row, input, k, this->nor );
			//get iterators
			size_t cit[ block_length_row ];
			for( size_t i = 0; i < block_length_row; ++i ) {
				cit[ i ] = 0;
			}
			//fill up c_ind and ds; loop until there are no more nonzeroes to add
			do {
				//default value
				depleted = true;
				//indicates when a nonzero has been added to the current block
				bool first = true;
				//remember top-left block index
				const size_t top_left = c_ind.size();
				//for each block row
				for( size_t i = 0; i < block_length_row; ++i ) {
					//add a maximum of block_length_col indices
					for( size_t j = 0; j < block_length_column; ++j ) {
						//check if this row is depleted
						if( cit[ i ] == indices[ i ].size() ) {
							//yes; add dummy zeroes
							c_ind.push_back( 0 ); //no column increment (so not to incur unnecessary cache misses)
							ds.push_back( 0 ); //zero value
							//DBG
							//std::cout << "Adding dummy zeros at block index " << i << ", " << j << std::endl;
						} else {
							//current nonzero index
							const size_t k = indices[ i ][ cit[ i ] ];
							//DBG
							//std::cout << "Column of nonzero " << k << " is " << input[ k ].j() << ", storing increment " << static_cast< _i_value >( input[ k ].j() - prev_col ) << " which is relative to " << prev_col << std::endl;
							//add this nonzero
							c_ind.push_back( input[ k ].j() - prev_col );
							ds   .push_back( input[ k ].value );
							//update prev_col if this was the first nonzero
							if( first ) {
								prev_col = input[ k ].j();
								first = false;
								//if this was not the very first block row, we need to signal an overflow to handle row changes
								if( first_row ) {
									//record overflow
									c_ind[ c_ind.size() - 1 ] += this->noc;
									//subsequent blocks will not trigger a row change
									first_row = false;
								}
								//if the first nonzero was not on the first row, we need to swap this increment with that of block index 0, 0
								if( i != 0 ) {
									std::swap( c_ind[ c_ind.size() - 1 ], c_ind[ top_left ]);
								}
							}
							//increment this row iterator
							++cit[ i ];
						}
					}
				}
				//check if depleted
				for( size_t i = 0; i < block_length_row; ++i ) {
					if( cit[ i ] < indices[ i ].size() ) {
						depleted = false;
						break;
					}
				}
			} while( !depleted );
			//clear all vectors
			for( size_t i = 0; i < block_length_row; ++i ) {
				indices[ i ].clear();
			}
			//the first nonzero that follows will be the first nonzero in this row of blocks
			first_row = true;
		} while( k < input.size() );

		//sanity check
		assert( c_ind.size() == ds.size() );

		//check allocsize is a multiple of block_length
		assert( c_ind.size() % block_length == 0 );

		//postprocess r_ind array
		postProcessRowIncrements( r_ind );

		//coda
		c_ind.push_back( this->noc ); //force overflow
		ds.   push_back( 0 );         //dummy nonzero

		//round up to block size
		for( size_t i = 1; i < block_length; ++i ) {
			c_ind.push_back( 0 ); //dummy increment
			ds   .push_back( 0 ); //dummy nonzero
		}

		//sanity check
		assert( c_ind.size() == ds.size() );

		//check allocsize is a multiple of block_length
		assert( c_ind.size() % block_length == 0 );

		//update allocsize
		allocsize = c_ind.size();

		//set the number of nonzeroes
		this->nnz = input.size();

		//calculate wasted space
		const size_t wasted_space = allocsize - this->nnz;
		fillIn = wasted_space;

		//DBG
		//std::cout << "Info: added " << wasted_space << " explicit nonzeroes (" << (100.0*((double)wasted_space)/((double)(this->nnz))) << "%) to enable full vectorisation.\n";

		//update jump size, add one more extra block that might be read during zax, but not used.
		jumpsize = r_ind.size() + block_length;

		//allocate arrays
		allocate();

		//copy row-jump vector
		for( size_t i = 0; i < r_ind.size(); ++i ) {
			//actual copy
			this->r_ind[ i ] = r_ind[ i ];
		}

		//pad row jump vector with zeroes
		for( size_t i = r_ind.size(); i < r_ind.size() + block_length; ++i ) {
			this->r_ind[ i ] = 0;
		}

		//copy column increments vector and nonzero values vector
		for( size_t i = 0; i < c_ind.size(); ++i ) {
			this->c_ind[ i ] = c_ind[ i ];
			this->ds   [ i ] = ds   [ i ];
		}

		//done
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = static_cast< ULI >( *( this->r_ind ) );
		col = static_cast< ULI >( *( this->c_ind ) );
	}

	/* Gets starting position (first nonzero location)
	void getStartingPos( ULI &row_start, ULI &column_start ) {
		row_start    = static_cast< ULI >( *( this->r_ind ) );
		column_start = static_cast< ULI >( *( this->c_ind ) );
	}*/

	/*
	 *  Sets starting position of matrix multiplication.
	 *  (Useful for example when the input/output vectors will be
	 *   shifted before passed on to this class' zax method.)
	 *
	 *  @param row_start    New row start location
	 *  @param column_start New column start location
	 /
	void setStartingPos( const ULI row_start, const ULI column_start ) {
		*( this->r_ind ) = row_start;
		*( this->c_ind ) = column_start;
	}*/

	/**
	 *  In-place z=xA function. Adapted from the master thesis of Joris Koster, Utrecht University.
	 *
	 *  @param pDataX Pointer to array x to multiply by the current matrix (Ax).
	 *  @param pDataZ Pointer to result array. Must be pre-allocated and its elements set to zero for correct results.
	 */
        virtual void zxa( const T*__restrict__ pDataX, T*__restrict__ pDataZ ) {
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) return;
		         T *__restrict__ pDataA    = ds;
		const    T *__restrict__ pDataAend = ds + this->allocsize - block_length;
		//const    T *__restrict__ const pDataXend = pDataX + this->nor;
		const T * const pDataZend = pDataZ + this->noc;
		//const    T *pDataZend = z + nor; //unused

		_i_value *__restrict__ pIncRow   = r_ind;
		_i_value *__restrict__ pIncCol   = c_ind;

		//define buffers
#ifdef __INTEL_COMPILER
		__declspec(align(64)) _i_value c_ind_buffer[ block_length ];
		__declspec(align(64)) T        input_buffer[ block_length ];
		__declspec(align(64)) T        outputbuffer[ block_length ];
#else
		_i_value c_ind_buffer[ block_length ] __attribute__((aligned));
		T        input_buffer[ block_length ] __attribute__((aligned));
		T        outputbuffer[ block_length ] __attribute__((aligned));
#endif

		//initialise kernel; fill c_ind_buffer
		for( size_t i = 0; i < block_length; ++i )
			c_ind_buffer[ i ] = *pIncCol++;
		//shift input vector
		pDataZ += c_ind_buffer[ 0 ];
		//reset column increment
		c_ind_buffer[ 0 ] = 0;
		//shift output vector
		pDataX += *pIncRow++;
		//start kernel
		while( pDataA < pDataAend ) {
			//while the row is not empty
			while( pDataZ < pDataZend ) {
				//fill input buffer
				for( size_t i = 0; i < block_length; ++i )
					input_buffer[ i ] = pDataX[ c_ind_buffer[ i ] ];
				//do vectorised multiplication
				for( size_t i = 0; i < block_length; ++i )
					outputbuffer[ i ] = pDataA[ i ] * input_buffer[ i ];
				//reduce into actual output vector
				for( size_t i = 0; i < block_length; ++i )
					*pDataZ += outputbuffer[ i ];
				//shift input data
				pDataA += block_length;
				//fill c_ind_buffer
				for( size_t i = 0; i < block_length; ++i )
					c_ind_buffer[ i ] = *pIncCol++;
				//shift input vector
				pDataZ += c_ind_buffer[ 0 ];
				//reset column increment
				c_ind_buffer[ 0 ] = 0;
			}
			//row jump, shift back input vector
			pDataZ -= this->noc;
			//shift output vector
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
		//boundary checks
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) return;
		      T *__restrict__ pDataA    = ds;

		//pointer aliases
		const T * const pDataAend      = ds     + this->allocsize - block_length;
		const T * const pDataXend      = pDataX + this->noc;
#ifndef NDEBUG
		const T * const pDataZst       = pDataZ;
		const T * const pDataXst       = pDataX;
		const T * const pDataZend      = pDataZ + this->nor;
#endif
		_i_value *__restrict__ pIncRow = r_ind;
		_i_value *__restrict__ pIncCol = c_ind;

		//sanity checks on start position
		assert( static_cast< ULI >( *r_ind ) < this->nor );
		assert( static_cast< ULI >( *c_ind ) < this->noc );

		//define buffers
#ifdef __INTEL_COMPILER
		__declspec(align(32)) _i_value c_ind_buffer[ block_length ];
		__declspec(align(32)) _i_value r_ind_buffer[ block_length ];
		__declspec(align(32)) T        input_buffer[ block_length ];
		__declspec(align(32)) T        outputbuffer[ block_length ];
		__declspec(align(32)) T        z_shadow    [ block_length ];
#else
		_i_value c_ind_buffer[ block_length ] __attribute__((aligned));
		_i_value r_ind_buffer[ block_length ] __attribute__((aligned));
		T        input_buffer[ block_length ] __attribute__((aligned));
		T        outputbuffer[ block_length ] __attribute__((aligned));
		T        z_shadow    [ block_length ] __attribute__((aligned));
#endif

		//initialise kernel; fill c_ind_buffer, fill r_ind_buffer, zero output buffer
		for( size_t i = 0; i < block_length; ++i ) {
			c_ind_buffer[ i ] = *pIncCol++;
			r_ind_buffer[ i ] = *pIncRow++;
			outputbuffer[ i ] = 0; //note: not this->zero_element, since those correspond to matrix `nonzeroes' only
		}
		//shift input vector, output vector
		pDataX += c_ind_buffer[ 0 ];
		pDataZ += r_ind_buffer[ 0 ];
		//reset start column index, start row index
		c_ind_buffer[ 0 ] = 0;
		r_ind_buffer[ 0 ] = 0;
		//keep track of z_shadow usage
		size_t blockrow = 0;
		//fill z_shadow
		for( size_t i = 0; i < block_length; ++i ) {
			z_shadow[ i ] = pDataZ[ r_ind_buffer[ i ] ];
		}
		//start kernel
		while( pDataA < pDataAend ) {
			//process row
			while( pDataX < pDataXend ) {
				//sanity checks
				assert( pDataA < pDataAend );
				assert( pDataX >= pDataXst );
				assert( pDataX < pDataXend );
				assert( pDataZ >= pDataZst );
				assert( pDataZ < pDataZend );
				assert( pDataX + this->noc >= pDataXend ); //otherwise pDataX is before the start of x!
				assert( pDataA + this->allocsize >= pDataAend );
				assert( pDataZ + this->nor >= pDataZend );
#ifdef _DEBUG
				std::cout << "Input  vector at " << ((pDataXend-pDataX)<0 ? (pDataX-pDataXend) : (this->noc+pDataX-pDataXend)) << std::endl;
				std::cout << "Output vector at " << (this->nor-(pDataZend-pDataZ)) << std::endl;
				std::cout << "Now at block nr. " << (this->allocsize-block_length-(pDataAend-pDataA))/block_length << std::endl;
#endif
				//DBG
				/*std::cout << "Input buffer { ";
				for( size_t j = 0; pIncCol != (c_ind+block_length) && j < block_length; ++j ) {
					std::cout << input_buffer[ j ] << " ";
				}
				std::cout << "} --> { ";*/

				//fill input buffer, z shadow
				for( size_t i = 0; i < block_length; ++i ) {
					input_buffer[ i ] = pDataX[ c_ind_buffer[ i ] ];
				}

				//DBG
				/*for( size_t j = 0; j < block_length; ++j ) {
					std::cout << input_buffer[ j ] << " ";
				}
				std::cout << "}" << std::endl;*/
				
				//do vectorised multiplication
				for( size_t i = 0; i < block_length; ++i ) {
					outputbuffer[ i ] += pDataA[ i ] * input_buffer[ i ];
				}
				//shift input data
				pDataA += block_length;
				//fill c_ind_buffer
				for( size_t i = 0; i < block_length; ++i ) {
					c_ind_buffer[ i ] = *pIncCol++;
				}
				//shift input vector
				pDataX += c_ind_buffer[ 0 ];
				//reset start column index
				c_ind_buffer[ 0 ] = 0;
			}
#if defined __INTEL_COMPILER && defined __KNC__
			switch( block_length_row ) {
			//check if we are in an all-reduce situation
			case 1:
				//yes, just reduce all elements of outputbuffer into z_shadow
				for( size_t i = 0; i < block_length; ++i ) {
					z_shadow[ blockrow ] += outputbuffer[ i ];
				}
				break;
			//check if we are in a no-reduce situation
			case block_length:
				//yes, just add all elements into z_shadow
				for( size_t i = 0; i < block_length; ++i ) {
					z_shadow[ i ] += outputbuffer[ i ];
				}
				break;
			//default case, partial reductions
			default:
				assert( block_length_row > 1 );
				//for each row in this vectorised block
				for( size_t b = 0; b < block_length_row; ++b ) {
					//get relative output vector index
					const size_t z_shadow_index = b * block_length_column + blockrow;
					//reduce into z_shadow
					for( size_t i = 0; i < block_length_column; ++i ) {
						//get outputbuffer index
						const size_t o_buffer_index = b * block_length_column + i;
						//reduce
						z_shadow[ z_shadow_index ] += outputbuffer[ o_buffer_index ];
					}
				}
			}
#else
			//for each row in this vectorised block 
			for( size_t b = 0; b < block_length_row; ++b ) {
				//reduce block rows into z_shadow
				for( size_t i = 0; i < block_length_column; ++i ) {
					z_shadow[ b + block_length_row * blockrow ] += outputbuffer[ b * block_length_column + i ];
					//DBG
					/*std::cout << "z_shadow at " << (b + block_length_row * blockrow) <<
						" is added with outputbuffer at " << (b * block_length_column + i) <<
						" = " << outputbuffer[ b * block_length_column + i ] << std::endl;
					std::cout << "z_shadow = { ";
					for( size_t j = 0; j < block_length; ++j ) {
						std::cout << z_shadow[ j ] << " ";
					}
					std::cout << "}" << std::endl;*/
				}
			}
#endif
			//we filled another block of rows of z_shadow
			++blockrow;

			//zero out output buffer
			for( size_t i = 0; i < block_length; ++i ) {
				outputbuffer[ i ] = 0;
			}
			//shift back input vector
			pDataX -= this->noc;
			//check if z_shadow is full
			if( blockrow == block_length_column ) {
				//write back z_shadow
				for( size_t i = 0; i < block_length; ++i ) {
					pDataZ[ r_ind_buffer[ i ] ] = z_shadow[ i ];
					//DBG
					/*std::cout << "Write back z_shadow[ " << i << " ] = "
						<< z_shadow[ i ] << " to output vector element at position "
						<< (this->nor-(pDataZend-pDataZ)+r_ind_buffer[i]) << std::endl;*/
				}
				//load new r_ind_buffer
				for( size_t i = 0; i < block_length; ++i ) {
					r_ind_buffer[ i ] = *pIncRow++;
				}
				//shift output vector
				pDataZ += r_ind_buffer[ 0 ];
				//reset start row index
				r_ind_buffer[ 0 ] = 0;
				//DBG
				/*std::cout << "r_ind_buffer = { ";
				for( size_t i = 0; i < block_length; ++i ) {
					std::cout << r_ind_buffer[ i ] << " ";
				}
				std::cout << "}" << std::endl;*/
				//read new z_shadow
				for( size_t i = 0; i < block_length; ++i ) {
					z_shadow[ i ] = pDataZ[ r_ind_buffer[ i ] ];
				}
				//reset block row
				blockrow = 0;
			}
		}
		//coda, write back z_shadow in case there is unwritten data there
		for( size_t i = 0; i < block_length; ++i ) {
			//write back
			pDataZ[ r_ind_buffer[ i ] ] = z_shadow[ i ];
			//DBG
			//std::cout << "Write back z_shadow[ " << i << " ] = "
			//	<< z_shadow[ i ] << " to output vector element at position "
			//	<< (this->nor-(pDataZend-pDataZ)+r_ind_buffer[i]) << " (coda)" << std::endl;
		}
		//DBG
		std::cout << "Goodbye from compiler-optimised " << block_length_row << "-by-" << block_length_column << " SpMV kernel." << std::endl;
        }

	/** Base deconstructor. */
	~oICRS() {
		deallocate();
	}

	virtual size_t bytesUsed() {
		return bytes;
	}

};

//derived from oICRS_MIC_template
#if defined __INTEL_COMPILER && defined __KNC__

//*** New implementation for gather-prefetch: ***//
#define _mm512_prefetch_i32gather_pd_0hint_init() \
register __mmask8 fullMask asm("sil") = 0xff;

#define _mm512_prefetch_i32gather_pd_0hint( index, mv ) \
__asm__ (                                            \
                       "vgatherpf0hintdpd (%1, %0, 8) {{%2}}\n"   \
                       : /* nothing */                            \
                       : "x" (index), "g" (mv), "r" (fullMask)    \
                       : /* nothing */                            \
                       )
//*** end implementation for gather-prefetch ***//

//*** Old implementation for gather-prefetch: ***//
//__mmask8 fullMask = 0xff;

//#define _mm512_prefetch_i32gather_pd_0hint( index, mv ) \
__asm__ (                                            \
                       "vgatherpf0hintdpd (%1, %0, 8) {{%2}}\n"   \
                       : /* nothing */                            \
                       : "x" (index), "g" (mv), "" (fullMask)     \
                       : /* nothing */                            \
                       )
//*** end implementation for gather-prefetch ***//

	/** union m512 data type for random readouts of vector registers (doubles). */
	union __m512d_union {
		__m512d internal;
		double  array[ 8 ];
	};

	/** union m512 data type for random readouts of vector registers (32-bit signed integers). */
	union __m512i_union {
		__m512i internal;
		int32_t array[ 16 ];
	};

	//16-byte unsigned integer indices specialisation for double data types, 8x1 blocking:
	void oICRS< double, 8, 1, int16_t >::allocate() {
		ds    = (double*)   _mm_malloc( this->allocsize * sizeof( double ), _SL_ALIGN_DOUBLE );
		c_ind = (int16_t*) _mm_malloc( this->allocsize * sizeof( int16_t ), _SL_ALIGN_INT16 );
		r_ind = (int16_t*) _mm_malloc( (this->jumpsize + 1) * sizeof( int16_t ), _SL_ALIGN_INT16 );
		bytes = sizeof( int16_t ) * (this->jumpsize + this->allocsize + 1) + sizeof( double ) * this->allocsize;
	}

	//16-byte unsigned integer indices specialisation for double data types, 2x4 blocking:
	void oICRS< double, 2, 4, int16_t >::allocate() {
		ds    = (double*)   _mm_malloc( this->allocsize * sizeof( double ), _SL_ALIGN_DOUBLE );
		c_ind = (int16_t*) _mm_malloc( this->allocsize * sizeof( int16_t ), _SL_ALIGN_INT16 );
		r_ind = (int16_t*) _mm_malloc( (this->jumpsize + 1) * sizeof( int16_t ), _SL_ALIGN_INT16 );
		bytes = sizeof( int16_t ) * (this->jumpsize + this->allocsize + 1) + sizeof( double ) * this->allocsize;
	}

	//16-byte unsigned integer indices specialisation for double data types, 4x2 blocking:
	void oICRS< double, 4, 2, int16_t >::allocate() {
		ds    = (double*)   _mm_malloc( this->allocsize * sizeof( double ), _SL_ALIGN_DOUBLE );
		c_ind = (int16_t*) _mm_malloc( this->allocsize * sizeof( int16_t ), _SL_ALIGN_INT16 );
		r_ind = (int16_t*) _mm_malloc( (this->jumpsize + 1) * sizeof( int16_t ), _SL_ALIGN_INT16 );
		bytes = sizeof( int16_t ) * (this->jumpsize + this->allocsize + 1) + sizeof( double ) * this->allocsize;
	}

	//16-byte unsigned integer indices specialisation for double data types:
	void oICRS< double, 1, 8, int16_t >::allocate() {
		ds    = (double*)   _mm_malloc( this->allocsize * sizeof( double ), _SL_ALIGN_DOUBLE );
		c_ind = (int16_t*) _mm_malloc( this->allocsize * sizeof( int16_t ), _SL_ALIGN_INT16 );
		r_ind = (int16_t*) _mm_malloc( (this->jumpsize + 1) * sizeof( int16_t ), _SL_ALIGN_INT16 );
		bytes = sizeof( int16_t ) * (this->jumpsize + this->allocsize + 1) + sizeof( double ) * this->allocsize;
	}

	//32-byte unsigned integer indices specialisation for double data types:
	void oICRS< double, 1, 8, int32_t >::allocate() {
		ds    = (double*)   _mm_malloc(  this->allocsize     * sizeof( double ),   _SL_ALIGN_DOUBLE );
		c_ind = (int32_t*)  _mm_malloc(  this->allocsize     * sizeof( int32_t ), _SL_ALIGN_INT32 );
		r_ind = (int32_t*)  _mm_malloc( (this->jumpsize + 1) * sizeof( int32_t ), _SL_ALIGN_INT32 );
		bytes = sizeof( int32_t ) * (this->jumpsize + this->allocsize + 1) + sizeof( double ) * this->allocsize;
	}

	//start 16-byte unsigned integer indices specialisation for double data types, 1x8 blocks.
        void oICRS< double, 1, 8, int16_t >::zax( const double *__restrict__ pDataX, double *__restrict__ pDataZ ) {

		//boundary checks
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) {
			return;
		}

		//initialise prefetch macro
		_mm512_prefetch_i32gather_pd_0hint_init();

		//pointer aliases
		const double * const pDataAend = ds     + this->allocsize - block_length;
		const double * const pDataXend = pDataX + this->noc;
		double  *__restrict__ pDataA   = ds;
		int16_t *__restrict__ pIncRow  = r_ind;
		int16_t *__restrict__ pIncCol  = c_ind;

		//define buffers
		__m512i c_ind_buffer;
		__m512d input_buffer;
		__m512d value_buffer;
		__m512d outputbuffer;
		__m512i zeroF = _mm512_set_epi32(
			1, 1, 1, 1, 1, 1, 1, 0,
			1, 1, 1, 1, 1, 1, 1, 0
		);

		//initialise kernel
		outputbuffer = _mm512_setzero_pd();

		//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
		c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );
		
		//shift input vector
		pDataX += *pIncCol;

		//move pIncCol up one block
		pIncCol += block_length;

		//reset start column index; alternatively, do this on 16-bit data earlier
		c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

		//prefetch input
		_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

		//shift output vector
		pDataZ += *pIncRow++;

		//start kernel
		while( pDataA < pDataAend ) {

			//process row
			while( pDataX < pDataXend ) {

				//manual unroll, stage 1: fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//fill nonzero buffer
				value_buffer = _mm512_load_pd( pDataA );

				//do vectorised multiply-add
				outputbuffer = _mm512_fmadd_pd( value_buffer, input_buffer, outputbuffer );

				//shift input data
				pDataA += block_length;

				//shift c_ind_buffer
				c_ind_buffer = _mm512_permute4f128_epi32( c_ind_buffer, _MM_PERM_DCDC );

				//shift input vector
				pDataX += *pIncCol;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

				//move pIncCol up
				pIncCol += block_length;
				
				//check for break
				if( pDataX >= pDataXend ) {

					//reduce row contribution
					*pDataZ +=_mm512_reduce_add_pd( outputbuffer );

					//reset sums buffer
					outputbuffer = _mm512_setzero_pd();

					//row jump, shift back input vector
					pDataX -= this->noc;

					//shift output vector
					pDataZ += *pIncRow++;

					//check for end of matrix
					if( pDataA >= pDataAend ) {
						break;
					}
				}

				//manual unroll, stage 2: fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//fill value buffer
				value_buffer = _mm512_load_pd( pDataA );

				//do vectorised multiply-add
				outputbuffer = _mm512_fmadd_pd( value_buffer, input_buffer, outputbuffer );

				//shift input data
				pDataA += block_length;

				//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
				c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_UINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

				//shift input vector
				pDataX += *pIncCol;

				//move pIncCol up
				pIncCol += block_length;

				//reset start column index; alternatively, do this on 16-bit data earlier
				c_ind_buffer =_mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

			}

			//reduce row contribution
			*pDataZ +=_mm512_reduce_add_pd( outputbuffer );

			//reset sums buffer
			outputbuffer = _mm512_setzero_pd();

			//row jump, shift back input vector
			pDataX -= this->noc;

			//shift output vector
			pDataZ += *pIncRow++;
		}
	}
	//end 16-byte unsigned integer indices specialisation for double data types, 1x8 blocks.
	
	//start 16-byte unsigned integer indices specialisation for double data types, 2x4 blocks.
        void oICRS< double, 2, 4, int16_t >::zax( const double *__restrict__ pDataX, double *__restrict__ pDataZ ) {
		//boundary checks
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) {
			return;
		}

		//initialise prefetch macro
		_mm512_prefetch_i32gather_pd_0hint_init();

		//DBG
		//const double * const pDataXst       = pDataX;
		//const double * const pDataZst       = pDataZ;

		//pointer aliases
		const double * const pDataAend = ds     + this->allocsize - block_length;
		const double * const pDataXend = pDataX + this->noc;
		double  *__restrict__ pDataA   = ds;
		int16_t *__restrict__ pIncRow  = r_ind;
		int16_t *__restrict__ pIncCol  = c_ind;

		//DBG
		//std::cout << "x start: " << pDataX << ", x end: " << pDataXend << std::endl;

		//define buffers
		__m512i c_ind_buffer;
		__m512d input_buffer;
		__m512d value_buffer;
		__m512d_union outputbuffer;
		__m512i zeroF = _mm512_set_epi32(
			1, 1, 1, 1, 1, 1, 1, 0,
			1, 1, 1, 1, 1, 1, 1, 0
		);

		//initialise kernel
		outputbuffer.internal = _mm512_setzero_pd();

		//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
		c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

		//DBG
		/*int16_t *dbgd = (int16_t*)pIncCol;
		std::cout << "c_ind_buffer= " << pDataX << ": "
			 << dbgd[ 0 ] << " "
			 << dbgd[ 1 ] << " "
			 << dbgd[ 2 ] << " "
			 << dbgd[ 3 ] << " "
			 << dbgd[ 4 ] << " "
			 << dbgd[ 5 ] << " "
			 << dbgd[ 6 ] << " "
			 << dbgd[ 7 ] << " "
			<< std::endl;*/
	
		//shift input vector
		pDataX += *pIncCol;

		//move pIncCol up one block
		pIncCol += block_length;

		//reset start column index; alternatively, do this on 16-bit data earlier
		c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

		//prefetch input
		_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

		//shift output vector
		pDataZ += *pIncRow++;

		//DBG
		//std::cout << "Entering loop..." << std::endl;

		//start kernel
		while( pDataA < pDataAend ) {

			//process row
			while( pDataX < pDataXend ) {

				//DBG
				/*std::cout << "In loop, stage 1, A-E: " << std::endl;
				__declspec(align(64)) int32_t dbg[ 16 ];
				_mm512_store_epi32( dbg, c_ind_buffer );
				std::cout << "Grabbing from " << pDataX << " ( " << (pDataX-pDataXst) << " ): "
						 << dbg[ 0 ] << " "
						 << dbg[ 1 ] << " "
						 << dbg[ 2 ] << " "
						 << dbg[ 3 ] << " "
						 << dbg[ 4 ] << " "
						 << dbg[ 5 ] << " "
						 << dbg[ 6 ] << " "
						 << dbg[ 7 ] << " "
						<< std::endl;*/

				//manual unroll, stage 1: fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//DBG
				//std::cout << "A" << std::endl;

				//fill nonzero buffer
				value_buffer = _mm512_load_pd( pDataA );

				//DBG
				//std::cout << "B" << std::endl;

				//do vectorised multiply-add
				outputbuffer.internal = _mm512_fmadd_pd( value_buffer, input_buffer, outputbuffer.internal );

				//DBG
				//std::cout << "C" << std::endl;

				//shift input data
				pDataA += block_length;

				//DBG
				//std::cout << "D" << std::endl;

				//shift c_ind_buffer
				c_ind_buffer = _mm512_permute4f128_epi32( c_ind_buffer, _MM_PERM_DCDC );

				//DBG
				//std::cout << "E" << std::endl;

				//shift input vector
				pDataX += *pIncCol;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

				//move pIncCol up
				pIncCol += block_length;
				
				//DBG
				//std::cout << "Unroll I done" << std::endl;

				//check for break
				if( pDataX >= pDataXend ) {

					//reduce half of row contribution
					const __m512d reduce_buffer1 = _mm512_swizzle_pd( outputbuffer.internal, _MM_SWIZ_REG_CDAB );
					outputbuffer.internal = _mm512_add_pd( outputbuffer.internal, reduce_buffer1 );
					//(output buffer now contains HG GH FE EF DC CD BA AB)
					const __m512d reduce_buffer2 = _mm512_swizzle_pd( outputbuffer.internal, _MM_SWIZ_REG_BADC );
					outputbuffer.internal = _mm512_add_pd( outputbuffer.internal, reduce_buffer2 );
					//(output buffer now contains HGFE GHEF FEHG EFGH DCBA CDAB BADC ABCD)

					//DBG
					//std::cout << "Writing out " << outputbuffer.array[ 0 ] << " to index " << (static_cast< size_t >(pDataZ - pDataZst)) << std::endl;

					//reduce part 1
					*pDataZ += outputbuffer.array[ 0 ];

					//row jump, shift back input vector
					pDataX -= this->noc;

					//shift output vector
					pDataZ += *pIncRow++;

					//DBG
					//std::cout << "Writing out " << outputbuffer.array[ 4 ] << " to index " << (static_cast< size_t >(pDataZ - pDataZst)) << std::endl;

					//reduce part 2
					*pDataZ += outputbuffer.array[ 4 ];

					//reset sums buffer
					outputbuffer.internal = _mm512_setzero_pd();

					//shift output vector
					pDataZ += *pIncRow++;

					//check for end of matrix
					if( pDataA >= pDataAend ) {
						break;
					}
				}

				//DBG
				/*std::cout << "In loop, stage 2, A-E: " << std::endl;
				_mm512_store_epi32( dbg, c_ind_buffer );
				std::cout << "Grabbing from " << pDataX << " ( " << (pDataX-pDataXst) << " ): "
						 << dbg[ 0 ] << " "
						 << dbg[ 1 ] << " "
						 << dbg[ 2 ] << " "
						 << dbg[ 3 ] << " "
						 << dbg[ 4 ] << " "
						 << dbg[ 5 ] << " "
						 << dbg[ 6 ] << " "
						 << dbg[ 7 ] << " "
						<< std::endl;*/

				//manual unroll, stage 2: fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//DBG
				//std::cout << "A" << std::endl;

				//fill value buffer
				value_buffer = _mm512_load_pd( pDataA );

				//DBG
				//std::cout << "B" << std::endl;

				//do vectorised multiply-add
				outputbuffer.internal = _mm512_fmadd_pd( value_buffer, input_buffer, outputbuffer.internal );

				//DBG
				//std::cout << "C" << std::endl;

				//shift input data
				pDataA += block_length;

				//DBG
				//std::cout << "D" << std::endl;

				//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
				c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

				//DBG
				//std::cout << "E" << std::endl;

				//shift input vector
				pDataX += *pIncCol;

				//move pIncCol up
				pIncCol += block_length;

				//DBG
				//std::cout << "F" << std::endl;

				//reset start column index; alternatively, do this on 16-bit data earlier
				c_ind_buffer =_mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

				//DBG
				//std::cout << "Unroll II done" << std::endl;

			}
					
			//reduce half of row contribution
			const __m512d reduce_buffer1 = _mm512_swizzle_pd( outputbuffer.internal, _MM_SWIZ_REG_CDAB );
			outputbuffer.internal = _mm512_add_pd( outputbuffer.internal, reduce_buffer1 );
			//(output buffer now contains HG GH FE EF DC CD BA AB)
			const __m512d reduce_buffer2 = _mm512_swizzle_pd( outputbuffer.internal, _MM_SWIZ_REG_BADC );
			outputbuffer.internal = _mm512_add_pd( outputbuffer.internal, reduce_buffer2 );
			//(output buffer now contains HGFE GHEF FEHG EFGH DCBA CDAB BADC ABCD)

			//DBG
			//std::cout << "Writing out " << outputbuffer.array[ 0 ] << " to index " << (static_cast< size_t >(pDataZ - pDataZst)) << std::endl;

			//reduce part 1
			*pDataZ += outputbuffer.array[ 0 ];

			//row jump, shift back input vector
			pDataX -= this->noc;

			//shift output vector
			pDataZ += *pIncRow++;

			//DBG
			//std::cout << "Writing out " << outputbuffer.array[ 4 ] << " to index " << (static_cast< size_t >(pDataZ - pDataZst)) << std::endl;

			//reduce part 2
			*pDataZ += outputbuffer.array[ 4 ];

			//reset sums buffer
			outputbuffer.internal = _mm512_setzero_pd();

			//shift output vector
			pDataZ += *pIncRow++;
		}
	}
	//end 16-byte unsigned integer indices specialisation for double data types, 2x4 blocks.
	
	//start 16-byte unsigned integer indices specialisation for double data types, 4x2 blocks.
	//NOTE: this code uses prefetching in both the input and output gathers. Enabling or disabling
	//      gather-prefetch on the output vector seems to have no discernable effect on performance.
	//      Presumably ICC filters out those prefetches based on internal heuristics. This is fine;
	//      the code is retained here for completeness.
        void oICRS< double, 4, 2, int16_t >::zax( const double *__restrict__ pDataX, double *__restrict__ pDataZ ) {
		//boundary checks
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) {
			return;
		}

		//initialise prefetch macro
		_mm512_prefetch_i32gather_pd_0hint_init();

		//DBG
		//const double * const pDataXst       = pDataX;
		//const double * const pDataZst       = pDataZ;

		//pointer aliases
		const double * const pDataAend = ds     + this->allocsize - block_length;
		const double * const pDataXend = pDataX + this->noc;
		double  *__restrict__ pDataA   = ds;
		int16_t *__restrict__ pIncRow  = r_ind;
		int16_t *__restrict__ pIncCol  = c_ind;

		//DBG
		//std::cout << "x start: " << pDataX << ", x end: " << pDataXend << std::endl;

		//define buffers
		__m512d z_shadow;
		__m512i c_ind_buffer;
		__m512i r_ind_buffer;
		__m512d input_buffer;
		__m512d value_buffer;
		__m512d_union outputbuffer;
		__m512i zeroF = _mm512_set_epi32(
			1, 1, 1, 1, 1, 1, 1, 0,
			1, 1, 1, 1, 1, 1, 1, 0
		);
		__mmask8 mask1 = 0x55;
		__mmask8 mask2 = 0xAA;

		//initialise kernel
		outputbuffer.internal = _mm512_setzero_pd();

		//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
		c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

		//DBG
		/*int16_t *dbgd = (int16_t*)pIncCol;
		std::cout << "c_ind_buffer= " << pDataX << ": "
			 << dbgd[ 0 ] << " "
			 << dbgd[ 1 ] << " "
			 << dbgd[ 2 ] << " "
			 << dbgd[ 3 ] << " "
			 << dbgd[ 4 ] << " "
			 << dbgd[ 5 ] << " "
			 << dbgd[ 6 ] << " "
			 << dbgd[ 7 ] << " "
			<< std::endl;*/
	
		//shift input vector
		pDataX += *pIncCol;

		//move pIncCol up one block
		pIncCol += block_length;

		//reset start column index; alternatively, do this on 16-bit data earlier
		c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

		//fill r_ind_buffer
		r_ind_buffer = _mm512_extload_epi32( pIncRow, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

		//nullify r_ind_buffer
		r_ind_buffer = _mm512_mullo_epi32( r_ind_buffer, zeroF );

		//prefetch index buffers
		_mm512_prefetch_i32gather_pd_0hint( r_ind_buffer, pDataZ );
		_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

		//shift output vector
		pDataZ += *pIncRow;

		//move pIncRow up one block
		pIncRow += block_length;

		//fill z shadow, this loads 8 output vector elements
		z_shadow = _mm512_i32loextgather_pd( r_ind_buffer, pDataZ, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

		//keep track of row/column switches (load vs register shuffle), plus mask shift for z_shadow
		bool colshift, rowshift, maskshift;
		colshift = rowshift = maskshift = true;

		//DBG
		//std::cout << "Entering loop..." << std::endl;

		//start kernel
		while( pDataA < pDataAend ) {

			//process row
			while( pDataX < pDataXend ) {

				//DBG
				/*std::cout << "In loop, stage 1, A-E: " << std::endl;
				__declspec(align(64)) int32_t dbg[ 16 ];
				_mm512_store_epi32( dbg, c_ind_buffer );
				std::cout << "Grabbing from " << pDataX << " ( " << (pDataX-pDataXst) << " ): "
						 << dbg[ 0 ] << " "
						 << dbg[ 1 ] << " "
						 << dbg[ 2 ] << " "
						 << dbg[ 3 ] << " "
						 << dbg[ 4 ] << " "
						 << dbg[ 5 ] << " "
						 << dbg[ 6 ] << " "
						 << dbg[ 7 ] << " "
						<< std::endl;*/

				//manual unroll, stage 1: fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//DBG
				//std::cout << "A" << std::endl;

				//fill nonzero buffer
				value_buffer = _mm512_load_pd( pDataA );

				//DBG
				//std::cout << "B" << std::endl;

				//do vectorised multiply-add
				outputbuffer.internal = _mm512_fmadd_pd( value_buffer, input_buffer, outputbuffer.internal );

				//DBG
				//std::cout << "C" << std::endl;

				//shift input data
				pDataA += block_length;

				//DBG
				//std::cout << "D" << std::endl;

				if( colshift ) {
					//DBG
					//std::cout << "D colshift A" << std::endl;

					//shuffle c_ind_buffer
					c_ind_buffer = _mm512_permute4f128_epi32( c_ind_buffer, _MM_PERM_DCDC );

					//DBG
					//std::cout << "D colshift B" << std::endl;
				} else {
					//DBG
					//std::cout << "D !colshift A" << std::endl;

					//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
					c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

					//DBG
					//std::cout << "D !colshift B" << std::endl;

					//reset main increment to 0
					c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

					//DBG
					//std::cout << "D !colshift C" << std::endl;
				}
				colshift = !colshift;

				//DBG
				//std::cout << "E" << std::endl;

				//shift input vector
				pDataX += *pIncCol;

				//DBG
				//std::cout << "F" << std::endl;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

				//move pIncCol up
				pIncCol += block_length;
				
				//DBG
				//std::cout << "end unroll I" << std::endl;

				//check for break
				if( pDataX >= pDataXend ) {
					//DBG
					//std::cout << "Row break in-between unroll I and unroll II" << std::endl;

					//reduce half of row contribution
					const __m512d reduce_buffer1 = _mm512_swizzle_pd( outputbuffer.internal, _MM_SWIZ_REG_CDAB );
					outputbuffer.internal = _mm512_add_pd( outputbuffer.internal, reduce_buffer1 );
					//(output buffer now contains HG GH FE EF DC CD BA AB)

					//DBG
					/*std::cout << "Writing out "
						<< outputbuffer.array[ 0 ] << ", "
						<< outputbuffer.array[ 2 ] << ", "
						<< outputbuffer.array[ 4 ] << ", "
						<< outputbuffer.array[ 6 ] << " to shadow of index "
						<< static_cast< size_t >(pDataZ - pDataZst)
						<< std::endl;*/

					//add to z_shadow
					if( maskshift ) {
						//DBG
						//std::cout << "Updating z_shadow..." << std::endl;

						z_shadow = _mm512_mask_add_pd( z_shadow, mask1, z_shadow, outputbuffer.internal );
					} else {
						//DBG
						//std::cout << "Updating z_shadow using mask2" << std::endl;

						z_shadow = _mm512_mask_add_pd( z_shadow, mask2, z_shadow, outputbuffer.internal );

						//DBG
						/*std::cout << "Writing out z_shadow { ";
						__declspec(align(64)) double dbgd[ 8 ];
						_mm512_store_pd( dbgd, z_shadow );
						std::cout << dbgd[ 0 ] << " "
								 << dbgd[ 1 ] << " "
								 << dbgd[ 2 ] << " "
								 << dbgd[ 3 ] << " "
								 << dbgd[ 4 ] << " "
								 << dbgd[ 5 ] << " "
								 << dbgd[ 6 ] << " "
								 << dbgd[ 7 ] << " "
								<< "}";
						std::cout << " to positions " << (pDataZ - pDataZst) << " + { ";
						__declspec(align(64)) int32_t dbg[ 16 ];
						_mm512_store_epi32( dbg, r_ind_buffer );
						std::cout << dbg[ 0 ] << " "
								 << dbg[ 1 ] << " "
								 << dbg[ 2 ] << " "
								 << dbg[ 3 ] << " "
								 << dbg[ 4 ] << " "
								 << dbg[ 5 ] << " "
								 << dbg[ 6 ] << " "
								 << dbg[ 7 ] << " "
								<< "}" << std::endl;*/

						//we are done with this z_shadow; first write out
						_mm512_i32loextscatter_pd( pDataZ, r_ind_buffer, z_shadow, _MM_DOWNCONV_PD_NONE, 8, _MM_HINT_NONE );

						//shift output vector
						pDataZ += *pIncRow;
				
						//DBG
						//std::cout << "Getting new r_ind_buffer" << std::endl;

						//get next r_ind_buffer
						if( rowshift ) {
							//shift r_ind_buffer
							r_ind_buffer = _mm512_permute4f128_epi32( r_ind_buffer, _MM_PERM_DCDC );
						} else {
							//fill r_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
							r_ind_buffer = _mm512_extload_epi32( pIncRow, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );
	
							//reset main increment to 0
							r_ind_buffer = _mm512_mullo_epi32( r_ind_buffer, zeroF ); //i.e., r_ind_buffer[ 0 ] = 0;
						}
						rowshift = !rowshift;

						//prefetch next z_shadow
						_mm512_prefetch_i32gather_pd_0hint( r_ind_buffer, pDataZ );

						//shift pIncRow
						pIncRow += block_length;

						//DBG
						//std::cout << "Getting new z_shadow..." << std::endl;

						//load next z_shadow (this loads 8 output vector elements)
						z_shadow = _mm512_i32loextgather_pd( r_ind_buffer, pDataZ, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );
					}
					maskshift = !maskshift;

					//row jump, shift back input vector
					pDataX -= this->noc;

					//reset sums buffer
					outputbuffer.internal = _mm512_setzero_pd();

					//check for end of matrix
					if( pDataA >= pDataAend ) {
						//write out z_shadow; there won't be a regular write-out after unroll II
						break;
					}

					//DBG
					//std::cout << "Exit row switch" << std::endl;
				}

				//DBG
				/*std::cout << "In loop, unroll II: " << std::endl;
				_mm512_store_epi32( dbg, c_ind_buffer );
				std::cout << "Grabbing from " << pDataX << " ( " << (pDataX-pDataXst) << " ): "
						 << dbg[ 0 ] << " "
						 << dbg[ 1 ] << " "
						 << dbg[ 2 ] << " "
						 << dbg[ 3 ] << " "
						 << dbg[ 4 ] << " "
						 << dbg[ 5 ] << " "
						 << dbg[ 6 ] << " "
						 << dbg[ 7 ] << " "
						<< std::endl;*/

				//manual unroll, stage 2: fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//DBG
				//std::cout << "A" << std::endl;

				//fill value buffer
				value_buffer = _mm512_load_pd( pDataA );

				//DBG
				//std::cout << "B" << std::endl;

				//do vectorised multiply-add
				outputbuffer.internal = _mm512_fmadd_pd( value_buffer, input_buffer, outputbuffer.internal );

				//DBG
				//std::cout << "C" << std::endl;

				//shift input data
				pDataA += block_length;

				//DBG
				//std::cout << "D" << std::endl;
				
				if( colshift ) {
					//DBG
					//std::cout << "D colshift A" << std::endl;

					//shuffle c_ind_buffer
					c_ind_buffer = _mm512_permute4f128_epi32( c_ind_buffer, _MM_PERM_DCDC );

					//DBG
					//std::cout << "D colshift B" << std::endl;
				} else {
					//DBG
					//std::cout << "D !colshift A" << std::endl;

					//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
					c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

					//DBG
					//std::cout << "D !colshift B" << std::endl;

					//reset main increment to 0
					c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

					//DBG
					//std::cout << "D !colshift C" << std::endl;
				}
				colshift = !colshift;

				//DBG
				//std::cout << "E" << std::endl;

				//shift input vector
				pDataX += *pIncCol;

				//DBG
				//std::cout << "F" << std::endl;

				//move pIncCol up
				pIncCol += block_length;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

				//DBG
				//std::cout << "Unroll II done" << std::endl;
			}
					
			//reduce half of row contribution
			const __m512d reduce_buffer1 = _mm512_swizzle_pd( outputbuffer.internal, _MM_SWIZ_REG_CDAB );
			outputbuffer.internal = _mm512_add_pd( outputbuffer.internal, reduce_buffer1 );
			//(output buffer now contains HG GH FE EF DC CD BA AB)
			
			//DBG
			/*std::cout << "Writing out "
				<< outputbuffer.array[ 0 ] << ", "
				<< outputbuffer.array[ 2 ] << ", "
				<< outputbuffer.array[ 4 ] << ", "
				<< outputbuffer.array[ 6 ] << " to shadow of index "
				<< (static_cast< size_t >(pDataZ - pDataZst))
				<< std::endl;*/

			//add to z_shadow
			if( maskshift ) {
				//DBG
				//std::cout << "Updating z_shadow..." << std::endl;

				z_shadow = _mm512_mask_add_pd( z_shadow, mask1, z_shadow, outputbuffer.internal );
			} else {
				//DBG
				//std::cout << "Updating z_shadow using mask2" << std::endl;

				z_shadow = _mm512_mask_add_pd( z_shadow, mask2, z_shadow, outputbuffer.internal );

				//DBG
				/*std::cout << "Writing out z_shadow { ";
				__declspec(align(64)) double dbgd[ 8 ];
				_mm512_store_pd( dbgd, z_shadow );
				std::cout << dbgd[ 0 ] << " "
						 << dbgd[ 1 ] << " "
						 << dbgd[ 2 ] << " "
						 << dbgd[ 3 ] << " "
						 << dbgd[ 4 ] << " "
						 << dbgd[ 5 ] << " "
						 << dbgd[ 6 ] << " "
						 << dbgd[ 7 ] << " "
						<< "}";
				std::cout << " to positions " << (pDataZ - pDataZst) << " + { ";
				__declspec(align(64)) int32_t dbg[ 16 ];
				_mm512_store_epi32( dbg, r_ind_buffer );
				std::cout << dbg[ 0 ] << " "
						 << dbg[ 1 ] << " "
						 << dbg[ 2 ] << " "
						 << dbg[ 3 ] << " "
						 << dbg[ 4 ] << " "
						 << dbg[ 5 ] << " "
						 << dbg[ 6 ] << " "
						 << dbg[ 7 ] << " "
						<< "}" << std::endl;*/

				//we are done with this z_shadow; first write out
				_mm512_i32loextscatter_pd( pDataZ, r_ind_buffer, z_shadow, _MM_DOWNCONV_PD_NONE, 8, _MM_HINT_NONE );

				//shift output vector
				pDataZ += *pIncRow;
				
				//DBG
				//std::cout << "Getting new r_ind_buffer" << std::endl;

				//get next r_ind_buffer
				if( rowshift ) {
					//shift r_ind_buffer
					r_ind_buffer = _mm512_permute4f128_epi32( r_ind_buffer, _MM_PERM_DCDC );
				} else {
					//fill r_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
					r_ind_buffer = _mm512_extload_epi32( pIncRow, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );
	
					//reset main increment to 0
					r_ind_buffer = _mm512_mullo_epi32( r_ind_buffer, zeroF ); //i.e., r_ind_buffer[ 0 ] = 0;
				}
				rowshift = !rowshift;

				//prefetch next z_shadow
				_mm512_prefetch_i32gather_pd_0hint( r_ind_buffer, pDataZ );

				//shift pIncRow
				pIncRow += block_length;

				//DBG
				//std::cout << "Getting new z_shadow..." << std::endl;

				//load next z_shadow (this loads 8 output vector elements)
				z_shadow = _mm512_i32loextgather_pd( r_ind_buffer, pDataZ, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );
			}
			maskshift = !maskshift;

			//reset sums buffer
			outputbuffer.internal = _mm512_setzero_pd();

			//row jump, shift back input vector
			pDataX -= this->noc;

			//DBG
			//std::cout << "End outer loop" << std::endl;
		}
		if( !maskshift ) {
			//there is a write-out of z_shadow still pending; do it
			_mm512_i32loextscatter_pd( pDataZ, r_ind_buffer, z_shadow, _MM_DOWNCONV_PD_NONE, 8, _MM_HINT_NONE );
		}
	}
	//end 16-byte unsigned integer indices specialisation for double data types, 4x2 blocks.
	
	//start 16-byte unsigned integer indices specialisation for double data types, 8x1 blocks.
	void oICRS< double, 8, 1, int16_t >::zax( const double *__restrict__ pDataX, double *__restrict__ pDataZ ) {
		//boundary checks
		if( this->nor == 0 || this->noc == 0 || this->nnz == 0 ) {
			return;
		}

		//initialise prefetch macro
		_mm512_prefetch_i32gather_pd_0hint_init();

		//DBG
		//const double * const pDataXst = pDataX;
		//const double * const pDataZst = pDataZ;

		//pointer aliases
		double *__restrict__ pDataA    = ds;
		const double * const pDataAend = ds     + this->allocsize - block_length;
		const double * const pDataXend = pDataX + this->noc;
		int16_t *__restrict__ pIncRow  = r_ind;
		int16_t *__restrict__ pIncCol  = c_ind;

		//define buffers
		__m512d z_shadow;
		__m512i c_ind_buffer;
		__m512i r_ind_buffer;
		__m512d input_buffer;
		__m512d value_buffer;
		__m512i zeroF = _mm512_set_epi32(
			1, 1, 1, 1, 1, 1, 1, 0,
			1, 1, 1, 1, 1, 1, 1, 0
		);

		//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
		c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

		//DBG
		/*__declspec(align(64)) int32_t dbg[ 16 ];
		_mm512_store_epi32( dbg, c_ind_buffer );
		std::cout << "c_ind_buffer equals " << dbg[ 0 ] << ", " << dbg[ 1 ] << ", " <<
			dbg[ 2 ] << ", " << dbg[ 3 ] << ", " << dbg[ 4 ] << ", " <<
			dbg[ 5 ] << ", " << dbg[ 6 ] << ", " << dbg[ 7 ] << "." << std::endl;*/

		//shift input vector
		pDataX += *pIncCol;

		//move pIncCol up one block
		pIncCol += block_length;

		//reset start column index; alternatively, do this on 16-bit data earlier
		c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;

		//fill r_ind_buffer
		r_ind_buffer = _mm512_extload_epi32( pIncRow, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

		//nullify r_ind_buffer
		r_ind_buffer = _mm512_mullo_epi32( r_ind_buffer, zeroF );

		//prefetch index buffers
		_mm512_prefetch_i32gather_pd_0hint( r_ind_buffer, pDataZ );
		_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

		//shift output vector
		pDataZ += *pIncRow;

		//move pIncRow up one block
		pIncRow += block_length;

		//fill z shadow, this loads 8 output vector elements
		z_shadow = _mm512_i32loextgather_pd( r_ind_buffer, pDataZ, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

		//keep track of row/column switches (load vs register shuffle), plus mask shift for z_shadow
		bool colshift, rowshift;
		colshift = rowshift = true;

		//DBG
		//std::cout << "Entering loop..." << std::endl;

		//start kernel
		while( pDataA < pDataAend ) {
			//for each block row
			while( pDataX < pDataXend ) {
				//DBG
				/*std::cout << "Looping with x @ " << (pDataX-pDataXst) << " and z @ " << (pDataZ-pDataZst) << std::endl;
				__declspec(align(64)) int32_t dbg[ 16 ];
				__declspec(align(64)) double dbgd[ 16 ];
				_mm512_store_epi32( dbg, c_ind_buffer );
				std::cout << "c_ind_buffer equals " << dbg[ 0 ] << ", " << dbg[ 1 ] << ", " <<
					dbg[ 2 ] << ", " << dbg[ 3 ] << ", " << dbg[ 4 ] << ", " <<
					dbg[ 5 ] << ", " << dbg[ 6 ] << ", " << dbg[ 7 ] << "." << std::endl;
				_mm512_store_epi32( dbg, r_ind_buffer );
				std::cout << "r_ind_buffer equals " << dbg[ 0 ] << ", " << dbg[ 1 ] << ", " <<
					dbg[ 2 ] << ", " << dbg[ 3 ] << ", " << dbg[ 4 ] << ", " <<
					dbg[ 5 ] << ", " << dbg[ 6 ] << ", " << dbg[ 7 ] << "." << std::endl;
				_mm512_store_pd( dbgd, z_shadow );
				std::cout << "z_shadow equals " << dbgd[ 0 ] << " " << dbgd[ 1 ] << " " 
						<< dbgd[ 2 ] << " " << dbgd[ 3 ] << " " << dbgd[ 4 ]<< " " 
						<< dbgd[ 5 ] << " " << dbgd[ 6 ] << " " << dbgd[ 7 ]
						<< std::endl;*/

				//fill input buffers
				input_buffer = _mm512_i32loextgather_pd( c_ind_buffer, pDataX, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

				//fill value buffer
				value_buffer = _mm512_load_pd( pDataA );

				//do vectorised multiply-add
				z_shadow = _mm512_fmadd_pd( value_buffer, input_buffer, z_shadow );

				//shift input data
				pDataA += block_length;

				//DBG
				//std::cout << "preshift" << std::endl;

				//check how to move to the next block
				if( colshift ) {
					//shift c_ind_buffer
					c_ind_buffer = _mm512_permute4f128_epi32( c_ind_buffer, _MM_PERM_DCDC );
				} else {
					//fill c_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
					c_ind_buffer = _mm512_extload_epi32( pIncCol, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

					//reset main increment to 0
					c_ind_buffer = _mm512_mullo_epi32( c_ind_buffer, zeroF ); //i.e., c_ind_buffer[ 0 ] = 0;
				}
				colshift = !colshift;

				//shift input vector
				pDataX += *pIncCol;

				//prefetch input
				_mm512_prefetch_i32gather_pd_0hint( c_ind_buffer, pDataX );

				//move pIncCol up
				pIncCol += block_length;
				
				//DBG
				//std::cout << "/inner loop" << std::endl;
			}

			//row jump, z_shadow is fully updated, write out z_shadow and load next. Scatter:
			_mm512_i32loextscatter_pd( pDataZ, r_ind_buffer, z_shadow, _MM_DOWNCONV_PD_NONE, 8, _MM_HINT_NONE );

			//get next r_ind_buffer
			if( rowshift ) {
				//shift r_ind_buffer
				r_ind_buffer = _mm512_permute4f128_epi32( r_ind_buffer, _MM_PERM_DCDC );
			} else {
				//fill r_ind_buffer NOTE: one load reads 16 increments; i.e., 2 blocks!
				r_ind_buffer = _mm512_extload_epi32( pIncRow, _MM_UPCONV_EPI32_SINT16, _MM_BROADCAST32_NONE, _MM_HINT_NT );

				//reset main increment to 0
				r_ind_buffer = _mm512_mullo_epi32( r_ind_buffer, zeroF ); //i.e., r_ind_buffer[ 0 ] = 0;
			}
			rowshift = !rowshift;

			//shift output vector
			pDataZ += *pIncRow;

			//prefetch next z_shadow
			_mm512_prefetch_i32gather_pd_0hint( r_ind_buffer, pDataZ );

			//shift pIncRow
			pIncRow += block_length;

			//load next z_shadow (this loads 8 output vector elements)
			z_shadow = _mm512_i32loextgather_pd( r_ind_buffer, pDataZ, _MM_UPCONV_PD_NONE, 8, _MM_HINT_NONE );

			//row jump, shift back input vector
			pDataX -= this->noc;
		}
	}
	//end 16-byte unsigned integer indices specialisation for double data types, 8x1 blocks.

#endif

#endif

