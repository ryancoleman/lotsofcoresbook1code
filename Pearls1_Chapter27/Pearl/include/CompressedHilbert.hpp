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
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2012.
 */


#ifndef _H_COMPRESSED_HILBERT
#define _H_COMPRESSED_HILBERT

#include <vector>

#include "Triplet.hpp"
#include "SparseMatrix.hpp"
#include "custom_quicksort.hpp"
#include "HilbertArray.hpp"

/** 
 * A Hilbert scheme backed by a specialised storage format.
 * This implementation does not use blocking, as that would require a
 * two-fold datastructure (COO or BICRS on the lower end). The block
 * size would then best be of the size of a char (256x256) or of a 
 * short int (65536x65536).
 */
template< typename T >
class CompressedHilbert: public SparseMatrix< T, ULI > {

    private:

    protected:

	/** Array of nonzero values */
	T* values;

	/** Pointer to indices in Hilbert-coordinate form */
	HilbertArrayInterface< T >* indices;

    public:

	/** Base deconstructor. */
	virtual ~CompressedHilbert() {
		delete [] values;
		delete indices;
	}

	/** Base constructor. */
	CompressedHilbert() {
		values  = NULL;
		indices = NULL;
	}

	/** Base constructor initialising from file. */
	CompressedHilbert( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );
	}

	/** Base constructor initialising from direct input. */
	CompressedHilbert( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero = 0 ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		this->zero_element = 0;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();

		//calculate hilbert coordinates
		std::vector< BigInt > hilbert_values;
		typename std::vector< Triplet< T > >::const_iterator it;
		unsigned long int h1, h2;
		for( it = input.begin(); it != input.end(); ++it ) {
			Matrix2HilbertCoordinates::IntegerToHilbert( it->i(), it->j(), h1, h2 );
			hilbert_values.push_back( BigInt( h1, h2 ) );
		}

		//sort triplets and hilbert coordinates according to hilbert coordinates
		custom_quicksort< BigInt, Triplet< T > >( hilbert_values, input, 0, input.size() );

		//check everything indeed is in ascending order
		std::vector< BigInt >::iterator b_it1 = hilbert_values.begin();
		std::vector< BigInt >::iterator b_it2 = hilbert_values.begin() + 1;
		for( ; b_it2 != hilbert_values.end(); ++b_it1, ++b_it2 ) {
			assert( *b_it1 < *b_it2 );
		}

		//now use the HilbertArray class to get a compressed version of hilbert_values
		indices = getDiffArray< T >( hilbert_values, m, n );

		//store nonzero values
		values   = new T[ this->nnz ];
		size_t c = 0;
		for( it  = input.begin(); it != input.end(); ++it )
			values[ c++ ] = it->value;

		//initialisation done
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = static_cast< ULI >( indices->getFirstRowIndex() );
		col = static_cast< ULI >( indices->getFirstColumnIndex() );
	}

	/** Calculates z=xA.
	 *  Note z is not implicitly zeroed before multiplication.
	 *
	 *  @param x The (initialised) input  vector.
	 *  @param z The (initialised) output vector.
	 */
	virtual void zxa( const T* x, T* z ) {
		//handle boundary case
		if( indices == NULL )
			return;
		//get pointers
		const T * const val_end = values + this->nnz;
		const T *__restrict__ p_v = values;
		//do zxa
		indices->zxa( x, z, p_v, val_end );
	}

	/** Calculates z=Ax.
	 *  Note z is not implicitly zeroed before multiplication.
	 *
	 *  @param x The (initialised) input  vector.
	 *  @param z The (initialised) output vector.
	 */
	virtual void zax( const T* x, T* z ) {
		//handle boundary case
		if( indices == NULL )
			return;
		//get pointers
		const T * const val_end = values + this->nnz;
		const T *__restrict__ p_v = values;
		//do zax
		indices->zax( x, z, p_v, val_end );

		/* Non-flat implementation:
		//get restricted pointers as well as the end position
		const T * const val_end = values + this->nnz;
		const T *__restrict__ const p_x = x;
		const T *__restrict__ p_v = values;
		T *__restrict__ const p_z = z;
		//initialise iterations through Hilbert indices
		indices->moveToStart( p_x, p_z, *p_v++ );
		//do iterations
		while( p_v < val_end )
			indices->moveToNext( p_x, p_z, *p_v++ );
		*/
	}

	/** @return The number of bytes used by this storage scheme. */
	virtual size_t bytesUsed() {
		if( indices == NULL )
			return 0;
		else
			return indices->bytesUsed();
	}

};

#endif

