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


#include <vector>
#include <assert.h>
#include "Triplet.hpp"
#include "SparseMatrix.hpp"

#ifndef _H_TS
#define _H_TS

/** The triplet scheme; a storage scheme for sparse matrices using triplets. */
template< typename T >
class TS: public SparseMatrix< T, ULI > {

   private:
	
   protected:

	/** The row indices of the nonzeros. */
	ULI* i;
	
	/** The column indices of the nonzeros. */
	ULI* j;

	/** The values of the nonzeros. */
	T* ds;

   public:

	/** Base constructor. */
	TS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	TS( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );
	}
	
	/** Base constructor.
	 *  @param input std::vector of triplets to be stored in this scheme.
	 *  @param m total number of rows.
	 *  @param n total number of columns.
	 *  @param zero what is to be considered the zero element.
	 */
	TS( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero = 0 ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		ULI offset = 0;

		this->zero_element = zero;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();
		ds = new T[ this->nnz ];
		i = new ULI[ this->nnz ];
		j = new ULI[ this->nnz ];
		for( ULI r=0; r<this->nnz; r++ )
			if( input[ r ].value != this->zero_element ) {
				//ds[ r ] = input[ r ];
				ds[ r - offset ] = input[ r ].value;
				i[ r - offset ] = input[ r ].i();
				j[ r - offset ] = input[ r ].j();
			} else {
				offset++;
			}
		this->nnz -= offset;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = this->i[ 0 ];
		col = this->j[ 0 ];
	}

	/**
	 *  In-place z=xA calculation algorithm.
	 *
	 *  @param x The vector x supplied for calculation of xA.
	 *  @param z The result vector z. Should be pre-allocated with entries set to zero.
	 */
        virtual void zxa( const T* x, T* z ) {
                for( ULI r=0; r<this->nnz; r++ ) {
                        assert( i[ r ] >= 0  );
                        assert( i[ r ] < this->nor );
                        assert( j[ r ] >= 0  );
                        assert( j[ r ] < this->noc );
                        z[ j[ r ] ] += ds[ r ] * x[ i[ r ] ];
                }
        }

	/**
	 *  In-place z=Ax calculation algorithm.
	 *
	 *  @param x The vector x supplied for calculation of Ax.
	 *  @param z The result vector z. Should be pre-allocated with entries set to zero.
	 */
        virtual void zax( const T* x, T* z ) {
                for( ULI r=0; r<this->nnz; r++ ) {
                        assert( i[ r ] >= 0  );
                        assert( i[ r ] < this->nor );
                        assert( j[ r ] >= 0  );
                        assert( j[ r ] < this->noc );
                        z[ i[ r ] ] += ds[ r ] * x[ j[ r ] ];
                }
        }

	virtual size_t bytesUsed() {
		return sizeof( ULI ) * 2 * this->nnz + sizeof( T ) * this->nnz;
	}

	/** Base destructor. */
	~TS() {
		delete [] i;
		delete [] j;
		delete [] ds;
	}

};

#endif

