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


#include <assert.h>
#include <vector>
#include "Triplet.hpp"
#include "SparseMatrix.hpp"

/** The dense diagonal matrix scheme; a storage scheme for sparse matrices consisting of only dense diagonals.
 *  @param T type of numerical values to store in the matrix (int, unsigned int, float, double, ...)
 *  @param number_of_diagonals number of (assumed) dense diagonals in the matrix
 *  @param diagonal_offsets must be defined global in calling code,
 *                          so that it is ensured constant at compiletime (necessary when template is used)
 */
template< typename T, int number_of_diagonals, int diagonal_offsets[] >
class DD_MATRIX: public SparseMatrix< T, unsigned long int > {

   private:
	
	/** Convenience typedef */
	typedef unsigned long int ULI;

   protected:

	/** The values of the nonzeros. */
	T** nzs;

	/** How many full length diagonals this (possible not square matrix) contains. */
	ULI full;

	/** What the main diagonal length is (longest diagonal). */
	ULI d;

	/** Whether or not nzs was allocated by this instance itself. */
	bool SELF_ALLOCATED;

   public:

	/** Base constructor. */
	DD_MATRIX() {
		load( std::vector< Triplet< T > >(), 0, 0, 0 );
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	DD_MATRIX( std::string file, T zero = 0 ) {
		loadFromFile( file, zero );
	}
	
	/** Base constructor.
	 *  @param input std::vector of triplets to be stored in this scheme.
	 *  @param m total number of rows.
	 *  @param n total number of columns.
	 *  @param zero what is to be considered the zero element.
	 */
	DD_MATRIX( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
		load( input, m, n, zero );
	}

	/** Dense diagonal matrix specific constructor.
	 *  @param nonzeroes The to this matrix belonging k times max(m,n)
	 *                   dense matrix representation.
	 *  @param m Number of matrix rows.
	 *  @param n Number of matrix columns.
	 *  @param zero What is considered to be the zero element.
	 *  @see DD_MATRIX( input, m, n, zero ) for specifics on the rest of the parameters.
	 */
	DD_MATRIX( T **nonzeroes, ULI m, ULI n, T zero ) {
		this->zero_element = zero;
		this->nor = m;
		this->noc = n;
		this->nnz = 0;

		for( ULI k=0; k<number_of_diagonals; k++ )
			if( diagonal_offsets[ k ] < 0 ) {
				this->nnz += m+diagonal_offsets[ k ] < n ? m+diagonal_offsets[ k ] : n;
			} else {
				this->nnz += n-diagonal_offsets[ k ] < m ? m-diagonal_offsets[ k ] : m;
			}

		d = m > n ? m : n;
		full = m > n ? m-n+1 : n-m+1;
		nzs = nonzeroes;

		SELF_ALLOCATED = false;
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {

		this->zero_element = zero;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();

		nzs = new T*[ number_of_diagonals ];
		d = m > n ? m : n;
		full = m > n ? m-n+1 : n-m+1;

		//first look which diagonals are populated
		ULI *diags = new ULI[m+n-1];
		for( ULI r=0; r<m+n-1; r++ )
			diags[ r ] = 0;
		for( ULI r=0; r<this->nnz; r++ )
			if( input[ r ].value != this->zero_element )
				diags[input[ r ].j()-input[ r ].i()+m-1] = 1;
		ULI nl1 = 0; //number of diagonals populated
		for( ULI r=0; r<m+n-1; r++ )
			if( diags[ r ] == 1 )
				nl1++;
		
		//assume perfect human input
		ULI c = 0;
		assert( number_of_diagonals == nl1 );
		for( ULI r=0; r<m+n-1; r++ )
			if( diags[ r ] == 1 )
				assert( static_cast< int >( r ) - static_cast< int >( m ) + 1 == diagonal_offsets[ c++ ] );

		//allocate nzs fully
		for( ULI r=0; r<number_of_diagonals; r++ ) {
			ULI curdiag = diagonal_offsets[ r ];
			nzs[r] = new T[ curdiag <= full ? d : d-(curdiag-full) ];
		}
	
		//build diagonal code to nzs-row array
		for( ULI r=0; r<m+n-1; r++ )
			diags[ r ] = static_cast< ULI >( -1 ); //signal overflow (hopefully)
		for( ULI r=0; r<number_of_diagonals; r++ )
			diags[ diagonal_offsets[ r ] + m - 1 ] = r;
	
		for( ULI r=0; r<this->nnz; r++ )
			if( input[ r ].value != this->zero_element ) {
				ULI cur = input[ r ].j() - input[ r ].i() + m-1; //diagonal code
				if( input[ r ].i() > input[ r ].j() )
					nzs[ diags[ cur ] ][ input[ r ].j() ] = input[ r ].value;
				else
					nzs[ diags[ cur ] ][ input[ r ].i() ] = input[ r ].value;
			}

		//set self allocation flag
		SELF_ALLOCATED = true;

		//done 
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( unsigned long int &row, unsigned long int &col ) {
		row = 0;
		col = 0;
	}

	/**
	 *  In-place z=xA calculation algorithm.
	 *
	 *  @param x The vector x supplied for calculation of xA.
	 *  @param z The result vector z. Should be pre-allocated with entries set to zero.
	 */
	virtual void zxa( T* x, T* z ) {
		for( ULI j=0; j<d; j++ ) {
			//theoretically, a compiler could unroll this inner loop perfectly
			//and optimise out the if-statements.
			for( ULI k=0; k<number_of_diagonals; k++ ) {
				const ULI i = j + diagonal_offsets[ k ];
				if( i >= this->nor ) continue;
				if( diagonal_offsets[ k ] < 0 )
					z[ j ] += nzs[ k ][ j ] * x[ i ];
				else
					z[ j ] += nzs[ k ][ i ] * x[ i ];
			}
		}
	}

	/**
	 *  In-place z=Ax calculation algorithm.
	 *
	 *  @param x The vector x supplied for calculation of Ax.
	 *  @param z The result vector z. Should be pre-allocated with entries set to zero.
	 */
        virtual void zax( T* x, T* z ) {
		//std::cout << "d= " << d << std::endl;
		for( ULI i=0; i<d; i++ ) {
			//theoretically, a compiler could unroll this inner loop perfectly
			//and optimise out the if-statements.
			for( ULI k=0; k<number_of_diagonals; k++ ) {
				const ULI j = i + diagonal_offsets[ k ];
				//std::cout << "j>=noc: " << j << " >= " << this->noc << std::endl;
				if( j >= this->noc ) continue;
				//std::cout << "i,j,k = " << i << "," << j << "," << k << std::endl;
				if( diagonal_offsets[ k ] < 0 )
					z[ i ] += nzs[ k ][ j ] * x[ j ];
				else
					z[ i ] += nzs[ k ][ i ] * x[ j ];
			}
		}
        }

	/** Base destructor. */
	~DD_MATRIX() {
		if( !SELF_ALLOCATED ) return;

		for( ULI k = 0; k<number_of_diagonals; k++ )
			delete [] nzs[ k ];
		delete [] nzs;
	}

};

