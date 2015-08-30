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


#ifndef _H_SM
#define _H_SM

#include <xmmintrin.h>

#include "Matrix.hpp"
#include "FileToVT.hpp"
#include "ULIdef.hpp"

/**
 *  Interface common to all sparse matrix storage schemes.
 */
template< typename T, typename IND >
class SparseMatrix: public Matrix< T > {
	private:
	
	protected:
		/** Number of rows. */
		IND nor;

		/** Number of columns */
		IND noc;

		/** Number of non-zeros. */
		IND nnz;

	public:

		/** The element considered to be zero. */
		T zero_element;

		/** Base constructor. */
		SparseMatrix() {}

		/** Base constructor. */
		SparseMatrix( const IND nzs, const IND nr, const IND nc, const T zero ):
			nnz( nzs ), nor( nr ), noc( nc ), zero_element( zero ) {}

		/** Base deconstructor. */
		virtual ~SparseMatrix() {}

		/** 
		 *  Function reading in from std::vector< Triplet< T > > format.
		 *  @param input The input matrix in triplet format.
		 *  @param m The number of rows of the input matrix.
		 *  @param n The number of columns of the input matrix.
		 *  @param zero Which element is to be considered zero.
		 */
		virtual void load( std::vector< Triplet< T > >& input, const IND m, const IND n, const T zero ) = 0;

		/**
		 *  Function which loads a matrix from a matrix market file.
		 *  @param file Filename (including path) to the matrix market file.
		 *  @param zero Which element is to be considered zero.
		 */
		void loadFromFile( const std::string file, const T zero=0 ) {
			ULI m, n;
			const size_t pos = file.find_last_of( '.' );
			const std::string ext = file.substr( pos + 1, file.length() );
			std::vector< Triplet< T > > VT;
			if( ext.compare( "trp" ) == 0 ) {
				VT = Triplet< T >::load( file, m, n );
			} else if( ext.compare( "crs" ) == 0 || ext.compare( "csr" ) == 0 ) {
				VT = Triplet< T >::loadCRS( file, m, n );
			} else //assume default matrix-market format
				VT = FileToVT::parse( file, m, n );
			this->load( VT, m, n, zero );
		}

		/**
                 * Queries the number of rows this matrix contains.
                 * @return The number of rows.
                 */
		virtual unsigned long int m() {
			return static_cast< unsigned long int >( nor );
		}

		/**
                 * Queries the number of columns this matrix contains.
                 * @return The number of columns.
                 */
		virtual unsigned long int n() {
			return static_cast< unsigned long int >( noc );
		}

		/**
                 * Queries the number of nonzeroes stored in this matrix.
                 * @return The number of nonzeros.
                 */
		virtual unsigned long int nzs() {
			return static_cast< unsigned long int >( nnz );
		}

		/** Returns the first nonzero index, per reference. */
		virtual void getFirstIndexPair( IND &row, IND &col ) = 0;

		/**
		 *  Calculates and returns z=Ax. The vectors x should have appropiately set length; this is not
		 *  enforced. Memory leaks, segmentation faults, the works; one or more of these will occur if dimensions
		 *  do not make sense.
		 *  The return array is allocated to length equal to the number of rows in this function.
		 *
		 *  *Warning* The output vector is aligned to 64-byte boundaries. Free the assigned memory using _mm_free!
		 *
		 *  @param x The input (dense) x-vector.
		 *  @return The matrix-vector multiplication Ax, where A corresponds to the currently stored matrix.
		 *  @see mv
		 */
		virtual T* mv( const T* x ) {
			T* ret = (T*) _mm_malloc( nor * sizeof( T ), 64 ); //instead of `new T[ nor ];'
			for( IND i=0; i<nor; i++ )
				ret[ i ] = zero_element;
			zax( x, ret );
			return ret;
		}

		/**
		 *  In-place z=Ax function.
		 *
		 *  @param x The x vector to multiply current matrix with.
		 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
		 */
		virtual void zax( const T*__restrict__ x, T*__restrict__ z ) = 0;

		/**
		 *  In-place z=xA function.
		 *
		 *  @param x The x vector to apply in left-multiplication to A
		 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
		 */
		virtual void zxa( const T*__restrict__x, T*__restrict__ z ) = 0;

};

#endif //_H_SM

