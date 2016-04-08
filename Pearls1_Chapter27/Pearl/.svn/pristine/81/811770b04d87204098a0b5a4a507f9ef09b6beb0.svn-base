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


#include "CRS.hpp"

extern "C" {
 #include <mkl_service.h>
 #include <mkl_spblas.h>
}

//#define _DEBUG

#ifndef _H_MKLCRS
#define _H_MKLCRS

/**
 *  The compressed row storage sparse matrix data structure.
 */
template< typename T >
class MKLCRS: public CRS< T > {

  private:

  protected:

	/** Required for call to MKL; a factor alpha=beta of 1. */
	double _one;

	/** Required for call to MKL; (no) transposition. */
	char trans;

	/** Required for call to MKL; matrix descriptor. */
	char *descr;

	/** Required for call to MKL; matrix row-size. */
	int _m;

	/** Required for call to MKL; matrix column-size. */
	int _n;

	/** Required for call to MKL; a plain-int version of col_ind. */
	int *_col_ind;

	/** Required for call to MKL; a plain-int version of row_start. */
	int *_row_start;

	/**
	 *  Does the required post-processing of Sparse Library's
	 *  CRS representation to one compatible with Intel MKL.
	 */
	void prepare() {
		_m = static_cast< int >( this->nor );
		_n = static_cast< int >( this->noc );
		_one = 1.0;

		_col_ind = new int[ this->nnz ];
		for( ULI i=0; i<this->nnz; ++i )
			_col_ind[ i ] = static_cast< int >( this->col_ind[ i ] );

		_row_start = new int[ _m + 1 ];
		for( int i=0; i<=_m; ++i )
			_row_start[ i ] = static_cast< int >( this->row_start[ i ] );

		trans = 'n';

		descr = new char[4];
		descr[ 0 ] = 'G';
		descr[ 1 ] = descr[ 2 ] = '*';
		descr[ 3 ] = 'C';
	}

  public:

	/** Base constructor. */
	MKLCRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	MKLCRS( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );

		prepare();
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid CRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns of the matrix.
	 *  @param zero The element considered to be zero.
	 */
	MKLCRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ) {
		this->nnz = number_of_nonzeros;
		this->nor = number_of_rows;
		this->noc = number_of_cols;
		this->zero_element = zero;
		this->row_start = new ULI[ this->nor + 1 ];
		this->ds = new T[ this->nnz ];
		this->col_ind = new ULI[ this->nnz ];

		prepare();
	}

	/** Copy constructor.
	 *  @param toCopy reference to the CRS datastructure to copy.
	 */
	MKLCRS( CRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->nor = toCopy.nor;
		this->row_start = new ULI[ this->nor + 1 ];
		this->ds = new T[ this->nnz ];
		this->col_ind = new ULI[ this->nnz ];
		for( ULI i=0; i<this->nnz; i++ ) {
			this->ds[ i ] = toCopy.ds[ i ];
			this->col_ind[ i ] = toCopy.col_ind[ i ];
		}
		for( ULI i=0; i<this->nor; i++ )
			this->row_start[ i ] = toCopy.row_start[ i ];

		prepare();
	}
	
	/**
	 *  Constructor which transforms a collection of input triplets to CRS format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *  @param input The input collection.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero The element considered to be zero.
	 */
	MKLCRS( std::vector< Triplet< T > > input, ULI m, ULI n, T zero ) {
		load( input, m, n, zero );

		prepare();
	}

#ifndef _NO_LIBNUMA
	/** Overloaded mv call; allocates output vector using numa_interleaved. */
	virtual T* mv( const T* x ) {
		T* ret = (T*) numa_alloc_interleaved( this->nor * sizeof( T ) );
		for( ULI i=0; i<this->nor; i++ )
			ret[ i ] = this->zero_element;
		this->zax( x, ret );
		return ret;
	}
#endif

	/** 
	 *  In-place z=xA function.
	 */
        virtual void zxa( const T*__restrict__ x, T*__restrict__ z ) {
		std::cerr << "MKLCRS does not implement the zxa, sorry!" << std::endl;
		exit( EXIT_FAILURE );
        }

	/** 
	 *  In-place z=Ax function.
	 *  
	 *  @param x The x vector to multiply current matrix with.
	 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 */
        virtual void zax( const T*__restrict__ x, T*__restrict__ z ) {
		//Warning: assumes T is a double
		mkl_dcsrmv( &trans, &_m, &_n, &_one, descr, this->ds, _col_ind, _row_start, &( _row_start[ 1 ] ), (T*__restrict__)x, &_one, z );
        }

	/** Base deconstructor. */
	virtual ~MKLCRS() {
		delete [] _col_ind;
		delete [] _row_start;
		delete [] descr;
	}

};

#endif

