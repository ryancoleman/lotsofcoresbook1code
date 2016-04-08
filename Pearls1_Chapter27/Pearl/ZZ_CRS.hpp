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


#include "SparseMatrix.hpp"
#include "Triplet.hpp"
#include <assert.h>
#include <vector>
#include <algorithm>
#include<iostream>

//#define _DEBUG

#ifndef _H_ZZ_CRS
#define _H_ZZ_CRS

/**
 *  The zig-zag compressed row storage sparse matrix data structure.
 */
template< typename T, typename _i_value=LI >
class ZZ_CRS: public SparseMatrix< T, LI > {

  private:

  protected:

	/** Array containing the actual this->nnz non-zeros. */
	T* ds;

	/** Array containing the column jumps. */
	LI* col_ind;

	/** Array containing the row jumps. */
	LI* row_start;

	/** Comparison function used for sorting input data. */
	static bool compareTripletsLTR( const Triplet< T >* one, const Triplet< T >* two ) {
		if( one->i() < two->i() )
			return true;
		if ( one->i() > two->i() )
			return false;
		return one->j() < two->j();
	}

	/** Comparison function used for sorting input data. */
	static bool compareTripletsRTL( const Triplet< T >* one, const Triplet< T >* two ) {
		if( one->i() < two->i() )
			return true;
		if ( one->i() > two->i() )
			return false;
		return one->j() > two->j();
	}

  public:

	/** Base constructor. */
	ZZ_CRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	ZZ_CRS( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid CRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns to be stored.
	 *  @param zero Which element is considered to be the zero element.
	 */
	ZZ_CRS( const long int number_of_nonzeros, const long int number_of_rows, const long int number_of_cols, T zero ):
		SparseMatrix< T, unsigned long int >( number_of_nonzeros, number_of_rows, number_of_cols, zero ) {
		ds = new T[ this->nnz ];
		col_ind = new long int[ this->nnz ];
		row_start = new long int[ this->nor + 1 ];
	}

	/**
	 *  Copy constructor.
	 *  @param toCopy reference to the CRS datastructure to copy.
	 */
	ZZ_CRS( ZZ_CRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->noc = toCopy.noc;
		this->nor = toCopy.nor;
		ds = new T[ this->nnz ];
		col_ind = new long int[ this->nnz ];
		row_start = new long int[ this->nor ];
		for( long int i=0; i<this->nnz; i = i + 1 ) {
			ds[ i ] = toCopy.ds[ i ];
			col_ind[ i ]= toCopy.col_ind[ i ];
		}
		for( long int i=0; i<=this->nor; i++ )
			row_start[ i ] = toCopy.row_start[ i ];
	}

	/**
	 *  Constructor which transforms a collection of input triplets to CRS format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *  @param input The input collection.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero Which element is considered zero.
	 */
	ZZ_CRS( std::vector< Triplet< T > >& input, const LI m, const LI n, const T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const LI m, const LI n, const T zero ) {
		this->zero_element = zero;
		this->nnz = input.size();
		this->nor = m;
		this->noc = n;

		std::vector< std::vector< Triplet< T >* > > ds( this->nor );
	
		//move input there
		typename std::vector< Triplet< T > >::iterator in_it = input.begin();
		for( ; in_it != input.end(); in_it++ ) {
			Triplet< T >* cur = &(*in_it);
			const long int currow = cur->i();
			const T value = cur->value;
			if( value == this->zero_element ) { this->nnz--; continue; }
			ds.at( currow ).push_back( cur );
		}

		//allocate arrays
		this->ds = new T[ this->nnz ];
		col_ind = new LI[ this->nnz ];
		row_start = new LI[ this->nor + 1 ];
		
		//make ZZ-CRS
		long int index = 0;
		bool LTR       = true;
		for( long int currow = 0; currow < this->nor; currow++ ) {
			row_start[ currow ] = index;
			if( ds.at( currow ).size() == 0 ) continue;
			if( LTR )
				sort( ds.at( currow ).begin(), ds.at( currow ).end(), compareTripletsLTR );
			else
				sort( ds.at( currow ).begin(), ds.at( currow ).end(), compareTripletsRTL );
			LTR = !LTR;
			typename std::vector< Triplet< T >* >::iterator row_it = ds.at( currow ).begin();
			for( ; row_it!=ds.at( currow ).end(); row_it++ ) {
				const Triplet< T > cur = *(*row_it);
				this->ds[ index ] = cur.value;
				col_ind[ index ] = cur.j();
				index++;
			}
		}
		assert( index == this->nnz );
		row_start[ this->nor ] = this->nnz;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( LI &row, LI &col ) {
		row = static_cast< unsigned long int >( this->row_start[ 0 ] );
		col = static_cast< unsigned long int >( this->col_ind[ 0 ] );
	}

	/**
	 *  In-place z=xA multiplication algorithm.
	 *  
	 *  @param x The vector x supplied for left-multiplication.
	 *  @param z The pre-allocated result vector. All elements should be set to zero in advance.
	 */
	virtual void zxa( const T* x, T* z ) {
		T*__restrict__ ds_p = ds;
		LI*__restrict__ col_ind_p = col_ind;
		LI row, index;
		for( row = 0; row < this->nor; row++, x++ ) {
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ )
				z[ *col_ind_p++ ] += *ds_p++ * *x;
		}
	}

	/**
	 *  In-place z=Ax multiplication algorithm.
	 *  
	 *  @param x The vector x supplied for multiplication.
	 *  @param z The pre-allocated result vector. All elements should be set to zero in advance.
	 */
	virtual void zax( const T*__restrict__ x, T*__restrict z ) {
		T*__restrict__ ds_p = ds;
		_i_value*__restrict__ col_ind_p = col_ind;
		long int index=0,row=0;
		for( ; row < this->nor; row++, z++ )
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ )
				*z += *ds_p++ * x[ *col_ind_p++ ];
	}

	virtual size_t bytesUsed() {
		return sizeof( LI ) * ( this->nnz + this->nor + 1 ) + sizeof( T ) * this->nnz;
	}

	/** Base deconstructor. */
	~ZZ_CRS() {
		delete [] ds;
		delete [] col_ind;
		delete [] row_start;
	}

};

#endif

