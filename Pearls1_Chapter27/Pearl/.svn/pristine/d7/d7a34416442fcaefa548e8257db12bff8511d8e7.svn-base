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
 *     A. N. Yzelman, Dept. of Mathematic, Utrecht University, 2008.
 */


#include "SparseMatrix.hpp"
#include <assert.h>
#include <iostream>

//Set when in debug mode
//#define _DEBUG

#ifndef _H_BICRS
#define _H_BICRS

/**
 *  Bi-directional Incremental Compressed Row Storage scheme.
 *  Supports jumping back and forward within columns.
 *  Supports jumping back and forward between rows.
 *  Main storage direction in column-wise.
 *  Storage requirements are 2nz plus the number of row jumps required.
 *  Many row jumps are disadvantageous to storage as well as speed.
 *  @param _t_value The type of the nonzeros in the matrix.
 *
 *  Warning: this class uses assertions! For optimal performance,
 *           define the NDEBUG flag (e.g., pass -DNDEBUG as a compiler
 *	     flag).
 */
template< typename _t_value, typename _i_value=LI >
class BICRS: public SparseMatrix< _t_value, ULI > {

     protected:

	/** Stores the row start position. */
	ULI r_start;

	/** Stores the column start position. */
	ULI c_start;

	/** Stores the row end position. */
	ULI r_end;

	/** Stores the column end position. */
	ULI c_end;

	/** Stores the number of row jumps. */
	ULI jumps;

	/** Stores the row jumps; size is _at maximum_ the number of nonzeros. */
	_i_value* r_inc;

	/** Stores the column jumps; size is exactly the number of nonzeros. */
	_i_value* c_inc;

	/** Stores the values of the individual nonzeros. Size is exactly the number of nonzeros. */
	_t_value* vals;

	/** Caches n times two. */
	_i_value ntt;
	
     public:

	/** Base deconstructor. */
	virtual ~BICRS() {
		delete [] r_inc;
		delete [] c_inc;
		delete [] vals;
	}

	/** Base constructor. */
	BICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	BICRS( std::string file, _t_value zero = 0 ) {
		this->loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor. Stores triplets in exactly the same order as passed
	 *  to this constructor.
	 *  @param row The row numbers of the individual nonzeros.
	 *  @param col The column numbers of the individual nonzeros.
	 *  @param val The values of the nonzeros.
	 *  @param m Number of matrix rows.
	 *  @param n Number of matrix columns.
	 *  @param nz Number of nonzeros.
	 *  @param zero Which value is to be regarded zero here.
	 */
	BICRS( _i_value* row, _i_value* col, _t_value* val, ULI m, ULI n, ULI nz, _t_value zero ) {
		this->load( row, col, val, m, n, nz, zero );
	}

	/** Base constructor.
	 *  @see SparseMatrix::SparseMatrix( input, m, n, zero )
	 */
	BICRS( std::vector< Triplet< _t_value > >& input, ULI m, ULI n, _t_value zero = 0 ) {
		this->load( input, m, n, zero );
	}

	/**
	 *  This function will rewrite the std::vector< Triplet > structure to one suitable
	 *  for the other load function.
	 *  @see load( row, col, val, m, n, nz ) 
	 *  @see SparseMatrix::load 
	 */
	virtual void load( std::vector< Triplet< _t_value > >& input, ULI m, ULI n, _t_value zero ) {
		ULI nz = input.size();
		_i_value* row = new _i_value[ nz ];
		_i_value* col = new _i_value[ nz ];
		_t_value* val = new _t_value[ nz ];
		unsigned long int c = 0;
		typename std::vector< Triplet< _t_value > >::iterator it = input.begin();
		for( ; it!=input.end(); it++, c++ ) {
			row[ c ] = (*it).i();
			col[ c ] = (*it).j();
			val[ c ] = (*it).value;
		}
		load( row, col, val, m, n, nz, zero );
		assert( vals != val );
		delete [] row;
		delete [] col;
		delete [] val;
	}

	/** @see BICRS( row, col, val, m, n, nz ) */
	void load( _i_value* row, _i_value* col, _t_value* val, ULI m, ULI n, ULI nz, _t_value zero ) {
#ifdef _DEBUG
		std::cerr << "Warning: _DEBUG flag set." << std::endl;
#endif
		this->zero_element = zero;
		this->nnz = nz;
		this->nor = m;
		this->noc = n;
		this->ntt = n;
		jumps = 0;
		if( nz == 0 ) {
			r_inc = c_inc = NULL;
			vals  = NULL;
			return; //done
		}
		_i_value prevrow = row[ 0 ];
		for( unsigned long int i=1; i<this->nnz; i++ ) {
			if( row[ i ] != prevrow )
				jumps++;
			prevrow = row[ i ];
		}
#ifdef _DEBUG
		std::cout << jumps << " row jumps found." << std::endl;
#endif
		r_inc = new _i_value[ jumps + 1 ];
		c_inc = new _i_value[ this->nnz ];
		vals  = new _t_value[ this->nnz ];
		for( unsigned long int i=0; i<this->nnz; ++i ) vals[i] = val[i];
	
		r_start = row[ 0 ];
		prevrow = row[ 0 ];	
		c_start = col[ 0 ];
		int prevcol = col[ 0 ];
		r_end = row[ nz - 1 ];
		c_end = col[ nz - 1 ];

#ifdef _DEBUG
		std::cout << "c_inc: " << prevcol << std::endl;
		std::cout << "r_inc: " << prevrow << std::endl;
#endif
		int c = 0;
		for( unsigned long int i=1; i<this->nnz; i++ ) {
			this->c_inc[ i-1 ] = col[ i ] - prevcol;
			if( row[ i ] != prevrow ) {
				this->c_inc[ i-1 ] += ntt;
				this->r_inc[ c++ ] = row[ i ] - prevrow;
#ifdef _DEBUG
				std::cout << "c_inc: " << ntt << std::endl;
				std::cout << "r_inc: " << row[ i ] - prevrow << std::endl;
#endif
				prevrow = row[ i ];
			}
#ifdef _DEBUG
			else
				std::cout << "c_inc: " << col[ i ] - prevcol << std::endl;
#endif
			prevcol = col[ i ];
		}
		//overflow so to signal end of matrix
		c_inc[ this->nnz - 1 ] = ntt;
		//initialise last row jump to zero (prevent undefined jump)
		r_inc[ c ] = 0;

#ifdef _DEBUG
		std::cout << "Construction done." << std::endl;
#endif
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = this->r_start;
		col = this->c_start;
	}

	/**
	 *  Calculates y=xA, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zxa( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		const _t_value * y		= y_p;
		_i_value *__restrict__ c_inc_p	= c_inc;
		_i_value *__restrict__ r_inc_p	= r_inc;
		_t_value *__restrict__ v_p	= vals;

#ifndef NDEBUG
		const _t_value * x				= x_p;
		const _t_value * const x_end			= x+this->nor;
		const _i_value *__restrict__ const c_inc_end	= c_inc+this->nnz+1;
#endif
		const _t_value * const y_end			= y+this->noc;
		const _t_value *__restrict__ const v_end	= vals+this->nnz;

		y_p += this->c_start;
		x_p += this->r_start;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p < v_end );
			assert( x_p >= x );
			assert( x_p < x_end );
			assert( c_inc_p >= c_inc );
			assert( c_inc_p <  c_inc_end );
			assert( r_inc_p >= r_inc );
			while( y_p < y_end ) {
#ifdef _DEBUG
				std::cout << (x_p-x) << "," << (y_p-y) << " next increment: " << (*(c_inc_p+1))<< std::endl;
#endif
				*y_p += *v_p++ * *x_p;
				y_p += *c_inc_p++;
			}
			y_p -= ntt;
			x_p += *r_inc_p++;
		}
	}

	/**
	 *  Calculates y=Ax, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zax( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		const _t_value * x		= x_p;
		_i_value *__restrict__ c_inc_p	= c_inc;
		_i_value *__restrict__ r_inc_p	= r_inc;
		_t_value *__restrict__ v_p	= vals;

#ifndef NDEBUG
		const _t_value * y				= y_p;
		const _t_value * const y_end			= y+this->nor;
		const _i_value *__restrict__ const c_inc_end	= c_inc+this->nnz;
#endif
		const _t_value * const x_end			= x+this->noc;
		const _t_value *__restrict__ const v_end	= vals+this->nnz;

		x_p += c_start;
		y_p += r_start;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p < v_end );
			assert( x_p >= x );
			assert( x_p < x_end );
			assert( c_inc_p >= c_inc );
			assert( c_inc_p <  c_inc_end );
			assert( r_inc_p >= r_inc );
			while( x_p < x_end ) {
#ifdef _DEBUG
				std::cout << (y_p-y) << "," << (x_p-x) << " next increment: " << (*(c_inc_p+1))<< std::endl;
#endif
				*y_p += *v_p++ * *x_p;
				 x_p += *c_inc_p++;
			}
			x_p -= ntt;
			y_p += *r_inc_p++;
		}
	}

	/**
	 *  Calculates y=Ax, but does not allocate y itself. Does a front-to-back Hilbert traversal.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zax_fb( _t_value*__restrict__ x_f, _t_value*__restrict__ y_f ) {
		const _t_value * x		= x_f;
		_i_value *__restrict__ c_inc_f	= c_inc;
		_i_value *__restrict__ r_inc_f	= r_inc;
		_i_value *__restrict__ c_inc_b	= c_inc+this->nnz - 1;
		_i_value *__restrict__ r_inc_b	= r_inc+this->jumps;
		_t_value *__restrict__ v_f	= vals;
		_t_value *__restrict__ v_b	= vals+this->nnz - 1;
		_t_value *__restrict__ x_b      = x_f + this->noc - 1;
		_t_value *__restrict__ y_b      = y_f + this->nor - 1;
#ifndef NDEBUG
		const _t_value * y				= y_f;
		const _t_value * const y_end			= y+this->nor;
		const _i_value *__restrict__ const c_inc_end	= c_inc+this->nnz;
#endif
		const _t_value * const x_end			= x+this->noc;
		const _t_value *__restrict__ const v_end	= vals+this->nnz;

		x_f += c_start;
		y_f += r_start;
		x_b += c_end;
		y_b += r_end;
		while( v_f < v_end && v_b >= vals ) {
			assert( y_f >= y );
			assert( y_f <  y_end );
			assert( y_b >= y );
			assert( y_b <  y_end );
			assert( v_f >= vals );
			assert( v_f < v_end );
			assert( v_b >= vals );
			assert( v_b < v_end );
			assert( x_f >= x );
			assert( x_f < x_end );
			assert( x_b >= x );
			assert( x_b < x_end );
			assert( c_inc_f >= c_inc );
			assert( c_inc_f <  c_inc_end );
			assert( c_inc_b >= c_inc );
			assert( c_inc_b <  c_inc_end );
			assert( r_inc_b >= r_inc );
			assert( r_inc_f >= r_inc );
#ifdef _DEBUG
			std::cout << (y_p-y) << "," << (x_p-x) << " next increment: " << (*(c_inc_p+1))<< std::endl;
#endif
			*y_f += *v_f++ * *x_f;
			 x_f += *c_inc_f++;
			if( x_f >= x_end ) {
				x_f -= ntt;
				y_f += *r_inc_f++;
			}
			if( v_b < v_f ) break;
			*y_b += *v_b-- * *x_b;
			 x_b -= *c_inc_b--;
			if( x_b < x ) {
				x_b += ntt;
				y_b -= *r_inc_b++;
			}
		}
	}

	virtual size_t bytesUsed() {
		return sizeof( ULI ) * 4 + sizeof( _i_value ) * ( this->nnz + jumps + 1 ) + sizeof( _t_value ) * this->nnz;
	}
};

#endif

