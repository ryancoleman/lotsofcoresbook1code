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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2011 (June).
 */


#include "Matrix.hpp"
#include "Triplet.hpp"
#include "FileToVT.hpp"
#include "BICRS.hpp"
#include "CCSWrapper.hpp"
#include <assert.h>
#include <iostream>
#include <vector>

// Set when in debug mode
//#define _DEBUG

//NOTE: this file defines two classes: CBICRS and CBICRS_FACTORY

#ifndef _H_CBICRS
#define _H_CBICRS

//#define TRACE

/**
 *  Compressed Bi-directional Incremental Compressed Row Storage (BICRS) scheme.
 *  Compression is done by using less large data types for storing increments;
 *  e.g., 8-bit signed chars instead of 64-bit signed long ints.
 *  The exact data type used are input parameters.
 *
 *  Refer to CBICRS_factory for an auto-tuned selection.
 *
 *  @param _t_value The type of the nonzeros in the matrix.
 *
 *  Warning: this class uses assertions! For optimal performance,
 *           define the NDEBUG flag (e.g., pass -DNDEBUG as a compiler
 *	     flag).
 */
template< typename _t_value,
	typename _master_i_value=signed long int, typename _master_j_value=signed long int,
	typename _i_value=LI, typename _j_value=LI >
class CBICRS: public SparseMatrix< _t_value, ULI > {

     protected:

	/** Stores the row chunk start increments; size is the number of nonzeros plus one. */
	_master_i_value* r_start;

	/** Stores the column chunk start increments; size is the number of nonzeros plus one.*/
	_master_i_value* c_start;

	/** Stores the row jumps; size is the number of nonzeros plus 2. */
	_i_value* r_inc;

	/** Stores the column jumps; size is exactly the number of nonzeros. */
	_i_value* c_inc;

	/** Bitmask used for switching between c_start and c_ind. */
	unsigned char* mask1;

	/** Bitmask used for switching between r_start and r_ind. */
	unsigned char* mask2;

	/** Stores the values of the individual nonzeros. Size is exactly the number of nonzeros. */
	_t_value* vals;

	/** Stores the number of bytes used for storage. */
	size_t bytes;

	/** Caches n times two. */
	_master_j_value ntt;
	
	/** 
 	 * Calculates the number of overflows given a triplet-form input.
 	 *
 	 *  @param ntt should equal the number of columns times two (n times two)
 	 */
	static void getNumberOfOverflows( const ULI nnz, ULI * const row, ULI * const col, const ULI ntt,
					 ULI &row_overflows, ULI &col_overflows, ULI &sim_overflows, ULI &jumps ) {
		unsigned long int row_max = (((unsigned long int)1) << (sizeof(_i_value)*8-1)) - 1; //2^{k-1}-1
		unsigned long int col_max = (((unsigned long int)1) << (sizeof(_j_value)*8-1)) - 1;
		ULI prevrow = row[ 0 ];
		ULI prevcol = col[ 0 ];
		for( unsigned long int i=1; i<nnz; i++ ) {
			bool overflow = false;
			if( row[ i ] > prevrow ) {
				if( row[ i ] - prevrow > row_max ) {
					row_overflows++;
					overflow = true;
				}
			} else if( prevrow > row[ i ] ) {
				if( prevrow - row[ i ] > row_max ) {
					row_overflows++;
					overflow = true;
				}
			}
			if( row[ i ] != prevrow ) {
				jumps++;
				if( col[ i ] - prevcol + ntt > col_max ) {
					if( overflow ) {
						row_overflows--;
						sim_overflows++;
					} else
						col_overflows++;
				}
			} else {
				if( col[ i ] > prevcol ) {
					if( col[ i ] - prevcol > col_max ) {
						col_overflows++;
					}
				} else {
					if( prevcol - col[ i ] > col_max ) {
						col_overflows++;
					}
				}
			}

			prevrow = row[ i ];
			prevcol = col[ i ];
		}
	}


	/** Estimates the number of bytes required by this data structure. */
	static unsigned long int memoryUsage( const ULI nnz, const ULI jumps, const ULI row_o, const ULI col_o, const ULI sim_o ) {
		return 	nnz/8 + (((nnz % 8) > 0) ? 2 : 1) +
			jumps/8 + (((jumps % 8) > 0) ? 1 : 0) +
			sizeof(_master_i_value) * (row_o + sim_o + 1) +
			sizeof(_master_j_value) * (col_o + sim_o + 2) +
			sizeof(_i_value)        * (jumps - row_o - sim_o) +
			sizeof(_j_value)        * (nnz   - col_o - sim_o) +
			sizeof(_t_value)        * nnz;
	}

     public:

	/** Calculates and returns the number of bytes used when employing this data structure. */
	static unsigned long int getMemoryUsage( ULI *row, ULI *col, const ULI nz, const ULI m, const ULI n ) {
		ULI jumps = 0;
		ULI row_overflows = 0; 
		ULI col_overflows = 0;
		ULI sim_overflows = 0; //simultaneous overflows
		getNumberOfOverflows( nz, row, col, 2l*n, row_overflows, col_overflows, sim_overflows, jumps );
		return memoryUsage( nz, jumps, row_overflows, col_overflows, sim_overflows );
	}

	/** Calculates and returns the number of bytes used when employing this data structure. */
	static unsigned long int getMemoryUsage( std::vector< Triplet< _t_value > > &input, const ULI m, const ULI n ) {
		ULI nz = input.size();
		ULI* row = new ULI[ nz ];
		ULI* col = new ULI[ nz ];
		unsigned long int c = 0;
		typename std::vector< Triplet< _t_value > >::iterator it = input.begin();
		for( ; it!=input.end(); it++, c++ ) {
			row[ c ] = (*it).i();
			col[ c ] = (*it).j();
		}
		unsigned long int ret = getMemoryUsage( row, col, nz, m, n );
		delete [] row;
		delete [] col;
		return ret;
	}

	/** Base deconstructor. */
	virtual ~CBICRS() {
		delete [] r_start;
		delete [] c_start;
		delete [] r_inc;
		delete [] c_inc;
		delete [] mask1;
		delete [] mask2;
		delete [] vals;
	}

	/** Base constructor. */
	CBICRS() {
		r_start = c_start = NULL;
		r_inc = c_inc = NULL;
		mask1 = mask2 = NULL;
		vals = NULL;
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	CBICRS( std::string file, _t_value zero = 0 ) {
		loadFromFile( file, zero );
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
	CBICRS( ULI* row, ULI* col, _t_value* val, ULI m, ULI n, ULI nz, _t_value zero ) {
		load( row, col, val, m, n, nz, zero );
	}

	/** Base constructor.
	 *  @see SparseMatrix::SparseMatrix( input, m, n, zero )
	 */
	CBICRS( std::vector< Triplet< _t_value > >& input, ULI m, ULI n, _t_value zero = 0 ) {
		load( input, m, n, zero );
	}

	/**
	 *  This function will rewrite the std::vector< Triplet > structure to one suitable
	 *  for the other load function.
	 *  @see load( row, col, val, m, n, nz ) 
	 *  @see SparseMatrix::load 
	 */
	virtual void load( std::vector< Triplet< _t_value > >& input, ULI m, ULI n, _t_value zero ) {
		ULI nz = input.size();
		ULI* row = new ULI[ nz ];
		ULI* col = new ULI[ nz ];
		_t_value* val = new _t_value[ nz ];
		unsigned long int c = 0;
		typename std::vector< Triplet< _t_value > >::iterator it = input.begin();
		for( ; it!=input.end(); it++, c++ ) {
			row[ c ] = (*it).i();
			col[ c ] = (*it).j();
			val[ c ] = (*it).value;
		}
		load( row, col, val, m, n, nz, zero );
		delete [] row;
		delete [] col;
		delete [] val;
	}

	/** @see CBICRS( row, col, val, m, n, nz ) */
	void load( ULI* row, ULI* col, _t_value* val, ULI m, ULI n, ULI nz, _t_value zero ) {
#ifdef _DEBUG
		std::cerr << "Warning: _DEBUG flag set." << std::endl;
#endif
		this->zero_element = zero;
		this->nnz = nz;
		this->nor = m;
		this->noc = n;
		this->ntt = 2l*n;
		unsigned long int row_max = (((unsigned long int)1) << (sizeof(_i_value)*8-1)) - 1; //2^{k-1}-1
		unsigned long int col_max = (((unsigned long int)1) << (sizeof(_j_value)*8-1)) - 1;
		ULI jumps = 0;
		ULI row_overflows = 0; 
		ULI col_overflows = 0;
		ULI sim_overflows = 0; //simultaneous overflows
		getNumberOfOverflows( nz, row, col, this->ntt, row_overflows, col_overflows, sim_overflows, jumps );
#ifndef NDEBUG
		std::cout << jumps << " row jumps found." << std::endl;
		std::cout << row_overflows << " exclusive row overflows found." << std::endl;
		std::cout << col_overflows << " exclusive column overflows found." << std::endl;
		std::cout << sim_overflows << " simultaneous row/column overflows found." << std::endl;
		std::cout << "Total array memory usage: "
				<< memoryUsage( this->nnz, jumps, row_overflows, col_overflows, sim_overflows )
				<< " bytes." << std::endl;
#endif

		//allocate arrays
		const size_t mask1size = (this->nnz+1)/8 + (((this->nnz+1)%8) > 0 ? 1 : 0) + 1;
		const size_t mask2size = jumps / 8 + ((jumps % 8) > 0 ? 1 : 0 ) + 1;
		mask1   = new unsigned char[ mask1size ];
		mask2   = new unsigned char[ mask2size ];
		r_start = new _master_i_value[ row_overflows + sim_overflows + 1 ];
		c_start = new _master_j_value[ col_overflows + sim_overflows + 2 ];
		r_inc = new _i_value[ jumps - row_overflows - sim_overflows ];
		c_inc = new _j_value[ this->nnz - col_overflows - sim_overflows - 1 ];
		vals  = new _t_value[ this->nnz ];
		
		//remember memory usage
		bytes  = sizeof( unsigned char ) * ( mask1size + mask2size );
		bytes += sizeof( _master_i_value ) * ( row_overflows + col_overflows + 2 * sim_overflows + 3 );
		bytes += sizeof( _i_value ) * ( jumps - row_overflows - sim_overflows );
		bytes += sizeof( _j_value ) * ( this->nnz - col_overflows - sim_overflows - 1 );
		bytes += sizeof( _t_value ) * this->nnz;

		//copy nonzero values
		for( unsigned long int i=0; i<this->nnz; ++i ) vals[i] = val[i];
	
		//fill initial values
		r_start[ 0 ] = (_master_i_value)row[ 0 ];
		ULI prevrow = row[ 0 ];	
		c_start[ 0 ] = (_master_j_value)col[ 0 ];
		ULI prevcol = col[ 0 ];

#ifdef _DEBUG
		std::cout << "r_start: " << r_start[0] << std::endl;
		std::cout << "c_start: " << c_start[0] << std::endl;
#endif

		unsigned long int cincc   = 0;
		unsigned long int rincc   = 0;
		unsigned long int cstartc = 1;
		unsigned long int rstartc = 1;
		unsigned long int mask2c  = 0;
		this->mask1[ 0 ] = 1;
		for( unsigned long int i=1; i<this->nnz; i++ ) {
			if( i%8 == 0 ) this->mask1[ i/8 ] = 0;
			if( mask2c%8 == 0 ) this->mask2[ mask2c/8 ] = 0;

			//different handling upon row changes
			if( row[ i ] != prevrow ) {
				if( static_cast< unsigned long int >( col[ i ] + ntt - prevcol ) > col_max ) {
					assert( cstartc < col_overflows + sim_overflows + 1 );
					this->c_start[ cstartc++ ] = col[ i ] - prevcol + ntt;
					this->mask1[ i/8 ] |= ((unsigned char)1)<<(i%8);
				} else {
					assert( cincc < this->nnz );
					this->c_inc[ cincc++ ] = col[ i ] - prevcol + ntt;
				}
				if( row[ i ] > prevrow ) {
					if( row[ i ] - prevrow > row_max ) {
						assert( rstartc < row_overflows + sim_overflows + 1 );
						this->r_start[ rstartc++ ] = row[ i ] - prevrow;
						this->mask2[ mask2c/8 ] |= ((unsigned char)1)<<(mask2c%8);mask2c++;
					} else {
						assert( rincc < jumps - row_overflows );
						this->r_inc[ rincc++ ] = row[ i ] - prevrow;
						mask2c++;
					}
				} else {
					if( prevrow - row[ i ] > row_max ) {
						assert( rstartc < row_overflows + sim_overflows + 1 );
						this->r_start[ rstartc++ ] = row[ i ] - prevrow;
						this->mask2[ mask2c/8 ] |= ((unsigned char)1)<<(mask2c%8);mask2c++;
					} else {
						assert( rincc < jumps - row_overflows );
						this->r_inc[ rincc++ ] = row[ i ] - prevrow;
						mask2c++;
					}
				}
				prevrow = row[ i ];
				prevcol = col[ i ];
#ifdef _DEBUG
				std::cout << i << ", (" << prevrow << "," << prevcol << ")), column increment = " <<
                	        ( this->mask1[ i/8 ] & ((unsigned char)1)<<(i%8) ? this->c_start[ cstartc-1 ] : this->c_inc[ cincc-1 ] ) <<
        	                ", row increment " << ( this->mask2[ (mask2c-1)/8 ] & ((unsigned char)1)<<((mask2c-1)%8) ?
	                        this->r_start[ rstartc-1 ] : this->r_inc[ rincc-1 ] ) << ", mask1 index " <<
				(i/8) << "(" << ( (this->mask1[ i/8 ] & ((unsigned char)1)<<(i%8)) > 0 ) << ")" << std::endl;
#endif
				continue;
			}

			//detect normal column overflow
			if( col[ i ] > prevcol ) {
				if( col[ i ] - prevcol > col_max ) {
					assert( cstartc < col_overflows + sim_overflows + 1 );
					this->c_start[ cstartc++ ] = col[ i ] - prevcol;
					this->mask1[ i/8 ] |= ((unsigned char)1)<<(i%8);
				} else {
					assert( cincc < this->nnz );
					this->c_inc[ cincc++ ] = col[ i ] - prevcol;
				}
			} else if( prevcol > col[ i ] ) {
				if( prevcol - col[ i ] > col_max ) {
					assert( cstartc < col_overflows + sim_overflows + 1 );
					this->c_start[ cstartc++ ] = col[ i ] - prevcol;
					this->mask1[ i/8 ] |= ((unsigned char)1)<<(i%8);
				} else {
					assert( cincc < this->nnz );
					this->c_inc[ cincc++ ] = col[ i ] - prevcol;
				}
			}
			prevcol = col[ i ];
#ifdef _DEBUG
			std::cout << i << ", (" << prevrow << "," << prevcol << "), column increment = " <<
			( this->mask1[ i/8 ] & ((unsigned char)1)<<(i%8) ? this->c_start[ cstartc-1 ] : this->c_inc[ cincc-1 ] ) <<
			", row increment " << ( this->mask2[ (mask2c-1)/8 ] & ((unsigned char)1)<<((mask2c-1)%8) ? 
			this->r_start[ rstartc-1 ] : this->r_inc[ rincc-1 ] ) << ", mask1 index = " << (i/8) <<
			"(" << ( (this->mask1[ i/8 ] & ((unsigned char)1)<<(i%8)) > 0 ) << ")" << std::endl;
#endif
		}
		assert( cincc == this->nnz - col_overflows - sim_overflows - 1 );
		assert( rincc == jumps - row_overflows - sim_overflows );
		assert( cstartc == col_overflows + sim_overflows + 1 ); //last one is added below
		assert( rstartc == row_overflows + sim_overflows + 1 );
		assert( mask2c == jumps );

		//force last overflow
		c_start[ col_overflows + sim_overflows + 1 ] = ntt;
		this->mask1[ this->nnz/8 ] |= ((unsigned char)1)<<(this->nnz%8);

		//check number of overflow masks
		assert( this->mask1[0] & 1 );
		assert( this->mask1[this->nnz/8] & ((unsigned char)1<<(this->nnz%8)) );
		unsigned long int mask1c = 0;
		for( unsigned long int k=0; k<=this->nnz; k++ )
			if( this->mask1[ k/8 ] & ((unsigned char)1<<(k%8)) )
				mask1c++;
		assert( mask1c == col_overflows + sim_overflows + 2 );
		assert( this->nnz+1-mask1c == this->nnz - col_overflows - sim_overflows - 1 );
#ifdef _DEBUG
		std::cout << "Construction done." << std::endl;
#endif
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = (ULI)(this->r_start[ 0 ]);
		col = (ULI)(this->c_start[ 0 ]);
	}

	/**
	 *  Calculates y=xA, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zxa( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		unsigned char *__restrict__ mask1_p	= this->mask1;
		unsigned char *__restrict__ mask2_p	= this->mask2;
		_master_i_value *__restrict__ r_start_p = r_start;
		_master_j_value *__restrict__ c_start_p = c_start;
		_i_value *__restrict__ r_inc_p		= r_inc;
		_j_value *__restrict__ c_inc_p		= c_inc;
		_t_value *__restrict__ v_p		= vals;
		char maskc1 = 0;
		char maskc2 = 0;

#ifndef NDEBUG
		const _t_value * const x	= x_p;
		const _t_value * const x_end	= x+this->nor;
#endif
		const _t_value * const y	= y_p;
		const _t_value * const y_end	= y+this->noc;
		const _t_value * const v_end	= vals+this->nnz;

		x_p += *r_start_p++;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p < v_end );
			assert( x_p >= x );
			assert( x_p < x_end );
			assert( c_inc_p >= c_inc );
			assert( r_inc_p >= r_inc );
			assert( mask1_p < this->mask1 + (this->nnz/8 + (this->nnz%8==0 ? 0 : 1)) );
			assert( mask1_p >= this->mask1 );
			assert( maskc1 == (v_p-vals) % 8 );
			assert( mask1_p == &(this->mask1[ (v_p-vals)/8 ]) );

			if( *mask1_p & ((unsigned char)1<<maskc1) ) {
#ifdef TRACE
				std::cout << "Overflowed column increment is " << *c_start_p << std::endl;
#endif
				y_p += *c_start_p++;
			} else {
#ifdef TRACE
				std::cout << "Compressed column increment is " << *c_inc_p << std::endl;
#endif
				y_p += *c_inc_p++;
			}
			if( ++maskc1 == 8 ) {
				maskc1 = 0;
				mask1_p++;
			}
			if( y_p >= y_end ) {
#ifdef TRACE
				std::cout << (y_p-y) << " > " << this->noc << " so performing a row jump." << std::endl;
#endif
				if( *mask2_p & ((unsigned char)1<<maskc2) ) {
					x_p += *r_start_p++;
				} else {
					x_p += *r_inc_p++;
				}
				if( ++maskc2 == 8 ) {
					maskc2 = 0;
					mask2_p++;
				}
				y_p -= ntt;
			}
#ifdef TRACE
			assert( mask1_p == &(this->mask1[ (v_p - vals + 1)/8 ]) );
			std::cout << (v_p-vals) << ", position: " << (x_p-x) << "(<=" << (this->noc) << ") by " << (y_p-y) << "(<=" << (this->nor) << "), mask1 index is " << (mask1_p-this->mask1) << ", mask1 was " << ((this->mask1[(v_p-vals)/8]&(unsigned char)1<<((v_p-vals)%8))>0) << std::endl;
#endif
			*y_p += *v_p++ * *x_p;
		}
	}

	/**
	 *  Calculates y=Ax, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zax( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		const _t_value * const x		= x_p;
		unsigned char *__restrict__ mask1_p	= this->mask1;
		unsigned char *__restrict__ mask2_p	= this->mask2;
		_master_i_value *__restrict__ r_start_p = r_start;
		_master_j_value *__restrict__ c_start_p = c_start;
		_i_value *__restrict__ r_inc_p		= r_inc;
		_j_value *__restrict__ c_inc_p		= c_inc;
		_t_value *__restrict__ v_p		= vals;
		/* CODE SWITCH, non-unrolled version: */
		char maskc1 = 0;
		char maskc2 = 0;
		unsigned char tmask1 = *mask1_p++;
		unsigned char tmask2 = *mask2_p++;
		//END, non-unrolled version */
		/* CODE SWITCH, unrolled version:
		char maskc2 = 1;
		unsigned char tmask1;
		unsigned char tmask2 = *mask2_p++;
		//END, unrolled version */
		
#ifndef NDEBUG
		const _t_value * const y	= y_p;
		const _t_value * const y_end	= y+this->nor;
#endif
		const _t_value * const x_end	= x+this->noc;
		const _t_value * const v_end	= vals+this->nnz;

		y_p += *r_start_p++;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p < v_end );
			assert( x_p >= x );
			assert( x_p < x_end );
			assert( c_inc_p >= c_inc );
			assert( r_inc_p >= r_inc );
			if ( mask1_p > this->mask1 + ((this->nnz+1)/8 + ((this->nnz+1)%8==0 ? 0 : 1)) ) {
				std::cout << "Mask1 is at start position for index " << (mask1_p-this->mask1)*8 << 
				" of " << this->nnz << ", maskc1=" << (int)maskc1 << std::endl;
			}
			assert( mask1_p <= this->mask1 + ((this->nnz+1)/8 + ((this->nnz+1)%8==0 ? 0 : 1)) );
			assert( mask1_p >= this->mask1 );

			/* CODE SWITCH, non-unrolled version: */
			if( tmask1 & 1 ) {
#ifdef TRACE
				if ( mask1_p >=this->mask1 + ((this->nnz+1)/8 + ((this->nnz+1)%8==0 ? 0 : 1)) )
				std::cout << "Overflowed column increment is " << (int)(*c_start_p) << std::endl;
#endif
				x_p += *c_start_p++;
			} else {
#ifdef TRACE
				if ( mask1_p >=this->mask1 + ((this->nnz+1)/8 + ((this->nnz+1)%8==0 ? 0 : 1)) )
				std::cout << "Compressed column increment is " << (int)(*c_inc_p) << std::endl;
#endif
				x_p += *c_inc_p++;
			}
			tmask1 >>= 1;
 			if( ++maskc1 == 8 ) {
				maskc1 = 0;
				tmask1 = *mask1_p++;
			}
			if( x_p >= x_end ) {
#ifdef TRACE
				std::cout << (x_p-x) << " > " << this->noc << " so performing a row jump." << std::endl;
#endif
				if( tmask2 & 1 ) {
					y_p += *r_start_p++;
				} else {
					y_p += *r_inc_p++;
				}
				tmask2 >>= 1;
				if( ++maskc2 == 8 ) {
					maskc2 = 0;
					tmask2 = *mask2_p++;
				}
				x_p -= ntt;
			}
#ifdef TRACE
			std::cout << (v_p-vals) << ", position: " << (y_p-y) << "(<=" << (this->nor) << ") by " << (x_p-x) << "(<=" << (this->noc) << "), mask1 index is " << (mask1_p-this->mask1) << std::endl;
#endif
			*y_p += *v_p++ * *x_p;
			// END, non-unrolled version. */
			/* CODE SWITCH, unrolled version: 
			tmask1 = *mask1_p++;
			// 1
			if( tmask1 & 1 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			std::cout << (v_p-vals) << ", position: " << (y_p-y) << "(<=" << (this->nor) << ") by " << (x_p-x) << "(<=" << (this->noc) << "), mask1 index is " << (mask1_p-this->mask1) << std::endl;
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 2
			if( tmask1 & 2 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 4
			if( tmask1 & 4 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 8
			if( tmask1 & 8 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 16
			if( tmask1 & 16 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 32
			if( tmask1 & 32 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 64
			if( tmask1 & 64 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			if( v_p >= v_end ) break;
			// 128
			if( tmask1 & 128 )
				x_p += *c_start_p++;
			else
				x_p += *c_inc_p++;
			if( x_p >= x_end ) {
				if( tmask2 & maskc2 )
					y_p += *r_start_p++;
				else
					y_p += *r_inc_p++;
				if( maskc2 == 128 ) {
					tmask2 = *mask2_p++;
					maskc2 = 1;
				} else
					maskc2 *= 2;
				x_p -= ntt;
			}
			*y_p += *v_p++ * *x_p;
			//END, unrolled version */
		}
	}

	virtual size_t bytesUsed() {
		return bytes;
	}
};


#endif


#ifndef _H_CBICRS_FACTORY
#define _H_CBICRS_FACTORY

/**
 *  Factory for the Compressed Bi-directional Incremental Compressed Row Storage scheme.
 *  Auto-tunes index type (signed char, short int, int or long int).
 *  May revert back to a plain BICRS implementation.
 *
 *  Warning: assertions may be used! For optimal performance,
 *           define the NDEBUG flag (pass -DNDEBUG as a compiler flag).
 */
template< typename _t_value >
class CBICRS_factory {

     protected:

	/** Used for auto-tunes the index type */
	template< typename _master_i_value, typename _master_j_value, typename _i_value, typename _j_value >
	static void investigateCCS( const std::string tn, std::vector< Triplet< _t_value > > input,
					unsigned long int m, unsigned long int n, unsigned long int &usage, unsigned long int &zle_usage ) {
		ULI nz   = input.size();
		ULI* row = new ULI[ nz ];
		ULI* col = new ULI[ nz ];
		unsigned long int c = 0;
		typename std::vector< Triplet< _t_value > >::iterator it = input.begin();
		for( ; it!=input.end(); it++, c++ ) {
			row[ c ] = (*it).i();
			col[ c ] = (*it).j();
		}
		usage = CBICRS< _t_value, _master_i_value, _master_j_value, _i_value, _j_value >::getMemoryUsage( col, row, nz, n, m ); //note the reversal
		delete [] row;
		delete [] col;

#ifndef NDEBUG
		std::cout << "Total array (" << tn << ") memory usage: " << usage << " bytes." << std::endl;
#endif
		zle_usage = usage + 1;
#ifndef NDEBUG
		std::cout << "Warning: no implementation for ZLE CBICRS yet, so no estimate is given!" << std::endl;
#endif
	}

	/** Used for auto-tuning of the index type */
	template< typename _master_i_value, typename _master_j_value, typename _i_value, typename _j_value >
	static void investigate( const std::string tn, std::vector< Triplet< _t_value > > triplets,
					ULI m, ULI n, unsigned long int &usage, unsigned long int &zle_usage ) {
		usage = CBICRS< _t_value, _master_i_value, _master_j_value, _i_value, _j_value >::getMemoryUsage( triplets, m, n );
#ifndef NDEBUG
		std::cout << "Total array (" << tn << ") memory usage: " << usage << " bytes." << std::endl;
#endif
		zle_usage = usage + 1;
#ifndef NDEBUG
		std::cout << "Warning: no implementation for ZLE CBICRS yet, so no estimate is given!" << std::endl;
#endif
	}

     public:

	static Matrix< _t_value >* getCBICRS( std::string file, _t_value zero = 0 ) {
		ULI m, n;
		std::vector< Triplet< double > > triplets = FileToVT::parse( file, m, n );
		return getCBICRS( triplets, m, n, zero ); 
	}

	static Matrix< _t_value >* getCBICCS( std::string file, _t_value zero = 0 ) {
		ULI m, n;
		std::vector< Triplet< double > > triplets = FileToVT::parse( file, m, n );
		return getCBICCS( triplets, m, n, zero ); 
	}

	static Matrix< _t_value >* getCBICRS( std::vector< Triplet< _t_value > > &triplets, unsigned long int m, unsigned long int n, _t_value zero = 0 ) {
		unsigned int rowbit = log2( m );
		if( ((unsigned long int)1)<<rowbit < m ) rowbit++;
		unsigned int colbit = log2( n );
		if( ((unsigned long int)1)<<colbit < n ) colbit++;
		std::cout << "Finding optimal expected index type." << std::endl;
		unsigned long int usage, zle_usage;
		/*if( rowbit > 31 ) {
			if( colbit > 31 )
				investigate< signed long int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigate< signed long int, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigate< signed long int, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigate< signed long int, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		} else if( rowbit > 15 ) {
			if( colbit > 31 )
				investigate< signed int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigate< signed int, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigate< signed int, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigate< signed int, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		} else if( rowbit > 7 ) {
			if( colbit > 31 )
				investigate< signed int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigate< signed int, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigate< signed int, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigate< signed int, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		} else {
			if( colbit > 31 )
				investigate< signed char, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigate< signed char, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigate< signed char, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigate< signed char, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		}*/
		investigate< signed long int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );

		char choice  = 0;
		bool zle     = zle_usage < usage;
		unsigned long int min = zle ? zle_usage : usage;
		investigate< signed long int, signed long int, signed short int, signed short int >( "short int", triplets, m, n, usage, zle_usage );
		unsigned long int curmin = zle_usage < usage ? zle_usage : usage;
		if( curmin < min ) {
			choice = 1;
			zle    = zle_usage < usage;
			min    = curmin;
		}
		investigate< signed long int, signed long int, signed int, signed int >( "int", triplets, m, n, usage, zle_usage );
		curmin = zle_usage < usage ? zle_usage : usage;
		if( curmin < min ) {
			choice = 2;
			zle    = zle_usage < usage;
			min    = curmin;
		}
		investigate<  signed long int, signed long int, signed long int, signed long int >( "long int", triplets, m, n, usage, zle_usage );
		curmin = zle_usage < usage ? zle_usage : usage;
		if( curmin < min ) {
			choice = 3;
			zle    = zle_usage < usage;
			min    = curmin;
		}
		switch( choice ) {
			case 0:
				std::cout << "Selecting `signed char' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed char)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed char)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new BICRS< _t_value, signed char >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CBICRS< _t_value, signed long int, signed long int, signed char, signed char >( triplets, m, n, zero );
				}
				break;
			case 1:
				std::cout << "Selecting `signed short int' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed short int)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed short int)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new BICRS< _t_value, signed short int >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CBICRS< _t_value, signed long int, signed long int, signed short int, signed short int >( triplets, m, n, zero );
				}
				break;
			case 2:
				std::cout << "Selecting `signed int' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed int)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed int)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new BICRS< _t_value, signed int >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CBICRS< _t_value, signed long int, signed long int, signed int, signed int >( triplets, m, n, zero );
				}
				break;
			case 3:
				std::cout << "Selecting `signed long int' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed long int)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed long int)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new BICRS< _t_value, signed long int >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CBICRS< _t_value, signed long int, signed long int, signed long int, signed long int >( triplets, m, n, zero );
				}
				break;
			default:
				std::cerr << "Error in tuning, invalid data type selected (" << choice << ")!" << std::endl;
				exit( EXIT_FAILURE );
		}
		std::cerr << "CBICRS not yet implemented!" << std::endl;
		exit( EXIT_FAILURE );
	}

	static Matrix< _t_value >* getCBICCS( std::vector< Triplet< _t_value > > &triplets, unsigned long int m, unsigned long int n, _t_value zero = 0 ) {
		unsigned int rowbit = log2( m );
		if( ((unsigned long int)1)<<rowbit < m ) rowbit++;
		unsigned int colbit = log2( n );
		if( ((unsigned long int)1)<<colbit < n ) colbit++;
		std::cout << "Finding optimal expected index type." << std::endl;
		unsigned long int usage, zle_usage;
		/*if( rowbit > 31 ) {
			if( colbit > 31 )
				investigateCCS< signed long int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigateCCS< signed long int, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigateCCS< signed long int, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigateCCS< signed long int, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		} else if( rowbit > 15 ) {
			if( colbit > 31 )
				investigateCCS< signed int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigateCCS< signed int, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigateCCS< signed int, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigateCCS< signed int, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		} else if( rowbit > 7 ) {
			if( colbit > 31 )
				investigateCCS< signed int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigateCCS< signed int, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigateCCS< signed int, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigateCCS< signed int, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		} else {
			if( colbit > 31 )
				investigateCCS< signed char, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 15 )
				investigateCCS< signed char, signed int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else if( colbit > 7 )
				investigateCCS< signed char, signed short int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
			else
				investigateCCS< signed char, signed char, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );
		}*/
		investigateCCS< signed long int, signed long int, signed char, signed char >( "signed char", triplets, m, n, usage, zle_usage );

		char choice  = 0;
		bool zle     = zle_usage < usage;
		unsigned long int min = zle ? zle_usage : usage;
		investigateCCS< signed long int, signed long int, signed short int, signed short int >( "short int", triplets, m, n, usage, zle_usage );
		unsigned long int curmin = zle_usage < usage ? zle_usage : usage;
		if( curmin < min ) {
			choice = 1;
			zle    = zle_usage < usage;
			min    = curmin;
		}
		investigateCCS< signed long int, signed long int, signed int, signed int >( "int", triplets, m, n, usage, zle_usage );
		curmin = zle_usage < usage ? zle_usage : usage;
		if( curmin < min ) {
			choice = 2;
			zle    = zle_usage < usage;
			min    = curmin;
		}
		investigateCCS<  signed long int, signed long int, signed long int, signed long int >( "long int", triplets, m, n, usage, zle_usage );
		curmin = zle_usage < usage ? zle_usage : usage;
		if( curmin < min ) {
			choice = 3;
			zle    = zle_usage < usage;
			min    = curmin;
		}
		switch( choice ) {
			case 0:
				std::cout << "Selecting `signed char' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed char)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed char)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new CCSWrapper< _t_value, BICRS< _t_value, signed char >, ULI >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CCSWrapper< _t_value, CBICRS< _t_value, signed long int, signed long int, signed char, signed char >, ULI >( triplets, m, n, zero );
				}
				break;
			case 1:
				std::cout << "Selecting `signed short int' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed short int)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed short int)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new CCSWrapper< _t_value, BICRS< _t_value, signed short int >, ULI >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CCSWrapper< _t_value, CBICRS< _t_value, signed long int, signed long int, signed short int, signed short int >, ULI >( triplets, m, n, zero );
				}
				break;
			case 2:
				std::cout << "Selecting `signed int' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed int)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed int)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new CCSWrapper< _t_value, BICRS< _t_value, signed int >, ULI >( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CCSWrapper< _t_value, CBICRS< _t_value, signed long int, signed long int, signed int, signed int >, ULI >( triplets, m, n, zero );
				}
				break;
			case 3:
				std::cout << "Selecting `signed long int' datatype";
				if( zle ) {
					std::cout << ", with ZLE" << std::endl;
				} else {
					std::cout << ", without ZLE";
					if( 2*m < (((unsigned long int)1)<<(sizeof(signed long int)*8-1))-1 &&
						n < (((unsigned long int)1)<<(sizeof(signed long int)*8-1))-1 ) {
						std::cout << "; matrix dimensions are small enough; reverting to plain BICRS." << std::endl;
						return new CCSWrapper< _t_value, BICRS< _t_value, signed long int > , ULI>( triplets, m, n, zero );
					} // else
					std::cout << std::endl;
					return new CCSWrapper< _t_value, CBICRS< _t_value, signed long int, signed long int, signed long int, signed long int >, ULI >( triplets, m, n, zero );
				}
				break;
			default:
				std::cerr << "Error in tuning, invalid data type selected (" << choice << ")!" << std::endl;
				exit( EXIT_FAILURE );
		}
		std::cerr << "CBICRS not yet implemented!" << std::endl;
		exit( EXIT_FAILURE );
	}

};

#endif

