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
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2011.
 */


#include "SparseMatrix.hpp"
#include "ICRS.hpp"
#include "CRS.hpp"
#include <assert.h>
#include <iostream>
#include <stdint.h>

#ifndef _H_FBICRS
#define _H_FBICRS

/**
 *  Hierarchical BICRS with fixed subblock size and distribution.
 *
 *  Uses compression on upper and lower level. Currently handles
 *  matrices with at most 39 bits indices (signed ints on upper
 *  level with unsigned chars on the lower level; it was either
 *  that or 72 bits, but the latter cannot efficiently be
 *  represented on 64-bit architectures).
 *
 *  @tparam _sub_ds Defaults chooses a vectorised BICRS structure
 *                  with unsigned short int index type as sub
 *                  data type. Other options:
 *         typedef  ICRS< _t_value, unsigned char > _sub_ds;
 *	   typedef  ICRS< _t_value, uint16_t >      _sub_ds;
 *	   typedef oICRS< _t_value, 8, 1, int16_t > _sub_ds; (default)
 *	   typedef  ICRS< _t_value, unsigned int >  _sub_ds;
 *	   typedef   CRS< _t_value > _sub_ds;
 *	   (using the latter with logBeta = 64 boils down to a
 *	    block row distributed CRS-based parallel SpMV)
 *  @tparam logBeta The log of the size of the maximum block size
 *                  represented by _sub_ds. Options:
 *	static const unsigned char logBeta = 7;  //=max ICRS/unsigned char
 *	static const unsigned char logBeta = 11; //=2k
 *	static const unsigned char logBeta = 12; //=4k
 *	static const unsigned char logBeta = 14; //=max oICRS/int16_t=16k
 *	static const unsigned char logBeta = 15; //=max ICRS/uint16_t=32k
 *	static const unsigned char logBeta = 17; //=128k
 *	static const unsigned char logBeta = 31; //=max ICRS/unsigned short
 *	static const unsigned char logBeta = 64; //=max CRS
 */
template< typename _t_value, typename _i_value=LI, typename _sub_ds = oICRS< _t_value, 8, 1, int16_t >, unsigned char logBeta = 14 >
class FBICRS: public SparseMatrix< _t_value, _i_value > {

    protected:

	ULI r_start;
	ULI c_start;
	ULI r_end;
	ULI c_end;
	ULI jumps;
	_i_value* r_inc;
	_i_value* c_inc;

   public:

	/** Maximum matrix size for _sub_ds with the above data type; row size */
	static const ULI beta_m = 1l << (logBeta-1);

	/** Maximum matrix size for _sub_ds with the above data type; column size */
	static const ULI beta_n = beta_m;

	/** Stores the total fillIn. */
	size_t fillIn;

	/** Stores the lower-level data structures */
	_sub_ds** dss;

	/** Caches the overflow forcing value */
	ULI n_overflow;

	virtual ~FBICRS() {
		if( r_inc != NULL ) delete [] r_inc;
		if( c_inc != NULL ) delete [] c_inc;
		for( size_t i = 0; i < static_cast< size_t >(this->nnz); i++ ) {
			if( dss[ i ] != NULL ) {
				delete dss[ i ];
			}
		}
		if( dss != NULL ) {
			delete [] dss;
		}
	}

	FBICRS(): r_start( -1 ), c_start( -1 ), r_inc( NULL ), c_inc( NULL ), dss( NULL ), n_overflow( 0 ) {}

	FBICRS( std::string file, _t_value zero = 0 ) {
		loadFromFile( file, zero );
	}

	FBICRS( _i_value* row, _i_value* col, _t_value* val, ULI m, ULI n, ULI nz,  _t_value zero ) {
		this->load( row, col, val, m, n, nz, zero );
	}

	FBICRS( std::vector< Triplet< _t_value > > &input, ULI m, ULI n, _t_value zero = 0 ) {
		this->load( input, m, n, zero );
	}

	FBICRS( std::vector< std::vector< Triplet< _t_value > > > &input, ULI m, ULI n, _t_value zero = 0 ) {
		this->load( input, m, n, zero );
	}

	virtual void load( std::vector< Triplet< _t_value > >& input, _i_value m, _i_value n, _t_value zero ) {
		//flat loader
		this->zero_element = zero;
		this->nnz = 1;
		this->nor = m;
		this->noc = n;
		this->n_overflow = n;
		r_start = c_start = 0;
		dss = new _sub_ds*[ 1 ];
		if( static_cast< ULI >( m ) > beta_m || static_cast< ULI >( n ) > beta_n ) {
			std::cerr << "Matrix is too large to fit into a flat FBICRS structure!" << std::endl;
			exit( 1 );
		}
		dss[ 0 ] = new _sub_ds( input, m, n, zero );
		r_inc = new _i_value[ 1 ];
		c_inc = new _i_value[ 1 ];
		c_inc[ 0 ] = this->n_overflow;
	}

	void load( _i_value* row, _i_value* col, _t_value* val, ULI m, ULI n, ULI nz, _t_value zero ) {
		//flat loader
		std::vector< Triplet< _t_value > > packed;
		for( unsigned long int i=0; i<nz; i++ ) {
			packed.push_back( Triplet< _t_value >( row[i], col[i], val[i] ) );
		}
		load( packed, m, n, zero );
	}

	void load( std::vector< std::vector< Triplet< _t_value > > > &input, ULI m, ULI n, _t_value zero ) {
		//full hierarchical loader using fixed blocks
		this->zero_element = zero;
		this->nor = m;
		this->noc = n;
		this->n_overflow = n;
		this->jumps = 0;

#ifndef NDEBUG
		unsigned long int cursum = 0;
		for( unsigned long int i=0; i<input.size(); i++ ) cursum += input[ i ].size();
		const unsigned long int local_nnz = cursum;
#endif

		//reconstruct input while ignoring empty subblocks
		std::vector< ULI > startpos;
		std::vector< std::vector< Triplet< _t_value > > > replace;
		typename std::vector< std::vector< Triplet< _t_value > > >::iterator it = input.begin();
		for( ; it != input.end(); ++it ) {
			if( it->size() > 0 ) {
				//put in CRS order
				qsort( &((*it)[0]), it->size(), sizeof( Triplet< _t_value > ), &(ICRS< _t_value >::compareTriplets) );
				//find lowest column index
				ULI lowest_col = (*it)[0].j(); //initial guess
				//full nonzero traversal
				for( size_t i = 1; i<it->size(); i++ )
					if( (*it)[ i ].j() < lowest_col )
						lowest_col = (*it)[i].j();
				//find lowest row index (note already CRS-ordered)
				const ULI lowest_row = (*it)[0].i();
				//get block start position
				const ULI startpos_row = beta_m * ( lowest_row / beta_m );
				const ULI startpos_col = beta_n * ( lowest_col / beta_n );
				startpos.push_back( startpos_row );
				startpos.push_back( startpos_col );

				//move blocks into a local view
				std::vector< Triplet< _t_value > > toPush;
				//reserve room
				toPush.reserve( it->size() );
				//copy
				for( unsigned long int c = 0; c<it->size(); c++ ) {
					const ULI new_i = (*it)[ c ].i() - startpos_row;
					const ULI new_j = (*it)[ c ].j() - startpos_col;
					assert( new_i < m );
					assert( new_j < n );
					toPush.push_back( Triplet< _t_value >( new_i, new_j, (*it)[ c ].value ) );
				}
				replace.push_back( toPush );
			}
		}

#ifndef NDEBUG
		//check for errors
		for( size_t r = 0, c = 0, q = 0; r < replace.size(); ++r, ++q ) {
			while( input[ q ].size() == 0 ) ++q;
			assert( q < input.size() );
			const ULI cur_block_i = startpos[ c++ ];
			const ULI cur_block_j = startpos[ c++ ];
			assert( replace[ r ].size() == input[ q ].size() );
			for( size_t k = 0; k < replace[ r ].size(); ++k ) {
				const ULI cur_i = cur_block_i + replace[ r ][ k ].i();
				const ULI cur_j = cur_block_j + replace[ r ][ k ].j();
				assert( cur_i == input[ q ][ k ].i() );
				assert( cur_j == input[ q ][ k ].j() );
			}
		}
#endif

		//check for empty matrix
		if( replace.size() == 0 ) {
			r_inc = c_inc = NULL;
			dss = NULL;
			this->nnz = 0;
			return;
		}

		//not empty, continue
		this->nnz = replace.size();

#ifndef NDEBUG
		cursum = 0;
		for( size_t i=0; i<replace.size(); i++ ) cursum += replace[i].size();
		assert( local_nnz == cursum );
#endif

		//initialise fillIn
		fillIn = 0;
		//fill dss
		dss = new _sub_ds*[ this->nnz ];
		for( size_t i = 0; i < static_cast< size_t >( this->nnz ); ++i ) {
			//note this implicitly assumes the input is properly cut s.t.
			//replace[i] indeed fits into a beta x beta matrix.
			ULI msize = beta_m;
			ULI nsize = beta_n;
			//if this is the last block, my size might be smaller
			//(the modulo check makes sure the last block is smaller than beta_m,
			// and not exactly equal to it)
			if( startpos[  2*i  ] == beta_m * ((m-1)/beta_m) && (m%beta_m) != 0 )
				msize = (m%beta_m);
			if( startpos[ 2*i+1 ] == beta_n * ((n-1)/beta_n) && (n%beta_n) != 0 )
				nsize = (n%beta_n);
#ifndef NDEBUG
			const ULI cur_block_i = startpos[  2*i  ];
			const ULI cur_block_j = startpos[ 2*i+1 ];
			for( size_t k = 0; k < replace[ i ].size(); ++k ) {
				const ULI cur_i = cur_block_i + replace[ i ][ k ].i();
				const ULI cur_j = cur_block_j + replace[ i ][ k ].j();
				assert( cur_i < m );
				assert( cur_j < n );
				assert( replace[ i ][ k ].i() < msize );
				assert( replace[ i ][ k ].j() < nsize );
			}
#endif
			dss[ i ] = new _sub_ds( replace[ i ], msize, nsize, zero );

			//update fillIn
			fillIn += dss[ i ]->fillIn;
		}
		
		//count row jumps
		size_t prevrow, prevcol;
		size_t currow, curcol;
		size_t walk = 0;
		prevrow = startpos[ walk++ ];
		prevcol = startpos[ walk++ ];
		for( size_t i = 1; i < replace.size(); i++ ) {
			currow = startpos[ walk++ ];
			curcol = startpos[ walk++ ];
			if( currow != prevrow ) jumps++;
			prevrow = currow;
		}
#ifdef _DEBUG
		std::cout << jumps << " upper-level jumps found." << std::endl;
#endif

		//allocate arrays
		r_inc = new _i_value[ jumps + 1 ];
		c_inc = new _i_value[ this->nnz ];

		//get starting position
		walk = 0;
		r_start = startpos[ walk++ ];
		c_start = startpos[ walk++ ];

		//get ending position
		c_end = startpos[ 2 * this->nnz - 1 ];
		r_end = startpos[ 2 * this->nnz - 2 ];

		//set prev variables		
		prevrow = r_start;
		prevcol = c_start;

		//and build the index arrays
		ULI c = 0;
		for( size_t i = 1; i < static_cast< size_t >(this->nnz); i++ ) {
			currow = startpos[ walk++ ];
			curcol = startpos[ walk++ ];
			assert( currow < static_cast< size_t >(this->nor) );
			assert( curcol < static_cast< size_t >(this->noc) );
			this->c_inc[ i-1 ] =  curcol - prevcol;
			if( currow != prevrow ) {
				this->c_inc[ i-1 ] += n_overflow;
				this->r_inc[ c++ ] = currow - prevrow;	
				assert( this->r_inc[c-1] + prevrow == currow );
				prevrow = currow;
			}
			assert( this->c_inc[i-1] + prevcol == curcol || this->c_inc[i-1] + prevcol - n_overflow == curcol );
			prevcol = curcol;
		}

		//set last c_inc value
		this->c_inc[ this->nnz - 1 ] = n_overflow;

#ifndef NDEBUG
		//check starting positions, first one is given by r_start and c_start
		ULI rowwalk = prevrow = r_start;
		ULI colwalk = c_start;
		walk = 0;
		currow = startpos[ walk++ ];
		curcol = startpos[ walk++ ];
		assert( currow == rowwalk );
		assert( curcol == colwalk );

		//now deduce and check other starting positions
		c = 0;
		for( size_t i = 1; i < static_cast< size_t >(this->nnz); ++i ) {
			colwalk += this->c_inc[ i - 1 ];
			if( colwalk >= n_overflow ) {
				colwalk -= n_overflow;
				rowwalk += this->r_inc[ c++ ];
			}
			assert( rowwalk < static_cast< size_t >(this->nor) );
			assert( colwalk < static_cast< size_t >(this->noc) );
			currow = startpos[ walk++ ];
			curcol = startpos[ walk++ ];
			assert( currow == rowwalk );
			assert( curcol == colwalk );
		}
		assert( c == jumps );

		//now do check in reverse direction
		rowwalk = prevrow = r_end;
		colwalk = c_end;
		walk = 2 * this->nnz - 1;
		curcol = startpos[ walk-- ];
		currow = startpos[ walk-- ];
		assert( currow == rowwalk );
		assert( curcol == colwalk );
		c = jumps - 1;
		for( size_t i = static_cast< size_t >(this->nnz) - 2; i < static_cast< size_t >(this->nnz) - 1; --i ) {
			colwalk -= this->c_inc[ i ];
			if( colwalk >= n_overflow ) {
				colwalk += n_overflow;
				rowwalk -= this->r_inc[ c-- ];
			}
			curcol = startpos[ walk-- ];
			currow = startpos[ walk-- ];
			assert( currow == rowwalk );
			assert( curcol == colwalk );
		}
#endif
	}

	virtual void getFirstIndexPair( _i_value &row, _i_value &col ) {
		row = this->r_start;
		col = this->c_start;
	}

	virtual void zxa( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		if( this->nnz == 0 ) return;

		//set parameters
		_sub_ds**__restrict__ ds_f = this->dss;
		_i_value*__restrict__ c_f = this->c_inc;
		_i_value*__restrict__ r_f = this->r_inc;
		const _t_value*__restrict__ x_f = x_p + this->r_start;
		_t_value*__restrict__ y_f = y_p + this->c_start;

		_t_value* const y_end = y_p + this->noc;
		_sub_ds**__restrict__ const ds_end = ds_f + this->nnz;

		//enter loop
		while( ds_f < ds_end ) {
			//forward direction
			(*ds_f++)->zxa( x_f, y_f );
			y_f += *c_f++;
			if( y_f >= y_end ) {
				y_f -= this->n_overflow;
				x_f += *r_f++;
			}
		}
	}

	virtual void zxa_fb( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		if( this->nnz == 0 ) return;

		//set parameters
		_sub_ds**__restrict__ ds_f = this->dss;
		_sub_ds**__restrict__ ds_b = ds_f + this->nnz - 1;
		_i_value*__restrict__ c_f = this->c_inc;
		_i_value*__restrict__ c_b = this->c_inc + this->nnz - 2;
		_i_value*__restrict__ r_f = this->r_inc;
		_i_value*__restrict__ r_b = this->r_inc + this->jumps - 1;
		const _t_value*__restrict__ x_f = x_p + this->r_start;
		const _t_value*__restrict__ x_b = x_p + this->r_end;
		_t_value*__restrict__ y_f = y_p + this->c_start;
		_t_value*__restrict__ y_b = y_p + this->c_end;

		const _t_value* y_end = y_p + this->noc;

		//enter loop
		while( ds_f <= ds_b ) {
			//forward direction
			(*ds_f++)->zxa( x_f, y_f );
			y_f += *c_f++;
			if( y_f >= y_end ) {
				y_f -= this->n_overflow;
				x_f += *r_f++;
			}

			//check whether we are done
			if( ds_f > ds_b ) break;

			//backward direction
			(*ds_b--)->zxa( x_b, y_b );
			y_b -= *c_b--;
			assert( c_b >= this->c_inc-1 );
			if( y_b < y_p ) {
				y_b += this->n_overflow;
				x_b -= *r_b--;
				assert( y_b >= y_p );
				assert( y_b <  y_end );
			}
		}
	}
	
	virtual void zax( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		if( this->nnz == 0 ) return;

		//set parameters
		_sub_ds**__restrict__ ds_f = this->dss;
		_i_value*__restrict__ c_f = this->c_inc;
		_i_value*__restrict__ r_f = this->r_inc;
		const _t_value*__restrict__ x_f = x_p + this->c_start;
		_t_value*__restrict__ y_f = y_p + this->r_start;

		const _t_value* const x_end = x_p + this->noc;
		_sub_ds**__restrict__ const ds_end = ds_f + this->nnz;

		//enter loop
		while( ds_f < ds_end ) {
			//forward direction
			assert( x_f >= x_p );
			assert( x_f < x_end );
			(*ds_f++)->zax( x_f, y_f );
			x_f += *c_f++;
			if( x_f < x_p || x_f >= x_end ) {
				x_f -= this->n_overflow;
				y_f += *r_f++;
			}
		}
	}

	virtual void zax_fb( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		if( this->nnz == 0 ) return;

		//set parameters
		_sub_ds**__restrict__ ds_f = this->dss;
		_sub_ds**__restrict__ ds_b = ds_f + this->nnz - 1;
		_i_value*__restrict__ c_f = this->c_inc;
		_i_value*__restrict__ c_b = this->c_inc + this->nnz - 2;
		_i_value*__restrict__ r_f = this->r_inc;
		_i_value*__restrict__ r_b = this->r_inc + this->jumps - 1;
		const _t_value*__restrict__ x_f = x_p + this->c_start;
		const _t_value*__restrict__ x_b = x_p + this->c_end;
		_t_value*__restrict__ y_f = y_p + this->r_start;
		_t_value*__restrict__ y_b = y_p + this->r_end;

		const _t_value* x_end = x_p + this->noc;

		//enter loop
		while( ds_f <= ds_b ) {
			//forward direction
			assert( x_f >= x_p );
			assert( x_f < x_end );
			(*ds_f++)->zax( x_f, y_f );
			x_f += *c_f++;
			if( x_f < x_p || x_f >= x_end ) {
				x_f -= this->n_overflow;
				y_f += *r_f++;
			}

			//check whether we are done
			if( ds_f > ds_b ) break;

			//backward direction
			assert( x_b >= x_p );
			assert( x_b < x_end );
			(*ds_b--)->zax( x_b, y_b );
			x_b -= *c_b--;
			assert( c_b >= this->c_inc-1 );
			if( x_b < x_p || x_b >= x_end ) {
				x_b += this->n_overflow;
				y_b -= *r_b--;
				assert( x_b >= x_p );
				assert( x_b <  x_end );
			}
		}
	}

	virtual size_t bytesUsed() {
		//upper-level bytes
		size_t ret = sizeof( ULI ) * 4 + sizeof( _i_value ) * ( this->nnz + jumps + 1 );
		//pointers to sub-level data structures
		ret += sizeof( void* ) * this->nnz;
		//recursively add lower-level data use
		for( size_t i = 0; i < static_cast< size_t >(this->nnz); ++i ) {
			ret += dss[ i ]->bytesUsed();
		}
		//return result
		return ret;
	}
};
#endif
