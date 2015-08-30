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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2010.
 */


#include "Matrix.hpp"
#include "SparseMatrix.hpp"
#include "TS.hpp"
#include "CRS.hpp"
#include "ICRS.hpp"
#include "ZZ_CRS.hpp"
#include "ZZ_ICRS.hpp"
#include "SVM.hpp"
#include "HTS.hpp"
#include "BICRS.hpp"
#include "CCSWrapper.hpp"


#include <assert.h>
#include <iostream>

// Set when in debug mode
//#define _DEBUG

#ifndef _H_HBICRS
#define _H_HBICRS

/**
 *  Hierarchical Bi-directional Incremental Compressed Row Storage scheme.
 *  Stores other SparseMatrix data structures, can include other HBICRS schemes,
 *  but this is not recommended with regards to efficiency.
 *  Supports jumping back and forward within columns.
 *  Supports jumping back and forward between rows.
 *  Main storage direction in column-wise.
 *  Storage requirements are 2nz plus the number of row jumps required, plus the
 *  storage requirements for each stored data structure, of course.
 *  Many row jumps are disadvantageous to storage as well as speed.
 *  @param _t_value The type of the nonzeros in the matrix.
 *
 *  This class is based on the BICRS as per revision 75.
 *  BICRS was chosen since it is an all-round improvement over Triplets, with
 *  the additional advantage of being pointer based (no jumps required when
 *  going recursive into deeper storage schemes).
 *
 *  Warning: this class uses assertions! For optimal performance,
 *           define the NDEBUG flag (e.g., pass -DNDEBUG as a compiler
 *	     flag).
 */
template< typename _t_value >
class HBICRS: public SparseMatrix< _t_value, LI > {

     protected:

	/** The number of row jumps. */
	size_t jumps;

	/** Stores the row jumps; size is _at maximum_ the number of nonzeros. */
	LI* r_inc;

	/** Stores the column jumps; size is exactly the number of nonzeros. */
	LI* c_inc;

	/** Stores the values individual storage schemes. */
	Matrix< _t_value > **vals;

	/** Caches n times two. */
	ULI ntt;
	
     public:

	/** Base deconstructor. */
	~HBICRS() {
		delete [] r_inc;
		delete [] c_inc;
		for( LI i=0; i<this->nnz; i++ )
			delete vals[ i ];
		delete [] vals;
	}

	/** Base constructor. */
	HBICRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	HBICRS( std::string file, _t_value zero = 0 ) {
		std::cerr << "Error: Hierarchical BICRS, loading in from MatrixMarket file not possible, since no way to infer a hierarchy!" << std::endl;
		exit( 1 );
	}

	/** Base constructor.
	 *  @see SparseMatrix::SparseMatrix( input, m, n, zero )
	 */
	HBICRS( std::vector< Triplet< _t_value > >& input, LI m, LI n, _t_value zero ) {
		load( input, m, n, zero );
	}

	/** Base constructor.
	 *  @see SparseMatrix::SparseMatrix( input, m, n, zero )
	 */
	HBICRS( std::vector< std::vector< Triplet< _t_value > > >& input, signed char* group_type, LI m, LI n, _t_value zero = 0 ) {
		load( input, group_type, m, n, zero );
	}

	/**
	 *  This function will rewrite the std::vector< Triplet > structure to one suitable
	 *  for the other load function.
	 *  @see load( row, col, val, m, n, nz ) 
	 *  @see SparseMatrix::load 
	 */
	virtual void load( std::vector< Triplet< _t_value > >& input, LI m, LI n, _t_value zero ) {
		std::cerr << "Error: Hierarchical BICRS, loading in from plain triplets not possible, since no way to infer a hierarchy!" << std::endl;
		exit( 1 );
	}

	/** Constructs the hierarchical part. Note that this function does not do anything with offsets, due to lack of support in the other datastructures used. */
	virtual void load( std::vector< std::vector< Triplet< _t_value > > >& input, signed char* group_type, LI m, LI n, _t_value zero ) {
		LI* rrow = new LI[ input.size() ];
		LI* ccol = new LI[ input.size() ];

		if( group_type == NULL ) {
			char *group_type = new char[ input.size() ];
			for( unsigned int i=0; i<input.size(); i++ ) group_type[i] = 2;
		}

		vals = new Matrix< _t_value >*[ input.size() ];

		unsigned long int g = 0;
		typename std::vector< std::vector< Triplet< _t_value > > >::iterator it = input.begin();
		for( ; it!=input.end(); it++, g++ ) {
			//skip empty blocks
			if( it->size() == 0 ) {
				g--;
				continue;
			}
			//std::cout << "\tHBICRS: Adding a group type " << static_cast< int >( group_type[g] ) << std::endl;
			switch( group_type[ g ] ) {
			case -7:
				vals[g] = new CCSWrapper< double, BICRS< _t_value >, ULI >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case -6:
				vals[g] = new CCSWrapper< double, HTS< _t_value >, ULI >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case -5:
				vals[g] = new CCSWrapper< double, SVM< _t_value >, ULI >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case -4:
				vals[g] = new CCSWrapper< _t_value, ZZ_ICRS< _t_value >, LI >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case -3:
				vals[g] = new CCSWrapper< double, ZZ_CRS< _t_value >, LI >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case -2:
				vals[g] = new CCSWrapper< double, ICRS< _t_value >, ULI >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case -1:
				vals[g] = new CCSWrapper< double, CRS< _t_value >, ULI >( *it, m, n, zero );
				rrow[g] = 0;	
				ccol[g] = 0;
				break;
			case 0:
				vals[g] = new TS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case 1:
				vals[g] = new CRS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;	
				ccol[g] = 0;
				break;
			case 2:
				vals[g] = new ICRS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case 3:
				vals[g] = new ZZ_CRS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case 4:
				vals[g] = new ZZ_ICRS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case 5:
				vals[g] = new SVM< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case 6:
				vals[g] = new HTS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			case 7:
				vals[g] = new BICRS< _t_value >( *it, m, n, zero );
				rrow[g] = 0;
				ccol[g] = 0;
				break;
			default:
				std::cerr << "Hierarchical sparse matrix: invalid subscheme ID (" << static_cast<int>(group_type[g]) << "), group " << g << "ignored!" << std::endl;
			}

		}
		this->zero_element = zero;
		load( rrow, ccol, m, n, g );
		delete [] rrow;
		delete [] ccol;
	}

	/** Builds the BICRS structure. */
	void load( LI* row, LI* col, LI m, LI n, LI nb ) {
#ifdef _DEBUG
		std::cerr << "Warning: _DEBUG flag set." << std::endl;
#endif

		this->nnz = nb;
		this->nor = m;
		this->noc = n;
		this->ntt=2*n;
		jumps = 0;
		int prevrow = row[ 0 ];
		for( LI i=1; i<this->nnz; i++ ) {
			if( row[ i ] != prevrow )
				jumps++;
			prevrow = row[ i ];
		}
#ifdef _DEBUG
		std::cout << jumps << " row jumps found." << std::endl;
#endif
		r_inc = new LI[ jumps + 2 ];
		c_inc = new LI[ this->nnz + 1 ];
	
		r_inc[ 0 ] = prevrow = row[ 0 ];	
		int prevcol = c_inc[ 0 ] = col[ 0 ];
#ifdef _DEBUG
		std::cout << "c_inc: " << prevcol << std::endl;
		std::cout << "r_inc: " << prevrow << std::endl;
#endif
		int c = 1;
		for( LI i=1; i<this->nnz; i++ ) {
			this->c_inc[ i ] = col[ i ] - prevcol;
			if( row[ i ] != prevrow ) {
				this->c_inc[ i ] += ntt;
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
		c_inc[ this->nnz ] = ntt;
		//initialise last row jump to zero (prevent undefined jump)
		r_inc[ c ] = 0;
#ifdef _DEBUG
		std::cout << "Construction done." << std::endl;
#endif
	}

	/** @see SparseMatrix::getFirstIndexPair */
	virtual void getFirstIndexPair( LI &i, LI &j ) {
		i = r_inc[ 0 ];
		j = c_inc[ 0 ];
	}

	/**
	 *  Calculates y=xA, but does not allocate y itself.
	 *  @param x The input vector should be initialised and of correct measurements.
	 *  @param y The output vector should be preallocated and of size m. Furthermore, y[i]=0 for all i, 0<=i<m.
	 */
	virtual void zxa( const _t_value*__restrict__ x_p, _t_value*__restrict__ y_p ) {
		const _t_value * y		= y_p;
		const _t_value * y_end		= y+this->noc;
		LI *__restrict__ c_inc_p	= c_inc;
		LI *__restrict__ r_inc_p	= r_inc;

#ifndef NDEBUG
		const _t_value * x		= x_p;
		const _t_value *__restrict__ x_end      = x+this->nor;
		const LI *__restrict__ c_inc_end       = c_inc+this->nnz+1;
#endif

		Matrix< _t_value >**__restrict__ v_end	= vals+this->nnz;
		Matrix< _t_value >**__restrict__ v_p	= vals;
#ifdef _DEBUG
		//simulation trace:
		int yc = c_inc[0];
		int xc = r_inc[0];
		std::cout << "( " << xc << " , " << yc << " )" << std::endl;
		int t=1;
		for( int i=1; i<this->nnz; i++ ) {
			yc += c_inc[i];
			if( yc > this->noc ) {
				yc -= ntt;
				xc += r_inc[t++];
			}
			assert( yc < this->noc );
			assert( xc < this->nor );
			assert( yc >=0 );
			assert( xc >=0 );
			std::cout << "( " << yc << " , " << xc << " )" << std::endl;
		}
#endif
		x_p += *r_inc_p++;
		y_p += *c_inc_p++;
		while( v_p < v_end ) {
			assert( y_p >= y );
			assert( y_p <  y_end );
			assert( v_p >= vals );
			assert( v_p <  v_end );
			assert( x_p >= x );
			assert( x_p <  x_end );
			assert( c_inc_p >= c_inc );
			assert( c_inc_p <  c_inc_end );
			assert( r_inc_p >= r_inc );
			while( y_p < y_end ) {
				(*(v_p++))->zxa( x_p, y_p );
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
		const _t_value * x_end		= x_p+this->noc;
		LI *__restrict__ c_inc_p	= c_inc;
		LI *__restrict__ r_inc_p	= r_inc;
#ifndef NDEBUG
		const _t_value * x			= x_p;
		const _t_value * const y		= y_p;
		const _t_value * y_end			= y+this->nor;
		const LI *__restrict__ c_inc_end	= c_inc+this->nnz+1;
#endif

		Matrix< _t_value >**__restrict__ v_end	= vals+this->nnz;
		Matrix< _t_value >**__restrict__ v_p	= vals;

#ifdef _DEBUG
		//simulation trace:
		int xc = c_inc[0];
		int yc = r_inc[0];
		std::cout << "( " << xc << " , " << yc << " )" << std::endl;
		int t=1;
		for( int i=1; i<this->nnz; i++ ) {
			xc += c_inc[i];
			if( xc > this->noc ) {
				xc -= ntt;
				yc += r_inc[t++];
			}
			assert( xc < this->nor );
			assert( yc < this->noc );
			assert( xc >=0 );
			assert( yc >=0 );
			std::cout << "( " << xc << " , " << yc << " )" << std::endl;
		}
#endif

		x_p += *c_inc_p++;
		y_p += *r_inc_p++;
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
				(*(v_p++))->zax( x_p, y_p );
				x_p += *c_inc_p++;
			}
			x_p -= ntt;
			y_p += *r_inc_p++;
		}
	}

	virtual size_t bytesUsed() {
		//base case; indices + pointers to substructures
		size_t ret = sizeof( LI ) * ( this->nnz + jumps + 2 ) + sizeof( void* ) * this->nnz;
		//get leaf data-structure memory use
		for( size_t i = 0; i < static_cast< size_t >(this->nnz); ++i )
			ret += vals[ i ]->bytesUsed();
		return ret;
	}

};

#endif

