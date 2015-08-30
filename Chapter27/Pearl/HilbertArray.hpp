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


#ifndef _H_HILBERTARRAY
#define _H_HILBERTARRAY

#include <vector>
#include <assert.h>

#include "BigInt.hpp"

/**
 * Class providing an interface to an efficient storage of a 1D array of
 * Hilbert coordinate increments. Includes functionality to translate
 * increments into 2D movements applied on two vectors of a template
 * type T. This interface is best viewed as an iterator.
 */
template< typename T >
class HilbertArrayInterface {

    protected:

    public:

	/** Does a zax-operation using the moveToStart and moveToNext functionalities. */
	virtual void zax( const T *__restrict__ const &x, T *__restrict__ const &y, const T *__restrict__ &v, const T *__restrict__ const &v_end ) = 0;

	/** Does a zxa-operation using the moveToStart and moveToNext functionalities. */
	virtual void zxa( const T *__restrict__ const &x, T *__restrict__ const &y, const T *__restrict__ &v, const T *__restrict__ const &v_end ) = 0;

	/** Gets the amount of storage used. */
	virtual size_t bytesUsed() = 0;

	/** Gets the first column index. */
	virtual ULI getFirstColumnIndex() = 0;

	/** Gets the first row index. */
	virtual ULI getFirstRowIndex() = 0;

	/** Moves this interface and the two given vectors to the start location. */
	virtual void moveToStart( const T *__restrict__ &x, T *__restrict__ &y, const T &v ) = 0;

	/** Moves this interface and the two given vectors to the next position. */
	virtual void moveToNext( const T *__restrict__ &x, T *__restrict__ &y, const T &v ) = 0;

	/** Base deconstructor. */
	virtual ~HilbertArrayInterface() {}

};

/**
 * Actual storage implementation 
 *
 * @tparam T  Nonzero data type.
 * @tparam I  Hilbert-coordinate increment type.
 * @tparam hI Full Hilbert-coordinate data type.
 * @tparam mI Full matrix coordinate data type.
 */
template< typename T, typename I, typename hI, typename mI >
class HilbertArray : public HilbertArrayInterface< T > {

    protected:

	/** Hilbert start coordinate. */
	hI start_coor;

	/** Index storage array. */
	I * array;

	/** Current position in the array. */
	I *__restrict__ curpos;

	/** Current Hilbert coordinate. */
	hI curcoor;

	/** Current row position. */
	mI currow;

	/** Current column position. */
	mI curcol;

	/** Keeps track of memory use. */
	size_t bytes;

	/**
	 * Translates a Hilbert coordinate into a row and column index.
	 *
	 *  @param coor Input Hilbert coordinate.
	 *  @param row  Where to decode the row index.
	 *  @param col  Where to decode the column index.
	 */
	void decode( const hI coor, mI &row, mI &col ) {
		//note: the decode must of course be in-tune with the encode
		//      operations in Matrix2HilbertCoordinates.cpp
		
		//max shift
		const unsigned char maxshift = 2 * sizeof( mI ) > sizeof( hI ) ? 4 * sizeof( hI ) : 8 * sizeof( mI ) ;

		//since maxshift might be smaller than the #bits for a mI type, zero them out explicitly
		row = col = 0;

#ifdef _DEBUG
		std::cout << "Decoding " << coor << "..." << std::endl;
#endif

		//set index mask
		mI mask = static_cast< mI >( 1 );
		mask <<= maxshift - 1;

		//set hilbert-index mask to double of matrix index mask
		hI Imask1 = static_cast< hI >( 1ul ) << ( 2 * maxshift - 2 );
		hI Imask2 = static_cast< hI >( 1ul ) << ( 2 * maxshift - 1 );

		//initialise inverse tables
		bool R[2][2] = {{ false, false }, { true, true }};
		bool C[2][2] = {{ false, true }, { true, false }};

		//loop over input bits, from most significant to least significant
		for( ; mask != 0; mask >>= 1, Imask1 >>= 2, Imask2 >>=2 ) {
			//calculate whether the current hilbert coordinate sets the
			//bit at mask
			const bool Imask1_true = (coor & Imask1) > 0;
			const bool Imask2_true = (coor & Imask2) > 0;
			const bool r = R[ Imask2_true ? 1 : 0 ][ Imask1_true ? 1 : 0 ];
			const bool c = C[ Imask2_true ? 1 : 0 ][ Imask1_true ? 1 : 0 ];

#ifdef _DEBUG
			std::cout << (&coor) << "@" << mask << ", " <<  (coor & Imask2) << (coor & Imask1) << " --> " << "(" << r << ", " << c << ")" << ", R={{"
				<< R[0][0] << "," << R[0][1] << "},{" << R[1][0] << "," << R[1][1] << "}}, C={{"
				<< C[0][0] << "," << C[0][1] << "},{" << C[1][0] << "," << C[1][1] << "}}" << std::endl;
#endif

			if( r )
				row |=  mask; //set row index bit to 1
			else
				row &= ~mask; //set row index bit to 0
			if( c )
				col |=  mask; //set column index bit to 1
			else
				col &= ~mask; //set column index bit to 0
			//update inverse tables
			if( Imask1_true == Imask2_true ) {
				if( Imask1_true ) { //permute
					std::swap( R[0][0], R[1][0] );
					std::swap( C[0][0], C[1][0] );
				} else { //swap
					std::swap( R, C );
				}
			} //otherwise no table update necessary
		}

#ifdef _DEBUG
		std::cout << "Decode complete: " << coor << " --> (" << row << ", " << col << ")" << std::endl;
#endif
	}

	/**
	 *  Translates a Hilbert coordinate into a row and column index.
	 *  Only updates the bits necessary to be updated, i.e.,
	 *  skips the first startPos bits of row and col.
	 *  Assumes the bits that will be skipped already have the correct
	 *  value (of course), i.e., this only works with the incremental
	 *  scheme!
	 *
	 *  @param cur  Current input Hilbert coordinate.
	 *  @param inc  Hilbert coordinate increment
	 *  @param row  Where to decode the row index.
	 *  @param col  Where to decode the column index.
	 */
	void decode( hI &cur, const I &inc, mI &row, mI &col ) {
		//note: the decode must of course be in-tune with the encode
		//      operations in Matrix2HilbertCoordinates.cpp

		//remember old value
		const hI old = cur;

		//get new value
		cur += inc;

		//max shift
		const unsigned char maxshift = 2 * sizeof( mI ) > sizeof( hI ) ? 4 * sizeof( hI ) : 8 * sizeof( mI ) ;

		//set index mask
		mI mask = static_cast< mI >( 1 );
		mask <<= maxshift - 1;

		//set hilbert-index mask to double of matrix index mask
		hI Imask1 = static_cast< hI >( 1ul ) << ( 2 * maxshift - 2 );
		hI Imask2 = static_cast< hI >( 1ul ) << ( 2 * maxshift - 1 );

		//initialise inverse tables
		bool R[2][2] = {{ false, false }, { true, true }};
		bool C[2][2] = {{ false, true }, { true, false }};

		//tracks whether changes are cascading
		bool cascade = false;

#ifdef _DEBUG
		std::cout << "Partial decode: " << cur << std::endl;
#endif

		//loop over input bits, from most significant to least significant
		for( ; mask != 0; mask >>= 1, Imask1 >>= 2, Imask2 >>=2 ) {
			//only update when old bits are out of date
			if( cascade || (cur & Imask1) != (old & Imask1) || (cur & Imask2) != (old & Imask2) ) {
				//calculate whether the current hilbert
				//coordinate sets the bit at mask
				const bool Imask1_true = (cur & Imask1) > 0;
				const bool Imask2_true = (cur & Imask2) > 0;
				const bool r = R[ Imask2_true ? 1 : 0 ][ Imask1_true ? 1 : 0 ];
				const bool c = C[ Imask2_true ? 1 : 0 ][ Imask1_true ? 1 : 0 ];
#ifdef _DEBUG
				std::cout << (cur & Imask2) << (cur & Imask1) << " --> " << "(" << r << ", " << c << ")" << ", R={{"
						<< R[0][0] << "," << R[0][1] << "},{" << R[1][0] << "," << R[1][1] << "}}, C={{"
						<< C[0][0] << "," << C[0][1] << "},{" << C[1][0] << "," << C[1][1] << "}}" << std::endl;
#endif
				if( r )
					row |=  mask; //set row index bit to 1
				else
					row &= ~mask; //set row index bit to 0
				if( c )
					col |=  mask; //set column index bit to 1
				else
					col &= ~mask; //set column index bit to 0
				//are we now at entry or exit points?
				if( Imask1_true == Imask2_true ) {
					//no matter where old was, it was definitely not at the same exit
					//point, so we are sure the Hilbert orientation has changed.
					//Therefore, changes will cascade down through the remaining bits.
					cascade = true; //the Hilbert orientation has changed in recursion!
					//update tables
					if( Imask1_true ) { //permute
						std::swap( R[0][0], R[1][0] );
						std::swap( C[0][0], C[1][0] );
					} else { //swap
						std::swap( R, C );
					}
				} else {
					//if old was at an entry or exit point, then again the curve
					//orientation has changed (with respect to old),and things
					//chancges can cascade.
					if( ((old & Imask1) > 0) == ((old & Imask2) > 0) )
						cascade = true;
				}
			} else { //even if we do not need to update row/col bits, we need to track the table
				if( ((cur & Imask1) > 0) == ((cur & Imask2) > 0) ) {
					if( (cur & Imask1) ) { //permute
						std::swap( R[0][0], R[1][0] );
						std::swap( C[0][0], C[1][0] );
					} else { //swap
						std::swap( R, C );
					}
				}
			}
		}
	}


    public:

	/** Base constructor. */
	HilbertArray( const std::vector< unsigned long int > &input ) {
		//set memory use to minimum case
		bytes = sizeof( hI ) * 2 + sizeof( mI ) * 2 + sizeof( I* );

		//check for empty input
		if( input.size() == 0 ) {
			start_coor = 0;
			curpos = NULL;
			array = NULL;
		} else {
			//store increment array in hilbert-Index (I) form
			array = new I[ input.size() - 1 ];
			//update memory use
			bytes += sizeof( I ) * ( input.size() - 1 );
			//store start location in matrix-Index form
			start_coor = input[ 0 ];
			//fill the array from input
			for( size_t i = 0; i < input.size() - 1; ++i )
				array[ i ] = static_cast< I >( input[ i + 1 ] );
			//initialise current position
			curpos = NULL;
		}
	}

	/**
	 * Gets the start-location of the first nonzero in this array.
	 *
	 * @return The row-position of this location.
	 */
	virtual ULI getFirstRowIndex() {
		decode( start_coor, currow, curcol );
		return static_cast< ULI >( currow );
	}

	/**
	 * Gets the start-location of the first nonzero in this array.
	 *
	 * @return The column-position of this location.
	 */
	virtual ULI getFirstColumnIndex() {
		decode( start_coor, currow, curcol );
		return static_cast< ULI >( curcol );
	}

	/**
	 * Flat implementation of the zax. Might perform better than using the moveToStart and
	 * moveToNext functions from an outside class.
	 *
	 * @param x Pointer to the input  vector.
	 * @param y Pointer to the output vector.
	 * @param v Pointer to the nonzero values array. Warning: the pointer position will be altered!
	 * @param v_end End-location of the nonzero values array (equal to v+nnz).
	 */
	virtual void zxa( const T *__restrict__ const &x, T *__restrict__ const &y, const T *__restrict__ &v, const T *__restrict__ const &v_end ) {
		//check for empty matrix
		if( v >= v_end )
			return;

		//translate first hilbert coordinate to 2D position
		decode( start_coor, currow, curcol );

		//prepare for iteration
		curpos  = array;
		curcoor = start_coor;
		
		//do initial muladd
		y[ curcol ] += *v++ * x[ currow ];

		//loop over all nonzeroes
		while( v < v_end ) {
			//update current hilbert coordinate and current coordinates
			decode( curcoor, *curpos++, currow, curcol );
			//do muladd
			y[ curcol ] += *v++ * x[ currow ];
		}
	}

	/**
	 * Flat implementation of the zax. Might perform better than using the moveToStart and
	 * moveToNext functions from an outside class.
	 *
	 * @param x Pointer to the input  vector.
	 * @param y Pointer to the output vector.
	 * @param v Pointer to the nonzero values array. Warning: the pointer position will be altered!
	 * @param v_end End-location of the nonzero values array (equal to v+nnz).
	 */
	virtual void zax( const T *__restrict__ const &x, T *__restrict__ const &y, const T *__restrict__ &v, const T *__restrict__ const &v_end ) {
		//check for empty matrix
		if( v >= v_end )
			return;

		//translate first hilbert coordinate to 2D position
#ifndef NDEBUG
		currow = curcol = 0; //otherwise valgrind may think currow and curcol are uninitialised after the below decode call.
#endif
		decode( start_coor, currow, curcol );
		
#ifndef NDEBUG
		mI checkR, checkC;
#endif		

		//prepare for iteration
		curpos  = array;
		curcoor = start_coor;
		
		//do initial muladd
		y[ currow ] += *v++ * x[ curcol ];

		//loop over all nonzeroes
		while( v != v_end ) {
			//update current hilbert coordinate and current coordinates
			decode( curcoor, *curpos++, currow, curcol );
#ifndef NDEBUG
			checkR = checkC = 0;
			decode( curcoor, checkR, checkC );
			assert( currow == checkR );
			assert( curcol == checkC );
			unsigned long int check1, check2;
			Matrix2HilbertCoordinates::IntegerToHilbert( currow, curcol, check1, check2 );
			assert( check2 == curcoor );
#endif

			//do muladd
			y[ currow ] += *v++ * x[ curcol ];
		}
	}

	/** @return The number of bytes used by this Hilbert array */
	virtual size_t bytesUsed() {
		return bytes;
	}

	/**
	 * Moves this interface to the start location and perform a
	 * single multiply-add there.
	 *
	 * @param x Points to the start of the input  vector.
	 * @param y Points to the start of the output vector.
	 * @param v Nonzero value to use.
	 */
	virtual void moveToStart( const T *__restrict__ &x, T *__restrict__ &y, const T &v ) {
		//translate first hilbert coordinate to 2D position
		decode( start_coor, currow, curcol );
		//do muladd
		y[ currow ] += v * x[ curcol ];
		//prepare for iteration
		curpos  = array;
		curcoor = start_coor;
	}

	/**
	 * Moves this interface and the two given vectors to the next position.
	 *
	 * @param x Points to the start of the input  vector.
	 * @param y Points to the start of the output vector.
	 * @param v Non-zero value to use.
	 */
	virtual void moveToNext( const T *__restrict__ &x, T *__restrict__ &y, const T &v ) {
		//update current hilbert coordinate and current coordinates
		decode( curcoor, *curpos++, currow, curcol );
		//do muladd
		y[ currow ] += v * x[ curcol ];
	}

	/** Base deconstructor. */
	virtual ~HilbertArray() {
		delete [] array;
	}

};

/**
 * @param  x An unsigned integer larger than 0.
 * @return floor(log_2(x))
 */
unsigned char binlog( unsigned long int x ) {
	assert( x != 0 );
	unsigned char ret = 0;
	bool thereIsARemainder = x > 1 ? x & 1 : false;
	while( (x >>= 1) > 0 ) {
		if( !thereIsARemainder && x > 1 )
			thereIsARemainder = x & 1;
		ret++;
	}
	return ret + ( thereIsARemainder ? 1 : 0 );
}

/**
 * @param x A 128-bit unsigned integer larger than 0.
 * @return floor(log_2(x))
 */
unsigned char binlog( BigInt x ) {
	assert( !(x.high == 0 && x.low == 0) );
	unsigned char ret = 0;
	bool thereIsARemainder = x.low > 1 ? x.low & 1 : false;
	if( !thereIsARemainder && x.high > 1 )
		thereIsARemainder = x.high & 1;
	while( (x.low >>= 1 ) > 0 ) {
		if( !thereIsARemainder && x.low > 1 )
			thereIsARemainder = x.low & 1;
		ret++;
	}
	while( (x.high >>= 1 ) > 0 ) {
		if( !thereIsARemainder && x.high > 1 )
			thereIsARemainder = x.high & 1;
		ret++;
	}
	return ret + ( thereIsARemainder ? 1 : 0 );
}

/** Factory method; gets an auto-tuned incremental Hilbert-coordinate array. */
template< typename T >
HilbertArrayInterface< T >* getDiffArray( const std::vector< BigInt > &values_in, const unsigned long int m, const unsigned long int n ) {
	//handle boundary case
	if( values_in.size() == 0 )
		return NULL;

	//get (log2 o floor) of max matrix dimension
	const unsigned char mSize = binlog( m ) > binlog( n ) ? binlog( m ) : binlog( n );

	//get (log2 o floor) of max Hilbert coordinate
	std::vector< BigInt >::const_iterator it = values_in.begin();
	BigInt maxH = *it;
	for( ++it; it != values_in.end(); ++it ) {
		if( maxH < *it )
			maxH = *it;
	}
	const unsigned char logMaxH = binlog( maxH );

	//make values array incremental (except first value)
	//simultaneously get (log2 o floor) of max Hilbert increment
	it = values_in.begin() + 1;
	unsigned long int curhigh = values_in[ 0 ].high;
	unsigned long int curlow  = values_in[ 0 ].low;
	unsigned long int maxDiff = it->low;
	std::vector< unsigned long int > values;
	values.push_back( curlow );
	for( ; it != values_in.end(); ++it ) {
		const unsigned long int diffhigh = it->high - curhigh;
		const unsigned long int difflow  = it->low  - curlow; 
		curhigh = it->high;
		curlow  = it->low;
		if( diffhigh > 0 ) {
			std::cerr << "Current implementation assumes increments of 64 bits or less, but a larger increment was found!" << std::endl;
			exit( EXIT_FAILURE );
		}
		values.push_back( difflow );
		if( maxDiff < difflow )
			maxDiff = difflow;
	}
	const unsigned char logMaxDiff = binlog( maxDiff );

	std::cout << "Matrix index bitlength:       " << static_cast< unsigned short int >( mSize ) << " (" << (m>n?m:n) << ")" << std::endl;
	std::cout << "Hilbert coordinate bitlength: " << static_cast< unsigned short int >( logMaxH ) << " (" << maxH.high << maxH.low << ")" << std::endl;
	std::cout << "Hilbert increment bitlength:  " << static_cast< unsigned short int >( logMaxDiff ) << " (" << maxDiff << ")" << std::endl;

	//construct and return tuned version
	if( mSize <= 8 ) { //mI is unsigned char
		if( logMaxDiff <= 8 ) { //I is unsigned char
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned char, unsigned char,      unsigned char >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned char, unsigned short int, unsigned char >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned char, unsigned int,       unsigned char >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned char, unsigned long int,  unsigned char >( values );
			}
		} else if( logMaxDiff <= 16 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned short int, unsigned char,      unsigned char >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned short int, unsigned short int, unsigned char >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned short int, unsigned int,       unsigned char >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned short int, unsigned long int,  unsigned char >( values );
			}
		} else if( logMaxDiff <= 32 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned int, unsigned char,      unsigned char >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned int, unsigned short int, unsigned char >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned int, unsigned int,       unsigned char >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned int, unsigned long int,  unsigned char >( values );
			}
		} else if( logMaxDiff <= 64 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned long int, unsigned char,      unsigned char >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned long int, unsigned short int, unsigned char >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned long int, unsigned int,       unsigned char >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned long int, unsigned long int,  unsigned char >( values );
			}
		}
	} else if( mSize <= 16 ) { //mI is unsigned short int
		if( logMaxDiff <= 8 ) { //I is unsigned char
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned char, unsigned char,      unsigned short int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned char, unsigned short int, unsigned short int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned char, unsigned int,       unsigned short int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned char, unsigned long int,  unsigned short int >( values );
			}
		} else if( logMaxDiff <= 16 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned short int, unsigned char,      unsigned short int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned short int, unsigned short int, unsigned short int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned short int, unsigned int,       unsigned short int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned short int, unsigned long int,  unsigned short int >( values );
			}
		} else if( logMaxDiff <= 32 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned int, unsigned char,      unsigned short int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned int, unsigned short int, unsigned short int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned int, unsigned int,       unsigned short int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned int, unsigned long int,  unsigned short int >( values );
			}
		} else if( logMaxDiff <= 64 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned long int, unsigned char,      unsigned short int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned long int, unsigned short int, unsigned short int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned long int, unsigned int,       unsigned short int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned long int, unsigned long int,  unsigned short int >( values );
			}
		}
	} else if( mSize <= 32 ) { //mI is unsigned int
		if( logMaxDiff <= 8 ) { //I is unsigned char
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned char, unsigned char,      unsigned int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned char, unsigned short int, unsigned int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned char, unsigned int,       unsigned int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned char, unsigned long int,  unsigned int >( values );
			}
		} else if( logMaxDiff <= 16 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned short int, unsigned char,      unsigned int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned short int, unsigned short int, unsigned int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned short int, unsigned int,       unsigned int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned short int, unsigned long int,  unsigned int >( values );
			}
		} else if( logMaxDiff <= 32 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned int, unsigned char,      unsigned int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned int, unsigned short int, unsigned int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned int, unsigned int,       unsigned int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned int, unsigned long int,  unsigned int >( values );
			}
		} else if( logMaxDiff <= 64 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned long int, unsigned char,      unsigned int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned long int, unsigned short int, unsigned int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned long int, unsigned int,       unsigned int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned long int, unsigned long int,  unsigned int >( values );
			}
		}
	} else if( mSize <= 64 ) { //mI is unsigned long int
		if( logMaxDiff <= 8 ) { //I is unsigned char
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned char, unsigned char,      unsigned long int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned char, unsigned short int, unsigned long int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned char, unsigned int,       unsigned long int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned char, unsigned long int,  unsigned long int >( values );
			}
		} else if( logMaxDiff <= 16 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned short int, unsigned char,      unsigned long int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned short int, unsigned short int, unsigned long int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned short int, unsigned int,       unsigned long int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned short int, unsigned long int,  unsigned long int >( values );
			}
		} else if( logMaxDiff <= 32 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned int, unsigned char,      unsigned long int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned int, unsigned short int, unsigned long int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned int, unsigned int,       unsigned long int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned int, unsigned long int,  unsigned long int >( values );
			}
		} else if( logMaxDiff <= 64 ) {
			if( logMaxH <= 8 ) { //hI is unsigned char
				return new HilbertArray< T, unsigned long int, unsigned char,      unsigned long int >( values );
			} else if( logMaxH <= 16 ) {
				return new HilbertArray< T, unsigned long int, unsigned short int, unsigned long int >( values );
			} else if( logMaxH <= 32 ) {
				return new HilbertArray< T, unsigned long int, unsigned int,       unsigned long int >( values );
			} else if( logMaxH <= 64 ) {
				return new HilbertArray< T, unsigned long int, unsigned long int,  unsigned long int >( values );
			}
		}
	} else {
		std::cerr << "Cannot handle matrix indices requiring " << mSize << " bits!" << std::endl;
		exit( EXIT_FAILURE );
	}
	return NULL;
}

#endif

