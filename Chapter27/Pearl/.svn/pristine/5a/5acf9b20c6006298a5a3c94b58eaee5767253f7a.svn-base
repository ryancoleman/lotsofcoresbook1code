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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2007.
 */


#include "Matrix2HilbertCoordinates.hpp"
#include<vector>
#include<string>
#include<fstream>

#ifndef _H_HILBERT_TRIPLET
#define _H_HILBERT_TRIPLET

/**
 *  Hilbert-coordinate-aware triplet. Can save vectors of HilbertTriplets to binary.
 */
template< typename T >
class HilbertTriplet {

   protected:

	/** The row coordinate of this triplet. */
	unsigned long int row;

	/** The column coordinate of this triplet. */
	unsigned long int column;

	/** Most significant part of a 128-bits Hilbert coordinate, for one-shot, non-iterative, calculation. */
	unsigned long int hilbert1;

	/** Least significant part of a 128-bits Hilbert coordinate, for one-shot, non-iterative, calculation. */
	unsigned long int hilbert2;

   public:

	/** @return Row index of this triplet. */	
	unsigned long int i() const { return row; }

	/** @return Column index of this triplet. */
	unsigned long int j() const { return column; }
	
	/** Value stored at this triplet. */
	T value;

	/** Base constructor.
	 *  @param i row index.
	 *  @param j column index.
	 *  @param val nonzero value.
	 */
	HilbertTriplet( unsigned long int i, unsigned long int j, T val ): row( i ), column( j ), hilbert1( 0 ), hilbert2( 0 ), value( val ) {}

	/** Base constructor. Sets all values to zero. */
	HilbertTriplet(): row( 0 ), column( 0 ), hilbert1( 0 ), hilbert2( 0 ), value( 0 ) {}

	/** Calculates the full Hilbert coordinate */
	void calculateHilbertCoordinate() {
		Matrix2HilbertCoordinates::IntegerToHilbert( row, column, hilbert1, hilbert2 );
	}

	/**
	 *	Gets the Hilbert coordinates.
	 *	Does not check if calculateHilbertCoordinate() was called first, otherwise (0,0) will be returned.
	 *	Note that the Hilbert coordinate is a 1D 128-bits unsigned integer.
	 *
	 *   @param h1 The first (most significant) 64-bits of the Hilbert coordinate
	 *   @param h2 the remainder of the Hilbert coordinate.
	 */
	void getHilbertCoordinate( unsigned long int &h1, unsigned long int &h2 ) {
		h1 = hilbert1;
		h2 = hilbert2;
	}

	/** @return h1 of HilbertTriplet<T>::getHilbertCoordinate() */
	unsigned long int getMostSignificantHilbertBits() {
		return hilbert1;
	}

	/** @return h2 of HilbertTriplet<T>::getHilbertCoordinate() */
	unsigned long int getLeastSignificantHilbertBits() {
		return hilbert2;
	}

	/**
	 *  Saves an array of Hilbert triplets to a file, in binary format.
	 *  Does NOT save the Hilbert coordinate!
	 *  (For reading-in the written file, use the regular Triplet scheme)
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Array of Hilbert triplets to write.
	 *  @param m Total number of rows in the matrix.
	 *  @param n Total number of columns in the matrix.
	 *  @param s Size of the array toWrite.
	 */
	static void save( std::string fn, HilbertTriplet< T >* toWrite, const unsigned long int m, const unsigned long int n, const unsigned long int s ) {
		std::fstream myFile ( fn.c_str(), std::ios::out | std::ios::binary);
		myFile.write( (char*) &m, sizeof( unsigned long int ) );
		myFile.write( (char*) &n, sizeof( unsigned long int ) );
		for( unsigned long int i = 0; i<s; i++ ) {
			const unsigned long int wi = toWrite[ i ].i();
			const unsigned long int wj = toWrite[ i ].j();
			const double wv = toWrite[ i ].value;
#ifdef _DEBUG
			std::cout << "Wrote: ( " << wi << " , " << wj << " , " << wv << " ) " << std::endl;
#endif
			myFile.write( (char*) &( wi ), sizeof( unsigned long int ) );
			myFile.write( (char*) &( wj ), sizeof( unsigned long int ) );
			myFile.write( (char*) &( wv ), sizeof( T ) );
		}
		myFile.close();
	}

	/**
	 *  Saves a std::vector of Hilbert triplets to a file, in binary format.
	 *  Does NOT save the Hilbert coordinate!
	 *  (For reading-in the written file, use the regular Triplet scheme)
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Vector of Hilbert triplets to write.
	 *  @param m Total number of rows in the matrix.
	 *  @param n Total number of columns in the matrix.
	 */
	static void save( std::string fn, std::vector< HilbertTriplet< T > > &toWrite, const unsigned long int m, const unsigned long int n ) {
		save( fn, &( toWrite[ 0 ] ), m, n, toWrite.size() );
	}

};

#endif

