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


#ifndef _H_MATRIX2HILBERTCOORDINATES
#define _H_MATRIX2HILBERTCOORDINATES

#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>

/**
 *  Class which maps coordinates to 1D Hilbert Coordinates.
 */
class Matrix2HilbertCoordinates {

   private:

	/** Base constructor. Does nothing. */
	Matrix2HilbertCoordinates() {}
   
   protected:

	static const unsigned char BITWIDTH = 64; //Warning: this should not exceed the bit length
						  //         of an unsigned long int on the system
						  //         used. Understand the algorithm before
						  //         changing this value.

   public:

	/** 
	 *  New method, October 2010. Maps any 2D coordinate (i,j),
	 *  with i and j 64-bits unsigned integers,
	 *  to a 1D 128-bits unsigned integer.
	 *
	 *  @param i  A 64-bits unsigned integer value in one dimension
	 *  @param j  A 64-bits unsigned integer value in the other dimension
	 *  @param h1 First part of the 128-bit Hilbert coordinate, unsigned integer format (most significant, first 64 bits)
	 *  @param h2 Second part of the 128-bit Hilbert coordinate (least significant, last 64 bits)
	 */
	static void IntegerToHilbert( const unsigned long int i, const unsigned long int j,
					unsigned long int &h1, unsigned long int &h2 );

};

#endif

