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

//#define _DEBUG

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
void Matrix2HilbertCoordinates::IntegerToHilbert( const unsigned long int i, const unsigned long int j,
					unsigned long int &h1, unsigned long int &h2 ) {
	//initialise hilbert coordinate to 0
	h1 = h2 = 0;
	//initialise translation mappings
	unsigned long int T1[2][2] = {{0,0},{1,1}};
	const unsigned long int T2[2][2] = {{0,1},{1,0}};
#ifdef _DEBUG
	std::cout << "T2: " << T2[0][0] << " " << T2[0][1] << std::endl;
	std::cout << "    " << T2[1][0] << " " << T2[1][1] << std::endl;
#endif
	//get half bit width for unsigned long ints (usually 32)
	const unsigned char HALFBITWIDTH = BITWIDTH / 2;
	//set which bit to write in h1
	unsigned char hpos = BITWIDTH-1;
	//process half of the bits in i and j to fill up h1
	for( unsigned char k=BITWIDTH; k>HALFBITWIDTH; k-- ) {
		const unsigned long int bitselect = static_cast<unsigned long int>(1)<<(k-1);
		const bool bi = ( i & bitselect ) > 0;
		const bool bj = ( j & bitselect ) > 0;
		const unsigned long int b1 = (T1[bi][bj]<<(hpos--));
		const unsigned long int b2 = (T2[bi][bj]<<(hpos--));
#ifdef _DEBUG
		std::cout << "Bit " << static_cast<unsigned int >(k) << ": " << bi << ", " << bj << std::endl;
		std::cout << "T1: " << T1[0][0] << " " << T1[0][1] << std::endl;
		std::cout << "    " << T1[1][0] << " " << T1[1][1] << std::endl;
		std::cout << "b1: " << b1 << std::endl;
		std::cout << "b2: " << b2 << std::endl;
#endif
		h1 |= b1;
		h1 |= b2;
		if( T1[bi][bj] == T2[bi][bj] ) {
			if( b2 == 0 ) {
				const unsigned long int temp = T1[1][0];
				T1[1][0] = T1[0][1];
				T1[0][1] = temp;
			} else {
				const unsigned long int temp = T1[0][0];
				T1[0][0] = T1[1][1];
				T1[1][1] = temp;
			}
		}
	}
	//reset hpos
	hpos = BITWIDTH - 1;
	//now continue and update h2 instead of h1
	for( unsigned char k=HALFBITWIDTH; k>0; k-- ) {
		const unsigned long int bitselect = static_cast<unsigned long int>(1)<<(k-1);
		const bool bi = ( i & bitselect ) > 0;
		const bool bj = ( j & bitselect ) > 0;
		const unsigned long int b1 = (T1[bi][bj]<<(hpos--));
		const unsigned long int b2 = (T2[bi][bj]<<(hpos--));
#ifdef _DEBUG
		std::cout << "Bit " << static_cast<unsigned int >(k) << ": " << bi << ", " << bj << std::endl;
		std::cout << "T1: " << T1[0][0] << " " << T1[0][1] << std::endl;
		std::cout << "    " << T1[1][0] << " " << T1[1][1] << std::endl;
		std::cout << "b1: " << b1 << std::endl;
		std::cout << "b2: " << b2 << std::endl;
#endif
		h2 |= b1;
		h2 |= b2;
		if( T1[bi][bj] == T2[bi][bj] ) {
			if( b2 == 0 ) {
				const unsigned long int temp = T1[1][0];
				T1[1][0] = T1[0][1];
				T1[0][1] = temp;
			} else {
				const unsigned long int temp = T1[0][0];
				T1[0][0] = T1[1][1];
				T1[1][1] = temp;
			}
		}
	}
}

