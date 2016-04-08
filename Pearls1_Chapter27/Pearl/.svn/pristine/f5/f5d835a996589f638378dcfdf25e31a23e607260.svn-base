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


#include "HilbertTriplet.hpp"

#ifndef _H_HILBERTTRIPLETCOMPARE
#define _H_HILBERTTRIPLETCOMPARE

/** Class-comparator of HilbertTriplet. */
template< typename T >
class HilbertTripletCompare {

	public:

	bool operator() ( HilbertTriplet< T > i, HilbertTriplet< T > j ) {
		const unsigned long int ihilbert1 = i.getMostSignificantHilbertBits();
		const unsigned long int jhilbert1 = j.getMostSignificantHilbertBits();
		if( ihilbert1 < jhilbert1 ) return true;
		if( ihilbert1 > jhilbert1 ) return false;
		const unsigned long int ihilbert2 = i.getLeastSignificantHilbertBits();
		const unsigned long int jhilbert2 = j.getLeastSignificantHilbertBits();
		if( ihilbert2 < jhilbert2 ) return true;
		return false;
	}
};

#endif

