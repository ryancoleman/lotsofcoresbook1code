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


#include "BlockHilbert.hpp"

#ifndef _H_BISECTIONHILBERT
#define _H_BISECTIONHILBERT

/** The Bisection Hilbert triplet scheme. In effect similar to (HierarchicalBICRS),
    but uses Hilbert coordinates to determine the order of the blocks,
    and a bisection algorithm to construct the individual blocks.
    Wraps around the BisectionHilbert class which already implements this scheme. */
template< typename T >
class BisectionHilbert: public BlockHilbert< T > {

   private:
	
   protected:

   public:

	/** Base deconstructor. */
	virtual ~BisectionHilbert() {}

	/** Base constructor. */
	BisectionHilbert() {
		this->bisection = 1;
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @param bisect Whether bisection-based blocking should be used.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	BisectionHilbert( std::string file, T zero = 0 ) {
		this->bisection = 1;
		this->loadFromFile( file, zero );
	}

	/** 
	 *  Base constructor.
	 *  Warning: the zero parameter is currently NOT USED!
	 *  @param input Raw input of normal triplets.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero What elements is considered to-be zero.
	 */
	BisectionHilbert( std::vector< Triplet< T > >& input, unsigned long int m, unsigned long int n, T zero ) {
		this->bisection = 1;
		this->load( input, m, n, zero );
	}

};

#endif

