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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2008.
 */


#ifndef _H_FTVT
#define _H_FTVT

#include <cstdlib>
#include <iostream>
#include <vector>

#include "mmio.h"
#include "Triplet.hpp"

// Flag indicating if we support cache-simulated triplets, as defined in the CACHE_SIM library.
#ifdef _SUPPORT_CS
	#include "CS_Triplet.hpp"
#endif

/** Class responsible for reading in matrix market files and converting them to vector< Triplet > format. */
class FileToVT {

   public:

	/** Parses a matrix-market input file */
	static std::vector< Triplet< double > > parse( std::string filename );
	/** Parses a matrix-market input file */
	static std::vector< Triplet< double > > parse( std::string filename, ULI &m, ULI &n );
	/** Parses a matrix-market input file */
	static std::vector< Triplet< double > > parse( std::string filename, ULI &m, ULI &n, unsigned long int &nnz );

#ifdef _SUPPORT_CS
	/** Parses a matrix-market input file */
	static std::vector< CS_Triplet< double > > cs_parse( std::string filename );
	/** Parses a matrix-market input file */
	static std::vector< CS_Triplet< double > > cs_parse( std::string filename, ULI &m, ULI &n );
	/** Parses a matrix-market input file */
	static std::vector< CS_Triplet< double > > cs_parse( std::string filename, ULI &m, ULI &n, unsigned long int &nnz );
#endif

};

#endif //_H_FTVT

