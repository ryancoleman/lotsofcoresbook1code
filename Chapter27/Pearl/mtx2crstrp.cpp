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


#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <time.h>

#include "SparseMatrix.hpp"
#include "TS.hpp"

bool compare( const Triplet< double > one, const Triplet< double > two ) {
	if ( one.i() < two.i() )
		return true;
	if ( one.i() > two.i() )
		return false;
	if ( one.j() < two.j() )
		return true;
	return false;
}

int main( int argc, char** argv ) {
	
	if( argc != 3 ) {
		std::cout << "Usage: " << argv[0] << " <input-matrix> <output-trp>" << std::endl << std::endl;
		std::cout << "sorts a matrix market or .trp file into CRS order and writes to output binary tiplets." << std::endl;
		std::cout << "Note: binary triplet format is machine dependent;";
		std::cout << "      take care when using the same binary files on different machine architectures." << std::endl;
		return EXIT_FAILURE;
	}

	std::string file = std::string( argv[1] );
	std::string out  = std::string( argv[2] );

	ULI m, n;
	std::vector< Triplet< double > > matrix = FileToVT::parse( file, m, n );

	std::cout << "Matrix dimensions: " << m << " times " << n << "." << std::endl;
	std::cout << "Sorting..." << std::endl;

	std::sort( matrix.begin(), matrix.end(), compare );

	std::cout << "Saving..." << std::endl;
	Triplet< double >::save( out, &(matrix[0]), m, n, matrix.size() );

	std::cout << "Done!" << std::endl;
	return EXIT_SUCCESS;
}

