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


#include "FileToVT.hpp"

std::vector< Triplet< double > > FileToVT::parse( std::string filename ) {
	ULI m, n;
	unsigned long int nnz;
	return parse( filename, m, n, nnz );
}

std::vector< Triplet< double > > FileToVT::parse( std::string filename, ULI &m, ULI &n ) {
	unsigned long int nnz;
	return parse( filename, m, n, nnz );
}

#ifdef _SUPPORT_CS
std::vector< CS_Triplet< double > > FileToVT::cs_parse( std::string filename ) {
	ULI m, n;
	unsigned long int nnz;
	return cs_parse( filename, m, n, nnz );
}

std::vector< CS_Triplet< double > > FileToVT::cs_parse( std::string filename, ULI &m, ULI &n ) {
	unsigned long int nnz;
	return cs_parse( filename, m, n, nnz );
}
#endif

std::vector< Triplet< double > > FileToVT::parse( std::string filename, ULI &m, ULI &n, unsigned long int &nnz ) {

	MM_typecode type;
	FILE *file;
	int M, N, NNZ;
	double temp_val;
	std::vector< Triplet< double > > ret;

	file = fopen( filename.c_str(), "r" );
        if (file == NULL) {
		std::cerr << "(FileToVT::parse) Could not open file: " << filename << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( mm_read_banner( file, &type ) != 0 ) {
		std::cerr << "Could not process Matrix Market banner (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( !mm_is_matrix( type ) ) {
		std::cerr << "This parser does not handle non-matrix types (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( !mm_is_sparse( type ) ) {
		std::cerr << "This parser will load matrices into a general Triplet format (i,j,v), and thus only applies to sparse matrices. (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( mm_is_complex( type ) ) {
		std::cerr << "Only non-complex matrices can be handled currently (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( mm_read_mtx_crd_size( file, &M, &N, &NNZ ) !=0 ) {
		std::cerr << "An error occured during matrix statistics retrieval (m,n,nnz) (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	} else {
		m = static_cast< ULI >( M );
		n = static_cast< ULI >( N );
		nnz = static_cast< unsigned long int >( NNZ );
	}

	for ( unsigned long int i=0; i<nnz; i++ ) {
		//fscanf( file, "%d %d %lg\n", &M, &N, &temp_val );

		double ignore;
		if( mm_read_mtx_crd_entry( file, &M, &N, &temp_val, &ignore, type ) != 0 )
			std::cerr << "Warning: Failed to read next matrix entry! (FileToVT::parse)" << std::endl;

    		if ( mm_is_pattern( type ) )
			temp_val = 1.0;

#ifdef _DEBUG		
		std::cout << "After read: (" << M << "," << N << "):\t\t" << temp_val << std::endl;
		if( M==8 || N==8 )
			std::cout << "After read: (" << M << "," << N << "):\t\t" << temp_val << std::endl;
#endif

		M--; //1--inf => 0--inf
		N--;
		if( static_cast< ULI >( M ) >= m ) {
			std::cerr << "Row index mismatch (FileToVT::parse)" << std::endl;
			exit( !EXIT_SUCCESS );
		}
		if( static_cast< ULI >( N ) >= n ) {
			std::cerr << "Column index mismatch (FileToVT::parse)" << std::endl;
			exit( !EXIT_SUCCESS );
		}

#ifdef _DEBUG		
		std::cout << "Adding   : (" << M << "," << N << "):\t\t" << temp_val << std::endl;
		if( M == 7 )
	                std::cout << "Adding   : (" << M << "," << N << "):\t\t" << temp_val << std::endl;
#endif
			
		ret.push_back( Triplet< double >( M, N, temp_val ) );
		if( mm_is_symmetric( type ) && ( M != N )) {
#ifdef _DEBUG
			if( N == 7 )
		                std::cout << "Adding   : (" << N << "," << M << "):\t\t" << temp_val << std::endl;
#endif
			ret.push_back( Triplet< double >( N, M, temp_val ) );
		}
	}

	fclose( file );

	return ret;
}

#ifdef _SUPPORT_CS
std::vector< CS_Triplet< double > > FileToVT::cs_parse( std::string filename, ULI &m, ULI &n, unsigned long int &nnz ) {

	MM_typecode type;
	FILE *file;
	int M, N, NNZ;
	double temp_val;
	std::vector< CS_Triplet< double > > ret;

	file = fopen( filename.c_str(), "r" );
        if (file == NULL) {
		std::cerr << "(FileToVT::parse) Could not open file: " << filename << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( mm_read_banner( file, &type ) != 0 ) {
		std::cerr << "Could not process Matrix Market banner (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( !mm_is_matrix( type ) ) {
		std::cerr << "This parser does not handle non-matrix types (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( !mm_is_sparse( type ) ) {
		std::cerr << "This parser will load matrices into a general Triplet format (i,j,v), and thus only applies to sparse matrices. (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( mm_is_complex( type ) ) {
		std::cerr << "Only non-complex matrices can be handled currently (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	}

	if ( mm_read_mtx_crd_size( file, &M, &N, &NNZ ) !=0 ) {
		std::cerr << "An error occured during matrix statistics retrieval (m,n,nnz) (FileToVT::parse)" << std::endl;
		exit( !EXIT_SUCCESS );
	} else {
		m = static_cast< ULI >( M );
		n = static_cast< ULI >( N );
		nnz = static_cast< unsigned long int >( NNZ );
	}

	for ( unsigned long int i=0; i<nnz; i++ ) {
		//fscanf( file, "%d %d %lg\n", &M, &N, &temp_val );

		double ignore;
		//NOTE matlab seems to read rows as columns and vice verca
		if( mm_read_mtx_crd_entry( file, &M, &N, &temp_val, &ignore, type ) != 0 )
			std::cerr << "Warning: Failed to read next matrix entry! (FileToVT::parse)" << std::endl;

#ifdef _DEBUG		
		std::cout << "After read: (" << M << "," << N << "):\t\t" << temp_val << std::endl;
		if( M==8 || N==8 )
			std::cout << "After read: (" << M << "," << N << "):\t\t" << temp_val << std::endl;
#endif

		M--; //1--inf => 0--inf
		N--;
		if( static_cast< ULI >( M ) >= m ) {
			std::cerr << "Row index mismatch (FileToVT::parse)" << std::endl;
			exit( !EXIT_SUCCESS );
		}
		if( static_cast< ULI >( N ) >= n ) {
			std::cerr << "Column index mismatch (FileToVT::parse)" << std::endl;
			exit( !EXIT_SUCCESS );
		}

#ifdef _DEBUG		
		std::cout << "Adding   : (" << M << "," << N << "):\t\t" << temp_val << std::endl;
		if( M == 7 )
	                std::cout << "Adding   : (" << M << "," << N << "):\t\t" << temp_val << std::endl;
#endif
			
		ret.push_back( CS_Triplet< double >( M, N, temp_val ) );
		if( mm_is_symmetric( type ) && ( M != N )) {
#ifdef _DEBUG
			if( N == 7 )
		                std::cout << "Adding   : (" << N << "," << M << "):\t\t" << temp_val << std::endl;
#endif
			ret.push_back( CS_Triplet< double >( N, M, temp_val ) );
		}
	}

	return ret;
}
#endif

