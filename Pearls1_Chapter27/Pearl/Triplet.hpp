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


#include<cstdlib>
#include<string>
#include<vector>
#include<fstream>
#include<iostream>

//#define _DEBUG

#include "ULIdef.hpp"

#ifndef _H_TRIPLET
#define _H_TRIPLET

/** A single triplet value. Stores, at minimum, a row
 *  position, a column position, and a value of the
 *  template type.
 */
template< typename T >
class Triplet {

   protected:
	
	/** The row coordinate of this triplet. */
	ULI row;
	
	/** The column coordinate of this triplet. */
	ULI column;

   public:

#ifdef TRIPLET_META
	/** Stores meta data */
	TRIPLET_META_TYPE meta;
#endif

	/** @return Row index of this triplet. */	
	ULI i() const { return row; }

	/** @return Column index of this triplet. */
	ULI j() const { return column; }

	/** Subtracts a given value from this nonzero row index. */
	void rowOffset( ULI offset ) { row -= offset; }

	/** Set the coordinates of this nonzero. */
	void setPosition( const unsigned long int i, const unsigned long int j ) {
		row = i;
		column = j;
	}

	/** Set the row coordinate of this nonzero. */
	void setRowPosition( const unsigned long int i ) {
		setPosition( i, column );
	}

	/** Set the column coordinate of this nonzero. */
	void setColumnPosition( const unsigned long int j ) {
		setPosition( row, j );
	}

	/** Value stored at this triplet. */
	T value;

	/** Base constructor.
	 *  @param i row index.
	 *  @param j column index.
	 *  @param val nonzero value.
	 */
	Triplet( ULI ii, ULI ij, T val ): row( ii ), column( ij ), value( val ) {}

	/** Base constructor. Sets all values to zero. */
	Triplet(): row( 0 ), column( 0 ), value( 0 ) {}

	/** Copy constructor. */
	Triplet( const Triplet< T > &toCopy ): row( toCopy.i() ), column( toCopy.j() ), value( toCopy.value ) {
#ifdef TRIPLET_META
		meta = toCopy.meta;
#endif
	}

	/** Loads an array of triplets from a binary file.
	 *  Warning: there is a difference between 32 and 64 bits files!
	 *  @param fn Filename of the file to load from.
	 *  @param m Reference to where the total number of rows is to-be stored.
	 *  @param n Reference to where the total number of columns is to-be stored.
	 *  @return A std::vector containing the triplets in the order they were loaded in.
	 */
	static std::vector< Triplet< T > > load( const std::string fn, ULI &m, ULI &n ) {
		std::fstream file( fn.c_str(), std::ios::in | std::ios::binary );
		ULI i; ULI j; double v;
		std::vector< Triplet< T > > ret;
		if( !file.is_open() ) {
			std::cerr << "Error while opening file" << std::endl;
			exit( 1 );
		}
		file.read( (char*) &i, sizeof( ULI ) );
		file.read( (char*) &j, sizeof( ULI ) );
		m = i;
		n = j;
#ifdef _DEBUG
		std::cout << "m: " << m << ", n: " << n << std::endl;
#endif
		while( true ) {
			file.read( (char*) &i, sizeof( ULI ) );
			if( !file ) break;
			file.read( (char*) &j, sizeof( ULI ) );
			file.read( (char*) &v, sizeof( T ) );
#ifdef _DEBUG
			std::cout << "Pushed back: ( " << i << " , " << j << " , " << v << " )" << std::endl;
#endif
			ret.push_back( Triplet< T >( i, j, v ) );
		}
		file.close();
		return ret;
	}

	/** 
	 * Loads a CRS text file and transforms it into a vector of Triplets.
	 * 
	 * @param fn Filename of the file to load from.
	 * @param m  Reference to where to store the number of rows of this matrix.
	 * @param n  Reference to where to store the number of columns of this matrix.
	 * @return A std::vector containing the read and converted triplets.
	 */
	static std::vector< Triplet< T > > loadCRS( const std::string fn, ULI &m, ULI &n ) {
		std::fstream file( fn.c_str(), std::ios::in );

		if ( !file.is_open() ) {
			std::cerr << "Error while opening file" << std::endl;
			exit( 1 );
		}

		ULI nonzeroes;
		file >> m;
		file >> n;
		file >> nonzeroes;
		
		std::vector< ULI > row_start( m+1 );
		for( size_t i = 0; i < m+1; ++i )
			file >> row_start[ i ];

		std::vector< ULI > cols( nonzeroes );
		for( size_t j = 0; j < nonzeroes; ++j )
			file >> cols[ j ];

		T   curval;
		ULI currow = 0;
		std::vector< Triplet< T > > ret;

		for( size_t k = 0; k < nonzeroes; ++k ) {
			while( k == row_start[ currow + 1 ] ) {
				currow++;
				if( currow > m + 1 ) {
					std::cerr << "Error in CRS file" << std::endl;
					exit( 1 );
				}
			}
			file >> curval;
			ret.push_back( Triplet< T >( currow, cols[ k ], curval ) );
		}

		file.close();

		return ret;		
	}

	/**
	 *  Saves an array of triplets to a file, in binary format.
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Array of triplets to write.
	 *  @param m Total number of rows in the matrix.+
	 *  @param n Total number of columns in the matrix.
	 *  @param s Size of the array toWrite.
	 */
        static void save( std::string fn, Triplet< T >* toWrite, const ULI m, const ULI n, const ULI s ) {
                std::fstream myFile ( fn.c_str(), std::ios::out | std::ios::binary);
                myFile.write( (char*) &m, sizeof( ULI ) );
                myFile.write( (char*) &n, sizeof( ULI ) );
                for( unsigned long int i = 0; i<s; i++ ) {
                        const ULI wi = toWrite[ i ].i();
                        const ULI wj = toWrite[ i ].j();
                        const double wv = toWrite[ i ].value;
#ifdef _DEBUG
                        std::cout << "Wrote: ( " << wi << " , " << wj << " , " << wv << " ) " << std::endl;
#endif
                        myFile.write( (char*) &( wi ), sizeof( ULI ) );
                        myFile.write( (char*) &( wj ), sizeof( ULI ) );
                        myFile.write( (char*) &( wv ), sizeof( T ) );
                }
                myFile.close();
        }

	/** Transposes this triplet, i.e., swapping the row and column value. */
	void transpose() { const ULI t = row; row = column; column = t; }

	/**
	 *  Saves a std::vector of triplets to a file, in binary format.
	 *  @param fn Filename to save to (overwrite mode).
	 *  @param toWrite Vector of triplets to write.
	 *  @param m Total number of rows in the matrix.
	 *  @param n Total number of columns in the matrix.
	 */
        static void save( std::string fn, std::vector< Triplet< T > > &toWrite, const ULI m, const ULI n ) {
                save( fn, &( toWrite[ 0 ] ), m, n, toWrite.size() );
        }

};

#endif

