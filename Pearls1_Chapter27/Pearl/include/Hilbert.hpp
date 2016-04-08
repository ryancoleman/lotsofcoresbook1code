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


#include <vector>
#include "SparseMatrix.hpp"
#include "HilbertTriplet.hpp"
#include "HilbertTripletCompare.hpp"
#include "Triplet.hpp"

#ifndef _H_HILBERT
#define _H_HILBERT

#define _HILBERT_COMPRESSED_BICRS

#ifdef _DEBUG
	#include <iostream>
#endif

#ifdef _HILBERT_COMPRESSED_BICRS
	#include "CBICRS.hpp"
#else
	#include "BICRS.hpp"
#endif

/** The Hilbert scheme backed by (C)BICRS.
 *  In effect similar to the Hilbert Triplet scheme (HTS),
 *   but uses BICRS to store the nonzeroes. */
template< typename T >
class Hilbert: public SparseMatrix< T, ULI > {

   private:

   protected:

	/** Minimum number of expansions */
	ULI minexp;

	/** Vector storing the non-zeros and their locations. */
	std::vector< HilbertTriplet< T > > ds;

	/** Actual data structure. */
	Matrix< T > *ads;

	/** Gets the data structure. Convience function, enables quick changes in Hilbert backing structure. */
	static Matrix< T >* getDataStructure( std::vector< Triplet< T > > &tds, const ULI m, const ULI n, T zero ) {
#ifdef _HILBERT_COMPRESSED_BICRS
		return CBICRS_factory< T >::getCBICRS( tds, m, n, zero );
#else
		return new BICRS< T >( tds, m, n, zero );
#endif
	}

   public:

	/** Base deconstructor. */
	virtual ~Hilbert() {
		if( ads != NULL )
			delete ads;
	}

	/** Base constructor. */
	Hilbert() {
		ads = NULL;
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	Hilbert( std::string file, T zero = 0 ) {
		ads = NULL;
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
	Hilbert( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
		ads = NULL;
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		this->zero_element = 0;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();
		for( ULI i=0; i<this->nnz; i++ ) {
			ds.push_back( HilbertTriplet< T >( input[ i ].i(), input[ i ].j(), input[ i ].value ) );
			ds[ i ].calculateHilbertCoordinate();
		}
		HilbertTripletCompare< T > compare;
		std::sort( ds.begin(), ds.end(), compare );
#ifdef _DEBUG
	        for( ULI i=0; i<this->nnz; i++ )
			std::cout << ds[i].getMostSignificantHilbertBits() << " " << ds[i].getLeastSignificantHilbertBits() << std::endl;
#endif
		//create temporary structure
		std::vector< Triplet< T > > tds;
		typename std::vector< HilbertTriplet< T > >::iterator it = ds.begin();
		for( ; it!=ds.end(); ++it )
			tds.push_back( Triplet< T >( it->i(), it->j(), it->value ) );
		//create actual structure
		ads = getDataStructure( tds, m, n, zero );
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = ds[ 0 ].i();
		col = ds[ 0 ].j();
	}

	/**
	 * Calculates z=xA.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zxa( const T* x, T* z ) {
		ads->zxa( x, z );
	}

	/**
	 * Calculates z=Ax.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zax( const T* x, T* z ) {
		ads->zax( x, z );
	}

	/** Gets the number of bytes used by this storage scheme. */
	virtual size_t bytesUsed() {
		return ads->bytesUsed();
	}

	/** Saves the current Hilbert structure in binary triplet form.
	 *  @param fn Filename to save to.
	 *  @see HilbertTriplet< T >::save
`	 */
	void saveBinary( const std::string fn ) {
		HilbertTriplet< T >::save( fn, ds, this->nor, this->noc );
	}

	/** Loads from binary triplets, assumes Hilbert ordering already done */
	void loadBinary( const std::string fn ) {
		std::cerr << "Warning: assuming binary file was saved by a HTS scheme, i.e., that it is pre-ordered using Hilbert coordinates." << std::endl;
		std::vector< Triplet< T > > tds = Triplet< T >::load( fn, this->nor, this->noc );
		ads = getDataStructure( tds, this->nor, this->noc, 0 );
	}

};

#endif

