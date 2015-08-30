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
 *     A. N. Yzelman, Dept. of Mathematics, Utrecht University, 2009.
 */


#include "SparseMatrix.hpp"
#include "Triplet.hpp"
#include <vector>

//#define _DEBUG

#ifndef _H_SVM
#define _H_SVM

#ifdef _DEBUG
#include<iostream>
#endif

/**
 *  The sparse vector matrix format.
 */
template< typename T >
class SVM: public SparseMatrix< T, ULI > {

  private:

  protected:

	/** Row major (0), column major (1) */
	char major;

	/** SVM structure */
	std::vector< std::vector< Triplet< T > > > ds;

 	/** Sorts 1D columnwise */
        static int compareTripletsR( const void * left, const void * right ) {
                const Triplet< T > one = *( (Triplet< T > *)left );
                const Triplet< T > two = *( (Triplet< T > *)right );
                if ( one.j() < two.j() )
                        return -1;
                if ( one.j() > two.j() )
                        return 1;
                return 0;
        }

	/** Sorts 1D rowwise */
	static int compareTripletsC( const void * left, const void * right ) {
		const Triplet< T > one = *( (Triplet< T > *)left );
		const Triplet< T > two = *( (Triplet< T > *)right );
		if ( one.i() < two.i() )
			return -1;
		if ( one.i() > two.i() )
			return 1;
		return 0;
	}

	/** 
         *  Helper function which finds a value with a given index.
         *  @param col_index The given column index.
         *  @param row_index The given row index.
	 *  @param ret Reference to the variable where the return triplet is stored (if found).
	 *  @return Whether or not a non-zero value was found.
	 */
	bool find( const ULI col_index, const ULI row_index, Triplet< double > &ret ) {
		if( major ) { //column major
			std::vector< Triplet< double > >::iterator it = ds.at( col_index ).begin();
			Triplet< double > cur;
			while( it!=ds.get( col_index ).end() ) {
				cur = *(it++);
				if( cur.i() == row_index ) {
					ret = cur;
					return true;
				}
			}
		} else {
			std::vector< Triplet< double > >::iterator it = ds.get( row_index ).begin();
			Triplet< double > cur;
			while( it!=ds.get( row_index ).end() ) {
				cur = *(it++);
				if( cur.j() == col_index ) {
					ret = cur;
					return true;
				}
			}
		}
		return false;
	}

  public:

	/** Base constructor. */
	SVM() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	SVM( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays.
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns to be stored.
	 *  @param zero Which element is considered to be the zero element.
	 */
	SVM( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, const T zero ):
		SparseMatrix< T, ULI >( number_of_nonzeros, number_of_rows, number_of_cols, zero ) {}

	/** Copy constructor.
	 *  @param toCopy reference to the CRS datastructure to copy.
	 */
	SVM( SVM< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->nor = toCopy.nor;
		this->noc = toCopy.noc;
		major = toCopy.major;
		for( ULI i=0; i<(major?this->noc:this->nor); i++ ) {
			std::vector< Triplet< double > > toAdd;
			std::vector< Triplet< double > >::iterator it = toCopy.ds.get( i ).begin();
			while( it!=toCopy.ds.get( i ).end() ) toAdd.push_back( *(it++) );
			ds.push_back( toAdd );
		}
	}
	
	/**
	 *  Constructor which transforms a collection of input triplets to SVM format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *  @param input The input collection.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero Which element is considered zero.
	 *  @param direction 0 for row-major, 1 for column major.
	 */
	SVM( const std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero, const char direction = 0 ) {
		load( input, m, n, zero, direction );
	}

	/** 
	 *  This will default to row-major SVM format.
	 *  @see SparseMatrix::load
	 */
	virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
		load( input, m, n, zero, 0 );
	}

	/**
	 *  Will load a collection of triplets into SVM format.
         *  The input collection is considered to have at most one triplet with unique
         *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
         *  met.
         *  @param input The input collection.
         *  @param m Total number of rows.
         *  @param n Total number of columns.
         *  @param zero Which element is considered zero.
         *  @param direction 0 for row-major, 1 for column major.
         */
	void load( const std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero, const char direction ) {
		this->zero_element = zero;
		this->nnz = input.size();
		this->nor = m;
		this->noc = n;
		major = direction;

		//build datastructure, move input there
		ds.resize(major?this->noc:this->nor);
		typename std::vector< Triplet< T > >::const_iterator in_it = input.begin();
		for( ; in_it != input.end(); in_it++ ) {
			Triplet< T > cur = *in_it;
			const ULI currow = major ? cur.j() : cur.i();
			const T value = cur.value;
			if( value == this->zero_element ) { this->nnz--; continue; }
			ds.at( currow ).push_back( cur );
		}

                //sort
                for( ULI currow = 0; currow < this->nor; currow++ ) {
                        if( ds.at( currow ).size() == 0 ) continue;
                        qsort( &( ds.at( currow )[ 0 ] ), ds.at( currow ).size(), sizeof( Triplet< T > ), 
				( major ? &compareTripletsC : &compareTripletsR ) );
                }
	}

	/**
	 *  @return Direct access to the SVM datastructure.
	 */
	std::vector< std::vector< Triplet< double > > >& getData() {
		return ds;
	}

	/**
	 *  Method which provides random matrix access to the stored sparse matrix.
	 *  @param i Row index.
	 *  @param j Column index.
	 *  @return Matrix valuei at (i,j).
	 */
	T& random_access( ULI i, ULI j ) {
		Triplet< T > ret;
		if ( find( i, j, ret ) ) {
			return ret.value;
		} else
			return this->zero_element;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = ds[0].at(0).i();
		col = ds[0].at(0).j();
	}

	/** 
	 *  In-place z=xA function.
	 *  
	 *  @param x The x vector to left-multiply current matrix with.
	 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 */
        virtual void zxa( const T*__restrict__ x, T*__restrict__ z ) {
		if( major ) {
			ULI col;
			for( col = 0; col < this->noc; col++, z++ ) {
				std::vector< Triplet< double > >::iterator it = ds.at( col ).begin();
				while( it!=ds.at(col).end() ) {
					const Triplet< double > cur = *(it++);
					*z += cur.value * x[ cur.i() ];
				}
			}
		} else {
			ULI row;
			for( row = 0; row < this->nor; row++, x++ ) {
				std::vector< Triplet< double > >::iterator it = ds.at( row ).begin();
				while( it!=ds.at( row ).end() ) {
					const Triplet< double > cur = *(it++);
					z[ cur.j() ] += cur.value * *x;
				}
			}
		}
      	}

	/** 
	 *  In-place z=Ax function.
	 *  
	 *  @param x The x vector to multiply current matrix with.
	 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 */
        virtual void zax( const T*__restrict__ x, T*__restrict__ z ) {
		if( major ) {
			ULI col;
			for( col = 0; col < this->noc; col++, x++ ) {
				std::vector< Triplet< double > >::iterator it = ds.at( col ).begin();
				while( it!=ds.at(col).end() ) {
					const Triplet< double > cur = *(it++);
					z[ cur.i() ] = cur.value * *x;
				}
			}
		} else {
			ULI row;
			for( row = 0; row < this->nor; row++, z++ ) {
				std::vector< Triplet< double > >::iterator it = ds.at( row ).begin();
				while( it!=ds.at( row ).end() ) {
					const Triplet< double > cur = *(it++);
					*z += cur.value * x[ cur.j() ];
				}
			}
		}
      	}

	virtual size_t bytesUsed() {
		//not optimal since Triplets are stored within vectors of vectors (the last term is unnecessary)
		return sizeof( ULI ) * 2 * this->nnz + sizeof( T ) * this->nnz + sizeof( void* ) * this->nor;
	}

	/** Base deconstructor. */
	~SVM() {}

};

#undef _DEBUG
#endif

