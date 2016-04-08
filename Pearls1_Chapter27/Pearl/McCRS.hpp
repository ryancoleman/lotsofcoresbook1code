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


#include "Triplet.hpp"
#include "SparseMatrix.hpp"
#include <vector>
#include <assert.h>

//#define _DEBUG

#ifndef _H_McCRS
#define _H_McCRS

#ifdef _DEBUG
#include<iostream>
#endif

//Uncomment the below to enable interleaved allocation of the matrix data structure
//(this can lead to performance gains,
// but on our testset average performance dropped with about 30 percent!)
#define INTERLEAVE_A

//no libnuma, no interleave
#ifdef _NO_LIBNUMA
 #undef INTERLEAVE_A
#endif

/**
 *  The compressed row storage sparse matrix data structure.
 */
template< typename T >
class McCRS: public SparseMatrix< T, ULI > {

  private:

  protected:

	/** Array keeping track of individual row starting indices. */
	ULI* row_start;

	/** Array containing the actual nnz non-zeros. */
	T* ds;

	/** Array containing the column indeces corresponding to the elements in ds. */
	ULI* col_ind;

 	/** Sorts 1D columnwise */
        static int compareTriplets( const void * left, const void * right ) {
                const Triplet< T > one = **( (Triplet< T > **)left );
                const Triplet< T > two = **( (Triplet< T > **)right );
                if ( one.j() < two.j() )
                        return -1;
                if ( one.j() > two.j() )
                        return 1;
                return 0;
        }     

	/** 
         *  Helper function which finds a value with a given column index on a given subrange of indices.
         *  @param col_index The given column index.
	 *  @param search_start The start index of the subrange (inclusive).
	 *  @param search_end The end index of the subrange (exlusive).
	 *  @param ret Reference to the variable where the return *index* is stored.
	 *  @return Whether or not a non-zero value should be returned.
	 */
	bool find( const ULI col_index, const ULI search_start, const ULI search_end, ULI &ret ) {
		for( ULI i=search_start; i<search_end; i++ )
			if( col_ind[ i ] == col_index ) {
				ret = i;
				return true;
			}
		return false;
	}

  public:

	/** Base constructor. */
	McCRS() {}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	McCRS( std::string file, T zero = 0 ) {
		this->loadFromFile( file, zero );
	}
	
	/**
	 *  Base constructor which only initialises the internal arrays. Note that to gain a valid McCRS structure,
	 *  these arrays have to be filled by some external mechanism (i.e., after calling this constructor, the
	 *  internal arrays contain garbage, resuling in invalid datastructure).
	 *  @param number_of_nonzeros The number of non-zeros to be stored.
	 *  @param number_of_rows The number of rows to be stored.
	 *  @param number_of_cols The number of columns of the matrix.
	 *  @param zero The element considered to be zero.
	 */
	McCRS( const ULI number_of_nonzeros, const ULI number_of_rows, const ULI number_of_cols, T zero ) {
		this->nnz = number_of_nonzeros;
		this->nor = number_of_rows;
		this->noc = number_of_cols;
		this->zero_element = zero;
#ifdef INTERLEAVE_A
		row_start = (ULI*) numa_alloc_interleaved( (this->nor + 1) * sizeof( ULI ) );
		ds = (T*) numa_alloc_interleaved( (this->nnz) * sizeof( T ) );
		col_ind = (ULI*) numa_alloc_interleaved( this->nnz * sizeof( ULI ) );
#else
		row_start = new ULI[ (this->nor + 1) ];
		ds        = new T[ this->nnz ];
		col_ind   = new ULI[ this->nnz ];
#endif
	}

	/** Copy constructor.
	 *  @param toCopy reference to the McCRS datastructure to copy.
	 */
	McCRS( McCRS< T >& toCopy ) {
		this->zero_element = toCopy.zero_element;
		this->nnz = toCopy.nnz;
		this->nor = toCopy.nor;
#ifdef INTERLEAVE_A
		row_start = (ULI*) numa_alloc_interleaved( (this->nor + 1) * sizeof( ULI ) );
		ds = (T*) numa_alloc_interleaved( (this->nnz) * sizeof( T ) );
		col_ind = (ULI*) numa_alloc_interleaved( (this->nnz) * sizeof( ULI ) );
#else
		row_start = new ULI[ (this->nor + 1) ];
		ds        = new T[ this->nnz ];
		col_ind   = new ULI[ this->nnz ];
#endif
		for( ULI i=0; i<this->nnz; i++ ) {
			ds[ i ] = toCopy.ds[ i ];
			col_ind[ i ] = toCopy.col_ind[ i ];
		}
		for( ULI i=0; i<this->nor; i++ )
			row_start[ i ] = toCopy.row_start[ i ];
	}
	
	/**
	 *  Constructor which transforms a collection of input triplets to McCRS format.
	 *  The input collection is considered to have at most one triplet with unique
	 *  pairs of indeces. Unspecified behaviour occurs when this assumption is not
	 *  met.
	 *  @param input The input collection.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param zero The element considered to be zero.
	 */
	McCRS( std::vector< Triplet< T > > input, ULI m, ULI n, T zero ) {
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
		std::cout << "\tLoading in a vector of " << input.size() << " triplets into McCRS..." << std::endl;

		this->zero_element = zero;
		//find nnz
		this->nnz = input.size();
	
		this->nor = m;
		this->noc = n;

		//build better datastructure
		std::vector< std::vector< Triplet< T >* > > ds( this->nor );
		
		//move input there
		typename std::vector< Triplet< T > >::iterator in_it = input.begin();
		for( ; in_it != input.end(); ++in_it ) {
			//Triplet< T >* cur = &(*in_it);
			const ULI currow = in_it->i();
			const T value = in_it->value;
			if( value == this->zero_element ) { 
				this->nnz--;
				continue;
			}
			ds.at( currow ).push_back( &(*in_it) );
		}

		//allocate arrays
#ifdef INTERLEAVE_A
		row_start = (ULI*) numa_alloc_interleaved( (this->nor + 1) * sizeof( ULI ) );
		this->ds = (T*) numa_alloc_interleaved( (this->nnz) * sizeof( T ) );
		col_ind = (ULI*) numa_alloc_interleaved( (this->nnz) * sizeof( ULI ) );
#else
		row_start = new ULI[ (this->nor + 1) ];
		this->ds  = new T[ this->nnz ];
		col_ind   = new ULI[ this->nnz ];
#endif

                //make McCRS
                ULI index = 0;
                for( ULI currow = 0; currow < this->nor; currow++ ) {
                        row_start[ currow ] = index;
                        if( ds.at( currow ).size() == 0 ) continue;
                        qsort( &( ds.at( currow )[ 0 ] ), ds.at( currow ).size(), sizeof( Triplet< T >* ), &compareTriplets );
                        typename std::vector< Triplet< T >* >::iterator row_it = ds.at( currow ).begin();
                        for( ; row_it!=ds.at( currow ).end(); row_it++ ) {
                                const Triplet< T > cur = *(*row_it);
                                this->ds[ index ] = cur.value;
                                col_ind[ index ] = cur.j();
                                index++;
                        }
                }
		row_start[ this->nor ] = this->nnz;
		std::cout << "\t" << index << " nonzeroes loaded into McCRS structure." << std::endl;
		assert( index == this->nnz );
	}

	/**
	 *  Method which provides random matrix access to the stored sparse matrix.
	 *  @param i Row index.
	 *  @param j Column index.
	 *  @return Matrix valuei at (i,j).
	 */
	T& random_access( ULI i, ULI j ) {
		ULI found_index;
		if ( find( j, row_start[ i ], row_start[ i+1 ], found_index ) ) {
#ifdef _DEBUG
			std::cout << "Searched col_ind between " << row_start[ i ] << " and " << row_start[ i + 1 ] << ", found: " << std::endl;
			std::cout << "Element (" << i << "," << j << ") found on index " << found_index << ", returning " << ds[ found_index ] << std::endl;
#endif
			return ds[ found_index ];
		} else
			return this->zero_element;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( ULI &row, ULI &col ) {
		row = this->row_start[ 0 ];
		col = this->col_ind[ 0 ];
	}

#ifndef _NO_LIBNUMA
	/** Overloaded mv call; allocates output vector using numa_interleaved. */
	virtual T* mv( const T* x ) {
		T* ret = (T*) numa_alloc_interleaved( this->nor * sizeof( T ) );
		for( ULI i=0; i<this->nor; i++ ) ret[ i ] = this->zero_element;
		zax( x, ret );
		return ret;
	}
#endif

	/** 
	 *  In-place z=xA function.
	 */
        virtual void zxa( const T*__restrict__ x, T*__restrict__ z ) {
		std::cerr << "CRS z=xA (left-sided SpMV multiplication) is not possible to do in parallel using CRS and OpenMP parallel fors" << std::endl;
		std::cerr << "Exiting..." << std::endl;
		exit( 1 );
        }

	/** 
	 *  In-place z=Ax function.
	 *  
	 *  @param x The x vector to multiply current matrix with.
	 *  @param z The result vector. Must be pre-allocated and its elements should be initialised to zero.
	 */
        virtual void zax( const T*__restrict__ x, T*__restrict__ z ) {
		const unsigned long int nor = this->nor;
		ULI * const col_ind = this->col_ind;
		ULI * const row_start = this->row_start;
		T * const ds = this->ds;
		#pragma omp parallel for shared( x, z ) schedule( guided )
		for( ULI row = 0; row < nor; row++ ) {
			ULI index, *j_p;
			T sum, *v_p, x_e;
			sum = 0.0;
			v_p = ds + row_start[ row ];
			j_p = col_ind + row_start[ row ];
			for( index = row_start[ row ]; index < row_start[ row + 1 ]; index++ ) {
				x_e = *(x + (*j_p++));
				sum += (*v_p++) * x_e;
			}
			z[ row ] = sum;
		}
        }

	virtual size_t bytesUsed() {
		return sizeof( ULI ) * ( this->nnz + this->nor ) + sizeof( T ) * this->nnz;
	}

	/** Returns pointer to the row_start vector. */
	ULI* rowJump() { return row_start; }

	/** Returns pointer to the column index vector. */
	ULI* columnIndices() { return col_ind; }

	/** Returns pointer to the matrix nonzeros vector. */
	T* values() { return ds; } 

	/** Base deconstructor. */
	virtual ~McCRS() {
#ifdef INTERLEAVE_A
		numa_free( row_start, (this->nor + 1) * sizeof( ULI ) );
		numa_free( ds,        (this->nnz) * sizeof( T ) );
		numa_free( col_ind,   (this->nnz) * sizeof( ULI ) );
#else
		delete [] row_start;
		delete [] ds;
		delete [] col_ind;
#endif
	}

};

#undef _DEBUG
#endif

