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
#include "HBICRS.hpp"

#ifdef _DEBUG
	#include <iostream>
#endif

#ifndef _H_BLOCKHILBERT
#define _H_BLOCKHILBERT

/** The Block Hilbert triplet scheme. In effect similar to (HierarchicalBICRS),
    but uses Hilbert coordinates to determine the order of the blocks,
    and a bisection algorithm to construct the individual blocks. */
template< typename T >
class BlockHilbert: public SparseMatrix< T, LI > {

   private:

   protected:

	/** Whether we use fixed block grid or a bisection-based grid */
	char bisection;

	/** Minimum number of nonzeroes per block */
	static const ULI BLOCK_M = 2048;
	static const ULI BLOCK_N = 2048;
	static const ULI MAX_BLOCKS = 1024*128;
	static const ULI MAX_DEPTH = 64; // NOTE: This CANNOT be larger than 64, unless unsigned long ints are larger than 64 bits!
	static const signed char DS = 2; // ICRS for use on the submatrices

	/** Minimum number of expansions */
	ULI minexp;

	/** Vector storing the block indices and their Hilbert coordinates */
	std::vector< HilbertTriplet< std::vector< Triplet < T > > > > hilbert;

	/** The actual data structure */
	HBICRS< T > *ds;

	char bisect( std::vector< HilbertTriplet< std::vector< Triplet< T > > > > &build,
			std::vector< Triplet< T > > &input, 
			const unsigned long int blocki, const unsigned long int blockj,
			const unsigned char depth ) {
		//find min/max row, column
		typename std::vector< Triplet< T > >::iterator it = input.begin();
		unsigned long int row, col; // number of rows/cols
		unsigned long int rl = it->i(); // initial min/max row/cols
		unsigned long int rh = it->i();
		unsigned long int cl = it->j();
		unsigned long int ch = it->j();
		++it;
		for( ; it != input.end(); ++it ) {
			if( rl > it->i() ) rl = it->i();
			if( rh < it->i() ) rh = it->i();
			if( cl > it->j() ) cl = it->j();
			if( ch < it->j() ) ch = it->j();
		}

		assert( rl <= rh );
		assert( cl <= ch );

		row = rh - rl + 1;
		col = ch - cl + 1;

		assert( row <= (unsigned long int)(this->nor) );
		assert( col <= (unsigned long int)(this->noc) );

		//remember bounds
		const unsigned long int row_lo = rl;
		const unsigned long int col_lo = cl;
		const unsigned long int row_hi = rh;
		const unsigned long int col_hi = ch;

		//build cumulative row/col arrays, count empty rows/cols
		unsigned long int re = 0;
		unsigned long int ce = 0;
		unsigned long int *rows = new unsigned long int[ row ];
		unsigned long int *cols = new unsigned long int[ col ];
		for( unsigned long int i=0; i<row; i++ ) rows[ i ] = 0;
		for( unsigned long int j=0; j<col; j++ ) cols[ j ] = 0;
		it = input.begin();
		for( ; it != input.end(); ++it ) {
			rows[ it->i()-rl ]++;
			cols[ it->j()-cl ]++;
		}
		for( unsigned long int i=1; i<row; i++ ) {
			if( rows[i] == 0 ) re++;
			rows[ i ] += rows[i-1];
		}
		for( unsigned long int j=1; j<col; j++ ) {
			if( cols[j] == 0 ) ce++;
			cols[ j ] += cols[j-1];
		}

		//check if not eligable to proceed, taking into account empty subrows, subcolumns
		if( row+col <= BLOCK_M + BLOCK_N + re + ce ) { //FIXME: empty count should take into account w (how many doubles in a cache line?)
		//if( row+col <= BLOCK_M + BLOCK_N ) {
			//already small enough block
			//so simply store current one again and exit
			build.push_back( HilbertTriplet< std::vector< Triplet< T > > >( blocki, blockj, input ) );
			delete [] rows; /*FIXME inefficient memory handling, but this is initialisation only */
			delete [] cols;
			return 2;
		}

		//set target
		const unsigned long int half = input.size() / 2;
		//prepare bisection
		rl = cl = 0; //reset to *local* column/row indices
		rh = row;
		ch = col;
		unsigned long int rs, cs;
		unsigned long int rm, cm;
		//do bisection search
		while( true ) {
			rm = (rl+rh)/2;
			cm = (cl+ch)/2;
			if( rows[ rm ] < half )      { rs = (rm - rl); rl = rm; }
			else if( rows[ rm ] > half ) { rs = (rh - rm); rh = rm; }
			else                           rs = 0;
			if( cols[ cm ] < half )      { cs = (cm - cl); cl = cm; }
			else if( cols[ cm ] > half ) { cs = (ch - cm); ch = cm; }
			else                           cs = 0;
			if( rs == 0 && cs == 0 ) break;
		}
		//note rh, ch, rl, ch, were used for search and now may not be the
		//values you may expect

		//prevent empty split if it all possible
		if( rm - row_lo == 0 && row_hi - rm > 1 ) rm++;
		if( row_hi - rm == 0 && rm - row_lo > 1 ) rm--;
		if( cm - col_lo == 0 && col_hi - cm > 1 ) cm++;
		if( col_hi - cm == 0 && cm - col_lo > 1 ) cm--;

		//choose a split direction
		char dir = 2; //0 for row, 1 for col
		if( rows[ rm ] > half ) rs = rows[ rm ] - half;
		else                    rs = half - rows[ rm ];
		if( cols[ cm ] > half ) cs = cols[ cm ] - half;
		else                    cs = half - cols[ cm ];
		if     ( rs > cs ) dir = 1;
		else if( cs > rs ) dir = 0;
		else               dir = rand() % 2; //tie break

		//don't split if one of the two parts will be empty
		if( dir ) {
			if( col_hi - cm == 0 || cm - col_lo == 0 ) {
				build.push_back( HilbertTriplet< std::vector< Triplet< T > > >( blocki, blockj, input ) );
				delete [] rows; /*FIXME inefficient memory handling, but this is initialisation only */
				delete [] cols;
				return 2;			
			}
		} else
			if( row_hi - rm == 0 || rm - row_lo == 0 ) {
				build.push_back( HilbertTriplet< std::vector< Triplet< T > > >( blocki, blockj, input ) );
				delete [] rows; /*FIXME inefficient memory handling, but this is initialisation only */
				delete [] cols;
				return 2;			
			}

		//split
		std::vector< Triplet< T > > t1, t2;
		it = input.begin();
		for( ; it != input.end(); ++it ) {
			if( dir ) { //column direction
				//remember to offset indices to go from local to global
				if( it->j() > cm+col_lo ) t2.push_back( *it );
				else t1.push_back( *it );
			} else {
				if( it->i() > rm+row_lo ) t2.push_back( *it );
				else t1.push_back( *it );
			}
		}

//		std::cout << "\t  storing block 1 of size " << (!dir?rm+1:row) << " by " << (dir?cm+1:col) << ", " << rows[rm] << " nnz" << std::endl;
//		std::cout << "\t  storing block 2 of size " << (!dir?row_hi-rm-row_lo+1:row) << " by " << (dir?col_hi-cm-col_lo+1:col) << ", " << (input.size()-rows[rm]) << " nnz" << std::endl;

		//add to return vector
		HilbertTriplet< std::vector< Triplet< T > > > p1( blocki, blockj, t1 );
		const unsigned long int new_i = dir ? blocki | (1<<(MAX_DEPTH-depth-1)) : blocki;
		const unsigned long int new_j =!dir ? blockj | (1<<(MAX_DEPTH-depth-1)) : blockj;
		HilbertTriplet< std::vector< Triplet< T > > > p2( new_i, new_j, t2 );
		if( t1.size() > 0 )
			build.push_back( p1 );
		if( t2.size() > 0 )
			build.push_back( p2 );

		//free memory
		delete [] rows;
		delete [] cols;

		//exit flag
		if( t1.size() == 0 || t2.size() == 0 ) return 0;
		return 1;
	}

	void bisect( std::vector< HilbertTriplet< std::vector< Triplet< T > > > > &build ) {
		char flag = 1;
		std::vector< char > go;
		go.reserve( build.size() );
		for( unsigned long int i=0; i<build.size(); ++i ) go.push_back( 1 );
		unsigned char depth = 0;
		while( flag && build.size() < MAX_BLOCKS && depth < MAX_DEPTH ) {
			std::vector< HilbertTriplet< std::vector< Triplet< T > > > > replace;
			std::vector< char > replace_go;
			replace_go.reserve( 2*build.size() );
			assert( build.size() == go.size() );
			for( unsigned long int i=0; i<build.size(); i++ )
				if( go[ i ] ) {
					char ret = bisect( replace, build[i].value, build[i].i(), build[i].j(), depth );
					if( ret == 1 ) {
						replace_go.push_back( 1 );
						replace_go.push_back( 1 );
					} else
						replace_go.push_back( 0 );
				} else {
					//just copy
					replace.push_back( build[i] );
					replace_go.push_back( 0 );
				}
			//check for stop condition
			if( replace.size() == build.size() ) flag = 0; //no new blocks created; done
			//replace
			build = replace;
			go    = replace_go;
//			std::cout << "\t  --------------------" << std::endl;
			depth++;
		}
	}
	
	/** Builds block array by bisection. resulting in variably-sized submatrices */
	std::vector< HilbertTriplet< std::vector< Triplet< T > > > > buildBisectionBlocks( std::vector< Triplet< T > > &input ) {
#ifndef NDEBUG
		const unsigned long int nnz = input.size();
#endif

		std::vector< HilbertTriplet< std::vector< Triplet< T > > > > ret;
		ret.push_back( HilbertTriplet< std::vector< Triplet< T > > >( 0, 0, input ) );
		bisect( ret );

		std::cout << "\t Bisection created " << ret.size() << " submatrices," << std::endl;

		unsigned long int sum = 0;
		typename std::vector< HilbertTriplet< std::vector< Triplet< T > > > >::iterator it = ret.begin();
		for( ; it != ret.end(); ++it ) {
			sum += it->value.size();
		}

		std::cout << "\t storing a total of " << sum << " nonzeroes." << std::endl;
	
		assert( nnz == sum );

		return ret;
	}

	/** Builds block array by ordering them in fixed-size sparse submatrices */
	std::vector< HilbertTriplet< std::vector< Triplet< T > > > > buildBlocks( std::vector< Triplet< T > > &input ) {
		const unsigned long int blockm = this->nor % BLOCK_M > 0 ? this->nor / BLOCK_M + 1 : this->nor / BLOCK_M;
		const unsigned long int blockn = this->noc % BLOCK_N > 0 ? this->nor / BLOCK_N + 1 : this->nor / BLOCK_N;
		const unsigned long int blocks = blockm * blockn;
		std::cout << "\tMaking " << blockm << " by " << blockn << " submatrices" << std::endl;
		std::cout << "\tyielding a total of " << blocks << " blocks. " << std::endl;
		std::vector< HilbertTriplet< std::vector< Triplet< T > > > > ret;
		for( unsigned long int i=0; i<blocks; i++ )
			ret.push_back( HilbertTriplet< std::vector< Triplet< T > > >( i / blockn, i % blockn, std::vector< Triplet< T > >() ) );
		typename std::vector< Triplet< T > >::iterator it = input.begin();
		for( ; it != input.end(); ++it ) {
			const unsigned long int blocki = it->i() / BLOCK_M;
			const unsigned long int blockj = it->j() / BLOCK_N;
			ret[ blockn * blocki + blockj ].value.push_back( *it );
		}
		return ret;
	}
	
	/** HilbertCoordinate comparison function. */
	bool cmp( HilbertTriplet< std::vector< Triplet< T > > > &left, HilbertTriplet< std::vector< Triplet< T > > > &right ) {

		if( left.i() == right.i() && left.j() == right.j() ) return true;

		const std::vector< unsigned long int > h_one = left.getHilbert();
		const std::vector< unsigned long int > h_two = right.getHilbert();

		unsigned long int max = h_one.size();
		bool max_is_one = true;
		if ( h_two.size() < max ) { max = h_two.size(); max_is_one = false; }
		for( unsigned long int i=0; i<max; i++ )
			if( h_one[ i ] != h_two[ i ] )
				return h_one[ i ] < h_two[ i ];
#ifdef _DEBUG		
		std::cout << "expand, ";
#endif
		if ( max_is_one )
			left.morePrecision( this->nor, this->noc );
		else
			right.morePrecision( this->nor, this->noc );

		return cmp( left, right );
	}

	/**
	 *  Binary search for finding a given triplet in a given range.
	 *  @param x triplet to-be found.
	 *  @param left Left bound of the range to search in.
	 *  @param right Right bound of the range to search in.
	 *  @return Index of the triplet searched.
	 */	
	unsigned long int find( HilbertTriplet< std::vector< Triplet< T > > > &x, ULI &left, ULI &right, ULI &minexp ) {
#ifdef _DEBUG
		std::cout << "left: " << left << ", right: " << right << std::endl;
#endif
		if( hilbert[left].getHilbert().size() > minexp ) minexp = hilbert[left].getHilbert().size();
		if( hilbert[right].getHilbert().size() > minexp ) minexp = hilbert[right].getHilbert().size();

		if ( left == right ) return left;
		if ( left+1 == right ) return right;
		if ( cmp( x, hilbert[ left ] ) ) return left;
		if ( cmp( hilbert[ right ], x ) ) return right+1;

		ULI mid = static_cast< unsigned long int >( ( left + right ) / 2 );
#ifdef _DEBUG
		std::cout << "mid: " << mid << std::endl;
#endif
		if ( cmp( x, hilbert[ mid ] ) )
			return find( x, left, mid, minexp );
		else
			return find( x, mid, right, minexp );
	}

   public:

	/** Base deconstructor. */
	virtual ~BlockHilbert() {
		if( ds != NULL )
			delete ds;
	}

	/** Base constructor. */
	BlockHilbert() {
		bisection = 0;
	}

	/** Base constructor.
	 *  Will read in from Matrix Market file.
	 *  @param bisect Whether bisection-based blocking should be used.
	 *  @see SparseMatrix::SparseMatrix( file, zero )
	 */
	BlockHilbert( std::string file, T zero = 0, char bisect = 0 ) {
		bisection = bisect;
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
	BlockHilbert( std::vector< Triplet< T > >& input, LI m, LI n, T zero ) {
		bisection = 0;
		load( input, m, n, zero );
	}


	/** 
	 *  Base constructor.
	 *  Warning: the zero parameter is currently NOT USED!
	 *  @param input Raw input of normal triplets.
	 *  @param m Total number of rows.
	 *  @param n Total number of columns.
	 *  @param bisect Whether bisection-based blocking should be used.
	 *  @param zero What elements is considered to-be zero.
	 */
	BlockHilbert( std::vector< Triplet< T > >& input, LI m, LI n, char bisect, T zero ) {
		bisection = bisect;
		load( input, m, n, zero );
	}

	/** @see SparseMatrix::load */
	virtual void load( std::vector< Triplet< T > >& input, const LI m, const LI n, const T zero ) {
		ds = NULL;
		this->zero_element = 0;
		this->nor = m;
		this->noc = n;
		this->nnz = input.size();

		std::cout << std::endl << "Loading into *experimental* block hilbert structure." << std::endl;
	
		if( bisection ) {
			std::cout << "\tUsing bisection to create blocks" << std::endl;
			hilbert = buildBisectionBlocks( input );
		} else {
			std::cout << "\tUsing fixed block sizes..." << std::endl;
			hilbert = buildBlocks( input );
		}

		std::cout << "\tGetting Hilbert coordinates" << std::endl;
		typename std::vector< HilbertTriplet< std::vector< Triplet< T > > > >::iterator it = hilbert.begin();
		for( ; it!=hilbert.end(); ++it )
			it->calculateHilbertCoordinate();
		std::cout << "\tUsing std::sort to get a Hilbert ordering on the sparse blocks" << std::endl;
		HilbertTripletCompare< std::vector< Triplet< T > > > compare;
		std::sort( hilbert.begin(), hilbert.end(), compare );

#ifdef _DEBUG
	        for( ULI i=0; i<this->nnz; i++ ) {
			for( ULI j=0; j<hilbert[i].getHilbert().size(); j++ )
				std::cout << hilbert[i].getHilbert()[j] << ",";
			std::cout << std::endl;
		}
#endif
		//load into HBICRS
		std::cout << "\tLoading into HBICRS..." << std::endl;
		std::vector< std::vector< Triplet< T > > > hierarchy;
		it = hilbert.begin();
		for( ; it!=hilbert.end(); ++it ) hierarchy.push_back( it->value );
		signed char *hierarchy_datatype = new signed char[ hierarchy.size() ];
		for( unsigned long int i=0; i<hierarchy.size(); i++ ) hierarchy_datatype[i] = DS; //ICRS default
		ds = new HBICRS< T >( hierarchy, hierarchy_datatype, this->nor, this->noc, 0 );
		delete [] hierarchy_datatype;

		std::cout << "BlockHilbert structure ready!" << std::endl;
	}

	/** Returns the first nonzero index, per reference. */
	virtual void getFirstIndexPair( LI &row, LI &col ) {
		ds->getFirstIndexPair( row, col );
	}

	/**
	 * Calculates z=xA.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zxa( const T* x, T* z ) {
		ds->zxa( x, z );
	}

	/**
	 * Calculates z=Ax.
	 * z is *not* set to 0 at the start of this method!
	 *
	 * @param x The (initialised) input vector.
	 * @param z The (initialised) output vector.
	 */ 
	virtual void zax( const T* x, T* z ) {
		ds->zax( x, z );
	}

	virtual size_t bytesUsed() {
		return ds->bytesUsed();
	}

};

#endif

