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
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2012.
 */


#ifndef _H_UPSCALER
#define _H_UPSCALER

#include "SBDTree.hpp"
#include "Triplet.hpp"

#include<algorithm>
#include<vector>
#include<map>

/** 
 *  Transforms SBD-structures over q parts into an SBD structure over p parts,
 *  with q>p, and q,p both powers of two.
 */
template< typename T >
class Upscaler : public SBDTree {

	private:

		/** ``Turn off'' standard constructor */
		Upscaler() {}

	protected:

		/** All nonzeroes, stored block-by-block. */
		std::vector< std::vector< Triplet< T > > > nonzeroes;

		/** Keeps track which processes are represented in which blocks. */
		std::vector< std::vector< bool > > containsPID;

		/** Pre-order tree iterator */
		class treeIterator {
			protected:

				/** Current position of the iterator. */
				unsigned long int walk;

				/** Starting position of the iterator. */
				unsigned long int ID;

				/**
				 * Processed[ i ] is true when this iterator has visited its i-th child. 
				 * May be unused if a pre-allocated vector is given to the iterator instead.
				 */
				std::vector< bool >  processed;

				/**
				 * Pointer to the `processed'-vector actually used.
				 * Can be a pointer to treeIterator::processed, or a user-defined vector instead.
				 */
				std::vector< bool > *p_processed;

				/** Class this iterator works on. */
				Upscaler< T > *p;

			public:

				/**
				 * Base constructor.
				 *  @param _p  Pointer to the base class
				 *  @param initial Initial position in the tree.
				 *  @param _pp Pointer to a vector where the visited subtrees may be stored.
				 */
				treeIterator( Upscaler< T > *_p, const unsigned long int initial, std::vector< bool > *_pp = NULL ):
					walk( initial ), ID( initial ), p_processed( _pp ), p( _p ){
					if( _pp == NULL ) {
						for( unsigned long int k=0; k<p->size(); ++k )
							processed.push_back( false );
						p_processed = &processed;
					} else {
						if( p_processed->size() < p->size() )
							p_processed->resize( p->size() );
						for( unsigned long int k=0; k<p->size(); ++k )
							p_processed->operator[]( k ) = false;
					}
					p_processed->operator[]( ID ) = true;
				}

				/** @return Returns the current position. */
				unsigned long int position() { return walk; }

				/** @return Returns false if there is no next element. */
				virtual bool next() {
					do {
						const unsigned long int P = p->parent[ walk ];
						const unsigned long int l = p->left_child[ walk ];
						const unsigned long int r = p->right_child[ walk ];
						if( l != NO_SUCH_NODE && !p_processed->operator[]( l ) ) {
							walk = l;
							break;
						} else if( r != NO_SUCH_NODE && !p_processed->operator[]( r ) ) {
							walk = r;
							break;
						} else
							walk = P;
					} while( walk != p->parent[ ID ] );
					
					if( walk == p->parent[ ID ] )
						return false; //done

					assert( !p_processed->operator[]( walk ) );
					//flag as visited
					p_processed->operator[]( walk ) = true;
					//return success
					return true;
				}
		};

		/** Same as treeIterator, but does in-order traversal instead of pre-order */
		class treeInOrderIterator : public treeIterator {
			public:

				/** 
				 * Base constructor. @see treeIterator::treeIterator.
				 *
				 * @param _p      Iterates on this object.
				 * @param initial Initial position.
				 * @param _pp     Buffer where the visited subtrees may be stored.
				 */
				treeInOrderIterator( Upscaler< T > *_p, const unsigned long int initial, std::vector< bool > *_pp = NULL ):
					treeIterator( _p, initial, _pp ) {
					this->p_processed->operator[]( this->ID ) = false; //undo pre-order traversal default visit of the top node
					next(); //assume at least one unprocessed node
				}

				/**
				 * Moves iterator to its next position.
				 *
				 * @return Whether there was a next position available.
				 */
				virtual bool next() {
					while( this->walk != this->p->parent[ this->ID ] ) {
						//run left as far as is possible
						while( this->p->left_child[ this->walk ] != NO_SUCH_NODE &&
							!this->p_processed->operator[]( this->p->left_child[ this->walk ] ) )
							this->walk = this->p->left_child[ this->walk ];

						//if the current position is already processed
						if( this->p_processed->operator[]( this->walk ) ) {
							//if I could not go left, try to go right and repeat
							if( this->p->right_child[ this->walk ] != NO_SUCH_NODE &&
								!this->p_processed->operator[]( this->p->right_child[ this->walk ] ) ) { 
								this->walk = this->p->right_child[ this->walk ];
							} else {
								//otherwise go up
								this->walk = this->p->parent[ this->walk ];
							}
						} else {
							//walk is the current one to be processed
							this->p_processed->operator[]( this->walk ) = true;
							return true;
						}
					}
					//I am back at the start location and the start location is processed
					return false;
				}
		};

		/** Same as treeIterator, but does post-order traversal instead of pre-order */
		class treePostOrderIterator : public treeIterator {
			public:

				/** 
				 * Base constructor. @see   treeIterator::treeIterator.
				 *
				 * @param _p      Iterates on this object.
				 * @param initial Initial position.
				 * @param _pp     Keeps track of which children are already visited.
				 */
				treePostOrderIterator( Upscaler< T > *_p, const unsigned long int initial, std::vector< bool > *_pp = NULL ):
					treeIterator( _p, initial, _pp ) {
					this->p_processed->operator[]( this->ID ) = false; //undo pre-order traversal default visit of the top node
					next(); //assume at least one unprocessed node
				}

				/**
				 * Moves iterator to next position.
				 *
				 * @return Returns false if there is no next element.
				 */
				virtual bool next() {
					while( this->walk != this->p->parent[ this->ID ] ) {
						//run left as far as is possible
						while( this->p->left_child[ this->walk ] != NO_SUCH_NODE &&
							!this->p_processed->operator[]( this->p->left_child[ this->walk ] ) )
							this->walk = this->p->left_child[ this->walk ];

						//if there is a right unprocessed child, go there and repeat
						if( this->p->right_child[ this->walk ] != NO_SUCH_NODE && 
							!this->p_processed->operator[]( this->p->right_child[ this->walk ] ) )
							this->walk = this->p->right_child[ this->walk ];
						else if( !this->p_processed->operator[]( this->walk ) ) {
							//I cannot go left nor right, but before I go up I should process current node
							this->p_processed->operator[]( this->walk ) = true;
							return true;
						} else {
							//I cannot go left nor right and current node is processed,
							//so go up and repeat
							this->walk = this->p->parent[ this->walk ];
						}
					}
					//I am back at the start location and the start location is processed
					return false;
				}
		};
		
		/**
		 * Determines, given a collection of nonzeroes, which rows and columns nonzeroes owned by s reside on,
		 * and flags these rows and columns used in the appropriate vectors.
		 *
		 * @param nonzeroes The collection of nonzeroes.
		 * @param s         The process ID.
		 * @param min_i	    Minimum row index to consider (other nonzeroes are ignored).
		 * @param min_j     Minimum column index to consider (other nonzeroes are ignored).
		 * @param emptyRows emptyRows[i] is true when no nonzeroes on row i are encountered.
		 * @param emptyCols emptyCols[j] is true when no nonzeroes on column j exist.
		 * @param empty     Will be true when process s has no nonzeroes in this collection.
		 */
		void determineEmpty( const std::vector< Triplet< T > > &nonzeroes, const unsigned long int s,
					const unsigned long int min_i, const unsigned long int min_j,
					std::vector< bool > &emptyRows, std::vector< bool > &emptyCols,
					bool &empty ) {
			empty = true;
			for( unsigned long int k = 0; k < nonzeroes.size(); ++k ) {
				if( nonzeroes[ k ].meta != s ) continue;
				const unsigned long int i = nonzeroes[ k ].i();
				const unsigned long int j = nonzeroes[ k ].j();
				assert( i - min_i < emptyRows.size() );
				assert( j - min_j < emptyCols.size() );
				if( emptyRows[ i - min_i ] ) emptyRows[ i - min_i ] = false;
				if( emptyCols[ j - min_j ] ) emptyCols[ j - min_j ] = false;
				if( empty ) empty = false;
			}
		}
	
		/** Determine min/max of nonzeroes owned by processor s inside node walk. */
		void updateMinMax( const unsigned long int walk, const unsigned long int s,
			unsigned long int &min_i, unsigned long int &max_i,
			unsigned long int &min_j, unsigned long int &max_j ) {
			for( unsigned long int k = 0; k < nonzeroes[ walk ].size(); k++ )
				if( nonzeroes[ walk ][ k ].meta == s ) {
					const unsigned long int i = nonzeroes[ walk ][ k ].i();
					const unsigned long int j = nonzeroes[ walk ][ k ].j();
					if( i < min_i ) min_i = i;
					if( i > max_i ) max_i = i;
					if( j < min_j ) min_j = j;
					if( j > max_j ) max_j = j;
				}
		}

		/** Determine min/max of nonzeroes owned by processor s,
		 *  contained in the subtree with root ID, as well as all separators
		 *  on the path from ID to the true root
		 */
		void determineMinMax( const unsigned long int ID, const unsigned long int s,
					unsigned long int &min_i, unsigned long int &max_i,
					unsigned long int &min_j, unsigned long int &max_j ) {
			//determine range of nonzeroes
			min_i = min_j = ULONG_MAX;
			max_i = max_j = 0;
			//go up
			unsigned long int walk = ID;
			for( ; walk != NO_SUCH_NODE; walk = parent[ walk ] ) {
				updateMinMax( walk, s, min_i, max_i, min_j, max_j );
			}
			//go down
			treeIterator it( this, ID );
			while( it.next() ) {
				//I am at a unprocessed node, so do processing
				updateMinMax( it.position(), s, min_i, max_i, min_j, max_j );
			}
		}

		/** adds nonzeroes from a given node into a given vector */
		void addNonzeroes( std::vector< Triplet< T > > &into,
					const unsigned long int from, const unsigned long int s,
					const std::map< unsigned long int, unsigned long int > &rowGlobalToLocal,
					const std::map< unsigned long int, unsigned long int > &colGlobalToLocal ) {
			for( unsigned long int i=0; i<nonzeroes[ from ].size(); i++ ) {
				if( nonzeroes[ from ][ i ].meta != s )
					continue;
				const unsigned long int old_i = nonzeroes[ from ][ i ].i();
				const unsigned long int old_j = nonzeroes[ from ][ i ].j();
				assert( rowGlobalToLocal.find( old_i ) != rowGlobalToLocal.end() );
				assert( colGlobalToLocal.find( old_j ) != colGlobalToLocal.end() );
				const unsigned long int new_i = rowGlobalToLocal.find( old_i )->second;
				const unsigned long int new_j = colGlobalToLocal.find( old_j )->second;
				into.push_back( Triplet< double >( new_i, new_j, nonzeroes[ from ][ i ].value ) );
				into.back().meta = s;
			}
		}

		/**
		 *  Reads out a subtree and returns the upscaled version.
		 *  Does global to local index translation. The non-binary
		 *  part of the tree is returned in remote_nonzeroes and
		 *  are not part of the upscaled_hierarchy structures.
		 *
		 *  @param ID root of the subtree to read out.
		 *  @param s  only output nonzeroes distributed to processor s.
		 *  @param local_nonzeroes nonzeroes corresponding to this new tree (flat vector).
		 *  @param remote_nonzeroes nonzeroes belonging to s contained in the path from ID to the root.
		 *  @param upscaled_hierarchy hierarchy corresponding to the upscaled nonzeroes.
		 *  @param upscaled_row_bounds row-wise boundaries corresponding to the upscaled nonzeroes.
		 *  @param upscaled_column_bounds column-wise boundaries corresponding to the upscaled nonzeroes.
		 *  @param rowLocalToGlobal maps local row indices to global indices.
		 *  @param columnLocalToGlobal maps local column indices to global indices.
		 *  @param rowGlobalToLocal maps global row indices to local indices.
		 *  @param colGlobalToLocal maps global column indices to local indices.
		 */
		void getSubTree( const unsigned long int ID,
				const unsigned long int s,
				std::vector< Triplet< T > > &local_nonzeroes,
				std::vector< Triplet< T > > &remote_nonzeroes,
				std::vector< unsigned long int > &upscaled_hierarchy,
				std::vector< unsigned long int > &upscaled_row_bounds,
				std::vector< unsigned long int > &upscaled_column_bounds,
				std::vector< unsigned long int > &rowLocalToGlobal,
				std::vector< unsigned long int > &columnLocalToGlobal,
				std::map< unsigned long int, unsigned long int > &rowGlobalToLocal,
				std::map< unsigned long int, unsigned long int > &colGlobalToLocal ) {

			//determine range of nonzeroes
			unsigned long int min_i, min_j, max_i, max_j;
	                determineMinMax( ID, s, min_i, max_i, min_j, max_j );

			//max should be an exclusive bound
			max_i++; max_j++;

			//report global ranges
			std::cout << "Processor " << s << " has global ranges \t( " << min_i << ", " << max_i << " ) x ( " << min_j << ", " << max_j << " )" << std::endl;

			//find empty rows and columns
			std::vector< bool > emptyRows, emptyCols, emptyNodes;
			for( unsigned long int i = min_i; i < max_i; ++i )
				emptyRows.push_back( true );
			for( unsigned long int j = min_j; j < max_j; ++j )
				emptyCols.push_back( true );
			for( unsigned long int k = 0; k < size(); k++ )
				emptyNodes.push_back( true );
			assert( emptyRows.size() == max_i - min_i );
			assert( emptyCols.size() == max_j - min_j );

			//down
			treeIterator it( this, ID );
			do {
				bool empty;
				determineEmpty( nonzeroes[ it.position() ], s, min_i, min_j, emptyRows, emptyCols, empty );
				emptyNodes[ it.position() ] = empty;
			} while( it.next() );

			//up
			unsigned long int walk = ID;
			while( walk != root ) {
				bool empty;
				walk = parent[ walk ];
				determineEmpty( nonzeroes[ walk ], s, min_i, min_j, emptyRows, emptyCols, empty );
				emptyNodes[ walk ] = empty;
			}

			//build local to global, and their inverses
			for( unsigned long int i = 0; i < max_i - min_i; ++i )
				if( !emptyRows[ i ] ) {
					rowGlobalToLocal[ i + min_i ] = rowLocalToGlobal.size();
					rowLocalToGlobal.push_back( i + min_i );
				}
			for( unsigned long int j = 0; j < max_j - min_j; ++j )
				if( !emptyCols[ j ] ) {
					colGlobalToLocal[ j + min_j ] = columnLocalToGlobal.size();
					columnLocalToGlobal.push_back( j + min_j );
				}

			//derive lowest and highest index by traversing the subtree
			unsigned long int lowest_index, highest_index;
			treeIterator it2( this, ID );
			lowest_index  = ULONG_MAX;
			highest_index = 0;
			do {
				const unsigned long int pos = it2.position();
				if( emptyNodes[ pos ] ) continue; //ignore empty nodes
				if( pos < lowest_index  ) lowest_index  = pos;
				if( pos > highest_index ) highest_index = pos;
			} while( it2.next() );

			//determine number of empty rows and column
			unsigned long int numEmptyRows, numEmptyCols;
			numEmptyRows = numEmptyCols = 0;
			for( unsigned long int i = 0; i < emptyRows.size(); i++ )
				if( emptyRows[ i ] ) numEmptyRows++;
			for( unsigned long int j = 0; j < emptyCols.size(); j++ )
				if( emptyCols[ j ] ) numEmptyCols++;

			//report local size
			std::cout << "Processor " << s << " has local size \t" << (max_i-numEmptyRows-min_i) << " x " << (max_j-numEmptyCols-min_j) << std::endl;

			//determine upscaled hierarchy and boundaries, write local nonzeroes
			treeInOrderIterator in_it( this, ID );
			std::map< unsigned long int, unsigned long int >::iterator glob2loc;
			do {
				const unsigned long int pos = in_it.position();

				//do not add out-of-range nodes
				if( pos < lowest_index || pos > highest_index ) continue;

				addNonzeroes( local_nonzeroes, pos, s, rowGlobalToLocal, colGlobalToLocal );

				glob2loc = rowGlobalToLocal.lower_bound( this->r_lo[ pos ] );
				assert( glob2loc->first >= this->r_lo[ pos ] );
				upscaled_row_bounds.push_back( glob2loc->second );

				glob2loc = colGlobalToLocal.lower_bound( this->c_lo[ pos ] );
				assert( glob2loc->first >= this->c_lo[ pos ] );
				upscaled_column_bounds.push_back( glob2loc->second );

				upscaled_hierarchy.push_back ( pos == ID ? 0 : parent[ pos ] - lowest_index );
			} while( in_it.next() );

			//insert final boundary points
			std::vector< unsigned long int >::iterator loc2glob;
			loc2glob = std::lower_bound( rowLocalToGlobal.begin(), rowLocalToGlobal.end(), this->r_hi[ highest_index ] );
			upscaled_row_bounds.push_back( loc2glob - rowLocalToGlobal.begin() );

			loc2glob = std::lower_bound( columnLocalToGlobal.begin(), columnLocalToGlobal.end(), this->c_hi[ highest_index ] );
			upscaled_column_bounds.push_back( loc2glob - columnLocalToGlobal.begin() );

			//write non-local nonzeroes
			walk = parent[ ID ];
			while( walk != NO_SUCH_NODE ) {
				addNonzeroes( remote_nonzeroes, walk, s, rowGlobalToLocal, colGlobalToLocal );
				walk = parent[ walk ];
			}
			//done!
		}

	public:

		/** Base constructor.*/
		Upscaler( const std::vector< Triplet< T > > &nonzeroes, const unsigned long int P,
				std::vector< unsigned long int > &r_hierarchy,
				std::vector< unsigned long int > &c_hierarchy,
				std::vector< unsigned long int > &r_bounds,
				std::vector< unsigned long int > &c_bounds,
				unsigned long int *row_perm  = NULL,
				unsigned long int *col_perm  = NULL,
				unsigned long int *proc2proc = NULL ) : SBDTree( r_hierarchy, c_hierarchy, r_bounds, c_bounds ) {
			//limitation check
			const unsigned long int trueP = r_bounds.size() / 2;
			if( floor(log2(trueP)) != log2(trueP) && trueP != P ) {
				std::cout << "P=" << trueP << " is not a power of 2. Current Upscaler cannot handle incomplete binary trees, exiting." << std::endl;
				exit( 1 );
			}

			//we will group nonzeroes and put them in the corresponding SBD tree nodes
			this->nonzeroes.resize( this->size() );
			containsPID.resize( this->size() );
			for( unsigned long int i=0; i<this->size(); ++i )
				for( unsigned long int s=0; s<P; ++s )
					containsPID[ i ].push_back( false );
			//to prevent inserting nonzeroes on at a time in a top-to-bottom way,
			//construct a boundary-to-node map to speed up finding the correct
			//internal node for each nonzero (note this is still of asymptotic
			//complexity O(nnz*log(hierarchy.size()))
			std::map< unsigned long int, unsigned long int > rowMap, colMap;
			//walk from back to front so that the correct node index is stored
			//for empty separators
			for( unsigned long int i=r_bounds.size()-1; i>0; --i ) {
				rowMap[ r_bounds[ i ] ] = i - 1;
			}
			//note that we identify blocks by their upper corner points!
			for( unsigned long int j=c_bounds.size()-1; j>0; --j ) {
				colMap[ c_bounds[ j ] ] = j - 1;
			}
			//for each nonzero, find the correct SBD tree spot
			for( unsigned long int k=0; k<nonzeroes.size(); ++k ) {
				const unsigned long int i = row_perm  == NULL ? nonzeroes[ k ].i()  : row_perm[ nonzeroes[ k ].i() ];
				const unsigned long int j = col_perm  == NULL ? nonzeroes[ k ].j()  : col_perm[ nonzeroes[ k ].j() ];
				const unsigned long int s = proc2proc == NULL ? nonzeroes[ k ].meta : proc2proc[ nonzeroes[ k ].meta ];
				const unsigned long int u = rowMap.upper_bound( i )->second;
				const unsigned long int v = colMap.upper_bound( j )->second;
				Triplet< T > newTriplet( i, j, nonzeroes[ k ].value );
				newTriplet.meta = s;
				assert( u < r_hierarchy.size() );
				assert( v < c_hierarchy.size() );
				assert( r_lo[ u ] <= i && i < r_hi[ u ] );
				assert( c_lo[ v ] <= j && j < c_hi[ v ] );
				if( u == v ) {
					containsPID[ u ][ s ] = true;
					this->nonzeroes[ u ].push_back( newTriplet );
				} else if( u == root || v == root ) {
					containsPID[ root ][ s ] = true;
					this->nonzeroes[ root ].push_back( newTriplet );
				} else if( (u && 1) == 1 || (v && 1 ) == 1)   {
					unsigned long int walk = (u && 1) == 1 ? parent[ u ] : parent[ v ];
					const unsigned long int compare = (u && 1) == 1 ? v : u;
					while( walk != root && walk != compare )
						walk = parent[ walk ];
					containsPID[ walk ][ s ] = true;
					this->nonzeroes[ walk ].push_back( newTriplet );
				} else {
					std::cerr << "Nonzero in two different pure SBD blocks, which is impossible." << std::endl <<
							"Quitting on account of a bad SBD hierarchy input." << std::endl;
					exit( EXIT_FAILURE );
				}
			}

			//check
			unsigned long int newnnz = 0;
			for( unsigned long int t = 0; t < size(); ++t )
				newnnz += this->nonzeroes[ t ].size();
			assert( newnnz = nonzeroes.size() );
		}

		/** Base deconstructor. */
		~Upscaler() {}

		/**
		 *  Reads out a subtree corresponding to only the nonzeroes owned by processor s,
		 *  and returns the upscaled version. Does global to local index translation.
		 *  The non-binary part of the tree is returned in remote_nonzeroes and are not
		 *  part of the upscaled_hierarchy structures.
		 *
		 *  @param s  only output nonzeroes distributed to processor s.
		 *  @param local_nonzeroes nonzeroes corresponding to this new tree (flat vector).
		 *  @param remote_nonzeroes nonzeroes belonging to s contained in the path from ID to the root.
		 *  @param upscaled_hierarchy hierarchy corresponding to the upscaled nonzeroes.
		 *  @param upscaled_row_bounds row-wise boundaries corresponding to the upscaled nonzeroes.
		 *  @param upscaled_column_bounds column-wise boundaries corresponding to the upscaled nonzeroes.
		 *  @param rowLocalToGlobal maps local row indices to global indices.
		 *  @param columnLocalToGlobal maps local column indices to global indices.
		 *  @param rowGlobalToLocal maps global row indices to local indices.
		 *  @param colGlobalToLocal maps global column indices to local indices.
		 *  @param P    The total number of parts to reduce the SBD tree to.
		 *  @param Pref The total number of parts (blocks) in this SBD tree.
		 */
		void getSubTree( const unsigned long int s,
				std::vector< Triplet< T > > &local_nonzeroes,
				std::vector< Triplet< T > > &remote_nonzeroes,
				std::vector< unsigned long int > &upscaled_hierarchy,
				std::vector< unsigned long int > &upscaled_row_bounds,
				std::vector< unsigned long int > &upscaled_column_bounds,
				std::vector< unsigned long int > &rowLocalToGlobal,
				std::vector< unsigned long int > &columnLocalToGlobal,
				std::map< unsigned long int, unsigned long int > &rowGlobalToLocal,
				std::map< unsigned long int, unsigned long int > &colGlobalToLocal,
				const unsigned long int P, const unsigned long int Pref ) {

			//determine which nodes contain nonzeroes belonging to s
			std::vector< bool > empty;
			for( unsigned long int i = 0; i < size(); ++i ) empty.push_back( true );

			//check the nonzeroes themselves
			for( unsigned long int i = 0; i < size(); ++i ) {
				for( unsigned long int k = 0; k < nonzeroes[ i ].size(); ++k ) {
					if( nonzeroes[ i ][ k ].meta == s ) {
						empty[ i ] = false;
						break;
					}
				}
			}

			//if one of my children is nonempty, I am nonempty too
			treePostOrderIterator post_it( this, this->root );
			do {
				const unsigned long int pos = post_it.position();
				if( this->left_child[ pos ] == NO_SUCH_NODE ) {
					//if this is a leaf node, do nothing
					assert( this->right_child[ pos ] == NO_SUCH_NODE );
				} else {
					//if this is an internal node, do the above check
					if( !empty[ this->left_child[ pos ] ] ) {
						if( empty[ pos ] ) empty[ pos ] = false;
					} else if( !empty[ this->right_child[ pos ] ] ) {
						if( empty[ pos ] ) empty[ pos ] = false;
					}
				}
			} while( post_it.next() );

			//the highest node with two nonempty children will be the new root
			unsigned long int root = ULONG_MAX;
			treeIterator pre_it( this, this->root );
			do {
				const unsigned long int pos = pre_it.position();
				if( !empty[ pos ] && 					//if there is no internal root,
					this->left_child[ pos ] == NO_SUCH_NODE &&	//take only the leaf node as
					this->right_child[ pos ] == NO_SUCH_NODE ) {	//root, assuming there is
						if( root == ULONG_MAX )			//exactly one
							root = pos; 				
						else {
							std::cout << "There seem to be multiple leaf nodes beloning to target processor, yet there are no shared separators?" << std::endl;
							exit( 1 );
						}
					}
				if( this->left_child[ pos ] != NO_SUCH_NODE ) {
					assert( this->right_child[ pos ] != NO_SUCH_NODE );
					if( !empty[ this->left_child[ pos ] ] &&
						!empty[ this->right_child[ pos ] ] ) {
						root = pos;
						break;
					}
				}
			} while( pre_it.next() );

			if( root == ULONG_MAX ) {
				std::cerr << "No root found! Error in hierarchy, exiting..." << std::endl;
				exit( 1 );
			}

			//continue with that node
			getSubTree( root, s, local_nonzeroes,
					remote_nonzeroes, upscaled_hierarchy,
					upscaled_row_bounds, upscaled_column_bounds,
					rowLocalToGlobal, columnLocalToGlobal,
					rowGlobalToLocal, colGlobalToLocal );
		}

		/** Reads out SBD data and prints to std::cout. Useful for debugging purposes. */
		void readout() {
			std::cout << "Index:     ";
			for( unsigned long int i=0; i<size(); ++i )
				std::cout << i << " ";
			std::cout << std::endl;
			std::cout << "Hierarchy: ";
			for( unsigned long int i=0; i<size(); ++i )
				std::cout << (this->parent[ i ]==NO_SUCH_NODE?0:this->parent[ i ]) << " ";
			std::cout << std::endl;
			std::cout << "Left:      ";
			for( unsigned long int i=0; i<size(); ++i )
				std::cout << (this->left_child[ i ]==NO_SUCH_NODE?0:this->left_child[ i ]) << " ";
			std::cout << std::endl;
			std::cout << "Right:     ";
			for( unsigned long int i=0; i<size(); ++i )
				std::cout << (this->right_child[ i ]==NO_SUCH_NODE?0:this->right_child[ i ]) << " ";
			std::cout << std::endl;
		}

};

#endif
