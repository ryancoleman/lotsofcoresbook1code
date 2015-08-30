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


#ifndef _H_DUCK
#define _H_DUCK

#include "BlockOrderer.hpp"

/** Codes the Minimal CCS block order */
template< typename T >
class Duck: public BlockOrderer< T > {
	protected:
		virtual void pre_readout ( const unsigned long int index ) {
			//switch between internal nodes and external nodes
			if( this->tree->isLeaf( index ) ) {
				//immediately add everything in range
				this->output->push_back( this->items[ index ] );
				if( this->datatype != NULL )
					this->datatype->push_back( NORMAL_DS );
			}
		}

		virtual void in_readout  ( const unsigned long int index ) {
			//skip leaf node
			if( this->tree->isLeaf( index ) ) return;
			//initialise
			typename std::vector< Triplet< T > >::iterator it = this->items[ index ].begin();
			std::vector< Triplet< T > > cur1;
			std::vector< Triplet< T > > cur2;
			std::vector< Triplet< T > > cur3;
			std::vector< Triplet< T > > cur4;
			std::vector< Triplet< T > > rep;
			//loop over this node's triplets
			for( ; it != this->items[ index ].end(); ++it ) {
				//now the 3 vertical separators (incl middle) in row order
				if( this->upper_vertical( index, *it ) )
					cur1.push_back( *it );
				else if( this->middle( index, *it ) )
					cur2.push_back( *it );
				else if( this->left_horizontal( index, *it ) )
					cur3.push_back( *it );
				else if( this->right_horizontal( index, *it ) )
					cur4.push_back( *it );
				else
					rep.push_back( *it );
			}
			if( cur1.size() + cur2.size() + cur3.size() + cur4.size() > 0 )
				this->items[ index ] = rep; //replace with smaller vector
			this->output->push_back( cur1 );
			this->output->push_back( cur2 );
			this->output->push_back( cur3 );
			this->output->push_back( cur4 );
			if( this->datatype != NULL ) {
				this->datatype->push_back( VERTICAL_DS );
				this->datatype->push_back( NORMAL_DS );
				this->datatype->push_back( HORIZONTAL_DS );
				this->datatype->push_back( HORIZONTAL_DS );
			}
		}

		virtual void post_readout( const unsigned long int index ) {
			//skip leaf node
			if( this->tree->isLeaf( index ) ) return;
			//initialise
			typename std::vector< Triplet< T > >::iterator it = this->items[ index ].begin();
			std::vector< Triplet< T > > cur;
			std::vector< Triplet< T > > rep;
			//loop over this node's triplets
			for( ; it != this->items[ index ].end(); ++it ) {
				//right horizontal separator last
				if( this->lower_vertical( index, *it ) )
					cur.push_back( *it );
				else
					rep.push_back( *it );
			}
			if( cur.size() > 0 )
				this->items[ index ] = rep; //replace with smaller vector
			assert( rep.size() == 0 );
			this->output->push_back( cur );
			if( this->datatype != NULL )
				this->datatype->push_back( VERTICAL_DS );
		}

	public:
		//(de)constructor of superclass gets called automagically
};

#endif

