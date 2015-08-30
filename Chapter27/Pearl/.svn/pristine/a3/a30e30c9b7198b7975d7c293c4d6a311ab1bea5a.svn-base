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


#include "SBDTree.hpp"

void SBDTree::build( std::vector< unsigned long int > &hierarchy,
			std::vector< unsigned long int > &r_bounds,
			std::vector< unsigned long int > &c_bounds ) {
	//keep track of root set (input check)
	root_is_set = false;
	assert( r_bounds.size() == hierarchy.size() + 1 );
	assert( c_bounds.size() == hierarchy.size() + 1 );

	//must have less than ULONG_MAX entries, as ULONG_MAX is code for root node or leaf nodes
	if( hierarchy.size() == ULONG_MAX ) {
		std::cerr << "SBDTree: Too many blocks (equals size of unsigned long)!" << std::endl;
		exit( EXIT_FAILURE );
	}
	
	parent = new unsigned long int[ hierarchy.size() ];
	left_child = new unsigned long int[ hierarchy.size() ];
	right_child = new unsigned long int[ hierarchy.size() ];
	r_lo = new unsigned long int[ hierarchy.size() ];
	r_hi = new unsigned long int[ hierarchy.size() ];
	c_lo = new unsigned long int[ hierarchy.size() ];
	c_hi = new unsigned long int[ hierarchy.size() ];
	nodes = hierarchy.size();

	//set left and right childs to NO_SUCH_NODE
	for( unsigned long int i=0; i<nodes; i++ )
		left_child[i] = right_child[i] = NO_SUCH_NODE;

	//read out hierarchy
	unsigned long int i = 0;
	std::vector< unsigned long int >::iterator it;
	for( it = hierarchy.begin(); it != hierarchy.end(); ++it, ++i ) {
		//node i has parent *it-1
		const unsigned long int parent = *it - 1;
		if( *it == 0 ) {
			if( root_is_set ) {
				std::cerr << "SBDTree: Root was already set, error in hierarchy array!" << std::endl;
				exit( EXIT_FAILURE );
			}
			this->parent[ i ] = NO_SUCH_NODE;
			root = i;
			root_is_set = true;
		} else
			this->parent[ i ] = parent;

		//set left or right child of the parent, unless this is the root node
		if( *it > 0 ) {
			if( i < parent )
				left_child[ parent ] = i;
			else
				right_child[ parent ] = i;
		}

		//set boundaries
		r_lo[ i ] = r_bounds[ i ];
		r_hi[ i ] = r_bounds[i+1];
		c_lo[ i ] = c_bounds[ i ];
		c_hi[ i ] = c_bounds[i+1];
	}

	//done
	if( !root_is_set ) {
		std::cerr << "SBDtree: No root set, error in hierarchy array!" << std::endl;
		exit( EXIT_FAILURE );
	}
}

/** Base constructor */
SBDTree::SBDTree( std::vector< unsigned long int > &r_hierarchy, std::vector< unsigned long int > &c_hierarchy,
			std::vector< unsigned long int > &r_bounds,
			std::vector< unsigned long int > &c_bounds ) {
	//r_hierarchy and c_hierarchy are supposed to be equal, check.
	assert( r_hierarchy.size() == c_hierarchy.size() );
	std::vector< unsigned long int >::iterator it1, it2;
	for( it1 = r_hierarchy.begin(), it2 = c_hierarchy.begin(); it1 != r_hierarchy.end() && it2 != c_hierarchy.end(); ++it1, ++it2 )
		assert( *it1 == *it2 );
	if( it1 != r_hierarchy.end() || it2 != c_hierarchy.end() )
		assert( false );

	build( r_hierarchy, r_bounds, c_bounds );
}

/** Base constructor. Warning: avoids some assertions! */
SBDTree::SBDTree( std::vector< unsigned long int > &hierarchy,
			std::vector< unsigned long int > &r_bounds,
			std::vector< unsigned long int > &c_bounds ) {
	build( hierarchy, r_bounds, c_bounds );
}

/** Base deconstructor. */
SBDTree::~SBDTree() {
	delete [] parent;
	delete [] left_child;
	delete [] right_child;
	delete [] r_lo;
	delete [] r_hi;
	delete [] c_lo;
	delete [] c_hi;
}

/** Gets, from a separator node, the bounding box of the nonzeroes contained in the separator.
    Note that this is *not* the r_lo/hi,c_lo/hi of the node itself; those are the bounds of the
    row-wise and column-wise separators.

    This is a logarithmic operation.
  */
void SBDTree::getSeparatorBB( const unsigned long int index,
				unsigned long int &r_lo, unsigned long int &r_hi,
				unsigned long int &c_lo, unsigned long int &c_hi ) {
	//first check if this indeed is a separator node
	if( index % 2 == 0 ) {
		std::cerr << "SBDTree::getSeparatorBB: passed node index does not correspond with a separator!" << std::endl;
		exit( EXIT_FAILURE );
	}
	//we have an SBD structure
	//so we need the lower row and lower column index of the upper-right block;
	//so recursively follow the left child
	unsigned long int cur = left_child[ index ];
	assert( cur != NO_SUCH_NODE ); //ensure left child exists
	while( left_child[ cur ] != NO_SUCH_NODE )
		cur = left_child[ cur ];
	assert( cur % 2 == 0 ); //this should not be a leaf node
	r_lo = this->r_lo[ cur ];
	c_lo = this->c_lo[ cur ];
	//for the upper indices, ask the rightmost child
	cur = right_child[ index ];
	assert( cur != NO_SUCH_NODE ); //ensure right child exists
	while( right_child[ cur ] != NO_SUCH_NODE )
		cur = right_child[ cur ];
	assert( cur % 2 == 0 ); //this should not be a leaf node
	r_hi = this->r_hi[ cur ];
	c_hi = this->c_hi[ cur ];
	//done!
}

/** Returns the parent of a given node. */
unsigned long int SBDTree::up( const unsigned long int index ) {
	if( parent[ index ] == NO_SUCH_NODE ) {
		std::cerr << "SBDTree::parent: root node has no parent!" << std::endl;
		exit( EXIT_FAILURE );
	}
	return parent[ index ];
}

/** Returns the left child of a given node. */
unsigned long int SBDTree::left( const unsigned long int index ) {
	if( left_child[ index ] == NO_SUCH_NODE ) {
		std::cerr << "SBDTree::left: node " << index << " has no left child!" << std::endl;
		exit( EXIT_FAILURE );
	}
	return left_child[ index ];
}

/** Returns the right child of a given node. */
unsigned long int SBDTree::right( const unsigned long int index ) {
	if( right_child[ index ] == NO_SUCH_NODE ) {
		std::cerr << "SBDTree::right: node " << index << " has no right child!" << std::endl;
		exit( EXIT_FAILURE );
	}
	return right_child[ index ];
}

/** Returns the row bounds corresponding to a given node. */
void SBDTree::rowBounds( const unsigned long int index,
		unsigned long int &r_lo,
		unsigned long int &r_hi ) {
	r_lo = this->r_lo[ index ];
	r_hi = this->r_hi[ index ];
}

/** Returns the column bounds corresponding to a given node. */
void SBDTree::columnBounds( const unsigned long int index,
			unsigned long int &c_lo,
			unsigned long int &c_hi ) {
	c_lo = this->c_lo[ index ];
	c_hi = this->c_hi[ index ];
}

/** Whether the given node is a leaf node */
char SBDTree::isLeaf( const unsigned long int index ) {
	if( left_child[ index ] == NO_SUCH_NODE ) {
		assert( left_child[ index ] == right_child[ index ] );
		return true;
	} else
		return false;
}

unsigned long int SBDTree::size() {
	return nodes;
}

unsigned long int SBDTree::getRoot() {
	return root;
}

