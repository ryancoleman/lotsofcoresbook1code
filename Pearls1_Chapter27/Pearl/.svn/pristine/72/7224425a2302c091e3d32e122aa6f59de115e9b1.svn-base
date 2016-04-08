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

#ifndef _H_CUSTOM_QUICKSORT
#define _H_CUSTOM_QUICKSORT

#include <vector>
#include <algorithm>

/** 
 * Determines a new pivot index for use in a customised quicksort implementation.
 * Note: an enhancement would be to look for the median value, if multiple same-valued keys can exist
 *       (currently they do not; but this becomes relevant when this scheme is combined with blocking)
 * 
 * @param key   The key value array
 * @param left  The left boundary of the search range (inclusive)
 * @param right The right boundary of the search range (exclusive)
 *
 * @return Index to the new pivot value (>= left and < right).
 */
template< typename KeyType >
size_t custom_quicksort_getpivot( const std::vector< KeyType > &key, const size_t left, const size_t right ) {
	size_t ret = right - left; //overflow-safe (as opposed to (right+left)/2)
	ret /= 2;
	ret += left;
	return ret;
}

/**
 * Custom quicksort inner-loop. Handles the swapping of the key array, as well as of the linked
 * data array.
 *
 * @param key   The key value array, along which the sort takes place.
 * @param data  The data value array, which will follow the same permutation as the to-be sorted key array.
 * @param left  Left boundary of the range within to sort (inclusive).
 * @param right Right boundary of the range within to sort (exclusive).
 * @param pivot Pivot index used in this inner-loop.
 *
 * @return Final pivot position after inner loop.
 */
template< typename KeyType, typename DataArrayType >
size_t custom_quicksort_inner( std::vector< KeyType > &key, DataArrayType &data,
				const size_t left, const size_t right, const size_t pivot ) {
	
	//move pivot to end
	std::swap( key [ pivot ], key [ right-1 ] );

	//data array mirrors every swap of the key array
	std::swap( data[ pivot ], data[ right-1 ] );

	//remember pivot value
	const BigInt pivotValue = key[ right-1 ];

	//screen from left to right
	size_t curIndex = left;
	for( size_t i = left; i < (right-1); ++i ) {

		//compare key to pivot value
		if( key[ i ] < pivotValue ) {
			//swap small entry to the left of the range
			if( i != curIndex ) {
				std::swap( key [ i ], key [ curIndex ] );
				std::swap( data[ i ], data[ curIndex ] );
			}
			//increment bound
			curIndex++;
		}
	}

	//move pivot value back to the middle
	std::swap( key [ curIndex ], key [ right-1 ] );
	std::swap( data[ curIndex ], data[ right-1 ] );

	//return new pivot position
	return curIndex;
}

/**
 * Custom quicksort outer-loop. Handles choosing the pivot, calling the inner-loop, and recursion.
 *
 * On entry, left must be less than right. This is not checked by this function.
 *
 * @param key   The key value array, along which the sort takes place.
 * @param data  The data value array, which will follow the same permutation as the to-be sorted key array.
 * @param left  Left boundary of the range within to sort (inclusive).
 * @param right Right boundary of the range within to sort (exclusive).
 */
template< typename KeyType, typename DataArrayType >
void custom_quicksort_outer( std::vector< KeyType > &key, DataArrayType &data, size_t left, size_t right ) {
	//determine pivot value
	const size_t pivot = custom_quicksort_getpivot( key, left, right );

	//perform inner-loop (swap-based sorting)
	const size_t newPivotIndex = custom_quicksort_inner< KeyType, DataArrayType >( key, data, left, right, pivot );
	
	//recurse
	if( newPivotIndex > left + 1 )
		custom_quicksort_outer< KeyType, DataArrayType >( key, data, left, newPivotIndex );
	if( newPivotIndex + 1 < right )
		custom_quicksort_outer< KeyType, DataArrayType >( key, data, newPivotIndex + 1, right );
}

/**
 * Custom quicksort entry-point. Switches between an indirect sort or a direct sort.
 *
 * @param key   The key value array, along which the sort takes place.
 * @param data  The data value array, which will follow the same permutation as the to-be sorted key array.
 * @param left  Left boundary of the range within to sort (inclusive).
 * @param right Right boundary of the range within to sort (exclusive).
 */
template< typename KeyType, typename DataType >
void custom_quicksort( std::vector< KeyType > &key, std::vector< DataType > &data, size_t left, size_t right ) {
	const bool implicit = sizeof( DataType ) > 2 * sizeof( size_t );

	//check bounds
	if( right <= left )
		return;

	//use a permutation array if the normal data type is too large
	if( implicit ) {

		//allocate permutation array
		size_t *permutation = new size_t[ right - left ];

		//initialise permutation array
		for( size_t i = left; i < right; ++i )
			permutation[ i ] = i;

		//perform quicksort using permutation array as the linked array
		custom_quicksort_outer< KeyType, size_t* >( key, permutation, left, right );

		//perform permutation on the actual data array
		std::vector< DataType > temp;

		//first create temporary array
		for( size_t i = left; i < right; ++i )
			temp.push_back( data[ permutation[ i ] ] );

		//copy back
		for( size_t i = left; i < right; ++i )
			data[ i ] = temp[ i - left ];

		//free permutation array
		delete [] permutation;

		//free temporary array
		temp.clear();

	} else {

		//otherwise do a direct linked quicksort
		custom_quicksort_outer< KeyType, std::vector< DataType > >( key, data, left, right );
	}
}

#endif

