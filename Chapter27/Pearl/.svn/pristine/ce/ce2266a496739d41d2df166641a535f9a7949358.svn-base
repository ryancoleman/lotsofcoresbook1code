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
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2013.
 * 
 * Defines default alignment variables. Currently AVX is standard.
 * Automatically adapts for Xeon Phi.
 */


#ifndef _H_SL_ALIGNMENT
#define _H_SL_ALIGNMENT

#ifdef __MIC__
 #define _SL_ALIGN_DOUBLE 64
 #define _SL_ALIGN_INT16 32
 #define _SL_ALIGN_INT32 64
 #define _SL_BLOCKSIZE 8
#else
 #define _SL_ALIGN_DOUBLE 64
 #define _SL_ALIGN_INT16 64
 #define _SL_ALIGN_INT32 64
 #define _SL_BLOCKSIZE 4
#endif

#endif

