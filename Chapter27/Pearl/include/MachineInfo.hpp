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
 *     A. N. Yzelman, Dept. of Computer Science, KU Leuven, 2011.
 */


#include <string>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <limits.h>

#ifndef _H_MACHINEINFO
#define _H_MACHINEINFO

/** Singleton class to get info on the current system */
class MachineInfo {

   private:
	/** Singleton class instance */
	static MachineInfo instance;

	/** Mutex to prevent concurrent accesses whenever appropriate */
	static pthread_mutex_t mutex;

	/** Whether the instance is initialised */
	static bool initialised;

	/** Base constructor (will load invalid values) */
	MachineInfo();

   protected:
	/** The number of available cores on this machine */
	unsigned long int P;

	/** Intialises the singleton instance */
	static void load();

   public:
	/** Gets a singleton instance */
	static const MachineInfo& getInstance();
	
	/** The number of available cores */
	unsigned long int cores() const;
	
};

#endif

