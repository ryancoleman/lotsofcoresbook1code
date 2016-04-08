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


#include "MachineInfo.hpp"

MachineInfo::MachineInfo() : P( ULONG_MAX ) {}

const MachineInfo& MachineInfo::getInstance() {
	//first fast check
	if( initialised ) return instance;
	//instance might not have been created, get lock
	pthread_mutex_lock( &mutex );
	//check again
	if( !initialised ) {
		//init instance
		load();
	}
	//return lock
	pthread_mutex_unlock( &mutex );
	//and return instance
	return instance;
}

void MachineInfo::load() {
	//create info stream
	std::ifstream info;
	//open info file
	info.open( "hardware.info" );
	if( !info.is_open() ) {
		std::cerr << "Warning: \tcould not open hardware.info, setting default values for MachineInfo" << std::endl;
		std::cerr << "         \tsetting the number of available cores to 1" << std::endl;
		instance.P = 1;
	} else {
		info >> instance.P; std::cout << "MachineInfo: number of available cores set to " << instance.P << std::endl;
		std::string line;
		getline( info, line );
		getline( info, line );
		while( info.good() ) {
			std::cerr << "MachineInfo: ignoring extra line " << line << std::endl;
			std::string line;
			getline( info, line );
		}
	}
	instance.initialised = true;
}

unsigned long int MachineInfo::cores() const { return P; }

MachineInfo MachineInfo::instance;

bool MachineInfo::initialised = false;

pthread_mutex_t MachineInfo::mutex = PTHREAD_MUTEX_INITIALIZER;

