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


#include <iostream>
#include <vector>
#include <map>
#include <pthread.h>
#include <numa.h>

#include "SparseMatrix.hpp"
#include "Matrix2HilbertCoordinates.hpp"

#include "csb/utility.h"
#include "csb/triple.h"
#include "csb/csc.h"
#include "csb/bicsb.h"

#ifndef _H_RDCSB
#define _H_RDCSB

/*
 * When defined, thread 0 will use the global y vector for its local
 * computations. This introduces extra work in the form of one sync,
 * and the amount of available processors for the parallel collect.
 * The advantage is less data replication of the output vector.
 *
 * The synchronisation cannot be prevented. Using all processors in
 * the collect code is possible but has not been programmed currently.
 */
//#define RDCSB_GLOBAL_Y

/*
 * When defined, RDCSB will not collect output results in the
 * global output vector passed to this library. Used for timing
 * the true SpMV speeds only.
 */
#define RDCSB_NO_COLLECT

/** Shared data for RDCSB threads. */
template< typename T >
class RDCSB_shared_data {

	public:

		unsigned long int id;

		unsigned long int P;

		/** 0 undef, 1 init, 2 zax, 3 zxa, 4 exit */
		unsigned char mode;

		/** how many times to repeat the operation set in `mode' (above, only for 2 and 3) */
		unsigned long int repeat;

		std::vector< Triplet< T > > *original;

		/** Will store rowsums */
		unsigned long int *nzb;

		/** Will store local timing */
		double time;

		pthread_mutex_t* mutex;

		pthread_cond_t*  cond;

		pthread_mutex_t* end_mutex;

		pthread_cond_t*  end_cond;

		unsigned long int *sync, *end_sync;

		unsigned long int output_vector_size;

		unsigned long int output_vector_offset;

		T *local_y;

		/** Base constructor. Will initialise all to invalid values or NULL. */
		RDCSB_shared_data(): id( -1 ), P( -1 ), mode( 0 ), time( 0 ), repeat( 0 ), original( NULL ), nzb( NULL ),
				mutex( NULL ), cond( NULL ), end_mutex( NULL ), end_cond( NULL ),
				sync( NULL ), end_sync( NULL ),
				output_vector_size( -1 ), output_vector_offset( -1 ) {}

		/** Recommended constructor */
		RDCSB_shared_data( unsigned long int _id, unsigned long int _P,
				std::vector< Triplet< double > > *_original,
				unsigned long int *_nzb, 
				pthread_mutex_t *_mutex, pthread_cond_t *_cond, pthread_mutex_t *_end_mutex, pthread_cond_t *_end_cond,
				unsigned long int *_sync, unsigned long int *_end_sync,
				unsigned long int _ovsize, unsigned long int _ovoffset ):
				id( _id ),  P( _P ), mode( 1 ), time( 0 ), repeat( 1 ), original( _original ), nzb( _nzb ),
				mutex( _mutex ), cond( _cond ), end_mutex( _end_mutex ), end_cond( _end_cond ),
				sync( _sync ), end_sync( _end_sync ),
				output_vector_size( _ovsize ), output_vector_offset( _ovoffset ), local_y( NULL ) {}
};

/** Full parallel row-distributed SpMV, based on CSB (BlockCRS + Morton curve + Cilk) and PThreads.
 *  Inspired by Aydin & Gilbert's CSB, and comments by Patrick Amestoy on the BICRS Hilbert scheme.
 *  May not compile due to PThreads/Cilk clashes. */
template< typename T >
class RDCSB: public SparseMatrix< T, ULI > {

	private:

	protected:

		/** Number of threads to fire up */
		static unsigned long int P;

		/** Input vector */
		static const T* input;

		/** Output vector */
		static T* output;

		/** p_threads associated to this data strcuture */
		pthread_t *threads;

		/** array of initial thread data */
		RDCSB_shared_data<T> *thread_data;

		/** Clock type used for thread-local timing */
		static clockid_t global_clock_id;

		/** Stop/continue mechanism: mutex */
		pthread_mutex_t mutex;

		/** Stop/continue mechanism: condition */
		pthread_cond_t cond;

		/** Wait for end mechanism: mutex */
		pthread_mutex_t end_mutex;

		/** Wait for end mechanism: condition */
		pthread_cond_t end_cond;

		/** Used for synchronising threads */
		unsigned long int sync;

		/** Used for construction end signal */
		unsigned long int end_sync;

	public:

		RDCSB( const std::string file, T zero ) {
			loadFromFile( file, zero );
		}

		RDCSB( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero ) {
			load( input, m, n, zero );
		}

		virtual ~RDCSB() {
			//set all daemon threads to exit mode
			for( unsigned long int i=0; i<P; i++ )
				thread_data[ i ].mode = 4;

			//wake up all daemon threads
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//allow threads to exit gracefully
			for( unsigned long int i=0; i<P; i++ )
				pthread_join( threads[ i ], NULL );

			//destroy data
			delete [] thread_data;
			delete [] threads;
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy(  &cond  );
		}

		void wait() {
			//wait for end signal
			pthread_cond_wait( &end_cond, &end_mutex );
			pthread_mutex_unlock( &end_mutex );
		}

		virtual void load( std::vector< Triplet< T > >& input, const ULI m, const ULI n, const T zero ) {
			//get number of cores available
			P = MachineInfo::getInstance().cores();

			//set kernel to local thread allocation if it wasn't already the case
			numa_set_localalloc();

			//base settings
			this->zero_element = zero;
			this->nor = m;
			this->noc = n;
			this->nnz = input.size();

			unsigned long int *nzb = new unsigned long int [ this->m() ];
			for( unsigned long int i=0; i<m; i++ ) nzb[ i ] = 0;

			//create P threads :)
			this->threads = new pthread_t[ P ];
			//initialize local initialisation data
			thread_data = new RDCSB_shared_data<T>[ P ];
			//initialize mutexes and conditions and a synchronisation counter
			pthread_mutex_init( &mutex, NULL );
			pthread_cond_init ( &cond,  NULL );
			pthread_mutex_init( &end_mutex, NULL );
			pthread_cond_init ( &end_cond,  NULL );
			sync     = 0;
			end_sync = 0;
			//lock end mutex (disallow threads that are created to signal for end
			//before this thread is done with spawning children)
			pthread_mutex_lock( &end_mutex );
			//go forth and multiply
			for( unsigned long int i=0; i<P; i++ ) {
				//build thread-local init data
				thread_data[ i ] = RDCSB_shared_data<T>( i, P, &input, nzb, &mutex, &cond, &end_mutex, &end_cond, &sync, &end_sync, -1, -1 );
				//set fixed affinity for threads
				cpu_set_t mask;
				CPU_ZERO( &mask );
				CPU_SET ( i, &mask );

				//TODO: use hwloc for better numa-aware pinning
				/*hwloc_topology_t topology;
				hwloc_topology_init ( &topology );
				hwloc_topology_load( topology );
				hwloc_bitmap_t cpuset;*/

				//prepare attributes
				pthread_attr_t attr;
				pthread_attr_init( &attr );
				//set fixed affinity in attribute, so that it starts binded immediately
				pthread_attr_setaffinity_np( &attr, sizeof( cpu_set_t ), &mask );
				//fire up thread
				pthread_create( &threads[i], &attr, &RDCSB::thread, (void*) &thread_data[i] );
				//free attr
				pthread_attr_destroy( &attr );
			}

			//wait for threads to finish initialisation
			wait();

			//delete temporary array
			delete [] nzb;
		}

		static void end( pthread_mutex_t* mutex, pthread_cond_t* cond, unsigned long int *sync, const unsigned long int P ) {
			pthread_mutex_lock( mutex );
			(*sync)++;
			if( *sync == P ) {
				pthread_cond_signal( cond );
				*sync = 0;
			}
			pthread_mutex_unlock( mutex );
		}

		static void synchronise( pthread_mutex_t* mutex, pthread_cond_t* cond, unsigned long int *sync, const unsigned long int P ) {
			pthread_mutex_lock( mutex );
			(*sync)++;
			if( *sync == P ) {
				*sync = 0;
				pthread_cond_broadcast( cond );
			} else
				pthread_cond_wait( cond, mutex );
			pthread_mutex_unlock( mutex );
		}

		static void* thread( void *data ) {
			//get short-hand notation
			RDCSB_shared_data<T>* shared  = (RDCSB_shared_data<T>*)data;
			const unsigned long int id  = shared->id;
			const unsigned long int P   = shared->P;
			const unsigned long int nnz = shared->original->size();
			pthread_mutex_t *mutex      = shared->mutex;
			pthread_cond_t  *cond       = shared->cond;

			cpu_set_t mask;
			CPU_ZERO( &mask );
			pthread_getaffinity_np( pthread_self(), sizeof( cpu_set_t ), &mask );
			if( !CPU_ISSET( id, &mask ) ) {
				std::cerr << "Incorrect pinning for thread " << id << "!" << std::endl;
				exit( 1 );
			}
			for( unsigned long int s=0; s<P; s++ ) {
				if( s==id ) continue;
				if( CPU_ISSET( s, &mask ) ) {
					std::cerr << "Thread " << id << " mask is larger than one core" << " (" << s << " is set)!" << std::endl;
					exit( 1 );
				}
			}

			//prepare to get global matrix dimensions
			ULI m, n;
			m = n = 0;
			//put rowsums in nzb
			const unsigned long int blocksize = (nnz % P) > 0 ? nnz / P + 1 : nnz / P;
			for( unsigned long int i=0; i<nnz; i++ ) {
				const unsigned long int currow = (*(shared->original))[ i ].i();
				const unsigned long int curcol = (*(shared->original))[ i ].j();
				if( currow >= id * blocksize && currow < (id + 1) * blocksize )
					shared->nzb[ currow ]++;
				if( currow > m ) m = currow;
				if( curcol > n ) n = curcol;
			}
			
			//dimensions are one higher than max indices
			m++;
			n++;

			//sync
			RDCSB::synchronise( mutex, cond, shared->sync, shared->P );

			//determine distribution
			const unsigned long int nnz_target = nnz / P;
			unsigned long int cursum = 0;

			//first sanity check
			for( unsigned long int i=0; i<m; i++ ) cursum += shared->nzb[ i ];
			assert( cursum == nnz );
			
			//continue
			cursum = 0;
			unsigned long int start, end, k = 0;
			start = end = -1;
			if( id == 0 ) start = 0;
			for( unsigned long int i=0; i<m; i++ ) {
				cursum += shared->nzb[ i ];	
				if( cursum >= nnz_target ) {
					if( k == id ) end   = i + 1;
					if(k+1== id ) start = i + 1;
					k++;
					cursum = 0;
				}
			}
			assert( k == P-1 );
			if( id == P-1 ) end = m;
			shared->output_vector_size   = end - start;
			shared->output_vector_offset = start;
			assert (shared->output_vector_size <= m );
			assert( shared->output_vector_offset + shared->output_vector_size <= m );

			//copy to local first			
			std::vector< Triplet< T > > local;
			for( unsigned long int i=0; i<nnz; i++ ) {
				const unsigned long int currow = (*(shared->original))[ i ].i();
				if( currow >= start && currow < end )
					local.push_back(
						Triplet< T >( (*(shared->original))[ i ].i() - start,
							(*(shared->original))[ i ].j(),
							(*(shared->original))[ i ].value )
					);
			}
			const unsigned long int local_nnz = local.size();
			m = shared->output_vector_size; //new matrix size is new m times old n

			//transform to CSB-compatible format
			unsigned int *row = new unsigned int[ local_nnz ];
			unsigned int *col = new unsigned int[ local_nnz ];
			double       *val = new double[ local_nnz ];
			for( unsigned long int i=0; i<local_nnz; i++ ) {
				row[ i ] = local[ i ].i();
				col[ i ] = local[ i ].j();
				val[ i ] = local[ i ].value;
			}
			local.clear();
			
			//load into CSB
			BiCsb< T, unsigned > dss( local_nnz, m, n, row, col, val, 0, 0 );

			//delete transform arrays
			delete [] row;
			delete [] col;
			delete [] val;

			//create local shadow of y to avoid write-contention
			T* y = NULL;
#ifdef RDCSB_GLOBAL_Y
			if( id > 0 ) {
#endif
				y = new T[ shared->output_vector_size ];
				for( unsigned long int i=0; i<shared->output_vector_size; i++ )
					y[ i ] = 0.0;
#ifdef RDCSB_GLOBAL_Y
			}
#endif
			shared->local_y = y;
	
			//exit construction mode
			shared->mode = 0;

			//signal end of construction
			pthread_mutex_lock( mutex );
			RDCSB::end( shared->end_mutex, shared->end_cond, shared->end_sync, shared->P );

			//enter daemon mode
			while( true ) {
				struct timespec clk_start, clk_stop; 
				pthread_cond_wait(  cond, mutex );
				pthread_mutex_unlock( mutex );

				if( shared->mode == 4 ) break;
	
				switch( shared->mode ) {
				case 3:
					assert( RDCSB<T>::input != NULL );
					assert( RDCSB<T>::output != NULL );
#ifdef RDCSB_GLOBAL_Y
					if( id == 0 ) {
						y = RDCSB::output;
						shared->local_y = y;
					}
#endif
					assert( y != NULL );
					std::cerr << "CSB interface to z=x*A has not been implemented!" << std::endl;
					exit( 1 );
#ifndef RDCSB_NO_COLLECT
					collectY( shared );
#endif
					break;
				case 2:
					assert( RDCSB<T>::input != NULL );
					assert( RDCSB<T>::output != NULL );
#ifdef RDCSB_GLOBAL_Y
					if( id == 0 ) {
						y = RDCSB::output;
						shared->local_y = y;
					}
#endif
					assert( y != NULL );

					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
					for( unsigned long int i=0; i<shared->repeat; ++i )
						bicsb_gaxpy( dss, RDCSB<T>::input, y );
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;

#ifndef RDCSB_NO_COLLECT
					collectY( shared );
#endif
					break;
				default:
					std::cout << "Thread " << id << ": Error, undefined operation (" << shared->mode << ")!" << std::endl;
					exit( -1 );
				}
				shared->mode = 0;

				//signal end of operation
				pthread_mutex_lock( mutex );
				RDCSB::end( shared->end_mutex, shared->end_cond, shared->sync, shared->P );
			}

			//done
#ifdef RDCSB_GLOBAL_Y
			if( id != 0 )
#endif
				delete [] y;
			return (NULL);
		}

		static void collectY( RDCSB_shared_data<T> *shared ) {

#ifdef RDCSB_GLOBAL_Y
			//FIXME It could be possible to distribute work over all processors
			//instead of p-1 processors, but this requires some extra balancing.
			const unsigned long int s = shared->id;
			if( s == 0 ) return;
#endif

			//do collect items of own block
			for( unsigned long int i = 0; i < shared->output_vector_size; i++ ) {
				assert( RDCSB<T>::output != NULL );
				assert( shared->local_y != NULL );
				RDCSB<T>::output[ shared->output_vector_offset + i ] += shared->local_y[ i ];
			}
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

		virtual void zxa( const T* x, T* z ) {
			zxa( x, z, 1 );
		}

		virtual void zxa( const T* x, T* z, const unsigned long int repeat ) {
			//set all daemon threads to do zxa
			for( unsigned long int i=0; i<P; i++ ) {
				thread_data[ i ].mode   = 3;
				thread_data[ i ].repeat = repeat;
			}

			//set input vector
			RDCSB<T>::input = x;

			//set output vector
			RDCSB<T>::output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//unset vectors
			RDCSB<T>::input  = NULL;
			RDCSB<T>::output = NULL;
		}

		virtual void zax( const T* x, T* z ) {
			zax( x, z, 1, 0, NULL );
		}

		virtual void zax( const T* x, T* z, const unsigned long int repeat, const clockid_t clock_id, double *elapsed_time ) {
			//set all daemon threads to do zax
			for( unsigned long int i=0; i<P; i++ ) {
				thread_data[ i ].mode   = 2;
				thread_data[ i ].repeat = repeat;
			}

			//set global clock ID
			global_clock_id = clock_id;

			//set input vector
			RDCSB<T>::input = x;

			//set output vector
			RDCSB<T>::output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//get elapsed time
			double maxtime = 0.0;
			for( unsigned long int i=0; i<P; i++ ) {
				const double curtime = thread_data[ i ].time;
				if( curtime > maxtime ) maxtime = curtime;
			}
			if( elapsed_time != NULL )
				*elapsed_time += maxtime;

			//unset vectors
			RDCSB<T>::input  = NULL;
			RDCSB<T>::output = NULL;
		}

		virtual void getFirstIndexPair( ULI &i, ULI &j ) {
			std::cerr << "Warning: RDCSB::getFirstIndexPair has no unique answer since it implements a parallel multiplication!\nIgnoring call..." << std::endl;
		}
};

template< typename T > unsigned long int RDCSB< T >::P = 0;

template< typename T > const T* RDCSB< T >::input  = NULL;

template< typename T > T* RDCSB< T >::output = NULL;

template< typename T > clockid_t RDCSB< T >::global_clock_id = 0;

#endif

