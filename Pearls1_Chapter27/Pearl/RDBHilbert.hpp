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

#ifndef _NO_LIBNUMA
 #include <numa.h>
#endif

//compiling with sniper will instrument a call to a repeated SpMV (zax)
#ifdef _WITH_SNIPER
 #include "sim_api.h"
#endif

//#include "hwloc.h"

#include "SparseMatrix.hpp"
#include "Matrix2HilbertCoordinates.hpp"
#include "FBICRS.hpp"
#include "MachineInfo.hpp"
#include "BigInt.hpp"

#ifndef _H_RDBHILBERT
#define _H_RDBHILBERT

/*
 * When defined, thread 0 will use the global y vector for its local
 * computations. This introduces extra work in the form of one sync,
 * and the amount of available processors for the parallel collect.
 * The advantage is less data replication of the output vector.
 *
 * The synchronisation cannot be prevented. Using all processors in
 * the collect code is possible but has not been programmed currently.
 */
#define RDBH_GLOBAL_Y

/*
 * When defined, RDBHilbert will not collect output results in the
 * global output vector passed to this library. Used for timing
 * the true SpMV speeds only.
 */
#define RDBH_NO_COLLECT

/** Shared data for RDBHilbert threads. */
template< typename T >
class RDB_shared_data {

	public:

		unsigned long int id;

		unsigned long int P;

		/** 0 undef, 1 init, 2 zax, 3 zxa, 4 exit, 5 reset */
		unsigned char mode;

		/** how many times to repeat the operation set in `mode' (above, only for 2 and 3) */
		unsigned long int repeat;

		std::vector< Triplet< T > > *original;

		/** Will store rowsums */
		unsigned long int *nzb;

		/** Will store local timing */
		double time;

		/** Will store the local amount of bytes used */
		size_t bytes;

		/** Will store the total fillIn at this thread */
		size_t fillIn;

		pthread_mutex_t* mutex;

		pthread_cond_t*  cond;

		pthread_mutex_t* end_mutex;

		pthread_cond_t*  end_cond;

		unsigned long int *sync, *end_sync;

		unsigned long int output_vector_size;

		unsigned long int output_vector_offset;

		T *local_y;

		const T ** input;

		T ** output;

		/** Base constructor. Will initialise all to invalid values or NULL. */
		RDB_shared_data(): id( -1 ), P( -1 ), mode( 0 ), repeat( 0 ), original( NULL ), nzb( NULL ), time( 0 ),
				mutex( NULL ), cond( NULL ), end_mutex( NULL ), end_cond( NULL ),
				sync( NULL ), end_sync( NULL ),
				output_vector_size( -1 ), output_vector_offset( -1 ),
				local_y( NULL ), input( NULL ), output( NULL ) {}

		/** Recommended constructor */
		RDB_shared_data( unsigned long int _id, unsigned long int _P,
				std::vector< Triplet< double > > *_original,
				unsigned long int *_nzb,
				pthread_mutex_t *_mutex, pthread_cond_t *_cond, pthread_mutex_t *_end_mutex, pthread_cond_t *_end_cond,
				unsigned long int *_sync, unsigned long int *_end_sync,
				unsigned long int _ovsize, unsigned long int _ovoffset,
				const T **_in, T **_out):
				id( _id ),  P( _P ), mode( 1 ), repeat( 1 ), original( _original ), nzb( _nzb ), time( 0 ),
				mutex( _mutex ), cond( _cond ), end_mutex( _end_mutex ), end_cond( _end_cond ),
				sync( _sync ), end_sync( _end_sync ),
				output_vector_size( _ovsize ), output_vector_offset( _ovoffset ),
				local_y( NULL ), input( _in ), output( _out ) {}
};

/** Compare function used for quicksort on an array of unsigned long ints */
int rdbh_uli_compare( const void *a, const void *b ) {
	return ( *(unsigned long int*)a - *(unsigned long int*)b );
}

/** The Beta Hilbert triplet scheme. Full parallel SpMV, based on Blocked Hilbert and PThreads.
 *  Inspired by Aydin & Gilbert's CSB, and comments by Patrick Amestoy on the BICRS Hilbert scheme. */
template< typename T, class MatrixType = FBICRS< T > >
class RDBHilbert: public SparseMatrix< T, ULI > {

	private:

	protected:

		/** Number of threads to fire up */
		static unsigned long int P;

		/** Which processors to pin threads to */
		std::vector< unsigned short int > p_translate;

		/** Input vector */
		const T* input;

		/** Output vector */
		T *output;

		/** p_threads associated to this data strcuture */
		pthread_t *threads;

		/** array of initial thread data */
		RDB_shared_data<T> *thread_data;

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

		/** Given FBICRS, the maximum value for columnwise matrix size */
		static const ULI max_n = FBICRS< T >::beta_n;

		/** Given FBICRS, the maximum value for the rowwise matrix size, assuming short ints on ICRS at the lower level */
		static const ULI max_m = FBICRS< T >::beta_m;

		/** Sets p_translate to 0..P-1 by default, or equal to the optionally supplied vector. */
		void set_p_translate( std::vector< unsigned short int > *_p_translate ) {
			if( _p_translate == NULL ) {
				//get number of cores available
				P = MachineInfo::getInstance().cores();
				for( unsigned short int i = 0; i < P; ++i ) {
#ifdef __MIC
					const unsigned short int num_threads_times_two = 120; //114 for the 57-core MIC
					//offset controls wether we are dividing over the first
					//2 HW threads or the latter two
					size_t offset = 0;
					//assume the number of cores is i div 2
					size_t core   = i / 2;
					//if i >= 120 we go for the latter two HW thread on each core
					if( i >= num_threads_times_two ) {
						//map to the same cores as for i<120
						core   = (i-num_threads_times_two) / 2;
						//but at the higher 2 HW threads
						offset = 2;
					}
					//assume the thread number on-core is i mod 2 (plus offset)
					const size_t hwthread = i % 2 + offset;
					//assume consecutively wrapped HW threads, 4 per core
					p_translate.push_back( 4 * core + hwthread );
#else
					//just use consecutive numbering otherwise
					p_translate.push_back( i );
#endif
				}
			} else {
				p_translate = *_p_translate;
				P = p_translate.size();
			}
		}

	public:

		RDBHilbert( const std::string file, T zero = 0, std::vector< unsigned short int > *_p_translate = NULL ) {
			set_p_translate( _p_translate );
			this->loadFromFile( file, zero );
		}

		RDBHilbert( std::vector< Triplet< T > >& input, ULI m, ULI n, T zero = 0, std::vector< unsigned short int > *_p_translate = NULL ) {
			set_p_translate( _p_translate );
			load( input, m, n, zero );
		}

		virtual ~RDBHilbert() {
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
#ifndef _NO_LIBNUMA
			//set kernel to local thread allocation if it wasn't already the case
			numa_set_localalloc();
#endif

			//base settings
			this->zero_element = zero;
			this->nor = m;
			this->noc = n;
			this->nnz = input.size();

			// Parallel distribution according to the Hilbert curve:
			// -get maximum Hilbert coordinate
			//  (Note: curve coordinates should be implented such that (0,0) maps to 0!
			//         Otherwise a minimum Hilbert coordinate should also be found, and
			//         below code must be adapted to work on a Hilbert coordinate range
			//         instead)
			// -simultaneously also get number of nonzeroes to locally store
			// -translate to a range
			// -call function to distribute this range over P threads
			//
			// Most operations are done in parallel (see thread function), basically
			// only shared array allocation and the actual forking is in this function.

			unsigned long int *nzb = new unsigned long int [ this->m() ];
			for( unsigned long int i=0; i<m; i++ ) nzb[ i ] = 0;

			//create P threads :)
			this->threads = new pthread_t[ P ];
			//initialize local initialisation data
			thread_data = new RDB_shared_data<T>[ P ];
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
				thread_data[ i ] = RDB_shared_data<T>( i, P, &input, nzb, &mutex, &cond, &end_mutex, &end_cond, &sync, &end_sync, -1, -1, &(this->input), &output );
				//set fixed affinity for threads
				cpu_set_t mask;
				CPU_ZERO( &mask );
				CPU_SET ( p_translate[ i ], &mask );

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
				pthread_create( &threads[i], &attr, &RDBHilbert::thread, (void*) &thread_data[i] );
				//free attr
				pthread_attr_destroy( &attr );
			}

			//wait for threads to finish initialisation
			wait();

			//report on fillIn
			size_t fillIn = 0;
			for( unsigned long int i=0; i<P; i++ ) {
				fillIn += thread_data[ i ].fillIn;
			}
			std::cout << "RDBHilbert: fill-in is " << fillIn << ", which is +" << 100.0*((double)fillIn)/((double)input.size()) << "%" << std::endl;

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
			RDB_shared_data<T>* shared  = (RDB_shared_data<T>*)data;
			const unsigned long int id  = shared->id;
			const unsigned long int P   = shared->P;
			const unsigned long int nnz = shared->original->size();
			pthread_mutex_t *mutex      = shared->mutex;
			pthread_cond_t  *cond       = shared->cond;

			//check if pinned to exactly one thread
			cpu_set_t mask;
			CPU_ZERO( &mask );
			pthread_getaffinity_np( pthread_self(), sizeof( cpu_set_t ), &mask );
			unsigned long int MyMask = ULONG_MAX;
			for( unsigned long int s=0; s<P; s++ ) {
				if( CPU_ISSET( s, &mask ) ) {
					if( MyMask == ULONG_MAX )
						MyMask = s;
					else {
						std::cerr << "Thread " << id << " mask is larger than one core" << " (" << MyMask << " and " << s << " are both set)!" << std::endl;
						exit( 1 );
					}
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
			RDBHilbert::synchronise( mutex, cond, shared->sync, shared->P );

			//determine distribution
			const unsigned long int nnz_target = nnz / P;
			unsigned long int cursum = 0;

			//first sanity check
			for( unsigned long int i=0; i<m; i++ ) cursum += shared->nzb[ i ];
			assert( cursum == nnz );
			
			//continue
			cursum = 0;
			unsigned long int start, end, k = 0;
			start = end = ULONG_MAX;
			if( id == 0 ) start = 0;
#ifndef NDEBUG
			if( id == 0 ) std::cout << "Breakpoints: ";
#endif
			for( unsigned long int i=0; i<m; i++ ) {
				cursum += shared->nzb[ i ];	
				if( cursum >= nnz_target ) {
#ifndef NDEBUG
					if( id == 0 ) std::cout << (i+1) << " ";
#endif
					if( k == id ) end   = i + 1;
					if(k+1== id ) start = i + 1;
					k++;
					cursum = 0;
				}
			}
			if( start == ULONG_MAX ) start = m;
			if( id == P-1 || end == ULONG_MAX ) end = m;
#ifndef NDEBUG
			if( id == 0 ) std::cout << std::endl;
			RDBHilbert::synchronise( mutex, cond, shared->sync, shared->P );
			std::cout << "Thread " << id << " has rows " << start << " to " << end << std::endl;
#endif
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

			//extract hilbert coordinate info
			unsigned long int h1, h2;
			BigInt maxh( 0, 0 ), minh( ULONG_MAX, ULONG_MAX );
			BigInt *t2h = new BigInt[ local_nnz ];
			for( unsigned long int i=0; i<local_nnz; i++ ) {
				Matrix2HilbertCoordinates::IntegerToHilbert( local[i].i() / max_m, local[i].j() / max_n, h1, h2 );
				const BigInt cur( h1, h2 );
				t2h[ i ] = cur;
				if( maxh < cur ) maxh = cur;
				if( minh > cur ) minh = cur;			
				assert( t2h[ i ].high == h1 );
				assert( t2h[ i ].low  == h2 );
			}

			//sequential quicksort O(nz*log(nz))
			qsort(  t2h, local_nnz, sizeof( BigInt ), bigint_compare );
			assert( local_nnz == 0 || t2h[ local_nnz - 1 ] == maxh );
			assert( local_nnz == 0 || t2h[ 0 ] == minh );
				
			//build local data structure from the triplets in shared->original
			//and beta block hilbert index to beta translation map
			std::map< BigInt, unsigned long int > h2b;
			std::vector< std::vector< Triplet< T > > > beta;

			//guard against invalid accesses due to overpartitioning
			if( local_nnz > 0 ) {
				//prepare beta and map
				//get first hilbert block coordinate
				BigInt cur = t2h[ 0 ];
				//remember index
				h2b[ cur ] = 0;
				//keep next index
				unsigned long int c = 0, size = 1;
				//remember previous index
				BigInt prev_h = cur;
				//do the same for the remainder of available hilbert coordinates
				for( unsigned long int i = 1; i < local_nnz; ++i, ++size ) {
					cur = t2h[ i ];
					if( cur != prev_h ) { //new coordinate
						//push back old vector with exact size
						beta.push_back( std::vector< Triplet< T > >() );
						beta.back().reserve( size );
						//store previous index
						h2b[ prev_h ] = c++;
						//reset
						prev_h = cur;
						size = 1;
					}
				}
				//push back last one
				beta.push_back( std::vector< Triplet< T > >() );
				beta.back().reserve( size );
				h2b[ cur ] = c;
			}
					
			//prepare to get matrix size (m, n) and some statistics too
			unsigned long int smin_m = ULONG_MAX;
			unsigned long int smin_n = ULONG_MAX;
			unsigned long int smax_m, smax_n;
			smax_m = smax_n = 0;
			unsigned long int *ms    = new unsigned long int[ h2b.size() ];
			unsigned long int *ns    = new unsigned long int[ h2b.size() ];
			unsigned long int *minms = new unsigned long int[ h2b.size() ];
			unsigned long int *minns = new unsigned long int[ h2b.size() ];
			for( unsigned long int i=0; i<h2b.size(); i++ ) {
				ms[ i ] = 0;
				ns[ i ] = 0;
				minms[ i ] = ULONG_MAX;
				minns[ i ] = ULONG_MAX;
			}

			//fill beta blocks in correct hilbert order O(local_nnz * log( nonzero_blocks ))
			for( unsigned long int i=0; i<local_nnz; i++ ) {
				const unsigned long int row = local[ i ].i();
				const unsigned long int col = local[ i ].j();
				Matrix2HilbertCoordinates::IntegerToHilbert( row / max_m, col / max_n, h1, h2 );
				const BigInt cur( h1, h2 );
				assert( cur >= minh );
				assert( cur <= maxh );
				assert( row <  m );
				const Triplet< T > toAdd = Triplet< T >( row, col, local[ i ].value );
#ifndef NDEBUG
				if( beta[ h2b[ cur ] ].size() > 1 ) {
					Matrix2HilbertCoordinates::IntegerToHilbert(
						beta[ h2b[ cur ] ].back().i() / max_m,
						beta[ h2b[ cur ] ].back().j() / max_n, h1, h2 );
					assert( cur == BigInt( h1, h2 ) );
					assert( beta[ h2b[ cur ] ].back().i() / max_m == row / max_m );
					assert( beta[ h2b[ cur ] ].back().j() / max_n == col / max_n );
				}
#endif
				beta[ h2b[ cur ] ].push_back( toAdd );
				if( row > ms[ h2b[ cur ] ] ) ms[ h2b[ cur ] ] = row;
				if( col > ns[ h2b[ cur ] ] ) ns[ h2b[ cur ] ] = col;
				if( row < minms[ h2b[ cur ] ] ) minms[ h2b[ cur ] ] = row;
				if( col < minns[ h2b[ cur ] ] ) minns[ h2b[ cur ] ] = col;
				if( row < smin_m ) smin_m = row;
				if( col < smin_n ) smin_n = col;
				if( row > smax_m ) smax_m = row;
				if( col > smax_n ) smax_n = col;
			}
			//size is max value + 1
			smax_m++; smax_n++;

#ifndef NDEBUG
			cursum = 0;
			for( unsigned long int i=0; i<beta.size(); i++ )
				cursum += beta[i].size();
			assert( cursum == local_nnz );
			std::cout << "Thread " << shared->id << ": " << smin_m << "," << smax_m << " times " << smin_n << "," << smax_n << " holding " << cursum << " nonzeroes." << std::endl;

			//sanity check: everything in one beta vector should share the
			//same block.
			for( unsigned long int i=0; i<beta.size(); ++i ) {
				const unsigned long int row_br = beta[i][0].i() / max_m;
				const unsigned long int col_br = beta[i][0].j() / max_n;
				for( unsigned long int k=1; k<beta[i].size(); ++k ) {
					assert( beta[i][k].i() / max_m == row_br );
					assert( beta[i][k].j() / max_n == col_br );
				}
			}
#endif

			//remove temporary values
			delete [] ms;
			delete [] ns;
			delete [] minms;
			delete [] minns;
			delete [] t2h;

			//load into FBICRS
			MatrixType dss( beta, smax_m - smin_m, n );
			beta.clear();

			//remember memory use
			shared->bytes  = dss.bytesUsed();
			shared->fillIn = dss.fillIn;

			//create local shadow of y to avoid write-contention
			T* y = NULL;
#ifdef RDBH_GLOBAL_Y
			if( id > 0 ) {
#endif
				y = new T[ shared->output_vector_size ];
				for( unsigned long int i=0; i<shared->output_vector_size; i++ )
					y[ i ] = 0.0;
#ifdef RDBH_GLOBAL_Y
			}
#endif
			shared->local_y = y;
	
			//exit construction mode
			shared->mode = 0;

			//signal end of construction
			pthread_mutex_lock( mutex );
			RDBHilbert::end( shared->end_mutex, shared->end_cond, shared->end_sync, shared->P );

			//enter daemon mode
			while( true ) {
				struct timespec clk_start, clk_stop; 
				pthread_cond_wait(  cond, mutex );
				pthread_mutex_unlock( mutex );

				if( shared->mode == 4 ) break;

				switch( shared->mode ) {
				case 5:
#ifdef RDBH_GLOBAL_Y
					if( shared->id != 0) {
#endif
						for( size_t i = 0; i < shared->output_vector_size; ++i ) {
							shared->local_y[ i ] = 0;
						}
#ifdef RDBH_GLOBAL_Y
					}
#endif
					break;
				case 3:
					assert( *(shared->input) != NULL );
					assert( *(shared->output) != NULL );
#ifdef RDBH_GLOBAL_Y
					if( id == 0 ) {
						y = *(shared->output);
						shared->local_y = y;
					}
#endif
					assert( y != NULL );

std::cout << "Here with " << shared->repeat<<std::endl;
					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
					for( unsigned long int i=0; i<shared->repeat; i++ ) {
						dss.zxa( *(shared->input), y );
					}
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;

#ifndef RDBH_NO_COLLECT
					collectY( shared );
#endif
					break;
				case 2:
					assert( *(shared->input) != NULL );
					assert( *(shared->output) != NULL );
#ifdef RDBH_GLOBAL_Y
					if( id == 0 ) {
						y = *(shared->output);
						shared->local_y = y;
					}
#endif
					assert( y != NULL );

					clock_gettime( global_clock_id, &clk_start);
					shared->time = 0.0;
#ifdef _WITH_SNIPER
					parmacs_roi_begin();
#endif
					for( unsigned long int i=0; i<shared->repeat; i++ ) {
						dss.zax( *(shared->input), y );
					}
#ifdef _WITH_SNIPER
					parmacs_roi_end();
#endif
					clock_gettime( global_clock_id, &clk_stop);
					shared->time  = (clk_stop.tv_sec-clk_start.tv_sec)*1000;
					shared->time += (clk_stop.tv_nsec-clk_start.tv_nsec)/1000000.0;

#ifndef RDBH_NO_COLLECT
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
				RDBHilbert::end( shared->end_mutex, shared->end_cond, shared->sync, shared->P );
			}

			//done
#ifdef RDBH_GLOBAL_Y
			if( id != 0 )
#endif
				delete [] y;
			return (NULL);
		}

		static void collectY( RDB_shared_data<T> *shared ) {

#ifdef RDBH_GLOBAL_Y
			//FIXME It could be possible to distribute work over all processors
			//instead of p-1 processors, but this requires some extra balancing.
			const unsigned long int s = shared->id;
			if( s == 0 ) return;
#endif

			//do collect items of own block
			for( unsigned long int i = 0; i < shared->output_vector_size; i++ ) {
				assert( *(shared->output) != NULL );
				assert( shared->local_y != NULL );
				(*(shared->output))[ shared->output_vector_offset + i ] += shared->local_y[ i ];
			}
		}

		/** Overloaded mv call; allocates output vector using numa_interleaved. */
		virtual T* mv( const T* x ) {
			//over-allocate in case we use matrices that are too small to be
			//vectorised internally (assuming we are using oICRS as sub_ds).
			size_t allocsize = (this->nor + 1) * sizeof( T );
			//allocate, either using libnuma (to interleave for i86)
			//or using dynamic aligned allocs (for Xeon Phi MICs)
#ifndef _NO_LIBNUMA
			T* ret = (T*) numa_alloc_interleaved( allocsize );
#else
			T* ret = (T*) _mm_malloc( allocsize, 64 ); //instead of `new T[ nor ];'
#endif
			//set to 0-vector
			for( ULI i=0; i<this->nor; i++ ) {
				ret[ i ] = this->zero_element;
			}

			//reset output vectors
			reset();

			//do in-place SpMV
			zax( x, ret );

			//return new output vector
			return ret;
		}

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
			input = x;

			//set output vector
			output = z;

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();

			//unset vectors
			input  = NULL;
			output = NULL;
		}

		void reset() {
			//reset all local output vectors
			for( unsigned long int i=0; i<P; i++ ) {
				thread_data[ i ].mode = 5;
			}

			//wake up all daemon threads
			pthread_mutex_lock( &end_mutex );
			pthread_mutex_lock( &mutex );
			pthread_cond_broadcast( &cond );
			pthread_mutex_unlock( &mutex );

			//wait for end of operation
			wait();
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
			input = x;

			//set output vector
			output = z;

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
			input  = NULL;
			output = NULL;
		}

		virtual void getFirstIndexPair( ULI &i, ULI &j ) {
			std::cerr << "Warning: RDBHilbert::getFirstIndexPair has no unique answer since it implements a parallel multiplication!\nIgnoring call..." << std::endl;
		}

		virtual size_t bytesUsed() {
			size_t ret = 0;
			for( size_t s = 0; s < P; ++s )
				ret += thread_data[ s ].bytes;
			return ret;
		}

};

template< typename T, class MatrixType > unsigned long int RDBHilbert< T, MatrixType >::P = 0;

template< typename T, class MatrixType > clockid_t RDBHilbert< T, MatrixType >::global_clock_id = 0;

#endif

