/* 
Copyright © 2014, Intel Corporation All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

Neither the name of Intel Corporation nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY Intel Corporation "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Intel Corporation BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
This N-Body direct kernel is based on the description from the paper "Test-driving Intel(r) Xeon Phi(tm) coprocessors with a basic N-body simulation" by Andrey Vladimirov and Vadim Karpusenko 
(http://research.colfaxinternational.com/post/2013/01/07/Nbody-Xeon-Phi.aspx)

Alejandro Duran, Intel Corporation
Larry Meadows, Intel Corporation
*/

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <mkl_vsl.h>
#include <omp.h>

#define ALIGN(x,alignment) ((x)+(alignment)-1)/(alignment)*(alignment)

template<class T>
struct PrecisionName
{
    static const char const * name;
};
const char * PrecisionName<float>::name="single";
const char * PrecisionName<double>::name="double";

struct Arch
{
    int l2_size;
    int cacheline_size;
} XeonPhiDesc = { 512*1024, 64 };

template<class T>
struct bf { static const int value; };

const int bf<float>::value = 4096;
const int bf<double>::value = bf<float>::value/2; 

template <class T>
struct ParticleSystem
{
   T *x,*y,*z;
   T *m;
   T *vx,*vy,*vz;   
   
   static const T softening = static_cast<T>(1e-9);
   static const T dt = static_cast<T>(0.01);
   const int nParticles;
   
   int BFI,BFJ;
   
   ParticleSystem(const int np) : nParticles(np) {};
   //~ParticleSystem() {};
   
   void allocate ( Arch &arch );
   void release ( );
      
   void computeTilingFactors(Arch &arch);
   double checksum ( );
   
   void initParticles ( );
};

template<class T>
double ParticleSystem<T>::checksum()
{
    double t = 0.0;
    for (int i = 0; i < nParticles; ++i)
    {
        double xsq = x[i] * x[i];
        double ysq = y[i] * y[i];
        double zsq = z[i] * z[i];
        t += xsq + ysq + zsq;
        double msq = m[i] * m[i];
        t += msq;
    }
    return std::sqrt(t);
}

template<>
void ParticleSystem<float>::initParticles()
{
    VSLStreamStatePtr rnStream; vslNewStream( &rnStream, VSL_BRNG_MT19937, 1 );
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, x, -1.0, 1.0);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, y, -1.0, 1.0);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, z, -1.0, 1.0);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, m, -1.0, 1.0);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, vx, -1.0, 1.0);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, vy, -1.0, 1.0);
    vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, vz, -1.0, 1.0);
}

template<>
void ParticleSystem<double>::initParticles()
{
    VSLStreamStatePtr rnStream; vslNewStream( &rnStream, VSL_BRNG_MT19937, 1 );
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, x, -1.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, y, -1.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, z, -1.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, m, -1.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, vx, -1.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, vy, -1.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnStream, nParticles, vz, -1.0, 1.0);
}


template<class T>
void ParticleSystem<T>::allocate(Arch &arch) 
{        
    int alignment = arch.cacheline_size/sizeof(T);
    int offset = ALIGN(nParticles,alignment);
    int bytes = offset * sizeof(T) * 7;        

    T *p = (T *) _mm_malloc(bytes, 64);
    
    x =  p + 0 * offset;
    y =  p + 1 * offset;
    z =  p + 2 * offset;
    m  = p + 3 * offset;
    vx = p + 4 * offset;
    vy = p + 5 * offset;
    vz = p + 6 * offset;

    initParticles();
}

template<class T >
void ParticleSystem<T>::release ()
{
    _mm_free(x);
}

template <class T>
void ParticleSystem<T>::computeTilingFactors(Arch &arch)
{
    int nthreads = omp_get_max_threads();
    int alignment = arch.cacheline_size/sizeof(T);
        
    int BFI_unaligned = (nParticles+nthreads-1)/nthreads;   
    BFI = ALIGN(BFI_unaligned, alignment);

    int iters = nParticles / BFI;
    if ( nParticles % BFI != 0 ) iters++;

    float balance = ((float) iters) / nthreads;
    if ( balance < 0.98 ) BFI = BFI_unaligned;

    BFJ = bf<T>::value;    

    float l2_ocupancy = ((float)nParticles)*sizeof(T)*4/ arch.l2_size;
    if ( l2_ocupancy < 4 ) BFJ *= 2;
    if ( l2_ocupancy < 2 ) BFJ *= 2;

    std::cout << "Block size for i-loop: " << BFI << std::endl
              << "Block size for j-loop: " << BFJ << std::endl;
}


template <class T>
void computeForces(ParticleSystem<T> p, const T dt) 
{
    const int n = p.nParticles;
    const int BFI = p.BFI;
    const int BFJ = p.BFJ;
        
    #pragma omp for
    for ( int ii = 0; ii < n ; ii += BFI ) {
        int imax = ii + BFI;
        if (imax > n ) imax = n;
        for ( int jj = 0; jj < n; jj += BFJ ) { 	   
            int jmax = jj + BFJ;
            if (jmax > n) jmax = n;
            for (int i = ii; i < imax ; i++ ) {
                T Fx = static_cast<T>(0.0); 
                T Fy = static_cast<T>(0.0); 
                T Fz = static_cast<T>(0.0);

                #pragma vector aligned
                for (int j = jj; j < jmax; j++) {
                    const T dx = p.x[j] - p.x[i];
                    const T dy = p.y[j] - p.y[i];
                    const T dz = p.z[j] - p.z[i];

                    const T drSquared = dx*dx + dy*dy + dz*dz + p.softening;
                    const T drPowerN12 = 1.0f / sqrtf(drSquared);
                    const T drPowerN32 = drPowerN12 * drPowerN12 * drPowerN12;        
                    const T s = p.m[j] * drPowerN32;

                    Fx += dx * s;
                    Fy += dy * s; 
                    Fz += dz * s;
                }
                p.vx[i] += dt*Fx;
                p.vy[i] += dt*Fy;
                p.vz[i] += dt*Fz;         
            }
        }
	}
}

template<class T>
void run ( ParticleSystem<T> p, const int nIters, const T dt ) 
{
    double start,end;
     
    #pragma omp parallel
    for (int iter = 0; iter < nIters; iter++)
    {
        computeForces(p,dt);        
        #pragma omp for
        for (int i = 0; i < p.nParticles; i++) {	    
            p.x[i] += p.vx[i] * dt;
            p.y[i] += p.vy[i] * dt;
            p.z[i] += p.vz[i] * dt;
	    }
    }
}

template< class T >
void run_experiment (int nParticles, int nIters)
{
    const char *precision_name = PrecisionName<T>::name;
    
    std::cout << "++++" << precision_name << " precision run ++++" << std::endl;    

    const T dt = static_cast<T>(0.01);
    ParticleSystem<T> p(nParticles);
    
    p.allocate(XeonPhiDesc);
    p.computeTilingFactors(XeonPhiDesc);
    
    /* warm-up run */
    run(p,1,dt);    
    
    /* real run */
    double start=omp_get_wtime();
    run(p,nIters,dt);    
    double end=omp_get_wtime();
    
    double avgTime = (end-start) / (double) (nIters);    
    
    const double flop_ratio = 20;
    const double g_interactions = 1e-9 * nParticles * nParticles / avgTime;
    
    std::cout << "Average iteration time: " << avgTime << " s." << std::endl
              << "Performance: " << g_interactions << " billion interactions per second" << std::endl
              << "GFLOPS: " << g_interactions * flop_ratio << std::endl
              << "Checksum: " << p.checksum() << std::endl;
    
    p.release();
    
    std::cout << "++++" << precision_name << " precision run ends ++++" << std::endl;    
}

void usage (char *pname)
{
    std::cout << pname << " [options]" << std::endl
             << "Options are:" << std::endl
             << "-s              Run only single precision version" << std::endl
             << "-d              Run only double precision version" << std::endl
             << "-i iters        Number of time steps to simulate" << std::endl
             << "-n particles    Number of particles to simulate" << std::endl
    ;
}

int main(const int argc, char * const *argv)
{
    bool run_both=true, run_single = false, run_double = false;    
    int nParticles = 30000;
    int nIters = 10;		

    char opt;
    while ((opt = getopt(argc, argv, "sdn:i:h")) != -1 ) {
        switch (opt) {
            case 's':   run_both = false; run_single = true;
                        break;
            case 'd':   run_both = false; run_double = true;
                        break;
            case 'n':   nParticles = atoi(optarg);
                        break;
            case 'i':   nIters = atoi(optarg);
                        break;
            default:
                        usage(argv[0]);
                        return 1;
        }
    }
    
    if ( optind < argc ) {
        usage(argv[0]);
        return 1;
    }
    
    std::cout << "Number of Particles: " << nParticles << std::endl
              << "Number of Iterations: " << nIters << std::endl
              << "OpenMP Threads: " << omp_get_max_threads() << std::endl;
        
    if (run_both || run_single) run_experiment<float>(nParticles,nIters);
    if (run_both || run_double) run_experiment<double>(nParticles,nIters);
}
