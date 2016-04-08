/***********************************************************************

 Copyright (c) 2013, James G. Dempsey

 All rights reserved.

 This code is a derivative work of Naoya Maruyama (his Copyright
 immediatly follows this Copyright)

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

*************************************************************************/
/*********************************************************

 Copyright (c) 2011-2012, Naoya Maruyama

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.

 * Neither the name of RIKEN AICS nor the names of its contributors may
   be used to endorse or promote products derived from this software
   without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

***********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include <immintrin.h>

#include "HyperThreadPhalanx.h"

#define REAL float
#if !defined(NX)
#define NX 256
#endif
#define NXP nx

// define VERIFY to verify consistency with serial program
// undefine VERIFY for release code
// #define VERIFY


#ifndef M_PI
#define M_PI (3.1415926535897932384626)
#endif

// global variables
int nThreads = -1;	// until initialized
int nCores = -1;	// until initialized
int nHTs = -1;		// until initialized

void init(REAL *buff, const int nx, const int ny, const int nz,
          const REAL kx, const REAL ky, const REAL kz,
          const REAL dx, const REAL dy, const REAL dz,
          const REAL kappa, const REAL time) {
  REAL ax, ay, az;
  int jz, jy, jx;
  ax = exp(-kappa*time*(kx*kx));
  ay = exp(-kappa*time*(ky*ky));
  az = exp(-kappa*time*(kz*kz));
  for (jz = 0; jz < nz; jz++) {
    for (jy = 0; jy < ny; jy++) {
      for (jx = 0; jx < nx; jx++) {
        int j = jz*NXP*ny + jy*NXP + jx;
        REAL x = dx*((REAL)(jx + 0.5));
        REAL y = dy*((REAL)(jy + 0.5));
        REAL z = dz*((REAL)(jz + 0.5));
        REAL f0 = (REAL)0.125
          *(1.0 - ax*cos(kx*x))
          *(1.0 - ay*cos(ky*y))
          *(1.0 - az*cos(kz*z));
        buff[j] = f0;
      }
    }
  }
}

REAL accuracy(const REAL *b1, REAL *b2, const int len, const int count) {
  REAL err = 0.0;
  int i;
  for (i = 0; i < len; i++) {
    err += (b1[i] - b2[i]) * (b1[i] - b2[i]);
  }
  return (REAL)sqrt(err/len) / count;
}


#if defined(VERIFY)
void
diffusion_baseline_verify(REAL *f1, REAL *f2, int nx, int ny, int nz,
                   REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
                   REAL cb, REAL cc) {
printf("Verifying\n");
  REAL largestErr = 0;
  int xErr, yErr, zErr;
  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        int c, w, e, n, s, b, t;
        c =  x + y * nx + z * nx * ny;
        w = (x == 0)    ? c : c - 1;
        e = (x == nx-1) ? c : c + 1;
        n = (y == 0)    ? c : c - nx;
        s = (y == ny-1) ? c : c + nx;
        b = (z == 0)    ? c : c - nx * ny;
        t = (z == nz-1) ? c : c + nx * ny;
        REAL check = cc * f1[c] + cw * f1[w] + ce * f1[e]
            + cs * f1[s] + cn * f1[n] + cb * f1[b] + ct * f1[t];
        REAL err = ((REAL)sqrt((f2[c] - check) * (check - f2[c]))) / check;
        if(err > largestErr)
        {
          largestErr = err;
          xErr = x;
          yErr = y;
          zErr = z;
        }
      }
    }
  }
  if(largestErr != 0.0)
  {
     printf("Largest error (%lf) at x=%d y=%d z=%d\n", (double)largestErr, xErr, yErr, zErr);
  }
  return;
}
#endif

diffusion_tiled(REAL *restrict f1, REAL *restrict f2, int nx, int ny, int nz,
              REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
              REAL cb, REAL cc, REAL dt, int count) {

    // assuming 3 threads/core, each core computing along x in adjacent (z) locations of the z/y plane
    //  C0  C1...
    // <--------- z ------------->
    // ^ abc
    // | abc
    // y abcabc
    // | ...abc
    // | ...abc
    // 
// .or.
    // assuming 4 threads/core, each core computing along x in adjacent (z) locations of the z/y plane
    //  C0  C1...
    // <--------- z ------------->
    // ^ abcd
    // | abcd
    // y abcdabcd
    // | ....abcd
    // | ....abcd
    // 

    // The code uses "Triad" to label the "abc" columns.
// .or.
    // The code uses "Quad" to label the "abcd" columns.
    // 
    // Each column, a for example, drills down 5 adjacent columns of x:
    //
    //   o
    //  oao  (the o's get pulled into cache while processing a)
    //   o
    //
    // By placing a core's Hyper Thread team's abc in adjacent columns
    // we can facillitate a higer L1 cache hit ratio:
    //
    //   ooo
    //  oabco  (the o's get pulled into cache while processing a)
    //   ooo
    //
// .or. 
    //
    //   oooo
    //  oabcdo  (the o's get pulled into cache while processing a)
    //   oooo
    // 
    // What was 5 cache loading columns per thread times 3 .or. 4 threads
    // used per core in original tiled code (15 columns)
    // now becomes 11 columns. Higher cache hit ratio, lower cache footprint.
    // Further, as each thread advanced in the original tiled code
    // 2 of the 5 columns experienced cache hits, 6 hits, 9 misses for three threads
    // The triad layout produces 12 hits, 5 misses for three threads
    // The quad layout produces  16 hits, 6 misses for four threads
  
  // Note, we are not using "parallel for"
  // each thread figures out the partitioning
  #pragma omp parallel
  {
    REAL *f1_t = f1;
    REAL *f2_t = f2;

    // number of Squads (singles/doublets/triplets/quadruples) across z dimension
    int nSquadsZ = (nz + nHTs - 1) / nHTs;
    // number of full (and partial) singles/doublets/triads/quads on z-y face
    int nSquadsZY = nSquadsZ * ny;
    int nSquadsZYPerCore = (nSquadsZY + nCores - 1) / nCores;
    // Determine this thread's triads/quads (TLS init setup myCore and myHT)
    int SquadBegin = nSquadsZYPerCore * myCore;
    int SquadEnd = SquadBegin + nSquadsZYPerCore; // 1 after last single/double/triad/quad for core
    if(SquadEnd > nSquadsZY)
	SquadEnd = nSquadsZY; // truncate if necessary

    // benchmark timing loop
    for (int i = 0; i < count; ++i) {
      // restrict current thread to it's subset of squads on the Z/Y face.
      for(int iSquad = SquadBegin; iSquad < SquadEnd; ++iSquad) {
        int z0 = (iSquad / ny) * nHTs;	// home z for 0'th team member for next squad
        int z = z0 + myHT;		// z for this team member
        int y = iSquad % ny;
        // last double/triad/quad along z may be partially filled
        // assure we are within z
        if(z < nz)
        {
            // determine the center cells and cells about the center
            int x = 0;
            int c, n, s, b, t;
            c =  x + y * nx + z * nx * ny;
            n = (y == 0)    ? c : c - nx;
            s = (y == ny-1) ? c : c + nx;
            b = (z == 0)    ? c : c - nx * ny;
            t = (z == nz-1) ? c : c + nx * ny;
            // c runs through x, n and s through y, b and t through z
            // x=0 special (no f1_t[c-1])
            f2_t[c] = cc * f1_t[c] + cw * f1_t[c] + ce * f1_t[c+1]
                + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
            // interior x's faster
#pragma noprefetch
#pragma simd  
            for (x = 1; x < nx-1; x++) {
              ++c;
              ++n;
              ++s;
              ++b;
              ++t;
              f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c+1]
                  + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
            } // for (x = 1; x < nx-1; x++)
            // final x special (f1_t[c+1])
            ++c;
            ++n;
            ++s;
            ++b;
            ++t;
            f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c]
                + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
        } // if(z < nz)
      } // for(int iSquad = SquadBegin; iSquad < SquadEnd; ++iSquad)
// barrier required because we removed implicit barrrier of #pragma omp for collapse(2)
          #pragma omp barrier
#if defined(VERIFY)
          #pragma omp master
          diffusion_baseline_verify(f1_t, f2_t, nx, ny, nz,
                   ce, cw, cn, cs, ct,
                   cb, cc);
          #pragma omp barrier
#endif

      REAL *t = f1_t;
      f1_t = f2_t;
      f2_t = t;
    } // count
  } // parallel
  return;
}

static double cur_second(void) {
  return omp_get_wtime();
}


void dump_result(REAL *f, int nx, int ny, int nz, char *out_path) {
  FILE *out = fopen(out_path, "w");
  assert(out);
  size_t nitems = nx * ny * nz;
  fwrite(f, sizeof(REAL), nitems, out);
  fclose(out);
}

int main(int argc, char *argv[]) 
{
  if(HyperThreadPhalanxInit())
    return -1;

  nThreads = HyperThreadPhalanx.nThreads;
  nCores = HyperThreadPhalanx.nCores;
  nHTs = HyperThreadPhalanx.nHTs;
  
  double time_begin, time_end;

  int    nx    = NX;
  int    ny    = NX;
  int    nz    = NX;

  REAL *f1 = (REAL *)malloc(sizeof(REAL)*NX*NX*NX);
  REAL *f2 = (REAL *)malloc(sizeof(REAL)*NX*NX*NX);
  REAL *answer = (REAL *)malloc(sizeof(REAL) * NXP*ny*nz);

  assert(f1 != NULL);
  assert(f2 != NULL);
  assert(answer != NULL);
  REAL *f_final = NULL;

  REAL   time  = 0.0;
  int    count = 0;  

  REAL l, dx, dy, dz, kx, ky, kz, kappa, dt;
  REAL ce, cw, cn, cs, ct, cb, cc;

  int    nthreads;
  #pragma omp parallel
  #pragma omp master
    nthreads = omp_get_num_threads();

  l = 1.0;
  kappa = 0.1;
  dx = dy = dz = l / nx;
  kx = ky = kz = 2.0 * M_PI;
  dt = 0.1*dx*dx / kappa;
  // original count computed with NX==256 and used
//  count = 0.1 / dt;
  // This produced a count of 6553
  // adjust count to provide for different value of NX
  // such that runtimes are approximately the same
#if defined(__MIC__)
  count = (6553. * (256.*256.*256.) / ((REAL)nx*(REAL)ny*(REAL)nz)) * nthreads / 240;
#else
  count = ((6553. * (256.*256.*256.) / ((REAL)nx*(REAL)ny*(REAL)nz)) * nthreads / 240) * 5;
#endif
  f_final = (count % 2)? f2 : f1;

  init(f1, nx, ny, nz, kx, ky, kz, dx, dy, dz, kappa, time);

  ce = cw = kappa*dt/(dx*dx);
  cn = cs = kappa*dt/(dy*dy);
  ct = cb = kappa*dt/(dz*dz);
  cc = 1.0 - (ce + cw + cn + cs + ct + cb);

  printf("Running diffusion kernel %d times with %d threads\n", count, nthreads); fflush(stdout);
  time_begin = cur_second();
  diffusion_tiled(f1, f2, nx, ny, nz, ce, cw, cn, cs, ct, cb, cc,
                 dt, count);
  time_end = cur_second();
  time = count * dt;
  dump_result(f_final, nx, ny, nz, "diffusion_result.dat");
  init(answer, nx, ny, nz, kx, ky, kz, dx, dy, dz, kappa, time);
// compute error per iteration
  REAL err = accuracy(f_final, answer, nx*ny*nz, count);
  double elapsed_time = time_end - time_begin;
  REAL mflops = (nx*ny*nz)*13.0*count/elapsed_time * 1.0e-06;
  double thput = (nx * ny * nz) * sizeof(REAL) * 3.0 * count
      / elapsed_time * 1.0e-09;

  fprintf(stderr, "Elapsed time : %.3f (s)\n", elapsed_time);
  fprintf(stderr, "FLOPS        : %.3f (MFlops)\n", mflops);
  fprintf(stderr, "Throughput   : %.3f (GB/s)\n", thput);  
  fprintf(stderr, "Accuracy     : %e\n", err);
  
  free(f1);
  free(f2);
  free(answer);
  return 0;
}

