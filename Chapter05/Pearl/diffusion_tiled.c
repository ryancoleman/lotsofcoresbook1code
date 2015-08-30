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
/*
Minor edits made to this program by Jim Dempsey of QuickThread Programming, LLC

	1) Change to use OpenMP timing function
	2) Change to permit NX to be optionally passed in as a compiler line -DNX=nnnn
	3) Add assert to assure array answer is allocated
	4) Add free(answer) at end of program
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <assert.h>

#define REAL float
#if !defined(NX)
#define NX 256
#endif

#ifndef M_PI
#define M_PI (3.1415926535897932384626)
#endif

// #define DIAGNOSTICS to enable diagnostic printouts
//#define DIAGNOSTICS

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
        int j = jz*nx*ny + jy*nx + jx;
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

diffusion_tiled(REAL *restrict f1, REAL *restrict f2, int nx, int ny, int nz,
              REAL ce, REAL cw, REAL cn, REAL cs, REAL ct,
              REAL cb, REAL cc, REAL dt, int count) {

  unsigned long tsc;
  int nthreads;

#pragma omp parallel
  {
    REAL *f1_t = f1;
    REAL *f2_t = f2;
    int mythread;

    for (int i = 0; i < count; ++i) {
#define YBF 16
#pragma omp for collapse(2)
      for (int yy = 0; yy < ny; yy += YBF) {
      for (int z = 0; z < nz; z++) {
        int ymax = yy + YBF;
        if (ymax >= ny) ymax = ny;
        for (int y = yy; y < ymax; y++) {
          int x;
          int c, n, s, b, t;
          x = 0;
          c =  x + y * nx + z * nx * ny;
          n = (y == 0)    ? c : c - nx;
          s = (y == ny-1) ? c : c + nx;
          b = (z == 0)    ? c : c - nx * ny;
          t = (z == nz-1) ? c : c + nx * ny;
          f2_t[c] = cc * f1_t[c] + cw * f1_t[c] + ce * f1_t[c+1]
              + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
#pragma simd  
          for (x = 1; x < nx-1; x++) {
            ++c;
            ++n;
            ++s;
            ++b;
            ++t;
            f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c+1]
                + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
          }
          ++c;
          ++n;
          ++s;
          ++b;
          ++t;
          f2_t[c] = cc * f1_t[c] + cw * f1_t[c-1] + ce * f1_t[c]
              + cs * f1_t[s] + cn * f1_t[n] + cb * f1_t[b] + ct * f1_t[t];
        } // tile ny
      } // tile nz
      } // block ny
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
  
  double time_begin, time_end;
  // requirement of cubic dimension
  int    nx    = NX;
  int    ny    = NX;
  int    nz    = NX;
  assert((ny==nx) && (nz==nx)); // requirement

  REAL *f1 = (REAL *)malloc(sizeof(REAL)*nx*ny*nz);
  REAL *f2 = (REAL *)malloc(sizeof(REAL)*nx*ny*nz);
  REAL *answer = (REAL *)malloc(sizeof(REAL) * nx*ny*nz);
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
