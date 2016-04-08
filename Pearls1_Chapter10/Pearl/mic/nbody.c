/*
Copyright (c) 2014, Rio Yokota, Mustafa AbdulJabbar
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <immintrin.h>

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)(tv.tv_sec+tv.tv_usec*1e-6);
}

int main() {
  // Initialize
  int N = 1 << 16;
  int NALIGN = 64;
  int i, j;
  float OPS = 20. * N * N * 1e-9;
  float EPS2 = 1e-6;
  double tic, toc;
  float * x = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * y = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * z = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * m = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * p = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * ax = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * ay = (float*) _mm_malloc(N * sizeof(float), NALIGN);
  float * az = (float*) _mm_malloc(N * sizeof(float), NALIGN);
#pragma omp parallel for
  for (i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    z[i] = drand48();
    m[i] = drand48() / N;
    p[i] = ax[i] = ay[i] = az[i] = 0;
  }
  printf("N : %d\n",N);

#pragma omp parallel private(j)
  {
#pragma omp single
    tic = get_time();
    // Vectorize target with intrinsics
#pragma omp for
    for (i=0; i<N; i+=16) {
      __m512 pi = _mm512_setzero_ps();
      __m512 axi = _mm512_setzero_ps();
      __m512 ayi = _mm512_setzero_ps();
      __m512 azi = _mm512_setzero_ps();
      __m512 xi = _mm512_load_ps(x+i);
      __m512 yi = _mm512_load_ps(y+i);
      __m512 zi = _mm512_load_ps(z+i);
      for (j=0; j<N; j++) {
	__m512 xj = _mm512_set1_ps(x[j]);
	xj = _mm512_sub_ps(xj, xi);
	__m512 yj = _mm512_set1_ps(y[j]);
	yj = _mm512_sub_ps(yj, yi);
	__m512 zj = _mm512_set1_ps(z[j]);
	zj = _mm512_sub_ps(zj, zi);
	__m512 R2 = _mm512_set1_ps(EPS2);
	R2 = _mm512_fmadd_ps(xj, xj, R2);
	R2 = _mm512_fmadd_ps(yj, yj, R2);
	R2 = _mm512_fmadd_ps(zj, zj, R2);
	__m512 mj = _mm512_set1_ps(m[j]);
	__m512 invR = _mm512_rsqrt23_ps(R2);
	mj = _mm512_mul_ps(mj, invR);
	pi = _mm512_add_ps(pi, mj);
	invR = _mm512_mul_ps(invR, invR);
	invR = _mm512_mul_ps(invR, mj);
	axi = _mm512_fmadd_ps(xj, invR, axi);
	ayi = _mm512_fmadd_ps(yj, invR, ayi);
	azi = _mm512_fmadd_ps(zj, invR, azi);
      }
      _mm512_store_ps(p+i, pi);
      _mm512_store_ps(ax+i, axi);
      _mm512_store_ps(ay+i, ayi);
      _mm512_store_ps(az+i, azi);
    }
#pragma omp single
    {
      toc = get_time();
      printf("Vectorize target with intrinsics : %e s : %lf GFlops\n",toc-tic, OPS/(toc-tic));

      // Vectorize source with intrinsics
      tic = get_time();
    }
#pragma omp for
    for (i=0; i<N; i++) {
      __m512 pi = _mm512_setzero_ps();
      __m512 axi = _mm512_setzero_ps();
      __m512 ayi = _mm512_setzero_ps();
      __m512 azi = _mm512_setzero_ps();
      __m512 xi = _mm512_set1_ps(x[i]);
      __m512 yi = _mm512_set1_ps(y[i]);
      __m512 zi = _mm512_set1_ps(z[i]);
      for (j=0; j<N; j+=16) {
	__m512 xj = _mm512_load_ps(x+j);
	xj = _mm512_sub_ps(xj, xi);
	__m512 yj = _mm512_load_ps(y+j);
	yj = _mm512_sub_ps(yj, yi);
	__m512 zj = _mm512_load_ps(z+j);
	zj = _mm512_sub_ps(zj, zi);
	__m512 R2 = _mm512_set1_ps(EPS2);
	R2 = _mm512_fmadd_ps(xj, xj, R2);
	R2 = _mm512_fmadd_ps(yj, yj, R2);
	R2 = _mm512_fmadd_ps(zj, zj, R2);
	__m512 mj = _mm512_load_ps(m+j);
	__m512 invR = _mm512_rsqrt23_ps(R2);
	mj = _mm512_mul_ps(mj, invR);
	pi = _mm512_add_ps(pi, mj);
	invR = _mm512_mul_ps(invR, invR);
	invR = _mm512_mul_ps(invR, mj);
	axi = _mm512_fmadd_ps(xj, invR, axi);
	ayi = _mm512_fmadd_ps(yj, invR, ayi);
	azi = _mm512_fmadd_ps(zj, invR, azi);
      }
      p[i] = _mm512_reduce_add_ps(pi);
      ax[i] = _mm512_reduce_add_ps(axi);
      ay[i] = _mm512_reduce_add_ps(ayi);
      az[i] = _mm512_reduce_add_ps(azi);
    }
#pragma omp single
    {
      toc = get_time();
      printf("Vectorize source with intrinsics : %e s : %lf GFlops\n",toc-tic, OPS/(toc-tic));

      // Vectorize target with pragma simd
      tic = get_time();
    }
#pragma simd
#pragma omp for
    for (i=0; i<N; i++) {
      float pi = 0;
      float axi = 0;
      float ayi = 0;
      float azi = 0;
      float xi = x[i];
      float yi = y[i];
      float zi = z[i];
      for (j=0; j<N; j++) {
	float dx = x[j] - xi;
	float dy = y[j] - yi;
	float dz = z[j] - zi;
	float R2 = dx * dx + dy * dy + dz * dz + EPS2;
	float invR = 1.0f / sqrtf(R2);
	float invR3 = m[j] * invR * invR * invR;
	pi += m[j] * invR;
	axi += dx * invR3;
	ayi += dy * invR3;
	azi += dz * invR3;
      }
      p[i] = pi;
      ax[i] = axi;
      ay[i] = ayi;
      az[i] = azi;
    }
#pragma omp single
    {
      toc = get_time();
      printf("Vectorize target with pragma simd: %e s : %lf GFlops\n",toc-tic, OPS/(toc-tic));

      // Vectorize source with pragma simd
      tic = get_time();
    }
#pragma omp for
    for (i=0; i<N; i++) {
      float pi = 0;
      float axi = 0;
      float ayi = 0;
      float azi = 0;
      float xi = x[i];
      float yi = y[i];
      float zi = z[i];
#pragma simd
      for (j=0; j<N; j++) {
	float dx = x[j] - xi;
	float dy = y[j] - yi;
	float dz = z[j] - zi;
	float R2 = dx * dx + dy * dy + dz * dz + EPS2;
	float invR = 1.0f / sqrtf(R2);
	float invR3 = m[j] * invR * invR * invR;
	pi += m[j] * invR;
	axi += dx * invR3;
	ayi += dy * invR3;
	azi += dz * invR3;
      }
      p[i] = pi;
      ax[i] = axi;
      ay[i] = ayi;
      az[i] = azi;
    }
#pragma omp single
    {
      toc = get_time();
      printf("Vectorize source with pragma simd: %e s : %lf GFlops\n",toc-tic, OPS/(toc-tic));
    }
  }

  _mm_free(x);
  _mm_free(y);
  _mm_free(z);
  _mm_free(m);
  _mm_free(p);
  _mm_free(ax);
  _mm_free(ay);
  _mm_free(az);
  return 0;
}
