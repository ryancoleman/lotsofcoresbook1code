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
  int NALIGN = 32;
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
    for (i=0; i<N; i+=8) {
      __m256 pi = _mm256_setzero_ps();
      __m256 axi = _mm256_setzero_ps();
      __m256 ayi = _mm256_setzero_ps();
      __m256 azi = _mm256_setzero_ps();
      __m256 xi = _mm256_load_ps(x+i);
      __m256 yi = _mm256_load_ps(y+i);
      __m256 zi = _mm256_load_ps(z+i);
      __m256 R2 = _mm256_set1_ps(EPS2);
      __m256 x2 = _mm256_set1_ps(x[0]);
      x2 = _mm256_sub_ps(x2, xi);
      __m256 y2 = _mm256_set1_ps(y[0]);
      y2 = _mm256_sub_ps(y2, yi);
      __m256 z2 = _mm256_set1_ps(z[0]);
      z2 = _mm256_sub_ps(z2, zi);
      __m256 mj = _mm256_set1_ps(m[0]);
      __m256 xj = x2;
      x2 = _mm256_mul_ps(x2, x2);
      R2 = _mm256_add_ps(R2, x2);
      __m256 yj = y2;
      y2 = _mm256_mul_ps(y2, y2);
      R2 = _mm256_add_ps(R2, y2);
      __m256 zj = z2;
      z2 = _mm256_mul_ps(z2, z2);
      R2 = _mm256_add_ps(R2, z2);
      __m256 invR;
      x2 = _mm256_set1_ps(x[1]);
      y2 = _mm256_set1_ps(y[1]);
      z2 = _mm256_set1_ps(z[1]);
      for (j=0; j<N-2; j++) {
	invR = _mm256_rsqrt_ps(R2);
	R2 = _mm256_set1_ps(EPS2);
	x2 = _mm256_sub_ps(x2, xi);
	y2 = _mm256_sub_ps(y2, yi);
	z2 = _mm256_sub_ps(z2, zi);
	mj = _mm256_mul_ps(mj, invR);
	pi = _mm256_add_ps(pi, mj);
	invR = _mm256_mul_ps(invR, invR);
	invR = _mm256_mul_ps(invR, mj);
	mj = _mm256_set1_ps(m[j+1]);
	xj = _mm256_mul_ps(xj, invR);
	axi = _mm256_add_ps(axi, xj);
	xj = x2;
	x2 = _mm256_mul_ps(x2, x2);
	R2 = _mm256_add_ps(R2, x2);
	x2 = _mm256_set1_ps(x[j+2]);
	yj = _mm256_mul_ps(yj, invR);
	ayi = _mm256_add_ps(ayi, yj);
	yj = y2;
	y2 = _mm256_mul_ps(y2, y2);
	R2 = _mm256_add_ps(R2, y2);
	y2 = _mm256_set1_ps(y[j+2]);
	zj = _mm256_mul_ps(zj, invR);
	azi = _mm256_add_ps(azi, zj);
	zj = z2;
	z2 = _mm256_mul_ps(z2, z2);
	R2 = _mm256_add_ps(R2, z2);
	z2 = _mm256_set1_ps(z[j+2]);
      }
      invR = _mm256_rsqrt_ps(R2);
      R2 = _mm256_set1_ps(EPS2);
      x2 = _mm256_sub_ps(x2, xi);
      y2 = _mm256_sub_ps(y2, yi);
      z2 = _mm256_sub_ps(z2, zi);
      mj = _mm256_mul_ps(mj, invR);
      pi = _mm256_add_ps(pi, mj);
      invR = _mm256_mul_ps(invR, invR);
      invR = _mm256_mul_ps(invR, mj);
      mj = _mm256_set1_ps(m[N-1]);
      xj = _mm256_mul_ps(xj, invR);
      axi = _mm256_add_ps(axi, xj);
      xj = x2;
      x2 = _mm256_mul_ps(x2, x2);
      R2 = _mm256_add_ps(R2, x2);
      yj = _mm256_mul_ps(yj, invR);
      ayi = _mm256_add_ps(ayi, yj);
      yj = y2;
      y2 = _mm256_mul_ps(y2, y2);
      R2 = _mm256_add_ps(R2, y2);
      zj = _mm256_mul_ps(zj, invR);
      azi = _mm256_add_ps(azi, zj);
      zj = z2;
      z2 = _mm256_mul_ps(z2, z2);
      R2 = _mm256_add_ps(R2, z2);
      invR = _mm256_rsqrt_ps(R2);
      mj = _mm256_mul_ps(mj, invR);
      pi = _mm256_add_ps(pi, mj);
      invR = _mm256_mul_ps(invR, invR);
      invR = _mm256_mul_ps(invR, mj);
      xj = _mm256_mul_ps(xj, invR);
      axi = _mm256_add_ps(axi, xj);
      yj = _mm256_mul_ps(yj, invR);
      ayi = _mm256_add_ps(ayi, yj);
      zj = _mm256_mul_ps(zj, invR);
      azi = _mm256_add_ps(azi, zj);
      _mm256_store_ps(p+i, pi);
      _mm256_store_ps(ax+i, axi);
      _mm256_store_ps(ay+i, ayi);
      _mm256_store_ps(az+i, azi);
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
      __m256 pi = _mm256_setzero_ps();
      __m256 axi = _mm256_setzero_ps();
      __m256 ayi = _mm256_setzero_ps();
      __m256 azi = _mm256_setzero_ps();
      __m256 xi = _mm256_set1_ps(x[i]);
      __m256 yi = _mm256_set1_ps(y[i]);
      __m256 zi = _mm256_set1_ps(z[i]);
      for (j=0; j<N; j+=8) {
	__m256 R2 = _mm256_set1_ps(EPS2);
	__m256 xj = _mm256_load_ps(x+j);
	xj = _mm256_sub_ps(xj, xi);
	__m256 yj = _mm256_load_ps(y+j);
	yj = _mm256_sub_ps(yj, yi);
	__m256 zj = _mm256_load_ps(z+j);
	zj = _mm256_sub_ps(zj, zi);
	__m256 mj = _mm256_load_ps(m+j);
	__m256 x2 = _mm256_mul_ps(xj, xj);
	R2 = _mm256_add_ps(R2, x2);
	__m256 y2 = _mm256_mul_ps(yj, yj);
	R2 = _mm256_add_ps(R2, y2);
	__m256 z2 = _mm256_mul_ps(zj, zj);
	R2 = _mm256_add_ps(R2, z2);
	__m256 invR = _mm256_rsqrt_ps(R2);
	mj = _mm256_mul_ps(mj, invR);
	pi = _mm256_add_ps(pi, mj);
	invR = _mm256_mul_ps(invR, invR);
	invR = _mm256_mul_ps(invR, mj);
	xj = _mm256_mul_ps(xj, invR);
	axi = _mm256_add_ps(axi, xj);
	yj = _mm256_mul_ps(yj, invR);
	ayi = _mm256_add_ps(ayi, yj);
	zj = _mm256_mul_ps(zj, invR);
	azi = _mm256_add_ps(azi, zj);
      }
      pi = _mm256_add_ps(_mm256_permute2f128_ps(pi,pi,1),pi);
      pi = _mm256_hadd_ps(pi,pi);
      pi = _mm256_hadd_ps(pi,pi);
      p[i] = ((float*)&pi)[0];
      axi = _mm256_add_ps(_mm256_permute2f128_ps(axi,axi,1),axi);
      axi = _mm256_hadd_ps(axi,axi);
      axi = _mm256_hadd_ps(axi,axi);
      ax[i] = ((float*)&axi)[0];
      ayi = _mm256_add_ps(_mm256_permute2f128_ps(ayi,ayi,1),ayi);
      ayi = _mm256_hadd_ps(ayi,ayi);
      ayi = _mm256_hadd_ps(ayi,ayi);
      ay[i] = ((float*)&ayi)[0];
      azi = _mm256_add_ps(_mm256_permute2f128_ps(azi,azi,1),azi);
      azi = _mm256_hadd_ps(azi,azi);
      azi = _mm256_hadd_ps(azi,azi);
      az[i] = ((float*)&azi)[0];
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
