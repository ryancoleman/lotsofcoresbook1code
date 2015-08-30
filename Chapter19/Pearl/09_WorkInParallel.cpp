/*
* Copyright (C) 2014 University of Nizhni Novgorod
* This file is part of the Black-Scholes Pricing Code, v 1.0
* 
* This file may be used under the terms of the GNU Lesser General Public 
* License version 2.1 as published by the Free Software Foundation 
* and appearing in the file LICENSE.LGPL included in the packaging of this 
* file. Please review the following information to ensure the GNU Lesser 
* General Public License version 2.1 requirements will be met: 
* http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.
*/

#include <math.h>
#include <mathimf.h>
#include <malloc.h>
#include <stdio.h>
#include "omp.h"

const float sig = 0.2f;      // Volatility (0.2 -> 20%)
const float r   = 0.05f;     // Interest rate (0.05 -> 5%)
const float T   = 3.0f;      // Maturity (3 -> 3 years)
const float S0  = 100.0f;    // Initial stock price
const float K   = 100.0f;    // Strike price

const float invsqrt2 = 0.707106781f;

float *pT, *pS0, *pK, *pC;
int N;
int numThreads = 1;

__declspec(noinline) float GetOptionPrice()
{
  float C;
  float d1, d2, p1, p2;

  d1 = (logf(S0 / K) + (r + sig * sig * 0.5) * T) / (sig * sqrtf(T));
  d2 = d1 - sig * sqrtf(T);
  p1 = cdfnorm(d1);
  p2 = cdfnorm(d2);
  C  = S0 * p1 - K * expf((-1.0) * r * T) * p2;

  return C;
}

__declspec(noinline) void GetOptionPrices(float *pT, float *pK, float *pS0, float *pC)
{
  int i;
  float d1, d2, erf1, erf2, invf;
  float sig2 = sig * sig;

#pragma simd
#pragma omp parallel for private(d1, d2, erf1, erf2, invf)
  for (i = 0; i < N; i++)
  {
    invf = invsqrtf(sig2 * pT[i]);
    d1 = (logf(pS0[i] / pK[i]) + (r + sig2 * 0.5f) * pT[i]) / invf;
    d2 = (logf(pS0[i] / pK[i]) + (r - sig2 * 0.5f) * pT[i]) / invf;
    erf1 = 0.5f + 0.5f * erff(d1 * invsqrt2);
    erf2 = 0.5f + 0.5f * erff(d2 * invsqrt2);
    pC[i]  = pS0[i] * erf1 - pK[i] * expf((-1.0f) * r * pT[i]) * erf2;
  }
}

double start, finish;
double t;

int main(int argc, char *argv[])
{
  int i;

  if (argc < 2)
  {
    printf("Usage: <executable> size [#of_threads]\n");
    return 1;
  }
  N = atoi(argv[1]);
  if (argc > 2)
    numThreads = atoi(argv[2]);

//  pT = (float *)memalign(32, 4 * N * sizeof(float));
  pT  = new float[4 * N];
  pK  = pT + N;
  pS0 = pT + 2 * N;
  pC  = pT + 3 * N;
  for (i = 0; i < N; i++)
  {
    pT[i] = T;
    pS0[i] = S0;
    pK[i] = K;
  }

  float res = GetOptionPrice();
  printf("%.8f;\n", res);

  omp_set_num_threads(numThreads);

  start = omp_get_wtime();
  GetOptionPrices(pT, pK, pS0, pC);
  finish = omp_get_wtime();
  t = finish - start;
  printf("v_09: %.8f; time = %lf\n", pC[0], t);

//  free(pT);
  delete [] pT;
  return 0;
}
