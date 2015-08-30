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
#include <stdio.h>
#include "omp.h"

const float sig = 0.2f;      // Volatility (0.2 -> 20%)
const float r   = 0.05f;     // Interest rate (0.05 -> 5%)
const float T   = 3.0f;      // Maturity (3 -> 3 years)
const float S0  = 100.0f;    // Initial stock price
const float K   = 100.0f;    // Strike price

float *pT, *pS0, *pK, *pC;
int N;

__declspec(noinline) float GetOptionPrice()
{
  float C;
  float d1, d2, p1, p2;

  d1 = (logf(S0 / K) + (r + sig * sig * 0.5) * T) / (sig * sqrtf(T));
  d2 = d1 - sig * sqrtf(T);
  p1 = cdfnormf(d1);
  p2 = cdfnormf(d2);
  C  = S0 * p1 - K * expf((-1.0) * r * T) * p2;

  return C;
}

__declspec(noinline) void GetOptionPrices(float *pT, float *pK, float *pS0, float *pC)
{
  int i;
  float d1, d2, p1, p2;

  for (i = 0; i < N; i++)
  {
    d1 = (log(pS0[i] / pK[i]) + (r + sig * sig * 0.5) * pT[i]) / (sig * sqrt(pT[i]));
    d2 = (log(pS0[i] / pK[i]) + (r - sig * sig * 0.5) * pT[i]) / (sig * sqrt(pT[i]));
    p1 = cdfnormf(d1);
    p2 = cdfnormf(d2);
    pC[i]  = pS0[i] * p1 - pK[i] * exp((-1.0) * r * pT[i]) * p2;
  }
}

double start, finish;
double t;

int main(int argc, char *argv[])
{
  int i;

  if (argc < 2)
  {
    printf("Usage: <executable> size\n");
    return 1;
  }
  N = atoi(argv[1]);

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

  start = omp_get_wtime();
  GetOptionPrices(pT, pK, pS0, pC);
  finish = omp_get_wtime();
  t = finish - start;
  printf("v_02: %.8f; time = %lf\n", pC[0], t);

  delete [] pT;
  return 0;
}
