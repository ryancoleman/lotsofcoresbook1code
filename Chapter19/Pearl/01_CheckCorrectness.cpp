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

const float sig = 0.2f;      // Volatility (0.2 -> 20%)
const float r   = 0.05f;     // Interest rate (0.05 -> 5%)
const float T   = 3.0f;      // Maturity (3 -> 3 years)
const float S0  = 100.0f;    // Initial stock price
const float K   = 100.0f;    // Strike price

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

int main(int argc, char *argv[])
{
  float res = GetOptionPrice();
  printf("%.8f;\n", res);

  return 0;
}
