//////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2014, FLORIAN WENDE, ZIB
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the ZIB nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL FLORIAN WENDE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
//////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <vector>

using namespace std;

int main( const int argc, const char **argv ) {

  vector<double>
    values_1;
  double
    temp_1,
    mean_1,
    err_1;
  uint32_t
    a,b;
  FILE
    *fp = NULL;
  char
    buffer[1024];

  values_1.clear();

  fp = fopen( argv[1], "r" );

  while( fgets( buffer, 1024, fp ) != NULL ) {

    sscanf( buffer, "%u%u%lf", &a, &b, &temp_1 );
    values_1.push_back( temp_1 );

  }

  mean_1 = 0.0;
  for( uint32_t i=0; i<values_1.size(); i++ )
    mean_1 += values_1[i];
  mean_1 /= values_1.size();

  err_1 = 0.0;
  for( uint32_t i=0; i<values_1.size(); i++ )
    err_1 += ( mean_1-values_1[i] )*( mean_1-values_1[i] );
  err_1 = sqrt( err_1/( values_1.size()*( values_1.size()-1 ) ) );

  printf( "%u\t%u\t%.6lf\t%.6lf\n", a, b, mean_1, err_1 );

  return 0;

}
