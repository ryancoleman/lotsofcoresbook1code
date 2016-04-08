// This code supplements the white paper
//    "Multithreaded Transposition of Square Matrices
//     with Common Code for 
//     Intel Xeon Processors and Intel Xeon Phi Coprocessors"
// available at the following URL:
//     http://research.colfaxinternational.com/post/2013/08/12/Trans-7110.aspx
// You are free to use, modify and distribute this code as long as you acknowledge
// the above mentioned publication.
// (c) Colfax International, 2013

#include "Transpose.h"
#include <cstdlib>

const int TILE = 16; // Empirically chosen tile size

void Transpose(FTYPE* const A, const int n) {

#pragma omp parallel for schedule(static)
  for (int ii = 0; ii < n; ii += TILE) {
    for (int jj = 0; jj <= ii; jj += TILE) {

      const int jMax = (jj+TILE < n ? jj+TILE : n);

      for (int j = jj; j < jMax; j++)  {

	const int iMin = (ii > j ? ii : j+1);
	const int iMax = (ii+TILE < n ? ii+TILE : n);
			    
	for (int i = iMin; i < iMax; i++) 
	  {
	    const FTYPE c = A[i*n + j];
	    A[i*n + j] = A[j*n + i];
	    A[j*n + i] = c;
	  }
      }
    }
  }
}
