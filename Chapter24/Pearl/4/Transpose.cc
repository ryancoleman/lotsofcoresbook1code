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

  // nEven is a multiple of TILE
  const int nEven = n - n%TILE;
  const int wTiles = nEven / TILE;                  // Complete tiles in each dimens.
  const int nTilesParallel = wTiles*(wTiles - 1)/2; // # of tiles in the matrix body


#pragma omp parallel
  { 

#pragma omp for schedule(static)
    for (int ii = TILE; ii < nEven; ii += TILE)
      for (int jj = 0; jj < ii; jj += TILE)

	// The main tile transposition microkernel
	for (int j = jj; j < jj+TILE; j++) 
	  for (int i = ii; i < ii+TILE; i++) 
	    {
	      const FTYPE c = A[i*n + j];
	      A[i*n + j] = A[j*n + i];
	      A[j*n + i] = c;
	    }



    // Transposing the tiles along the main diagonal
#pragma omp for schedule(static)
    for (int ii = 0; ii < nEven; ii += TILE) {
      const int jj = ii;

      for (int j = jj; j < jj+TILE; j++) 
	for (int i = ii; i < j; i++) {
	  const FTYPE c = A[i*n + j];
	  A[i*n + j] = A[j*n + i];
	  A[j*n + i] = c;
	}
    }

    // Transposing ``peel'' around the right and bottom perimeter
#pragma omp for schedule(static)
    for (int j = 0; j < nEven; j++) {
      for (int i = nEven; i < n; i++) {
	const FTYPE c = A[i*n + j];
	A[i*n + j] = A[j*n + i];
	A[j*n + i] = c;
      }
    }

  }

  // Transposing the bottom-right corner
  for (int j = nEven; j < n; j++) {
    for (int i = nEven; i < j; i++) {
      const FTYPE c = A[i*n + j];
      A[i*n + j] = A[j*n + i];
      A[j*n + i] = c;
    }
  }

}

