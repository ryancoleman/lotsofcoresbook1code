// (c) 2013-2014, Colfax International

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "Transpose.h"

void InitMatrix(FTYPE* A, const int n) {

#pragma omp parallel for schedule(static)
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      A[i*n+j] = (FTYPE)(i*n+j);

}

FTYPE VerifyTransposed(FTYPE* A, const int n) {
  // Calculate the norm of the deviation of the transposed matrix elements
  // from the pre-determined pattern
  FTYPE err = 0;
#pragma omp parallel for reduction(+:err)
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      const FTYPE diff = (A[i*n+j] - (FTYPE)(j*n+i));
      err += diff*diff;
    }
  return sqrt(err);
}

void FlushCaches(char* const cacheFlush, const int nCacheFlush) {
  // To avoid the retention of small matrices in cache, read/write a large dummy array
#pragma omp parallel for schedule(guided,1)
  for (int i = 0; i < nCacheFlush; i += 64)
    cacheFlush[i] = (cacheFlush[i] + cacheFlush[i+1])/2;
}

int main(int argc, char** argv){

  if (argc < 2) {
    printf("Format:\n        %s [ size [ trials ] ]\n\n", argv[0]);
    printf("where:\n  size is the matrix size, n, for benchmarking\n");
    printf("  trials is the number of trials (default is 100)\n\n");
    exit(1);
  }

  const int n = atoi(argv[1]);

  int nTrials=100;

  if (argc >= 3)  nTrials = atoi(argv[2]);

  printf("# Benchmarking in-place transposition of [%d x %d] matrix\n", n, n);
  printf("# Platform: %s, threads: %d, trials: %d, method: %d\n", 
	 #ifdef __MIC__ 
	 "MIC",
	 #else
	 "CPU",
	 #endif
	 omp_get_max_threads(), nTrials, METHOD);

  const int skipTrials = 2; // Skip the first two trials (initialization, verification in trial 0)
  const int nSizes = 7;    // How many matrix sizes to test

  // Matrix data container, aligned on a 64-byte boundary
  FTYPE* A = (FTYPE*)_mm_malloc(n*n*sizeof(FTYPE), 64);

  const size_t nCacheFlush = ((size_t)n*(size_t)n*sizeof(FTYPE) < (1L<<27L) ? 1L << 27L : 1L);
  char* cacheFlush = (char*)malloc(nCacheFlush);

  // Container for benchmark results
  double t[nTrials];

  for (int iTrial = 0; iTrial < nTrials; iTrial++) {

    // For the first trial only, initialize the matrix
    if (iTrial == 0) InitMatrix(A, n);

    // Perform and time the transposition operation
    const double t0 = omp_get_wtime();
    Transpose(A, n);
    const double t1 = omp_get_wtime();

    if (iTrial == 0) {
      // For the first trial only, verify the result of the transposition
      if (VerifyTransposed(A,n) > 1e-6) {
	fprintf(stderr,"Result of transposition is incorrect!\n");
	exit(1);
      }
    }

    // To simulate a realistic situation where the matrix at the beginning
    // of a calculation is not in the cache, clear the previously transposed
    // matrix from caches
#ifdef FLUSH_CACHES
    FlushCaches(cacheFlush, nCacheFlush);
#endif
	
    // Record the benchmark result
    t[iTrial] = t1-t0;

  }

  // Calculating transposition rate statistics
  double Ravg = 0;
  double Rmin = 1e100;
  double Rmax = 0;
  double dR = 0;
  for (int i = skipTrials; i < nTrials; i++) {
    // Transposition rate in GB/s:
    const double R = sizeof(FTYPE)*(FTYPE)2*(FTYPE)n*(FTYPE)n/t[i]*1.0e-9;
    Ravg += R;
    dR += R*R;
    Rmin = (R < Rmin ? R : Rmin); // Minimum observed value
    Rmax = (R > Rmax ? R : Rmax); // Maximum observed value
  }
  Ravg /= (nTrials - skipTrials); // Mean
  dR = sqrt(dR/(nTrials-skipTrials) - Ravg*Ravg); // Mean square deviation

  printf("n=%6d  rate= %6.1f +- %4.1f GB/s    range=[ %6.1f ... %6.1f ]\n",
	 n, Ravg, dR, Rmin, Rmax);
  fflush(0);

  free(cacheFlush);

  _mm_free(A);

}
