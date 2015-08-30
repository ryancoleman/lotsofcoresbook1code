/* ----------------------------------------------------------------------
   miniMD is a simple, parallel molecular dynamics (MD) code.   miniMD is
   an MD microapplication in the Mantevo project at Sandia National
   Laboratories ( http://www.mantevo.org ). The primary
   authors of miniMD are Steve Plimpton (sjplimp@sandia.gov) , Paul Crozier
   (pscrozi@sandia.gov) and Christian Trott (crtrott@sandia.gov).

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This library is free software; you
   can redistribute it and/or modify it under the terms of the GNU Lesser
   General Public License as published by the Free Software Foundation;
   either version 3 of the License, or (at your option) any later
   version.

   This library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this software; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA.  See also: http://www.gnu.org/licenses/lgpl.txt .

   For questions, contact Paul S. Crozier (pscrozi@sandia.gov) or
   Christian Trott (crtrott@sandia.gov).

   Please read the accompanying README and LICENSE files.
---------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include <algorithm>
#include <immintrin.h>

#include "neighbor.h"
#include "openmp.h"

#define FACTOR 0.999
#define SMALL 1.0e-6
#define AVX

#ifdef AVX
#define _mm256_cmpneq_ps(x,y) _mm256_cmp_ps(x,y,_CMP_NEQ_UQ)            
#define _mm256_cmpge_ps(x,y)  _mm256_cmp_ps(x,y,_CMP_GE_OS)
#define _mm256_cmple_ps(x,y)  _mm256_cmp_ps(x,y,_CMP_LE_OS)      
#endif

Neighbor::Neighbor()
{
  ncalls = 0;
  max_totalneigh = 0;
  numneigh = NULL;
  numrem = NULL;
  neighbors = NULL;
  maxneighs = 100;
  nmax = 0;
  bincount = NULL;
  bins = NULL;
  atoms_per_bin = 8;
  stencil = NULL;
  threads = NULL;
  halfneigh = 0;
  ghost_newton = 1;
}

Neighbor::~Neighbor()
{
#ifdef ALIGNMALLOC
  if(numneigh) _mm_free(numneigh);
  if(numrem) _mm_free(numrem);
  if(neighbors) _mm_free(neighbors);
#else 
  if(numneigh) free(numneigh);
  if(numrem) free(numrem);
  if(neighbors) free(neighbors);
#endif
  
  if(bincount) free(bincount);

  if(bins) free(bins);
}

// HPPG: Helper struct for std::sort that separates local/remote neighbors.
struct neighbor_sort
{
  int start;
  int end;
  bool operator() (int i, int j)
  {
    bool i_local = (i >= start && i < end);
    bool j_local = (j >= start && j < end);
    if (i_local != j_local) return i_local;
    else return (i < j);
  }  
};

// HPPG: Helper struct for std::stable_partition that separates local/remote neighbors.
struct neighbor_partition
{
  int start;
  int end;
  inline bool operator() (int i)
  {
    return (i >= start && i < end);
  }
};

/* function to pack neighbors into a list using SSE */
__m128i SHUFFLE_TABLE[16];
unsigned inline pack_neighbors(int* p, __m128i ids, __m128 mask) {
    unsigned int _mask = _mm_movemask_ps(mask);
    __m128i result = _mm_shuffle_epi8(ids, SHUFFLE_TABLE[_mask]);        
    _mm_storeu_si128((__m128i*) p, result);
    return _mm_popcnt_u32(_mask);
}

/* binned neighbor list construction with full Newton's 3rd law
   every pair stored exactly once by some processor
   each owned atom i checks its own bin and other bins in Newton stencil */
void Neighbor::build(Atom &atom)
{
  ncalls++;
  const int nlocal = atom.nlocal;
  const int nall = atom.nlocal + atom.nghost;

  /* extend atom arrays if necessary */
  #pragma omp master
  if(nall > nmax) {
    nmax = nall;
#ifdef ALIGNMALLOC
    if(numneigh) _mm_free(numneigh);
    if(numrem) _mm_free(numrem);
    numneigh = (int*) _mm_malloc(nmax * sizeof(int) + ALIGNMALLOC, ALIGNMALLOC);
    numrem = (int*) _mm_malloc(nmax * sizeof(int) + ALIGNMALLOC, ALIGNMALLOC);
    if(neighbors) _mm_free(neighbors);	
    neighbors = (int*) _mm_malloc(nmax * maxneighs * sizeof(int*) + ALIGNMALLOC, ALIGNMALLOC);	
#else

    if(numneigh) free(numneigh);
    if(numrem) free(numrem);

    if(neighbors) free(neighbors);

    numneigh = (int*) malloc(nmax * sizeof(int));
    numrem = (int*) malloc(nmax * sizeof(int));
    neighbors = (int*) malloc(nmax * maxneighs * sizeof(int*));
#endif
  }

  int omp_me = omp_get_thread_num();
  int num_omp_threads = threads->omp_num_threads;
  int master = -1;

  #pragma omp master
  {
    master = omp_me;

    #ifdef AVX
    SHUFFLE_TABLE[ 0] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);     // 0000
    SHUFFLE_TABLE[ 1] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 0);     // 0001
    SHUFFLE_TABLE[ 2] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 6, 5, 4);     // 0010
    SHUFFLE_TABLE[ 3] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 7, 6, 5, 4, 3, 2, 1, 0);     // 0011
    SHUFFLE_TABLE[ 4] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 10, 9, 8);   // 0100
    SHUFFLE_TABLE[ 5] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 11, 10, 9, 8, 3, 2, 1, 0);   // 0101
    SHUFFLE_TABLE[ 6] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 11, 10, 9, 8, 7, 6, 5, 4);   // 0110
    SHUFFLE_TABLE[ 7] = _mm_set_epi8(0, 0, 0, 0, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);   // 0111
    SHUFFLE_TABLE[ 8] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 15, 14, 13, 12); // 1000
    SHUFFLE_TABLE[ 9] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 14, 13, 12, 3, 2, 1, 0); // 1001
    SHUFFLE_TABLE[10] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 14, 13, 12, 7, 6, 5, 4); // 1010
    SHUFFLE_TABLE[11] = _mm_set_epi8(0, 0, 0, 0, 15, 14, 13, 12, 7, 6, 5, 4, 3, 2, 1, 0); // 1011
    SHUFFLE_TABLE[12] = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 15, 14, 13, 12, 11, 10, 9, 8); // 1100
    SHUFFLE_TABLE[13] = _mm_set_epi8(0, 0, 0, 0, 15, 14, 13, 12, 11, 10, 9, 8, 3, 2, 1, 0); // 1101
    SHUFFLE_TABLE[14] = _mm_set_epi8(0, 0, 0, 0, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4); // 1110
    SHUFFLE_TABLE[15] = _mm_set_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0); // 1111
    #endif
  }

  // HPPG: Set up decomposition across threads.
  int tid = omp_get_thread_num();
  int perthread = nlocal / threads->omp_num_threads;
  int leftover = nlocal % threads->omp_num_threads;
  int start_atom = tid * perthread;
  int end_atom = start_atom + perthread;
  if (tid == threads->omp_num_threads-1) end_atom = nlocal;
  threads->data[tid].start_atom = start_atom;
  threads->data[tid].end_atom = end_atom;
  
  struct neighbor_sort sorter;
  sorter.start = start_atom;
  sorter.end = end_atom;
  
  struct neighbor_partition partitioner;
  partitioner.start = start_atom;
  partitioner.end = end_atom;  

  #pragma omp barrier
  /* bin local & ghost atoms */

  binatoms(atom);
  count = 0;
  /* loop over each atom, storing neighbors */

  const MMD_float* x = &atom.x[0][0];

#ifdef AVX
  __m256 mcutneighsq = _mm256_set1_ps(cutneighsq);
#endif

  resize = 1;
  #pragma omp barrier

  while(resize) {
    #pragma omp barrier
    int new_maxneighs = maxneighs;
    resize = 0;
    #pragma omp barrier

    int local_tot = 0;
    int remote_tot = 0;
    int cached_bin = -1; // we cache the stencil for this bin
    int ncache = 0;
    int CACHE_SIZE = 1024; // hard-code the size of our cache
    if (!threads->data[tid].stencil_cache) {
      threads->data[tid].stencil_cache =  (MMD_float*) _mm_malloc(4 * CACHE_SIZE * sizeof(MMD_float), 64);
    }
    MMD_float* stencil_cache = threads->data[tid].stencil_cache;
    for (int i = start_atom; i < end_atom; i++) {

      int* neighptr = &neighbors[i * maxneighs];

      int n = 0;
      int local = 0;
      int remote = 0;

      const MMD_float xtmp = x[i * PAD + 0];
      const MMD_float ytmp = x[i * PAD + 1];
      const MMD_float ztmp = x[i * PAD + 2];

#ifdef AVX
      __m128i mi = _mm_set1_epi32(i);
      __m128i mstart = _mm_set1_epi32(start_atom);
      __m128i mend = _mm_set1_epi32(end_atom-1); // need to get around lack of cmpge/cmple
      __m256 mxtmp = _mm256_set1_ps(xtmp);
      __m256 mytmp = _mm256_set1_ps(ytmp);
      __m256 mztmp = _mm256_set1_ps(ztmp);
#endif

      // If we encounter a new bin, cache its contents.
      const int ibin = coord2bin(xtmp, ytmp, ztmp);
      if (ibin != cached_bin)
      {
        ncache = 0;
        for (int k = 0; k < nstencil; k++)
        {
          const int jbin = ibin + stencil[k];
          int* loc_bin = &bins[jbin * atoms_per_bin];
          for (int m = 0; m < bincount[jbin]; m++)
          {
            const int j = loc_bin[m];
            *((int*)&stencil_cache[0*CACHE_SIZE + ncache]) = j;
            stencil_cache[1*CACHE_SIZE + ncache] = x[j * PAD + 0];
            stencil_cache[2*CACHE_SIZE + ncache] = x[j * PAD + 1];
            stencil_cache[3*CACHE_SIZE + ncache] = x[j * PAD + 2];
            ncache++;
          }
        }
        if (ncache >= CACHE_SIZE) printf("ERROR: Too many atoms in the stencil - %d > %d\n", ncache, CACHE_SIZE);
        cached_bin = ibin;
      }

      // Otherwise, we just look at the neighbors in the cache.
      int c = 0;
#ifdef AVX
      for (; c < (ncache/8)*8; c += 8)
      {
        const __m256 delx = _mm256_sub_ps(mxtmp, _mm256_load_ps(&stencil_cache[1*CACHE_SIZE + c]));
        const __m256 dely = _mm256_sub_ps(mytmp, _mm256_load_ps(&stencil_cache[2*CACHE_SIZE + c]));
        const __m256 delz = _mm256_sub_ps(mztmp, _mm256_load_ps(&stencil_cache[3*CACHE_SIZE + c]));
        const __m256 rsq = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(delx, delx), _mm256_mul_ps(dely, dely)), _mm256_mul_ps(delz, delz));
        __m256 mask = _mm256_cmple_ps(rsq, mcutneighsq);

        __m128i j1 = _mm_load_si128((__m128i*) &stencil_cache[c + 0]);
        __m128 cmask1 = _mm256_castps256_ps128(mask);
        __m128 nmask1 = _mm_castsi128_ps(_mm_cmpgt_epi32(j1, mi));
        __m128 rmask1 = _mm_castsi128_ps(_mm_or_si128(_mm_cmplt_epi32(j1, mstart), _mm_cmpgt_epi32(j1, mend)));
        __m128 lmask1 = _mm_andnot_ps(rmask1, nmask1);
        rmask1 = _mm_and_ps(cmask1, rmask1);
        lmask1 = _mm_and_ps(cmask1, lmask1);
        n += pack_neighbors(&neighptr[n], j1, _mm_or_ps(lmask1, rmask1));
        local += _mm_popcnt_u32(_mm_movemask_ps(lmask1));
        remote += _mm_popcnt_u32(_mm_movemask_ps(rmask1));
        
        __m128i j2 = _mm_load_si128((__m128i*) &stencil_cache[c + 4]);
        __m128 cmask2 = _mm256_extractf128_ps(mask, 1);
        __m128 nmask2 = _mm_castsi128_ps(_mm_cmpgt_epi32(j2, mi));
        __m128 rmask2 = _mm_castsi128_ps(_mm_or_si128(_mm_cmplt_epi32(j2, mstart), _mm_cmpgt_epi32(j2, mend)));
        __m128 lmask2 = _mm_andnot_ps(rmask2, nmask2);
        rmask2 = _mm_and_ps(cmask2, rmask2);
        lmask2 = _mm_and_ps(cmask2, lmask2);
        n += pack_neighbors(&neighptr[n], j2, _mm_or_ps(lmask2, rmask2));
        local += _mm_popcnt_u32(_mm_movemask_ps(lmask2));
        remote += _mm_popcnt_u32(_mm_movemask_ps(rmask2));
      }
#endif

      for (; c < ncache; c++)
      {

        const int j = *((int*)&stencil_cache[0*CACHE_SIZE +c]);
        if (j <= i && j >= start_atom && j < end_atom) continue;

        const MMD_float delx = xtmp - stencil_cache[1*CACHE_SIZE + c];
        const MMD_float dely = ytmp - stencil_cache[2*CACHE_SIZE + c];
        const MMD_float delz = ztmp - stencil_cache[3*CACHE_SIZE + c];
        const MMD_float rsq = delx * delx + dely * dely + delz * delz;

        if (rsq <= cutneighsq)
        {
           neighptr[n++] = j;
           if (j >= start_atom && j < end_atom) local++;
           else remote++;
        }

      }

      numneigh[i] = local;
      numrem[i] = remote;
      std::stable_partition(neighptr, neighptr + n, partitioner);
      local_tot += local;
      remote_tot += remote;

      if(n >= maxneighs) {
        resize = 1;
        if(n >= new_maxneighs) new_maxneighs = n;
      }
    }
    double tot = local_tot + remote_tot;
    //printf("%d: local = %d (%g%%), remote = %d (%g%%)\n", tid, local_tot, local_tot / tot * 100, remote_tot, remote_tot / tot * 100);

    // #pragma omp barrier

    if(resize) {
      #pragma omp master
      {
        maxneighs = new_maxneighs * 1.2;
#ifdef ALIGNMALLOC
  		_mm_free(neighbors);
  		neighbors = (int*) _mm_malloc(nmax* maxneighs * sizeof(int) + ALIGNMALLOC, ALIGNMALLOC);
#else
  		free(neighbors);
        neighbors = (int*) malloc(nmax* maxneighs * sizeof(int));
#endif
      }
      #pragma omp barrier
    }
  }

  #pragma omp barrier

}

void Neighbor::binatoms(Atom &atom, int count)
{
  const int omp_me = omp_get_thread_num();
  const int num_omp_threads = threads->omp_num_threads;

  const int nlocal = atom.nlocal;
  const int nall = count<0?atom.nlocal + atom.nghost:count;
  const MMD_float* x = &atom.x[0][0];

  xprd = atom.box.xprd;
  yprd = atom.box.yprd;
  zprd = atom.box.zprd;

  resize = 1;

  #pragma omp barrier

  while(resize > 0) {
    #pragma omp barrier
    resize = 0;
    #pragma omp barrier
    #pragma omp for schedule(static)
    for(int i = 0; i < mbins; i++) bincount[i] = 0;


    OMPFORSCHEDULE
    for(int i = 0; i < nall; i++) {
      const int ibin = coord2bin(x[i * PAD + 0], x[i * PAD + 1], x[i * PAD + 2]);

      if(bincount[ibin] < atoms_per_bin) {
        int ac;
#ifdef OpenMP31
        #pragma omp atomic capture
        ac = bincount[ibin]++;
#else
        ac = __sync_fetch_and_add(bincount + ibin, 1);
#endif
        bins[ibin * atoms_per_bin + ac] = i;
      } else resize = 1;
    }

    // #pragma omp barrier

    #pragma omp master

    if(resize) {
      free(bins);
      atoms_per_bin *= 2;
      bins = (int*) malloc(mbins * atoms_per_bin * sizeof(int));
    }

    // #pragma omp barrier
  }

  #pragma omp barrier

}

/* convert xyz atom coords into local bin #
   take special care to insure ghost atoms with
   coord >= prd or coord < 0.0 are put in correct bins */

inline int Neighbor::coord2bin(MMD_float x, MMD_float y, MMD_float z)
{
  int ix, iy, iz;

  if(x >= xprd)
    ix = (int)((x - xprd) * bininvx) + nbinx - mbinxlo;
  else if(x >= 0.0)
    ix = (int)(x * bininvx) - mbinxlo;
  else
    ix = (int)(x * bininvx) - mbinxlo - 1;

  if(y >= yprd)
    iy = (int)((y - yprd) * bininvy) + nbiny - mbinylo;
  else if(y >= 0.0)
    iy = (int)(y * bininvy) - mbinylo;
  else
    iy = (int)(y * bininvy) - mbinylo - 1;

  if(z >= zprd)
    iz = (int)((z - zprd) * bininvz) + nbinz - mbinzlo;
  else if(z >= 0.0)
    iz = (int)(z * bininvz) - mbinzlo;
  else
    iz = (int)(z * bininvz) - mbinzlo - 1;

  return (iz * mbiny * mbinx + iy * mbinx + ix + 1);
}


/*
setup neighbor binning parameters
bin numbering is global: 0 = 0.0 to binsize
                         1 = binsize to 2*binsize
                         nbin-1 = prd-binsize to binsize
                         nbin = prd to prd+binsize
                         -1 = -binsize to 0.0
coord = lowest and highest values of ghost atom coords I will have
        add in "small" for round-off safety
mbinlo = lowest global bin any of my ghost atoms could fall into
mbinhi = highest global bin any of my ghost atoms could fall into
mbin = number of bins I need in a dimension
stencil() = bin offsets in 1-d sense for stencil of surrounding bins
*/

int Neighbor::setup(Atom &atom)
{
  int i, j, k, nmax;
  MMD_float coord;
  int mbinxhi, mbinyhi, mbinzhi;
  int nextx, nexty, nextz;
  int num_omp_threads = threads->omp_num_threads;

  cutneighsq = cutneigh * cutneigh;

  xprd = atom.box.xprd;
  yprd = atom.box.yprd;
  zprd = atom.box.zprd;

  /*
  c bins must evenly divide into box size,
  c   becoming larger than cutneigh if necessary
  c binsize = 1/2 of cutoff is near optimal

  if (flag == 0) {
    nbinx = 2.0 * xprd / cutneigh;
    nbiny = 2.0 * yprd / cutneigh;
    nbinz = 2.0 * zprd / cutneigh;
    if (nbinx == 0) nbinx = 1;
    if (nbiny == 0) nbiny = 1;
    if (nbinz == 0) nbinz = 1;
  }
  */

  binsizex = xprd / nbinx;
  binsizey = yprd / nbiny;
  binsizez = zprd / nbinz;
  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;

  coord = atom.box.xlo - cutneigh - SMALL * xprd;
  mbinxlo = static_cast<int>(coord * bininvx);

  if(coord < 0.0) mbinxlo = mbinxlo - 1;

  coord = atom.box.xhi + cutneigh + SMALL * xprd;
  mbinxhi = static_cast<int>(coord * bininvx);

  coord = atom.box.ylo - cutneigh - SMALL * yprd;
  mbinylo = static_cast<int>(coord * bininvy);

  if(coord < 0.0) mbinylo = mbinylo - 1;

  coord = atom.box.yhi + cutneigh + SMALL * yprd;
  mbinyhi = static_cast<int>(coord * bininvy);

  coord = atom.box.zlo - cutneigh - SMALL * zprd;
  mbinzlo = static_cast<int>(coord * bininvz);

  if(coord < 0.0) mbinzlo = mbinzlo - 1;

  coord = atom.box.zhi + cutneigh + SMALL * zprd;
  mbinzhi = static_cast<int>(coord * bininvz);

  /* extend bins by 1 in each direction to insure stencil coverage */

  mbinxlo = mbinxlo - 1;
  mbinxhi = mbinxhi + 1;
  mbinx = mbinxhi - mbinxlo + 1;

  mbinylo = mbinylo - 1;
  mbinyhi = mbinyhi + 1;
  mbiny = mbinyhi - mbinylo + 1;

  mbinzlo = mbinzlo - 1;
  mbinzhi = mbinzhi + 1;
  mbinz = mbinzhi - mbinzlo + 1;

  /*
  compute bin stencil of all bins whose closest corner to central bin
  is within neighbor cutoff
  for partial Newton (newton = 0),
  stencil is all surrounding bins including self
  for full Newton (newton = 1),
  stencil is bins to the "upper right" of central bin, does NOT include self
  next(xyz) = how far the stencil could possibly extend
  factor < 1.0 for special case of LJ benchmark so code will create
  correct-size stencil when there are 3 bins for every 5 lattice spacings
  */

  nextx = static_cast<int>(cutneigh * bininvx);

  if(nextx * binsizex < FACTOR * cutneigh) nextx++;

  nexty = static_cast<int>(cutneigh * bininvy);

  if(nexty * binsizey < FACTOR * cutneigh) nexty++;

  nextz = static_cast<int>(cutneigh * bininvz);

  if(nextz * binsizez < FACTOR * cutneigh) nextz++;

  nmax = (2 * nextz + 1) * (2 * nexty + 1) * (2 * nextx + 1);

  if(stencil) free(stencil);

  stencil = (int*) malloc(nmax * sizeof(int));

  nstencil = 0;
  int kstart = -nextz;

  // HPPG: Always look in all directions for neighbors.
  /*if(halfneigh && ghost_newton) {
    kstart = 0;
    stencil[nstencil++] = 0;
  }*/

  // HPPG: Always look in all directions for neighbors.
  for(k = kstart; k <= nextz; k++) {
    for(j = -nexty; j <= nexty; j++) {
      for(i = -nextx; i <= nextx; i++) {
        //if(!ghost_newton || !halfneigh || (k > 0 || j > 0 || (j == 0 && i > 0)))
          if(bindist(i, j, k) < cutneighsq) {
            stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
          }
      }
    }
  }

  mbins = mbinx * mbiny * mbinz;

  if(bincount) free(bincount);

  bincount = (int*) malloc(mbins * num_omp_threads * sizeof(int));

  if(bins) free(bins);

  bins = (int*) malloc(mbins * num_omp_threads * atoms_per_bin * sizeof(int));
  return 0;
}

/* compute closest distance between central bin (0,0,0) and bin (i,j,k) */

MMD_float Neighbor::bindist(int i, int j, int k)
{
  MMD_float delx, dely, delz;

  if(i > 0)
    delx = (i - 1) * binsizex;
  else if(i == 0)
    delx = 0.0;
  else
    delx = (i + 1) * binsizex;

  if(j > 0)
    dely = (j - 1) * binsizey;
  else if(j == 0)
    dely = 0.0;
  else
    dely = (j + 1) * binsizey;

  if(k > 0)
    delz = (k - 1) * binsizez;
  else if(k == 0)
    delz = 0.0;
  else
    delz = (k + 1) * binsizez;

  return (delx * delx + dely * dely + delz * delz);
}
