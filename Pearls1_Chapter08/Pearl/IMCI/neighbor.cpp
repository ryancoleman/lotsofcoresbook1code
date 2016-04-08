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
#include <assert.h>
#include <immintrin.h>

#include "neighbor.h"
#include "openmp.h"

#define FACTOR 0.999
#define SMALL 1.0e-6
#define IMCI 

#define PRINT_IMCI(v) \
{ \
  __declspec(align(64)) float tmp[16]; \
  _mm512_store_ps(tmp, v); \
  printf("%s: %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", #v, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], tmp[9], tmp[10], tmp[11], tmp[12], tmp[13], tmp[14], tmp[15]); \
} 

#define PRINT_IMCI_EPI32(v) \
{ \
  __declspec(align(64)) int tmp[16]; \
  _mm512_store_epi32(tmp, v); \
  printf("%s: %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", #v, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], tmp[9], tmp[10], tmp[11], tmp[12], tmp[13], tmp[14], tmp[15]); \
} 

Neighbor::Neighbor()
{
  ncalls = 0;
  max_totalneigh = 0;
  numneigh = NULL;
  numrem = NULL;
  neighbors = NULL;
  maxneighs = 90;
  maxneighs = (maxneighs + 15) & ~15;
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
  inline bool operator() (int i, int j)
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
  partitioner.start = start_atom*4;
  partitioner.end = end_atom*4;

  #pragma omp barrier
  /* bin local & ghost atoms */

  binatoms(atom);
  count = 0;
  /* loop over each atom, storing neighbors */

  const MMD_float* x = &atom.x[0][0];

#ifdef IMCI 
  __m512 mcutneighsq = _mm512_set_1to16_ps(cutneighsq);
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
#define TRANSPOSE_CACHE
#ifdef TRANSPOSE_CACHE	
    int cached_bin = -1; // we cache the stencil for this bin
    int ncache = 0;
    int CACHE_SIZE = 1024; // hard-code the size of our cache
    if (!threads->data[tid].stencil_cache) {
      threads->data[tid].stencil_cache =  (MMD_float*) _mm_malloc(4 * CACHE_SIZE * sizeof(MMD_float), 64);
    }
    MMD_float* stencil_cache = threads->data[tid].stencil_cache;
#endif
    for (int i = start_atom; i < end_atom; i++) {

      int* neighptr = &neighbors[i * maxneighs];

      int n = 0;
      int local = 0;
      int remote = 0;

      const MMD_float xtmp = x[i * PAD + 0];
      const MMD_float ytmp = x[i * PAD + 1];
      const MMD_float ztmp = x[i * PAD + 2];

#ifdef IMCI 
      __m512i mi = _mm512_set_1to16_epi32(i);
      __m512i mstart = _mm512_set_1to16_epi32(start_atom);
      __m512i mend = _mm512_set_1to16_epi32(end_atom-1); // need to get around lack of cmpge/cmple
      __m512 mxtmp = _mm512_set_1to16_ps(xtmp);
      __m512 mytmp = _mm512_set_1to16_ps(ytmp);
      __m512 mztmp = _mm512_set_1to16_ps(ztmp);
	  __m512 zeroes = _mm512_set_1to16_ps(0.0f);
#endif

      const int ibin = coord2bin(xtmp, ytmp, ztmp);

#ifdef TRANSPOSE_CACHE
      // If we encounter a new bin, cache its contents.
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
#ifdef IMCI 
      for (; c < (ncache/16)*16; c += 16)
      {
        const __m512 delx = _mm512_sub_ps(mxtmp, _mm512_load_ps(&stencil_cache[1*CACHE_SIZE + c]));
        const __m512 dely = _mm512_sub_ps(mytmp, _mm512_load_ps(&stencil_cache[2*CACHE_SIZE + c]));
        const __m512 delz = _mm512_sub_ps(mztmp, _mm512_load_ps(&stencil_cache[3*CACHE_SIZE + c]));
        const __m512 rsq = _mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(delx, delx), _mm512_mul_ps(dely, dely)), _mm512_mul_ps(delz, delz));
        __mmask16 cmask = _mm512_cmple_ps_mask(rsq, mcutneighsq);

        __m512i mj = _mm512_load_epi32(&stencil_cache[c]);
        __mmask16 nmask = _mm512_cmpgt_epi32_mask(mj, mi);
        __mmask16 rmask = _mm512_kor(_mm512_cmplt_epi32_mask(mj, mstart), _mm512_cmpgt_epi32_mask(mj, mend));
        __mmask16 lmask = _mm512_kandn(rmask, nmask);
        rmask = _mm512_kand(cmask, rmask);
        lmask = _mm512_kand(cmask, lmask);
        _mm512_mask_extpackstorelo_epi32(&neighptr[n], _mm512_kor(lmask, rmask), _mm512_slli_epi32(mj, 2), _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
        _mm512_mask_extpackstorehi_epi32(((char*)&neighptr[n])+64, _mm512_kor(lmask,rmask), _mm512_slli_epi32(mj, 2), _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
        local += _mm_countbits_32(lmask);
        remote += _mm_countbits_32(rmask);
        n += _mm_countbits_32(lmask) + _mm_countbits_32(rmask);
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
           neighptr[n++] = j*4;
           if (j >= start_atom && j < end_atom) local++;
           else remote++;
        }

      }
#else
      __declspec(align(64)) int ttmp[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
      int sixteen = 16;
      __m512i all_16 = _mm512_extload_epi32(&sixteen, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
      for(int k = 0; k < nstencil; k++) {
        const int jbin = ibin + stencil[k];

        int* loc_bin = &bins[jbin * atoms_per_bin];

		int loopbound = bincount[jbin];
	    __m512i iter_count = _mm512_load_epi32(&ttmp);
	    __m512i limit = _mm512_extload_epi32(&loopbound, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);			
        for (int m = 0; m < loopbound; m += 16) {
			
			// Check if we have less than 16 neighbors remaining.
			__mmask16 valid = _mm512_cmplt_epi32_mask(iter_count, limit);
			__m512i mj = _mm512_mask_extloadunpacklo_epi32(_mm512_undefined_epi32(), valid, &loc_bin[m], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
			mj = _mm512_mask_extloadunpackhi_epi32(mj, valid, &loc_bin[m], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
			iter_count = _mm512_add_epi32(iter_count, all_16);
			
			// Gather the x/y/z data for the neighbors.
			__m512 xjtmp = _mm512_mask_i32gather_ps(zeroes, valid, _mm512_slli_epi32(mj, 2), &x[0], _MM_SCALE_4);
			__m512 yjtmp = _mm512_mask_i32gather_ps(zeroes, valid, _mm512_slli_epi32(mj, 2), &x[1], _MM_SCALE_4);
			__m512 zjtmp = _mm512_mask_i32gather_ps(zeroes, valid, _mm512_slli_epi32(mj, 2), &x[2], _MM_SCALE_4);
			
			// Compute squared distance.
			__m512 delx = _mm512_sub_ps(mxtmp, xjtmp);
			__m512 dely = _mm512_sub_ps(mytmp, yjtmp);
			__m512 delz = _mm512_sub_ps(mztmp, zjtmp);
			__m512 delxsq = _mm512_mul_ps(delx, delx);
			__m512 delysq = _mm512_mul_ps(dely, dely);
			__m512 delzsq = _mm512_mul_ps(delz, delz);
			__m512 rsq = _mm512_add_ps(delxsq, _mm512_add_ps(delysq, delzsq));
			__mmask16 cmask = _mm512_kand(valid, _mm512_cmple_ps_mask(rsq, mcutneighsq));
			
			// Append to neighbor list.
			__mmask16 nmask = _mm512_cmpgt_epi32_mask(mj, mi);
			__mmask16 rmask = _mm512_kor(_mm512_cmplt_epi32_mask(mj, mstart), _mm512_cmpgt_epi32_mask(mj, mend));
			__mmask16 lmask = _mm512_kandn(rmask, nmask);
			rmask = _mm512_kand(cmask, rmask);
			lmask = _mm512_kand(cmask, lmask);
			_mm512_mask_extpackstorelo_epi32(&neighptr[n], _mm512_kor(lmask, rmask), _mm512_slli_epi32(mj, 2), _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
			_mm512_mask_extpackstorehi_epi32(((char*)&neighptr[n])+64, _mm512_kor(lmask,rmask), _mm512_slli_epi32(mj, 2), _MM_DOWNCONV_EPI32_NONE, _MM_HINT_NONE);
			local += _mm_countbits_32(lmask);
			remote += _mm_countbits_32(rmask);
			n += _mm_countbits_32(lmask) + _mm_countbits_32(rmask);
			
        }
      }
#endif

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
		maxneighs = (maxneighs + 15) & ~15;
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
