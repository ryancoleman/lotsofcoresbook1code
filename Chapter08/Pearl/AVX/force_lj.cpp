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
#include "math.h"
#include "force_lj.h"
#include "openmp.h"

#define AVX
#include <immintrin.h>

#ifndef VECTORLENGTH
#define VECTORLENGTH 4
#endif

ForceLJ::ForceLJ()
{
  cutforce = 0.0;
  cutforcesq = 0.0;
  use_oldcompute = 0;
  reneigh = 1;
  style = FORCELJ;

  epsilon = 1.0;
  sigma6 = 1.0;
  sigma = 1.0;

}
ForceLJ::~ForceLJ() {}

void ForceLJ::setup()
{
  cutforcesq = cutforce * cutforce;
}


void ForceLJ::compute(Atom &atom, Neighbor &neighbor, Comm &comm, int me)
{
  eng_vdwl = 0;
  virial = 0;

  if(evflag) {
    if(use_oldcompute)
      return compute_original<1>(atom, neighbor, me);

    if(neighbor.halfneigh) {
      if(neighbor.ghost_newton) {
          return compute_halfneigh_threaded<1, 1>(atom, neighbor, me);
      } else {
          return compute_halfneigh_threaded<1, 0>(atom, neighbor, me);
      }
    } else return compute_fullneigh<1>(atom, neighbor, me);
  } else {
    if(use_oldcompute)
      return compute_original<0>(atom, neighbor, me);

    if(neighbor.halfneigh) {
      if(neighbor.ghost_newton) {
          return compute_halfneigh_threaded<0, 1>(atom, neighbor, me);
      } else {
          return compute_halfneigh_threaded<0, 0>(atom, neighbor, me);
      }
    } else return compute_fullneigh<0>(atom, neighbor, me);

  }
}

//original version of force compute in miniMD
//  -MPI only
//  -not vectorizable
template<int EVFLAG>
void ForceLJ::compute_original(Atom &atom, Neighbor &neighbor, int me)
{
  int i, j, k, nlocal, nall, numneigh;
  MMD_float xtmp, ytmp, ztmp, delx, dely, delz, rsq;
  MMD_float sr2, sr6, force;
  int* neighs;
  MMD_float** x, **f;

  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;
  x = atom.x;
  f = atom.f;

  eng_vdwl = 0;
  virial = 0;
  // clear force on own and ghost atoms

  for(i = 0; i < nall; i++) {
    f[i][0] = 0.0;
    f[i][1] = 0.0;
    f[i][2] = 0.0;
  }

  // loop over all neighbors of my atoms
  // store force on both atoms i and j

  for(i = 0; i < nlocal; i++) {
    neighs = &neighbor.neighbors[i * neighbor.maxneighs];
    numneigh = neighbor.numneigh[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    for(k = 0; k < numneigh; k++) {
      j = neighs[k];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if(rsq < cutforcesq) {
        sr2 = 1.0 / rsq;
        sr6 = sr2 * sr2 * sr2 * sigma6;
        force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
        f[i][0] += delx * force;
        f[i][1] += dely * force;
        f[i][2] += delz * force;
        f[j][0] -= delx * force;
        f[j][1] -= dely * force;
        f[j][2] -= delz * force;

        if(EVFLAG) {
          eng_vdwl += (4.0 * sr6 * (sr6 - 1.0)) * epsilon;
          virial += (delx * delx + dely * dely + delz * delz) * force;
        }
      }
    }
  }
}


//optimised version of compute
//  -MPI only
//  -use temporary variable for summing up fi
//  -enables vectorization by:
//     -getting rid of 2d pointers
//     -use pragma simd to force vectorization of inner loop
template<int EVFLAG, int GHOST_NEWTON>
void ForceLJ::compute_halfneigh(Atom &atom, Neighbor &neighbor, int me)
{
  int* neighs;
  int tid = omp_get_thread_num();

  const int nlocal = atom.nlocal;
  const int nall = atom.nlocal + atom.nghost;
  MMD_float* x = &atom.x[0][0];
  MMD_float* f = &atom.f[0][0];

  // clear force on own and ghost atoms
  for(int i = 0; i < nall; i++) {
    f[i * PAD + 0] = 0.0;
    f[i * PAD + 1] = 0.0;
    f[i * PAD + 2] = 0.0;
  }

  // loop over all neighbors of my atoms
  // store force on both atoms i and j
  MMD_float t_energy = 0;
  MMD_float t_virial = 0;

  for(int i = 0; i < nlocal; i++) {
    neighs = &neighbor.neighbors[i * neighbor.maxneighs];
    const int numneighs = neighbor.numneigh[i];
    const MMD_float xtmp = x[i * PAD + 0];
    const MMD_float ytmp = x[i * PAD + 1];
    const MMD_float ztmp = x[i * PAD + 2];

    MMD_float fix = 0.0;
    MMD_float fiy = 0.0;
    MMD_float fiz = 0.0;

#ifdef USE_SIMD
    #pragma simd reduction (+: fix,fiy,fiz)
#endif
    for(int k = 0; k < numneighs; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - x[j * PAD + 0];
      const MMD_float dely = ytmp - x[j * PAD + 1];
      const MMD_float delz = ztmp - x[j * PAD + 2];
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;

      if(rsq < cutforcesq) {
        const MMD_float sr2 = 1.0 / rsq;
        const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
        const MMD_float force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;

        fix += delx * force;
        fiy += dely * force;
        fiz += delz * force;

        if(GHOST_NEWTON || j < nlocal) {
          f[j * PAD + 0] -= delx * force;
          f[j * PAD + 1] -= dely * force;
          f[j * PAD + 2] -= delz * force;
        }

        if(EVFLAG) {
          const MMD_float scale = (GHOST_NEWTON || j < nlocal) ? 1.0 : 0.5;
          t_energy += scale * (4.0 * sr6 * (sr6 - 1.0)) * epsilon;
          t_virial += scale * (delx * delx + dely * dely + delz * delz) * force;
        }

      }
    }

    f[i * PAD + 0] += fix;
    f[i * PAD + 1] += fiy;
    f[i * PAD + 2] += fiz;

  }

  eng_vdwl += t_energy;
  virial += t_virial;

}

inline void AOS_TO_SOA_GATHER(float4* restrict array, __m256* xout, __m256* yout, __m256* zout, int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8) {
    
    __m256 in15 = _mm256_castps128_ps256(_mm_load_ps((float*) &array[j1]));
    __m256 in26 = _mm256_castps128_ps256(_mm_load_ps((float*) &array[j2]));
    __m256 in37 = _mm256_castps128_ps256(_mm_load_ps((float*) &array[j3]));
    __m256 in48 = _mm256_castps128_ps256(_mm_load_ps((float*) &array[j4]));
    
    in15 = _mm256_insertf128_ps(in15, _mm_load_ps((float*) &array[j5]), 1);
    in26 = _mm256_insertf128_ps(in26, _mm_load_ps((float*) &array[j6]), 1);
    in37 = _mm256_insertf128_ps(in37, _mm_load_ps((float*) &array[j7]), 1);
    in48 = _mm256_insertf128_ps(in48, _mm_load_ps((float*) &array[j8]), 1);
    
    __m256 xy1256 = _mm256_shuffle_ps(in15, in26, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 xy3478 = _mm256_shuffle_ps(in37, in48, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 z01256 = _mm256_shuffle_ps(in15, in26, _MM_SHUFFLE(3, 2, 3, 2));
    __m256 z03478 = _mm256_shuffle_ps(in37, in48, _MM_SHUFFLE(3, 2, 3, 2));
    
    *xout = _mm256_shuffle_ps(xy1256, xy3478, _MM_SHUFFLE(2, 0, 2, 0));
    *yout = _mm256_shuffle_ps(xy1256, xy3478, _MM_SHUFFLE(3, 1, 3, 1));
    *zout = _mm256_shuffle_ps(z01256, z03478, _MM_SHUFFLE(2, 0, 2, 0));    
    
}

inline void SOA_TO_AOS_SCATTER(float4* restrict array, __m256 xin, __m256 yin, __m256 zin, int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8) {

    const __m256 zeroes = _mm256_setzero_ps();
    
    __m256 xy1256 = _mm256_shuffle_ps(xin, yin, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 xy3478 = _mm256_shuffle_ps(xin, yin, _MM_SHUFFLE(3, 2, 3, 2));
    __m256 z01256 = _mm256_shuffle_ps(zin, zeroes, _MM_SHUFFLE(1, 0, 1, 0));
    __m256 z03478 = _mm256_shuffle_ps(zin, zeroes, _MM_SHUFFLE(3, 2, 3, 2));
    
    __m256 out15 = _mm256_shuffle_ps(xy1256, z01256, _MM_SHUFFLE(2, 0, 2, 0));
    __m256 out26 = _mm256_shuffle_ps(xy1256, z01256, _MM_SHUFFLE(3, 1, 3, 1));
    __m256 out37 = _mm256_shuffle_ps(xy3478, z03478, _MM_SHUFFLE(2, 0, 2, 0));
    __m256 out48 = _mm256_shuffle_ps(xy3478, z03478, _MM_SHUFFLE(3, 1, 3, 1));
    
    _mm_store_ps((float*) &array[j1], _mm256_castps256_ps128(out15));
    _mm_store_ps((float*) &array[j2], _mm256_castps256_ps128(out26));
    _mm_store_ps((float*) &array[j3], _mm256_castps256_ps128(out37));
    _mm_store_ps((float*) &array[j4], _mm256_castps256_ps128(out48));
    
    _mm_store_ps((float*) &array[j5], _mm256_extractf128_ps(out15, 1));
    _mm_store_ps((float*) &array[j6], _mm256_extractf128_ps(out26, 1));
    _mm_store_ps((float*) &array[j7], _mm256_extractf128_ps(out37, 1));
    _mm_store_ps((float*) &array[j8], _mm256_extractf128_ps(out48, 1));   
        
}

inline float REDUCE(__m256 val) {

    float retval;
    __m128 xmm0, xmm1;
    xmm0 = _mm256_extractf128_ps(val, 1);
    xmm1 = _mm_add_ps(_mm256_castps256_ps128(val), xmm0);
    xmm0 = _mm_movehl_ps(xmm0, xmm1);
    xmm1 = _mm_add_ps(xmm0, xmm1);
    xmm0 = _mm_movehdup_ps(xmm1);
    xmm1 = _mm_add_ps(xmm0, xmm1);
    _mm_store_ss(&retval, xmm1);
    return retval;
  
}

// HPPG optimized version of compute
//  -MPI + OpenMP (separates local/remote neighbors to avoid conflicts)
//  -use temporary variable for summing up fi
//  -enables vectorization by:
//    -getting rid of 2d pointers
//    -use intrinsics to force vectorization of inner loop
template<int EVFLAG, int GHOST_NEWTON>
void ForceLJ::compute_halfneigh_threaded(Atom &atom, Neighbor &neighbor, int me)
{
  int nlocal, nall;
  int* neighs;
  MMD_float* x, *f;
  int tid = omp_get_thread_num();

  MMD_float t_eng_vdwl = 0;
  MMD_float t_virial = 0;

  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;
  x = &atom.x[0][0];
  f = &atom.f[0][0];

  float4* xs = (float4*) x;
  float4* fs = (float4*) f;

#ifdef AVX
  // Useful constants.
  __m256 mcutforcesq = _mm256_set1_ps(cutforcesq);
  __m256 zeroes = _mm256_set1_ps(0.0f);
  __m256 ones = _mm256_set1_ps(1.0f);
  __m256 halfs = _mm256_set1_ps(0.5f);
  __m256 fours = _mm256_set1_ps(4.0f);
  __m256 m48 = _mm256_set1_ps(48.0f);
  __m256 msigma6 = _mm256_set1_ps(sigma6);
  __m256 mepsilon = _mm256_set1_ps(epsilon);
#endif  

  #pragma omp barrier
  // clear force on own and ghost atoms

  OMPFORSCHEDULE
  for(int i = 0; i < nall; i++) {
    f[i * PAD + 0] = 0.0;
    f[i * PAD + 1] = 0.0;
    f[i * PAD + 2] = 0.0;
  }

  // loop over all neighbors of my atoms
  // store force on both atoms i and j
  int start_atom = threads->data[tid].start_atom;
  int end_atom = threads->data[tid].end_atom;
#ifdef AVX
  __m256 m_eng_vdwl = zeroes;
  __m256 m_virial = zeroes;
#endif  
  for(int i = start_atom; i < end_atom; i++) {
    neighs = &neighbor.neighbors[i * neighbor.maxneighs];
    const int numloc = neighbor.numneigh[i];
    const int numrem = neighbor.numrem[i];
    const MMD_float xtmp = x[i * PAD + 0];
    const MMD_float ytmp = x[i * PAD + 1];
    const MMD_float ztmp = x[i * PAD + 2];
    MMD_float fix = 0.0;
    MMD_float fiy = 0.0;
    MMD_float fiz = 0.0;
	
#ifdef AVX	
	// Broadcast x/y/z to all 8 AVX lanes.
	__m256 xitmp = _mm256_broadcast_ss(&xs[i].x);
	__m256 yitmp = _mm256_broadcast_ss(&xs[i].y);
	__m256 zitmp = _mm256_broadcast_ss(&xs[i].z);
	__m256 fxitmp = zeroes;
	__m256 fyitmp = zeroes;
	__m256 fzitmp = zeroes;	
#endif

    int k = 0;
#ifdef AVX
	int loopbound = (numloc/8)*8;
	for(; k < loopbound; k += 8) {
	
		// Load neighbor indices.
        unsigned long long j1 = neighs[k+0];
        unsigned long long j2 = neighs[k+1];
        unsigned long long j3 = neighs[k+2];
        unsigned long long j4 = neighs[k+3];
        unsigned long long j5 = neighs[k+4];
        unsigned long long j6 = neighs[k+5];
        unsigned long long j7 = neighs[k+6];
        unsigned long long j8 = neighs[k+7];	
		
		// Load 8 neighbors in AoS, and transpose their positions to SoA.
		__m256 xjtmp, yjtmp, zjtmp;
		AOS_TO_SOA_GATHER(xs, &xjtmp, &yjtmp, &zjtmp, j1, j2, j3, j4, j5, j6, j7, j8);
		
		// Compute squared distance.
        __m256 delx = _mm256_sub_ps(xitmp, xjtmp);
        __m256 dely = _mm256_sub_ps(yitmp, yjtmp);
        __m256 delz = _mm256_sub_ps(zitmp, zjtmp);
        __m256 delxsq = _mm256_mul_ps(delx, delx);
        __m256 delysq = _mm256_mul_ps(dely, dely);
        __m256 delzsq = _mm256_mul_ps(delz, delz);
        __m256 rsq = _mm256_add_ps(delxsq, _mm256_add_ps(delysq, delzsq));
		
		// Compute force.
        __m256 sr2 = _mm256_div_ps(ones, rsq);     
        __m256 sr6 = _mm256_mul_ps(_mm256_mul_ps(sr2, _mm256_mul_ps(sr2, sr2)), msigma6);
        __m256 F = _mm256_mul_ps(m48, _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(sr6, _mm256_sub_ps(sr6, halfs)), sr2), mepsilon));
        F = _mm256_blendv_ps(zeroes, F, _mm256_cmp_ps(rsq, mcutforcesq, _CMP_LT_OS));
		
		// Update force[i]
		__m256 Fx = _mm256_mul_ps(F, delx);
        __m256 Fy = _mm256_mul_ps(F, dely);
        __m256 Fz = _mm256_mul_ps(F, delz);
        fxitmp = _mm256_add_ps(fxitmp, Fx);
        fyitmp = _mm256_add_ps(fyitmp, Fy);
        fzitmp = _mm256_add_ps(fzitmp, Fz);        
		
		// Gather force[j], update force[j] and scatter force[j].
		// Requires an AoS->SoA transpose and an SoA->AoS transpose.   
		__m256 fxjtmp, fyjtmp, fzjtmp;
		AOS_TO_SOA_GATHER(fs, &fxjtmp, &fyjtmp, &fzjtmp, j1, j2, j3, j4, j5, j6, j7, j8);
		fxjtmp = _mm256_sub_ps(fxjtmp, Fx);
        fyjtmp = _mm256_sub_ps(fyjtmp, Fy);
        fzjtmp = _mm256_sub_ps(fzjtmp, Fz);    
		SOA_TO_AOS_SCATTER(fs, fxjtmp, fyjtmp, fzjtmp, j1, j2, j3, j4, j5, j6, j7, j8);
		
		// Update energy/virial.
		if (EVFLAG) {
			m_eng_vdwl = _mm256_add_ps(m_eng_vdwl, _mm256_mul_ps(_mm256_mul_ps(fours, _mm256_mul_ps(sr6, _mm256_sub_ps(sr6, ones))),mepsilon));
			m_virial = _mm256_add_ps(m_virial, _mm256_mul_ps(rsq, F));
		}
		
	}
#endif
#ifdef USE_SIMD
    #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for (; k < numloc; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - xs[j].x;
      const MMD_float dely = ytmp - xs[j].y;
      const MMD_float delz = ztmp - xs[j].z;
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;

	  const MMD_float sr2 = 1.0f / rsq;
	  const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
	  const MMD_float force = (rsq < cutforcesq) ? 48.0f * sr6 * (sr6 - 0.5f) * sr2 * epsilon : 0.0f;

	  fix += delx * force;
	  fiy += dely * force;
	  fiz += delz * force;

	  fs[j].x -= delx * force;
	  fs[j].y -= dely * force;
	  fs[j].z -= delz * force;

	  if (EVFLAG) {
  	    t_eng_vdwl += (4.0f * sr6 * (sr6 - 1.0f)) * epsilon;
	    t_virial += (delx * delx + dely * dely + delz * delz) * force;
	  }
    }
#ifdef AVX
	loopbound = numloc + (numrem/8)*8;
	for(; k < loopbound; k += 8) {
	
		// Load neighbor indices.
        unsigned long long j1 = neighs[k+0];
        unsigned long long j2 = neighs[k+1];
        unsigned long long j3 = neighs[k+2];
        unsigned long long j4 = neighs[k+3];
        unsigned long long j5 = neighs[k+4];
        unsigned long long j6 = neighs[k+5];
        unsigned long long j7 = neighs[k+6];
        unsigned long long j8 = neighs[k+7];	
		
		// Load 8 neighbors in AoS, and transpose their positions to SoA.
		__m256 xjtmp, yjtmp, zjtmp;
		AOS_TO_SOA_GATHER(xs, &xjtmp, &yjtmp, &zjtmp, j1, j2, j3, j4, j5, j6, j7, j8);
		
		// Compute squared distance.
        __m256 delx = _mm256_sub_ps(xitmp, xjtmp);
        __m256 dely = _mm256_sub_ps(yitmp, yjtmp);
        __m256 delz = _mm256_sub_ps(zitmp, zjtmp);
        __m256 delxsq = _mm256_mul_ps(delx, delx);
        __m256 delysq = _mm256_mul_ps(dely, dely);
        __m256 delzsq = _mm256_mul_ps(delz, delz);
        __m256 rsq = _mm256_add_ps(delxsq, _mm256_add_ps(delysq, delzsq));
		
		// Compute force.
        __m256 sr2 = _mm256_div_ps(ones, rsq);     
        __m256 sr6 = _mm256_mul_ps(_mm256_mul_ps(sr2, _mm256_mul_ps(sr2, sr2)), msigma6);
        __m256 F = _mm256_mul_ps(m48, _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(sr6, _mm256_sub_ps(sr6, halfs)), sr2), mepsilon));
        F = _mm256_blendv_ps(zeroes, F, _mm256_cmp_ps(rsq, mcutforcesq, _CMP_LT_OS));
		
		// Update force[i]
		__m256 Fx = _mm256_mul_ps(F, delx);
        __m256 Fy = _mm256_mul_ps(F, dely);
        __m256 Fz = _mm256_mul_ps(F, delz);
        fxitmp = _mm256_add_ps(fxitmp, Fx);
        fyitmp = _mm256_add_ps(fyitmp, Fy);
        fzitmp = _mm256_add_ps(fzitmp, Fz);
		
		// Update energy/virial.
		if (EVFLAG) {
			m_eng_vdwl = _mm256_add_ps(m_eng_vdwl, _mm256_mul_ps(halfs, _mm256_mul_ps(_mm256_mul_ps(fours, _mm256_mul_ps(sr6, _mm256_sub_ps(sr6, ones))),mepsilon)));
			m_virial = _mm256_add_ps(m_virial, _mm256_mul_ps(halfs, _mm256_mul_ps(rsq, F)));
		}
		
	}
#endif	
#ifdef USE_SIMD
      #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for (; k < numloc + numrem; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - xs[j].x;
      const MMD_float dely = ytmp - xs[j].y;
      const MMD_float delz = ztmp - xs[j].z;
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;

      const MMD_float sr2 = 1.0f / rsq;
      const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
      const MMD_float force = (rsq < cutforcesq) ? 48.0f * sr6 * (sr6 - 0.5f) * sr2 * epsilon : 0.0f;

      fix += delx * force;
      fiy += dely * force;
      fiz += delz * force;

      if (EVFLAG) {
        t_eng_vdwl += 0.5f * (4.0f * sr6 * (sr6 - 1.0f)) * epsilon;
        t_virial += 0.5f * (delx * delx + dely * dely + delz * delz) * force;
      }
    }
	
#ifdef AVX
	// Include reduced force contribution from all SIMD lanes.
	fix += REDUCE(fxitmp);
	fiy += REDUCE(fyitmp);
	fiz += REDUCE(fzitmp);
#endif
    fs[i].x += fix;
    fs[i].y += fiy;
    fs[i].z += fiz;

  }
  
  // Reduction and accumulate to eng/virial.
  if (EVFLAG)
  {
#ifdef AVX
	t_eng_vdwl += REDUCE(m_eng_vdwl);
	t_virial += REDUCE(m_virial);
#endif
    #pragma omp atomic
    eng_vdwl += t_eng_vdwl;
    #pragma omp atomic
    virial += t_virial;
  }

  #pragma omp barrier
}

//optimised version of compute
//  -MPI + OpenMP (using full neighborlists)
//  -gets rid of fj update (read/write to memory)
//  -use temporary variable for summing up fi
//  -enables vectorization by:
//    -get rid of 2d pointers
//    -use pragma simd to force vectorization of inner loop
template<int EVFLAG>
void ForceLJ::compute_fullneigh(Atom &atom, Neighbor &neighbor, int me)
{
  int nlocal, nall;
  int* neighs;
  MMD_float* x, *f;
  int tid = omp_get_thread_num();

  MMD_float t_eng_vdwl = 0;
  MMD_float t_virial = 0;
  nlocal = atom.nlocal;
  nall = atom.nlocal + atom.nghost;
  x = &atom.x[0][0];
  f = &atom.f[0][0];

  #pragma omp barrier
  // clear force on own and ghost atoms

  OMPFORSCHEDULE
  for(int i = 0; i < nlocal; i++) {
    f[i * PAD + 0] = 0.0;
    f[i * PAD + 1] = 0.0;
    f[i * PAD + 2] = 0.0;
  }

  // loop over all neighbors of my atoms
  // store force on atom i

  OMPFORSCHEDULE
  for(int i = 0; i < nlocal; i++) {
    neighs = &neighbor.neighbors[i * neighbor.maxneighs];
    const int numneighs = neighbor.numneigh[i];
    const MMD_float xtmp = x[i * PAD + 0];
    const MMD_float ytmp = x[i * PAD + 1];
    const MMD_float ztmp = x[i * PAD + 2];
    MMD_float fix = 0;
    MMD_float fiy = 0;
    MMD_float fiz = 0;

    //pragma simd forces vectorization (ignoring the performance objections of the compiler)
    //also give hint to use certain vectorlength for MIC, Sandy Bridge and WESTMERE this should be be 8 here
    //give hint to compiler that fix, fiy and fiz are used for reduction only

#ifdef USE_SIMD
    #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for(int k = 0; k < numneighs; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - x[j * PAD + 0];
      const MMD_float dely = ytmp - x[j * PAD + 1];
      const MMD_float delz = ztmp - x[j * PAD + 2];
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;
      if(rsq < cutforcesq) {
        const MMD_float sr2 = 1.0 / rsq;
        const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
        const MMD_float force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
        fix += delx * force;
        fiy += dely * force;
        fiz += delz * force;

        if(EVFLAG) {
          t_eng_vdwl += sr6 * (sr6 - 1.0) * epsilon;
          t_virial += (delx * delx + dely * dely + delz * delz) * force;
        }
      }
      
    }

    f[i * PAD + 0] += fix;
    f[i * PAD + 1] += fiy;
    f[i * PAD + 2] += fiz;

  }

  t_eng_vdwl *= 4.0;
  t_virial *= 0.5;

  #pragma omp atomic
  eng_vdwl += t_eng_vdwl;
  #pragma omp atomic
  virial += t_virial;
  #pragma omp barrier
}


