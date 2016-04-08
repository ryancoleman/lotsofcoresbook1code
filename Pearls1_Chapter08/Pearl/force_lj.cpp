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
        if(threads->omp_num_threads > 1)
          return compute_halfneigh_threaded<0, 1>(atom, neighbor, me);
        else
          return compute_halfneigh<0, 1>(atom, neighbor, me);
      } else {
        if(threads->omp_num_threads > 1)
          return compute_halfneigh_threaded<0, 0>(atom, neighbor, me);
        else
          return compute_halfneigh<0, 0>(atom, neighbor, me);
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

#define _MM_BCAST_PS(a) _mm512_extload_ps(a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define _MM_BCAST4_PS(a) _mm512_extload_ps(a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)
#define _MM_MASK_BCAST4_PS(v, m, a) _mm512_mask_extload_ps(v, m, a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)
#define _MM_STORE4_PS(a, m, v) _mm512_mask_store_ps(a, m, v)
#define _MM_PACKSTORE4_PS(a, m, v) _mm512_mask_extpackstorelo_ps(a, m, v, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE)
#define _MM_STORE_SS(a, v) _mm512_mask_extpackstorelo_ps(a, 0x1, v, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE)

float inline REDUCE(__m512 val)
{
	float retval;
	__m512 reduce_1 = _mm512_castsi512_ps(_mm512_permute4f128_epi32(_mm512_castps_si512(val), _MM_PERM_CDAB));
	reduce_1 = _mm512_add_ps(reduce_1, val);
	__m512 reduce_2 = _mm512_castsi512_ps(_mm512_permute4f128_epi32(_mm512_castps_si512(reduce_1), _MM_PERM_AACC));
	reduce_1 = _mm512_add_ps(reduce_1, reduce_2);
	reduce_1 = _mm512_add_ps(reduce_1, _mm512_swizzle_ps(reduce_1, _MM_SWIZ_REG_CDAB));
	reduce_1 = _mm512_add_ps(reduce_1, _mm512_swizzle_ps(reduce_1, _MM_SWIZ_REG_BADC));
	_MM_STORE_SS(&retval, reduce_1);
	return retval;
}

// HPPG
//optimised version of compute
//  -MPI + OpenMP (separates local/remote neighbors to avoid conflicts)
//  -use temporary variable for summing up fi
//  -enables vectorization by:
//    -getting rid of 2d pointers
//    -use pragma simd to force vectorization of inner loop (not currently supported due to OpenMP atomics
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

    int k = 0;
#ifdef USE_SIMD
    #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for (; k < numloc; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - xs[j].x;
      const MMD_float dely = ytmp - xs[j].y;
      const MMD_float delz = ztmp - xs[j].z;
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutforcesq) {
        const MMD_float sr2 = 1.0f / rsq;
        const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
        const MMD_float force = 48.0f * sr6 * (sr6 - 0.5f) * sr2 * epsilon;

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
    }
#ifdef USE_SIMD
      #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for (; k < numloc + numrem; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - xs[j].x;
      const MMD_float dely = ytmp - xs[j].y;
      const MMD_float delz = ztmp - xs[j].z;
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutforcesq) {
        const MMD_float sr2 = 1.0f / rsq;
        const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
        const MMD_float force = 48.0f * sr6 * (sr6 - 0.5f) * sr2 * epsilon;

        fix += delx * force;
        fiy += dely * force;
        fiz += delz * force;

        if (EVFLAG) {
          t_eng_vdwl += 0.5f * (4.0f * sr6 * (sr6 - 1.0f)) * epsilon;
          t_virial += 0.5f * (delx * delx + dely * dely + delz * delz) * force;
        }
      }
    }

    fs[i].x += fix;
    fs[i].y += fiy;
    fs[i].z += fiz;

  }
  
  // Reduction and accumulate to eng/virial.
  if (EVFLAG)
  {
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


