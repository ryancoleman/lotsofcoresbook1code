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

#define IMCI
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

#define _MM_BCAST_PS(a) _mm512_extload_ps(a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE)
#define _MM_BCAST4_PS(a) _mm512_extload_ps(a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)
#define _MM_MASK_BCAST4_PS(v, m, a) _mm512_mask_extload_ps(v, m, a, _MM_UPCONV_PS_NONE, _MM_BROADCAST_4X16, _MM_HINT_NONE)
#define _MM_STORE4_PS(a, m, v) _mm512_mask_store_ps(a, m, v)
#define _MM_PACKSTORE4_PS(a, m, v) _mm512_mask_extpackstorelo_ps(a, m, v, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE)
#define _MM_STORE_SS(a, v) _mm512_mask_extpackstorelo_ps(a, 0x1, v, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE)

#define PRINT_IMCI(v) \
{ \
  __declspec(align(64)) float tmp[16]; \
  _mm512_store_ps(tmp, v); \
  printf("%s: %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", #v, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[6], tmp[7], tmp[8], tmp[9], tmp[10], tmp[11], tmp[12], tmp[13], tmp[14], tmp[15]); \
} 

#define SWIZZLES
#ifdef SWIZZLES
#define AOS_TO_SOA_GATHER(array, xout, yout, zout, indices) \
{ \
	__m512 in15 = _MM_BCAST4_PS(&array[indices[0]]); \
	__m512 in26 = _MM_BCAST4_PS(&array[indices[1]]); \
	__m512 in37 = _MM_BCAST4_PS(&array[indices[2]]); \
	__m512 in48 = _MM_BCAST4_PS(&array[indices[3]]); \
	if (valid > 0x000F) {\
	in15 = _MM_MASK_BCAST4_PS(in15, mask_00F0, &array[indices[4]]); \
	in26 = _MM_MASK_BCAST4_PS(in26, mask_00F0, &array[indices[5]]); \
	in37 = _MM_MASK_BCAST4_PS(in37, mask_00F0, &array[indices[6]]); \
	in48 = _MM_MASK_BCAST4_PS(in48, mask_00F0, &array[indices[7]]); \
	} if (valid > 0x00FF) {\
	in15 = _MM_MASK_BCAST4_PS(in15, mask_0F00, &array[indices[8]]); \
	in26 = _MM_MASK_BCAST4_PS(in26, mask_0F00, &array[indices[9]]); \
	in37 = _MM_MASK_BCAST4_PS(in37, mask_0F00, &array[indices[10]]); \
	in48 = _MM_MASK_BCAST4_PS(in48, mask_0F00, &array[indices[11]]); \
	} if (valid > 0x0FFF) {\
	in15 = _MM_MASK_BCAST4_PS(in15, mask_F000, &array[indices[12]]); \
	in26 = _MM_MASK_BCAST4_PS(in26, mask_F000, &array[indices[13]]); \
	in37 = _MM_MASK_BCAST4_PS(in37, mask_F000, &array[indices[14]]); \
	in48 = _MM_MASK_BCAST4_PS(in48, mask_F000, &array[indices[15]]); \
	} \
	__m512 xz01 = _mm512_mask_swizzle_ps(in15, _mm512_int2mask(0xAAAA), in26, _MM_SWIZ_REG_CDAB); \
	__m512 yw01 = _mm512_mask_swizzle_ps(in26, _mm512_int2mask(0x5555), in15, _MM_SWIZ_REG_CDAB); \
	__m512 xz23 = _mm512_mask_swizzle_ps(in37, _mm512_int2mask(0xAAAA), in48, _MM_SWIZ_REG_CDAB); \
	__m512 yw23 = _mm512_mask_swizzle_ps(in48, _mm512_int2mask(0x5555), in37, _MM_SWIZ_REG_CDAB); \
	\
	xout = _mm512_mask_swizzle_ps(xz01, _mm512_int2mask(0xCCCC), xz23, _MM_SWIZ_REG_BADC); \
	yout = _mm512_mask_swizzle_ps(yw01, _mm512_int2mask(0xCCCC), yw23, _MM_SWIZ_REG_BADC); \
	zout = _mm512_mask_swizzle_ps(xz23, _mm512_int2mask(0x3333), xz01, _MM_SWIZ_REG_BADC); \
	\
}

/*#define AOS_TO_SOA_GATHER(array, xout, yout, zout, indices) \
{ \
	__m512 in15, in26, in37, in48; \
	switch (valid) \
	{ \
		case 0xFFFF: in48 = _MM_MASK_BCAST4_PS(in48, mask_F000, &array[indices[15]]); \
		case 0x7FFF: in37 = _MM_MASK_BCAST4_PS(in37, mask_F000, &array[indices[14]]); \
		case 0x3FFF: in26 = _MM_MASK_BCAST4_PS(in26, mask_F000, &array[indices[13]]); \
		case 0x1FFF: in15 = _MM_MASK_BCAST4_PS(in15, mask_F000, &array[indices[12]]); \
		case 0x0FFF: in48 = _MM_MASK_BCAST4_PS(in48, mask_0F00, &array[indices[11]]); \
		case 0x07FF: in37 = _MM_MASK_BCAST4_PS(in37, mask_0F00, &array[indices[10]]); \
		case 0x03FF: in26 = _MM_MASK_BCAST4_PS(in26, mask_0F00, &array[indices[9]]); \
		case 0x01FF: in15 = _MM_MASK_BCAST4_PS(in15, mask_0F00, &array[indices[8]]); \
		case 0x00FF: in48 = _MM_MASK_BCAST4_PS(in48, mask_00F0, &array[indices[7]]); \
		case 0x007F: in37 = _MM_MASK_BCAST4_PS(in37, mask_00F0, &array[indices[6]]); \
		case 0x003F: in26 = _MM_MASK_BCAST4_PS(in26, mask_00F0, &array[indices[5]]); \
		case 0x001F: in15 = _MM_MASK_BCAST4_PS(in15, mask_00F0, &array[indices[4]]); \
		case 0x000F: in48 = _MM_MASK_BCAST4_PS(in48, mask_000F, &array[indices[3]]); \
		case 0x0007: in37 = _MM_MASK_BCAST4_PS(in37, mask_000F, &array[indices[2]]); \
		case 0x0003: in26 = _MM_MASK_BCAST4_PS(in26, mask_000F, &array[indices[1]]); \
		case 0x0001: in15 = _MM_MASK_BCAST4_PS(in15, mask_000F, &array[indices[0]]); \
	} \
        __m512 xz01 = _mm512_mask_swizzle_ps(in15, _mm512_int2mask(0xAAAA), in26, _MM_SWIZ_REG_CDAB); \
        __m512 yw01 = _mm512_mask_swizzle_ps(in26, _mm512_int2mask(0x5555), in15, _MM_SWIZ_REG_CDAB); \
        __m512 xz23 = _mm512_mask_swizzle_ps(in37, _mm512_int2mask(0xAAAA), in48, _MM_SWIZ_REG_CDAB); \
        __m512 yw23 = _mm512_mask_swizzle_ps(in48, _mm512_int2mask(0x5555), in37, _MM_SWIZ_REG_CDAB); \
        \
        xout = _mm512_mask_swizzle_ps(xz01, _mm512_int2mask(0xCCCC), xz23, _MM_SWIZ_REG_BADC); \
        yout = _mm512_mask_swizzle_ps(yw01, _mm512_int2mask(0xCCCC), yw23, _MM_SWIZ_REG_BADC); \
        zout = _mm512_mask_swizzle_ps(xz23, _mm512_int2mask(0x3333), xz01, _MM_SWIZ_REG_BADC); \
        \
}*/

#define SOA_TO_AOS_SCATTER(array, xin, yin, zin, indices) \
{ \
	\
	__m512 xz01 = _mm512_mask_swizzle_ps(xin, _mm512_int2mask(0xCCCC), zin, _MM_SWIZ_REG_BADC); \
	__m512 yw01 = _mm512_mask_swizzle_ps(yin, _mm512_int2mask(0xCCCC), zeroes, _MM_SWIZ_REG_BADC); \
	__m512 xz23 = _mm512_mask_swizzle_ps(zin, _mm512_int2mask(0x3333), xin, _MM_SWIZ_REG_BADC); \
	__m512 yw23 = _mm512_mask_swizzle_ps(zeroes, _mm512_int2mask(0x3333), yin, _MM_SWIZ_REG_BADC); \
	\
	__m512 out15 = _mm512_mask_swizzle_ps(xz01, _mm512_int2mask(0xAAAA), yw01, _MM_SWIZ_REG_CDAB); \
	__m512 out26 = _mm512_mask_swizzle_ps(yw01, _mm512_int2mask(0x5555), xz01, _MM_SWIZ_REG_CDAB); \
	__m512 out37 = _mm512_mask_swizzle_ps(xz23, _mm512_int2mask(0xAAAA), yw23, _MM_SWIZ_REG_CDAB); \
	__m512 out48 = _mm512_mask_swizzle_ps(yw23, _mm512_int2mask(0x5555), xz23, _MM_SWIZ_REG_CDAB); \
	\
	_MM_PACKSTORE4_PS(&array[indices[0]], mask_000F, out15); \
	_MM_PACKSTORE4_PS(&array[indices[1]], mask_000F, out26); \
	_MM_PACKSTORE4_PS(&array[indices[2]], mask_000F, out37); \
	_MM_PACKSTORE4_PS(&array[indices[3]], mask_000F, out48); \
	if (valid > 0x000F) {\
	_MM_PACKSTORE4_PS(&array[indices[4]], mask_00F0, out15); \
	_MM_PACKSTORE4_PS(&array[indices[5]], mask_00F0, out26); \
	_MM_PACKSTORE4_PS(&array[indices[6]], mask_00F0, out37); \
	_MM_PACKSTORE4_PS(&array[indices[7]], mask_00F0, out48); \
	} if (valid > 0x00FF) {\
	_MM_PACKSTORE4_PS(&array[indices[8]], mask_0F00, out15); \
	_MM_PACKSTORE4_PS(&array[indices[9]], mask_0F00, out26); \
	_MM_PACKSTORE4_PS(&array[indices[10]], mask_0F00, out37); \
	_MM_PACKSTORE4_PS(&array[indices[11]], mask_0F00, out48); \
	} if (valid > 0x0FFF) {\
	_MM_PACKSTORE4_PS(&array[indices[12]], mask_F000, out15); \
	_MM_PACKSTORE4_PS(&array[indices[13]], mask_F000, out26); \
	_MM_PACKSTORE4_PS(&array[indices[14]], mask_F000, out37); \
	_MM_PACKSTORE4_PS(&array[indices[15]], mask_F000, out48); \
	} \
}

/*#define SOA_TO_AOS_SCATTER(array, xin, yin, zin, indices) \
{ \
        \
        __m512 xz01 = _mm512_mask_swizzle_ps(xin, _mm512_int2mask(0xCCCC), zin, _MM_SWIZ_REG_BADC); \
        __m512 yw01 = _mm512_mask_swizzle_ps(yin, _mm512_int2mask(0xCCCC), zeroes, _MM_SWIZ_REG_BADC); \
        __m512 xz23 = _mm512_mask_swizzle_ps(zin, _mm512_int2mask(0x3333), xin, _MM_SWIZ_REG_BADC); \
        __m512 yw23 = _mm512_mask_swizzle_ps(zeroes, _mm512_int2mask(0x3333), yin, _MM_SWIZ_REG_BADC); \
        \
        __m512 out15 = _mm512_mask_swizzle_ps(xz01, _mm512_int2mask(0xAAAA), yw01, _MM_SWIZ_REG_CDAB); \
        __m512 out26 = _mm512_mask_swizzle_ps(yw01, _mm512_int2mask(0x5555), xz01, _MM_SWIZ_REG_CDAB); \
        __m512 out37 = _mm512_mask_swizzle_ps(xz23, _mm512_int2mask(0xAAAA), yw23, _MM_SWIZ_REG_CDAB); \
        __m512 out48 = _mm512_mask_swizzle_ps(yw23, _mm512_int2mask(0x5555), xz23, _MM_SWIZ_REG_CDAB); \
        \
	switch (valid) \
	{ \
		case 0xFFFF: _MM_PACKSTORE4_PS(&array[indices[15]], mask_F000, out48); \
		case 0x7FFF: _MM_PACKSTORE4_PS(&array[indices[14]], mask_F000, out37); \
		case 0x3FFF: _MM_PACKSTORE4_PS(&array[indices[13]], mask_F000, out26); \
		case 0x1FFF: _MM_PACKSTORE4_PS(&array[indices[12]], mask_F000, out15); \
		case 0x0FFF: _MM_PACKSTORE4_PS(&array[indices[11]], mask_0F00, out48); \
		case 0x07FF: _MM_PACKSTORE4_PS(&array[indices[10]], mask_0F00, out37); \
		case 0x03FF: _MM_PACKSTORE4_PS(&array[indices[9]], mask_0F00, out26); \
		case 0x01FF: _MM_PACKSTORE4_PS(&array[indices[8]], mask_0F00, out15); \
		case 0x00FF: _MM_PACKSTORE4_PS(&array[indices[7]], mask_00F0, out48); \
		case 0x007F: _MM_PACKSTORE4_PS(&array[indices[6]], mask_00F0, out37); \
		case 0x003F: _MM_PACKSTORE4_PS(&array[indices[5]], mask_00F0, out26); \
		case 0x001F: _MM_PACKSTORE4_PS(&array[indices[4]], mask_00F0, out15); \
		case 0x000F: _MM_PACKSTORE4_PS(&array[indices[3]], mask_000F, out48); \
		case 0x0007: _MM_PACKSTORE4_PS(&array[indices[2]], mask_000F, out37); \
		case 0x0003: _MM_PACKSTORE4_PS(&array[indices[1]], mask_000F, out26); \
		case 0x0001: _MM_PACKSTORE4_PS(&array[indices[0]], mask_000F, out15); \
	} \
}*/

#else

#define AOS_TO_SOA_GATHER(array, xout, yout, zout, indices) \
{ \
	__m512 tmp0 = _MM_BCAST4_PS(&array[indices[0]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_00F0, &array[indices[1]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_0F00, &array[indices[2]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_F000, &array[indices[3]]); \
	_mm512_mask_extpackstorelo_ps(&tmp[0], mask_AAAA, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[16], mask_BBBB, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[32], mask_CCCC, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	\
	tmp0 = _MM_BCAST4_PS(&array[indices[4]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_00F0, &array[indices[5]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_0F00, &array[indices[6]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_F000, &array[indices[7]]); \
	_mm512_mask_extpackstorelo_ps(&tmp[4], mask_AAAA, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[20], mask_BBBB, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[36], mask_CCCC, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	\
	tmp0 = _MM_BCAST4_PS(&array[indices[8]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_00F0, &array[indices[9]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_0F00, &array[indices[10]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_F000, &array[indices[11]]); \
	_mm512_mask_extpackstorelo_ps(&tmp[8], mask_AAAA, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[24], mask_BBBB, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[40], mask_CCCC, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	\
	tmp0 = _MM_BCAST4_PS(&array[indices[12]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_00F0, &array[indices[13]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_0F00, &array[indices[14]]); \
	tmp0 = _MM_MASK_BCAST4_PS(tmp0, mask_F000, &array[indices[15]]); \
	_mm512_mask_extpackstorelo_ps(&tmp[12], mask_AAAA, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[28], mask_BBBB, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	_mm512_mask_extpackstorelo_ps(&tmp[44], mask_CCCC, tmp0, _MM_DOWNCONV_PS_NONE, _MM_HINT_NONE); \
	\
	xout = _mm512_load_ps(&tmp[0]); \
	yout = _mm512_load_ps(&tmp[16]); \
	zout = _mm512_load_ps(&tmp[32]); \
}

#define SOA_TO_AOS_SCATTER(array, xin, yin, zin, indices) \
{ \
	_mm512_store_ps(&tmp[0], xin); \
	_mm512_store_ps(&tmp[16], yin); \
	_mm512_store_ps(&tmp[32], zin); \
	\
	__m512 tmp0 = _mm512_mask_extloadunpacklo_ps(zeroes, mask_AAAA, &tmp[0], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_BBBB, &tmp[16], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_CCCC, &tmp[32], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	_MM_PACKSTORE4_PS(&array[indices[0]], mask_000F, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[1]], mask_00F0, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[2]], mask_0F00, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[3]], mask_F000, tmp0); \
	\
	tmp0 = _mm512_mask_extloadunpacklo_ps(zeroes, mask_AAAA, &tmp[4], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_BBBB, &tmp[20], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_CCCC, &tmp[36], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	_MM_PACKSTORE4_PS(&array[indices[4]], mask_000F, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[5]], mask_00F0, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[6]], mask_0F00, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[7]], mask_F000, tmp0); \
	\
	tmp0 = _mm512_mask_extloadunpacklo_ps(zeroes, mask_AAAA, &tmp[8], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_BBBB, &tmp[24], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_CCCC, &tmp[40], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	_MM_PACKSTORE4_PS(&array[indices[8]], mask_000F, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[9]], mask_00F0, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[10]], mask_0F00, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[11]], mask_F000, tmp0); \
	\
	tmp0 = _mm512_mask_extloadunpacklo_ps(zeroes, mask_AAAA, &tmp[12], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_BBBB, &tmp[28], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	tmp0 = _mm512_mask_extloadunpacklo_ps(tmp0, mask_CCCC, &tmp[44], _MM_UPCONV_PS_NONE, _MM_HINT_NONE); \
	_MM_PACKSTORE4_PS(&array[indices[12]], mask_000F, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[13]], mask_00F0, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[14]], mask_0F00, tmp0); \
	_MM_PACKSTORE4_PS(&array[indices[15]], mask_F000, tmp0); \
}

#endif

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

  float* xs = (float*) x;
  float* fs = (float*) f;

#ifdef IMCI
  // Useful constants.
  __declspec(align(64)) float tmp[64];
  const __mmask mask_AAAA = 0x1111;
  const __mmask mask_BBBB = 0x2222;
  const __mmask mask_CCCC = 0x4444;
  const __mmask mask_DDDD = 0x8888;
  const __mmask mask_000F = 0x000F;
  const __mmask mask_00F0 = 0x00F0;
  const __mmask mask_0F00 = 0x0F00;
  const __mmask mask_F000 = 0xF000;  
  __m512 mcutforcesq = _mm512_set_1to16_ps(cutforcesq);
  __m512 zeroes = _mm512_set_1to16_ps(0.0f);
  __m512 ones = _mm512_set_1to16_ps(1.0f);
  __m512 halfs = _mm512_set_1to16_ps(0.5f);
  __m512 fours = _mm512_set_1to16_ps(4.0f);
  __m512 m48 = _mm512_set_1to16_ps(48.0f);
  __m512 msigma6 = _mm512_set_1to16_ps(sigma6);
  __m512 mepsilon = _mm512_set_1to16_ps(epsilon);
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
#ifdef IMCI
  __m512 m_eng_vdwl = zeroes;
  __m512 m_virial = zeroes;
  __declspec(align(64)) int ttmp[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  int sixteen = 16;
  __m512i all_16 = _mm512_extload_epi32(&sixteen, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);  
#endif  
  int i4 = start_atom*4;
  for(int i = start_atom; i < end_atom; i++,i4+=4) {
    neighs = &neighbor.neighbors[i * neighbor.maxneighs];
    const int numloc = neighbor.numneigh[i];
    const int numrem = neighbor.numrem[i];
    const MMD_float xtmp = x[i*4 + 0];
    const MMD_float ytmp = x[i*4 + 1];
    const MMD_float ztmp = x[i*4 + 2];
    MMD_float fix = 0.0;
    MMD_float fiy = 0.0;
    MMD_float fiz = 0.0;
	
#ifdef IMCI	
	// Broadcast x/y/z to all 16 IMCI lanes.
	__m512i mi = _mm512_slli_epi32(_mm512_extload_epi32(&i, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE), 2);
	__m512 xitmp = _MM_BCAST_PS(&xs[i*4+0]);
	__m512 yitmp = _MM_BCAST_PS(&xs[i*4+1]);
	__m512 zitmp = _MM_BCAST_PS(&xs[i*4+2]);
	__m512 fxitmp = zeroes;
	__m512 fyitmp = zeroes;
	__m512 fzitmp = zeroes;	
#endif

    int k = 0;
#ifdef IMCI
	int loopbound = numloc;
	__m512i iter_count = _mm512_add_epi32(_mm512_load_epi32(&ttmp), _mm512_set1_epi32(k));
	__m512i limit = _mm512_extload_epi32(&loopbound, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);	
	for(; k < loopbound; k += 16) {
		// Check if we have less than 16 neighbors remaining.
		// If we do, use "i" as a dummy index.
		__declspec(align(64)) int indices[16];
		__mmask16 valid = _mm512_cmplt_epi32_mask(iter_count, limit);
		__m512i mj = _mm512_mask_extloadunpacklo_epi32(mi, valid, &neighs[k], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
		mj = _mm512_mask_extloadunpackhi_epi32(mj, valid, &neighs[k], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
		iter_count = _mm512_add_epi32(iter_count, all_16);
		_mm512_store_epi32(indices, mj);
			
		// Load 16 neighbors in AoS, and transpose their positions to SoA.
		__m512 xjtmp, yjtmp, zjtmp;
		AOS_TO_SOA_GATHER(xs, xjtmp, yjtmp, zjtmp, (&neighs[k]));
		
		// Compute squared distance.
        __m512 delx = _mm512_sub_ps(xitmp, xjtmp);
        __m512 dely = _mm512_sub_ps(yitmp, yjtmp);
        __m512 delz = _mm512_sub_ps(zitmp, zjtmp);
        __m512 delxsq = _mm512_mul_ps(delx, delx);
        __m512 delysq = _mm512_mul_ps(dely, dely);
        __m512 delzsq = _mm512_mul_ps(delz, delz);
        __m512 rsq = _mm512_add_ps(delxsq, _mm512_add_ps(delysq, delzsq));
		
		// Compute force.
        //__m512 sr2 = _mm512_div_ps(ones, rsq);     
        __m512 sr2 = _mm512_rcp23_ps(rsq);
        __m512 sr6 = _mm512_mul_ps(_mm512_mul_ps(sr2, _mm512_mul_ps(sr2, sr2)), msigma6);
		__mmask16 cutoff_mask = _mm512_kand(valid, _mm512_cmplt_ps_mask(rsq, mcutforcesq));
        __m512 F = _mm512_mask_mul_ps(zeroes, cutoff_mask, m48, _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(sr6, _mm512_sub_ps(sr6, halfs)), sr2), mepsilon));
		
		// Update force[i]
		__m512 Fx = _mm512_mul_ps(F, delx);
        __m512 Fy = _mm512_mul_ps(F, dely);
        __m512 Fz = _mm512_mul_ps(F, delz);
        fxitmp = _mm512_add_ps(fxitmp, Fx);
        fyitmp = _mm512_add_ps(fyitmp, Fy);
        fzitmp = _mm512_add_ps(fzitmp, Fz);        
		
		// Gather force[j], update force[j] and scatter force[j].
		// Requires an AoS->SoA transpose and an SoA->AoS transpose.   
		__m512 fxjtmp, fyjtmp, fzjtmp;
		AOS_TO_SOA_GATHER(fs, fxjtmp, fyjtmp, fzjtmp, (&neighs[k]));
		fxjtmp = _mm512_sub_ps(fxjtmp, Fx);
        fyjtmp = _mm512_sub_ps(fyjtmp, Fy);
        fzjtmp = _mm512_sub_ps(fzjtmp, Fz);    
		SOA_TO_AOS_SCATTER(fs, fxjtmp, fyjtmp, fzjtmp, (&neighs[k]));
		
		// Update energy/virial.
		if (EVFLAG) {
			m_eng_vdwl = _mm512_mask_add_ps(m_eng_vdwl, cutoff_mask, m_eng_vdwl, _mm512_mul_ps(_mm512_mul_ps(fours, _mm512_mul_ps(sr6, _mm512_sub_ps(sr6, ones))), mepsilon));
			m_virial = _mm512_mask_add_ps(m_virial, cutoff_mask, m_virial, _mm512_mul_ps(rsq, F));
		}
		
	}
#else
#ifdef USE_SIMD
    #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for (; k < numloc; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - xs[j+0];
      const MMD_float dely = ytmp - xs[j+1];
      const MMD_float delz = ztmp - xs[j+2];
      const MMD_float rsq = delx * delx + dely * dely + delz * delz;

	  const MMD_float sr2 = 1.0f / rsq;
	  const MMD_float sr6 = sr2 * sr2 * sr2 * sigma6;
	  const MMD_float force = (rsq < cutforcesq) ? 48.0f * sr6 * (sr6 - 0.5f) * sr2 * epsilon : 0.0f;

	  fix += delx * force;
	  fiy += dely * force;
	  fiz += delz * force;

	  fs[j+0] -= delx * force;
	  fs[j+1] -= dely * force;
	  fs[j+2] -= delz * force;

	  if (EVFLAG) {
  	    t_eng_vdwl += (4.0f * sr6 * (sr6 - 1.0f)) * epsilon;
	    t_virial += (delx * delx + dely * dely + delz * delz) * force;
	  }
    }
#endif
#ifdef IMCI
	k = numloc;
    loopbound = numloc + numrem;
	iter_count = _mm512_add_epi32(_mm512_load_epi32(&ttmp), _mm512_set1_epi32(k));
	limit = _mm512_extload_epi32(&loopbound, _MM_UPCONV_EPI32_NONE, _MM_BROADCAST_1X16, _MM_HINT_NONE);
	for(; k < loopbound; k += 16) {
	
		// Check if we have less than 16 neighbors remaining.
		// If we do, use "i" as a dummy index.
		__declspec(align(64)) int indices[16];
		__mmask16 valid = _mm512_cmplt_epi32_mask(iter_count, limit);
		__m512i mj = _mm512_mask_extloadunpacklo_epi32(mi, valid, &neighs[k], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
		mj = _mm512_mask_extloadunpackhi_epi32(mj, valid, &neighs[k], _MM_UPCONV_EPI32_NONE, _MM_HINT_NONE);
		iter_count = _mm512_add_epi32(iter_count, all_16);
		_mm512_store_epi32(indices, mj);
		
		// Load 16 neighbors in AoS, and transpose their positions to SoA.
		__m512 xjtmp, yjtmp, zjtmp;
		AOS_TO_SOA_GATHER(xs, xjtmp, yjtmp, zjtmp, (&neighs[k]));
		
		// Compute squared distance.
        __m512 delx = _mm512_sub_ps(xitmp, xjtmp);
        __m512 dely = _mm512_sub_ps(yitmp, yjtmp);
        __m512 delz = _mm512_sub_ps(zitmp, zjtmp);
        __m512 delxsq = _mm512_mul_ps(delx, delx);
        __m512 delysq = _mm512_mul_ps(dely, dely);
        __m512 delzsq = _mm512_mul_ps(delz, delz);
        __m512 rsq = _mm512_add_ps(delxsq, _mm512_add_ps(delysq, delzsq));
		
		// Compute force.
        //__m512 sr2 = _mm512_div_ps(ones, rsq);     
        __m512 sr2 = _mm512_rcp23_ps(rsq);
        __m512 sr6 = _mm512_mul_ps(_mm512_mul_ps(sr2, _mm512_mul_ps(sr2, sr2)), msigma6);
        __mmask16 cutoff_mask = _mm512_kand(valid, _mm512_cmplt_ps_mask(rsq, mcutforcesq));
        __m512 F = _mm512_mask_mul_ps(zeroes, cutoff_mask, m48, _mm512_mul_ps(_mm512_mul_ps(_mm512_mul_ps(sr6, _mm512_sub_ps(sr6, halfs)), sr2), mepsilon));
		
		// Update force[i]
		__m512 Fx = _mm512_mul_ps(F, delx);
        __m512 Fy = _mm512_mul_ps(F, dely);
        __m512 Fz = _mm512_mul_ps(F, delz);
        fxitmp = _mm512_add_ps(fxitmp, Fx);
        fyitmp = _mm512_add_ps(fyitmp, Fy);
        fzitmp = _mm512_add_ps(fzitmp, Fz);
		
		// Update energy/virial.
		if (EVFLAG) {
			m_eng_vdwl = _mm512_mask_add_ps(m_eng_vdwl, cutoff_mask, m_eng_vdwl, _mm512_mul_ps(halfs, _mm512_mul_ps(_mm512_mul_ps(fours, _mm512_mul_ps(sr6, _mm512_sub_ps(sr6, ones))),mepsilon)));
			m_virial = _mm512_mask_add_ps(m_virial, cutoff_mask, m_virial, _mm512_mul_ps(halfs, _mm512_mul_ps(rsq, F)));
		}
		
	}
#else	
#ifdef USE_SIMD
      #pragma simd reduction (+: fix,fiy,fiz,t_eng_vdwl,t_virial)
#endif
    for (; k < numloc + numrem; k++) {
      const int j = neighs[k];
      const MMD_float delx = xtmp - xs[j+0];
      const MMD_float dely = ytmp - xs[j+1];
      const MMD_float delz = ztmp - xs[j+2];
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
#endif
	
#ifdef IMCI
	// Include reduced force contribution from all SIMD lanes.
	fix += REDUCE(fxitmp);
	fiy += REDUCE(fyitmp);
	fiz += REDUCE(fzitmp);
#endif
    fs[i*4+0] += fix;
    fs[i*4+1] += fiy;
    fs[i*4+2] += fiz;

  }
  
  // Reduction and accumulate to eng/virial.
  if (EVFLAG)
  {
#ifdef IMCI
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


