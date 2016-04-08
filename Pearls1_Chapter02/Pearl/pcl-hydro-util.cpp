/*
  A 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
  (C) Jason Sewall : Intel -- for the 'pcl-hydro' C++ version
*/
/*

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/
#include "pcl-hydro.hpp"
#include <cstring>

void init_hydro(hydro *h)
{
    if(h->nxystep == -1)
        h->nxystep = std::max(h->global_n[0], h->global_n[1]);
    h->ystride   = h->global_n[0] + 2*2;
    h->varstride = h->ystride * (h->global_n[1] + 2*2);
    h->q         = (REAL_T *) malloc(sizeof(REAL_T) * h->varstride * 4);

    h->inv_slope_type = 1.0/h->slope_type;

    for(int i = 0; i < h->varstride; ++i)
        h->q[i + 0*h->varstride] = (REAL_T) 1.0;
    for(int i = 0; i < h->varstride; ++i)
        h->q[i + 1*h->varstride] = (REAL_T) 0.0;
    for(int i = 0; i < h->varstride; ++i)
        h->q[i + 2*h->varstride] = (REAL_T) 0.0;
    for(int i = 0; i < h->varstride; ++i)
        h->q[i + 3*h->varstride] = (REAL_T) 1e-5;

    switch(h->testcase)
    {
    case 0:
        {
            const int x                             = h->global_n[0] / 2 + 2;
            const int y                             = h->global_n[1] / 2 + 2;
            h->q[h->ystride*y + x + 3*h->varstride] = ((REAL_T) 1.0) / h->dx / h->dx;
        }
        break;
    case 1:
        {
            const int x                             = 2;
            const int y                             = 2;
            h->q[h->ystride*y + x + 3*h->varstride] = ((REAL_T) 1.0) / h->dx / h->dx;
        }
        break;
    case 2:
        {
            const int x = 2;
            for(int j = 0; j < h->global_n[1]; ++j)
            {
                const int y                             = j + 2;
                h->q[h->ystride*y + x + 3*h->varstride] = ((REAL_T) 1.0) / h->dx / h->dx;
            }
        }
        break;
    default:
        die("Test case %d not implemented!\n", h->testcase);
    }
}

void write_hydro_ts(timeseries_writer *tw, const hydro *h)
{
    const int xw = h->global_n[0] + 2*2;
    const int yw = h->global_n[1] + 2*2;

    tw->new_frame(h->t, xw * yw * 4 * sizeof(REAL_T));
    for(int v = 0; v < 4; ++v)
        for(int j = 0; j < yw; ++j)
            tw->append(h->q + v*h->varstride + j * h->ystride, xw * sizeof(REAL_T));
}

static void set_boundary(      REAL_T *restrict dest,
                  const REAL_T           sign,
                  const int              width,
                  const int              stride)
{
    for(int i = 0; i < width; ++i)
        dest[i*stride] = sign*dest[(2*width - 1 - i)*stride];
}

void set_boundaries(      REAL_T *restrict dest_base,
                    const REAL_T *restrict signs,
                    const int              width,
                    const int              stride,
                    const int              nv,
                    const int              vstride)
{
    for(int v = 0; v < nv; ++v)
        set_boundary(dest_base + v*vstride, signs[v], width, stride);
}

REAL_T compute_timestep(const hydro *h)
{
    REAL_T courantv = SMALLC;
    for(int j = 2; j < h->global_n[1] + 2; ++j)
        for(int i = 2; i < h->global_n[0] + 2; ++i)
        {
            const int offs = j * h->ystride + i;
            REAL_T    prim_rho;
            REAL_T    prim_inv_rho;
            REAL_T    prim_u;
            REAL_T    prim_v;
            REAL_T    E_internal;

            conservative_to_primitive(&prim_rho,                  &prim_inv_rho, &prim_u,                     &prim_v,                     &E_internal,
                                      h->q[offs + 0*h->varstride],               h->q[offs + 1*h->varstride], h->q[offs + 2*h->varstride], h->q[offs + 3*h->varstride]);
            const REAL_T prim_p = equation_of_state(prim_rho, E_internal);
            const REAL_T prim_c = speed_of_sound   (prim_inv_rho,     prim_p);

            courant(&courantv, prim_u, prim_v, prim_c);
        }

    return h->courant_number * h->dx / courantv;
}

bool set_scheme(hydro::hscheme s, const REAL_T dt_dx)
{
    switch(s)
    {
    case hydro::MUSCL:
        ZEROL   = -((REAL_T) 100.0)/dt_dx;
        ZEROR   =  ((REAL_T) 100.0)/dt_dx;
        PROJECT =  (REAL_T) 1.0;
        break;
    case hydro::PLMDE:
        ZEROL   =  (REAL_T) 0;
        ZEROR   =  (REAL_T) 0;
        PROJECT =  (REAL_T) 1.0;
        break;
    case hydro::COLLELA:
        ZEROL   = (REAL_T) 0.0;
        ZEROR   = (REAL_T) 0.0;
        PROJECT = (REAL_T) 0.0;
        break;
    default:
        return false;
    }
    return true;
}
