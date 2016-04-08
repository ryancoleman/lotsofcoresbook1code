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

#include <cstdio>
#include <unistd.h>
#include <ctime>
#include <limits>
#include <getopt.h>
#include <omp.h>
#include <vector>

struct tile
{
    int   position[2];
    tile *neighbors[4]; // w, e, s, n

    int offset[2];
    int n[2];

    int ystride;
    int varstride;

    int     x_e_ystride;
    int     x_e_varstride;
    REAL_T *x_edges[2];
    int     y_e_ystride;
    int     y_e_varstride;
    REAL_T *y_edges[2];

    REAL_T *q;
    REAL_T *rho;
    REAL_T *rhou;
    REAL_T *rhov;
    REAL_T *E;
};

struct hydro_decomp
{
    int   decomp[2];
    int   ntiles;
    tile *tiles;
};

static void gather(hydro *h, const hydro_decomp *hd)
{
    for(int t = 0; t < hd->ntiles; ++t)
    {
        const tile *tl = hd->tiles + t;
        const int base = (tl->offset[1] + 2)*h->ystride + tl->offset[0] + 2;
        for(int j = 0; j < tl->n[1]; ++j)
        {
            for(int v = 0; v < 4; ++v)
                std::memcpy( h->q + base + j    * h->ystride + v*h->varstride,
                            tl->q +        (j+2)*tl->ystride + v*tl->varstride + 2,
                            sizeof(REAL_T)*tl->n[0]);
        }
    }
}

static void scatter(hydro_decomp *hd, const hydro *h)
{
    for(int t = 0; t < hd->ntiles; ++t)
    {
        const tile *tl = hd->tiles + t;
        const int base = (tl->offset[1] + 2)*h->ystride + tl->offset[0] + 2;
        for(int j = 0; j < tl->n[1]; ++j)
        {
            for(int v = 0; v < 4; ++v)
                std::memcpy(tl->q +        (j+2)*tl->ystride + v*tl->varstride + 2,
                             h->q + base + j    * h->ystride + v*h->varstride,
                            sizeof(REAL_T)*tl->n[0]);
        }
    }
}

static std::vector<int> get_factors(int num)
{
    std::vector<int> res;
    while(num > 1)
    {
        if((num % 2) == 0)
        {
            num /= 2;
            res.push_back(2);
        }
        else
        {
            bool divisor = false;
            for(int i = 3; i <= std::sqrt(num); i+=2)
                if((num % i) == 0)
                {
                    num /= i;
                    res.push_back(i);
                    divisor = true;
                    break;
                }
            if(!divisor)
            {
                res.push_back(num);
                num = 1;
            }
        }
    }
    return res;
}

static double nth_subset_sum(unsigned long long num, const std::vector<double> *set)
{
    double res = 0.0;
    for(int i = 0; i < set->size(); ++i)
        if(num & (1ULL << i))
            res += (*set)[i];
    return res;
}

static int nth_subset_prod(unsigned long long num, const std::vector<int> *set)
{
    int res = 1;
    for(int i = 0; i < set->size(); ++i)
        if(num & (1ULL << i))
            res *= (*set)[i];
    return res;
}

static unsigned long long best_subset(const double split, const std::vector<double> *ln_factors)
{
    double best = std::numeric_limits<double>::max();
    unsigned long long argbest = -1;
    for( unsigned long long i = 1; i < (1 << ln_factors->size())-1; ++i)
    {
        const double prod = std::abs(nth_subset_sum(i, ln_factors) - split);
        if(prod < best)
        {
            best = prod;
            argbest = i;
        }
    }
    return argbest;
}

static std::vector<int> best_decomp(int num, int groups)
{
    const double split = std::log((double)num)/groups;
    std::vector<int> res(groups);

    int todiv = num;
    for(int i = 0; i < groups-1; ++i)
    {
        std::vector<int>    factors = get_factors(todiv);
        std::vector<double> ln_factors(factors.size());
        for(int j = 0; j < ln_factors.size(); ++j)
            ln_factors[j] = std::log(factors[j]);

        const unsigned long long set_idx = best_subset(split, &ln_factors);
        res[i] = nth_subset_prod(set_idx, &factors);
        todiv /= res[i];
    }
    res[groups-1] = todiv;
    std::sort(res.begin(), res.end());
    return res;
}

static void init_tile(tile *tl, int xstart, int xend, int ystart, int yend)
{
    static const int target_alignment = 64;
    static const int arith_alignment  = target_alignment/sizeof(REAL_T);

    tl->offset[0] = xstart;
    tl->offset[1] = ystart;

    tl->n[0] = xend - xstart;
    tl->n[1] = yend - ystart;

    const int min_stride = tl->n[0] + 2*2;
    tl->ystride           = round_to_alignment(min_stride, arith_alignment);
    tl->varstride         = tl->ystride * (tl->n[1] + 2*2);

    const int alloc_offset = arith_alignment - 2;
    tl->q    = ((REAL_T *)_mm_malloc( sizeof(REAL_T) * (4 * tl->varstride + arith_alignment), target_alignment)) + alloc_offset;
    tl->rho  = tl->q + 0*tl->varstride;
    tl->rhou = tl->q + 1*tl->varstride;
    tl->rhov = tl->q + 2*tl->varstride;
    tl->E    = tl->q + 3*tl->varstride;

    const int min_x_e_varstride = tl->n[1] * 2;
    tl->x_e_ystride   = 2;
    tl->x_e_varstride = round_to_alignment(min_x_e_varstride, arith_alignment);
    tl->x_edges[0] = (REAL_T *) _mm_malloc( sizeof(REAL_T) * 4 * tl->x_e_varstride, target_alignment);
    tl->x_edges[1] = (REAL_T *) _mm_malloc( sizeof(REAL_T) * 4 * tl->x_e_varstride, target_alignment);

    const int min_y_e_varstride = tl->n[0] * 2;
    tl->y_e_ystride   = tl->n[0];
    tl->y_e_varstride = round_to_alignment(min_y_e_varstride, arith_alignment);
    tl->y_edges[0] = (REAL_T *) _mm_malloc( sizeof(REAL_T) * 4 * tl->y_e_varstride, target_alignment);
    tl->y_edges[1] = (REAL_T *) _mm_malloc( sizeof(REAL_T) * 4 * tl->y_e_varstride, target_alignment);
}

static void init_hydro_decomp(hydro_decomp *hd, const hydro *h, int ntiles, int quiet)
{
    const int x_work = (int) std::ceil(h->global_n[0]/(double)SIMD_WIDTH);
    const int y_work = h->global_n[1];
    const int total_cells = h->global_n[0]*h->global_n[1];

    std::vector<int> decomp(best_decomp(ntiles, 2));

    hd->decomp[0] = decomp[0];
    hd->decomp[1] = decomp[1];

    if(x_work/hd->decomp[0] < 1)
        die("Decomposition %d, %d doesn't allow for at least SIMD_WIDTH columns per thread!\n", hd->decomp[0], hd->decomp[1]);
    if(SIMD_WIDTH == 1 && x_work/hd->decomp[0] < 2)
        die("Decomposition %d, %d doesn't allow for at least 2 columns per thread!\n", hd->decomp[0], hd->decomp[1]);
    if(h->global_n[1]/hd->decomp[1] < 2)
        die("Decomposition %d, %d doesn't allow for at least 2 rows per thread!\n", hd->decomp[0], hd->decomp[1]);

    hd->ntiles    = hd->decomp[0] * hd->decomp[1];
    hd->tiles     = (tile*)malloc(sizeof(tile)*hd->ntiles);

    for(int j = 0; j < hd->decomp[1]; ++j)
    {
        u64 ystart, yend;
        divvy(&ystart, &yend, y_work, j, hd->decomp[1]);

        for(int i = 0; i < hd->decomp[0]; ++i)
        {
            int   tileno = hd->decomp[0] * j + i;
            tile *tl     = hd->tiles + tileno;

            u64 xstart, xend;
            divvy(&xstart, &xend, x_work, i, hd->decomp[0]);
            xstart *= SIMD_WIDTH;
            xend    = std::min(xend*SIMD_WIDTH, (u64)h->global_n[0]);

            init_tile(tl, xstart, xend, ystart, yend);

            tl->position[0] = i;
            tl->position[1] = j;

            tl->neighbors[0] = i > 0               ? hd->tiles + tileno - 1            : 0;
            tl->neighbors[1] = i < hd->decomp[0]-1 ? hd->tiles + tileno + 1            : 0;
            tl->neighbors[2] = j > 0               ? hd->tiles + tileno - hd->decomp[0] : 0;
            tl->neighbors[3] = j < hd->decomp[1]-1 ? hd->tiles + tileno + hd->decomp[0] : 0;

            const int tile_cells = (yend-ystart) * (xend-xstart);
            if(quiet < 1)
                printf("[config] Tile % 4d: (% 3d, % 3d) [% 5llu, % 5llu]x[% 5llu, % 5llu] (%4.2lf%% total/%4.2lf%% ideal) n: %d %d %d %d\n", tileno, i, j, xstart, xend, ystart, yend, 100.0*((double)tile_cells)/total_cells, 100.0*1.0/hd->ntiles, tl->neighbors[0] ? (int)(tl->neighbors[0] - hd->tiles) : -1  , tl->neighbors[1] ? (int)(tl->neighbors[1] - hd->tiles) : -1, tl->neighbors[2] ? (int)(tl->neighbors[2] - hd->tiles) : -1, tl->neighbors[3] ? (int)(tl->neighbors[3] - hd->tiles) : -1);
        }
    }

    scatter(hd, h);
}

static REAL_T tile_x_step(tile *restrict tl, const REAL_T dtdx, const hydro *restrict h, int jstart, int jend, bool do_courant)
{
    static const REAL_T x_signs[4] = {1.0, -1.0, 1.0, 1.0};

    VMASK_T alltrue;
    mask_true(&alltrue);
    VINT_T linear;
    linear_offset(&linear);

    VREAL_T courantv       = (VREAL_T) 0.0;
    REAL_T  final_courantv = 0.0;
    for(int j = jstart; j < jend; ++j)
    {
        const int   ob = (j + 2) * tl->ystride + 0;
        if(tl->neighbors[0])
        {
            for(int v = 0; v < 4; ++v)
            {
                tl->q[v*tl->varstride + ob + 0] = tl->x_edges[0][v*tl->x_e_varstride + j * tl->x_e_ystride + 0];
                tl->q[v*tl->varstride + ob + 1] = tl->x_edges[0][v*tl->x_e_varstride + j * tl->x_e_ystride + 1];
            }
        }
        else
            set_boundaries(tl->q + ob,                x_signs, 2,  1, 4, tl->varstride);

        if(tl->neighbors[1])
        {
            for(int v = 0; v < 4; ++v)
            {
                tl->q[v*tl->varstride + ob + tl->n[0] + 2 + 0] = tl->x_edges[1][v*tl->x_e_varstride + j * tl->x_e_ystride + 0];
                tl->q[v*tl->varstride + ob + tl->n[0] + 2 + 1] = tl->x_edges[1][v*tl->x_e_varstride + j * tl->x_e_ystride + 1];
            }
        }
        else
            set_boundaries(tl->q + ob + tl->n[0] + 3, x_signs, 2, -1, 4, tl->varstride);

        const int   o  = ob + 2;
        strip_work sw;
        strip_prime(&sw, tl->rho + o, tl->rhou + o, tl->rhov + o, tl->E + o, h,                             1, dtdx);

#ifndef HAVE_SIMD_TYPE
        int i = 0;
        for(; i + SIMD_WIDTH <= tl->n[0]; i+=SIMD_WIDTH)
        {
            const VREAL_T cv = strip_stable(h,  tl->rho + o, tl->rhou + o, tl->rhov + o, tl->E + o, &sw, i, 1, (VREAL_T) dtdx, do_courant);
            courantv         = std::max(courantv, cv);
        }
#else
        vstrip_work vsw;
        for(int v = 0; v < 4; ++ v)
        {
            vsw.flux      [v][0][SIMD_WIDTH-1] = sw.flux      [v][0];
            vsw.left_flux [v][0][SIMD_WIDTH-1] = sw.left_flux [v][0];
        }
        for(int v = 0; v < 5; ++ v)
        {
            vsw.prim      [v][0][SIMD_WIDTH-2] = sw.prim      [v][0];
            vsw.prim      [v][0][SIMD_WIDTH-1] = sw.prim      [v][1];
            vsw.prim      [v][1][SIMD_WIDTH-1] = sw.prim      [v][1];
        }

        int i = 0;
        for(; i + SIMD_WIDTH <= tl->n[0]; i+=SIMD_WIDTH)
        {
            const VREAL_T cv = hstrip_stable(h,  tl->rho + o, tl->rhou + o, tl->rhov + o, tl->E + o, &vsw, i, 1, (VREAL_T) dtdx, alltrue, do_courant);
            courantv         = std::max(courantv, cv);
        }

        for(; i < tl->n[0]; i+=SIMD_WIDTH)
        {
            const VMASK_T in_bounds = mask_lt(linear + (VINT_T) (double)i, (VINT_T) (double)tl->n[0]);
            const VREAL_T cv        = hstrip_stable(h,  tl->rho + o, tl->rhou + o, tl->rhov + o, tl->E + o, &vsw, i, 1, (VREAL_T) dtdx, in_bounds, do_courant);
            courantv                = std::max(courantv, cv);
        }
#endif // HAVE_SIMD_TYPE
    }

#ifdef HAVE_SIMD_TYPE
    for(int i = 0; i < SIMD_WIDTH; ++i)
        final_courantv = std::max(final_courantv, courantv[i]);
#else
    final_courantv = courantv;
#endif // HAVE_SIMD_TYPE
    return final_courantv;
}

static void tile_y_step(tile *restrict tl, const REAL_T dtdx, const hydro *restrict h, int istart, int iend)
{
    static const REAL_T y_signs[4] = {1.0, 1.0, -1.0, 1.0};

    VMASK_T alltrue;
    mask_true(&alltrue);
    VINT_T  linear;
    linear_offset(&linear);

    int i = istart;
    for(; i + SIMD_WIDTH <= iend; i+=SIMD_WIDTH)
    {
        const int    ob = i + 2;
        for(int k = 0; k < SIMD_WIDTH; ++k)
        {
            if(tl->neighbors[2])
            {
                for(int v = 0; v < 4; ++v)
                {
                    tl->q[v*tl->varstride +                 ob + k] = tl->y_edges[0][v*tl->y_e_varstride                     + k + i];
                    tl->q[v*tl->varstride + 1*tl->ystride + ob + k] = tl->y_edges[0][v*tl->y_e_varstride + 1*tl->y_e_ystride + k + i];
                }
            }
            else
                set_boundaries(tl->q + ob + k,                             y_signs, 2,  tl->ystride, 4, tl->varstride);

            if(tl->neighbors[3])
            {
                for(int v = 0; v < 4; ++v)
                {
                    tl->q[v*tl->varstride + (2 + tl->n[1])*tl->ystride + ob + k] = tl->y_edges[1][v*tl->y_e_varstride                     + k + i];
                    tl->q[v*tl->varstride + (3 + tl->n[1])*tl->ystride + ob + k] = tl->y_edges[1][v*tl->y_e_varstride + 1*tl->y_e_ystride + k + i];
                }
            }
            else
                set_boundaries(tl->q + ob + (tl->n[1] + 3)*tl->ystride +k, y_signs, 2, -tl->ystride, 4, tl->varstride);
        }
        const int    o  = ob + 2 * tl->ystride;
        vstrip_work   sw;
        vstrip_prime(&sw, tl->rho + o, tl->rhov + o, tl->rhou + o, tl->E + o, h,          tl->ystride, (VREAL_T) dtdx);
        for(int j = 0; j < tl->n[1]; ++j)
            vstrip_stable(h,  tl->rho + o, tl->rhov + o, tl->rhou + o, tl->E + o, &sw, j, tl->ystride, (VREAL_T) dtdx, alltrue, false);
    }
    for(; i < iend; i+=SIMD_WIDTH)
    {
        const VMASK_T in_bounds = mask_lt(linear + (VINT_T) (double) i, (VINT_T)  (double) iend);

        const int    ob = i + 2;
        for(int k = 0; k < SIMD_WIDTH && k + i < iend; ++k)
        {
            if(tl->neighbors[2])
            {
                for(int v = 0; v < 4; ++v)
                {
                    tl->q[v*tl->varstride +                 ob + k] = tl->y_edges[0][v*tl->y_e_varstride                     + k + i];
                    tl->q[v*tl->varstride + 1*tl->ystride + ob + k] = tl->y_edges[0][v*tl->y_e_varstride + 1*tl->y_e_ystride + k + i];
                }
            }
            else
                set_boundaries(tl->q + ob + k, y_signs, 2,  tl->ystride, 4, tl->varstride);

            if(tl->neighbors[3])
            {
                for(int v = 0; v < 4; ++v)
                {
                    tl->q[v*tl->varstride + (2 + tl->n[1])*tl->ystride + ob + k] = tl->y_edges[1][v*tl->y_e_varstride +                     k + i];
                    tl->q[v*tl->varstride + (3 + tl->n[1])*tl->ystride + ob + k] = tl->y_edges[1][v*tl->y_e_varstride + 1*tl->y_e_ystride + k + i];
                }
            }
            else
                set_boundaries(tl->q + ob + (tl->n[1] + 3)*tl->ystride +k, y_signs, 2, -tl->ystride, 4, tl->varstride);
        }
        const int    o  = ob + 2 * tl->ystride;
        vstrip_work   sw;
        vstrip_prime(&sw, tl->rho + o, tl->rhov + o, tl->rhou + o, tl->E + o, h,          tl->ystride, (VREAL_T) dtdx);
        for(int j = 0; j < tl->n[1]; ++j)
            vstrip_stable(h,  tl->rho + o, tl->rhov + o, tl->rhou + o, tl->E + o, &sw, j, tl->ystride, (VREAL_T) dtdx, in_bounds, false);
    }
}

static REAL_T x_step(hydro_decomp *hd, const hydro *restrict h, const REAL_T dtdx, int tid, bool do_courant)
{
    return tile_x_step(hd->tiles + tid, dtdx, h, 0, hd->tiles[tid].n[1], do_courant);
}

static void   y_step(hydro_decomp *hd, const hydro *restrict h, const REAL_T dtdx, int tid)
{
    tile_y_step(hd->tiles + tid, dtdx, h, 0, hd->tiles[tid].n[0]);
}

static void send_low_x_edge(tile *tl)
{
    if(tl->neighbors[0])
    {
        tile      *neighbor       = tl->neighbors[0];
        const int  dest_e_ystride = neighbor->x_e_ystride;
        for(int v = 0; v < 4; ++v)
        {
            REAL_T *dst_base = neighbor->x_edges[1] + v*neighbor->x_e_varstride;
            REAL_T *src_base = tl->q                + v*tl->varstride           + 2*tl->ystride + 2;
            for(int j = 0; j < tl->n[1]; ++j)
            {
                dst_base[j * dest_e_ystride + 0] = src_base[j * tl->ystride + 0];
                dst_base[j * dest_e_ystride + 1] = src_base[j * tl->ystride + 1];
            }
        }
    }
}

static void send_high_x_edge(tile *tl)
{
    if(tl->neighbors[1])
    {
        tile      *neighbor       = tl->neighbors[1];
        const int  dest_e_ystride = neighbor->x_e_ystride;
        for(int v = 0; v < 4; ++v)
        {
            REAL_T *dst_base = neighbor->x_edges[0] + v*neighbor->x_e_varstride;
            REAL_T *src_base = tl->q                + v*tl->varstride            + 2*tl->ystride + tl->n[0];
            for(int j = 0; j < tl->n[1]; ++j)
            {
                dst_base[j * dest_e_ystride + 0] = src_base[j * tl->ystride + 0];
                dst_base[j * dest_e_ystride + 1] = src_base[j * tl->ystride + 1];
            }
        }
    }
}

static void send_low_y_edge(tile *tl)
{
    if(tl->neighbors[2])
    {
        tile      *neighbor       = tl->neighbors[2];
        const int  dest_e_ystride = neighbor->y_e_ystride;
        for(int v = 0; v < 4; ++v)
        {
            REAL_T *dst_base = neighbor->y_edges[1] + v*neighbor->y_e_varstride;
            REAL_T *src_base = tl->q                + v*tl->varstride            + 2*tl->ystride + 2;
            memcpy(dst_base + 0 * dest_e_ystride, src_base + 0 * tl->ystride, sizeof(REAL_T)*tl->n[0]);
            memcpy(dst_base + 1 * dest_e_ystride, src_base + 1 * tl->ystride, sizeof(REAL_T)*tl->n[0]);
        }
    }
}

static void send_high_y_edge(tile *tl)
{
    if(tl->neighbors[3])
    {
        tile      *neighbor       = tl->neighbors[3];
        const int  dest_e_ystride = neighbor->y_e_ystride;
        for(int v = 0; v < 4; ++v)
        {
            REAL_T *dst_base = neighbor->y_edges[0] + v*neighbor->y_e_varstride;
            REAL_T *src_base = tl->q                + v*tl->varstride           + tl->n[1]*tl->ystride + 2;
            memcpy(dst_base + 0 * dest_e_ystride, src_base + 0 * tl->ystride, sizeof(REAL_T)*tl->n[0]);
            memcpy(dst_base + 1 * dest_e_ystride, src_base + 1 * tl->ystride, sizeof(REAL_T)*tl->n[0]);
        }
    }
}

static void print_config(FILE *fp)
{
    fprintf(fp, "GIT_VERSION         = %s\n", GIT_VERSION);
    fprintf(fp, "COMPILER_VERSION    = %s\n", COMPILER_VERSION);
    fprintf(fp, "BUILD_DATE          = %s %s\n", __DATE__, __TIME__);
    char hostbuff[256];
    gethostname(hostbuff, 256);
    fprintf(fp, "HOSTNAME            = %s\n", hostbuff);
    fprintf(fp, "NTHREADS            = %d\n",   omp_get_max_threads());
    fprintf(fp, "KMP_AFFINITY        = %s\n",   getenv("KMP_AFFINITY"));
    fprintf(fp, "VARIANT             = %s\n",   "PCL-HYDRO-TILE");
    fprintf(fp, "SIMD_ARCH           = %s\n",   SIMD_STR);
    fprintf(fp, "SIMD_LEN            = %d\n",   SIMD_WIDTH);
    fprintf(fp, "SIZEOF(REAL_T)      = %zu\n",  sizeof(REAL_T));
}

static void print_params(FILE *fp, const hydro *h)
{
    time_t thetime;
    time(&thetime);
    char timebuff[256];
    char *timestr = ctime_r(&thetime, timebuff);
    fprintf(fp, "RUN_START_DATE      = %s",     timestr);
    fprintf(fp, "NX                  = %d\n",   h->global_n[0]);
    fprintf(fp, "NY                  = %d\n",   h->global_n[1]);
    fprintf(fp, "TESTCASE            = %d\n",   h->testcase);
    fprintf(fp, "SCHEME              = %d\n",   h->scheme);
    fprintf(fp, "NSTEPMAX            = %u\n",   h->nstepmax);
    fprintf(fp, "IORDER              = %d\n",   h->iorder);
    fprintf(fp, "COURANT_NUMBER      = %lf\n",  h->courant_number);
    fprintf(fp, "DX                  = %lf\n",  h->dx);
    fprintf(fp, "TEND                = %lf\n",  h->tend);
}

static const char usage_str[] = "USAGE:\t%s <input file> [-h] [-v]\n";

static void usage(const char *name)
{
    die(usage_str, basename(name));
}

static void help(const char *name)
{
    fprintf(stderr, usage_str, name);
    fprintf(stderr, "DESCRIPTION\n"
            "\t Compute FVM solution to the Euler equations in 2d.\n"
            "\t Reimplementation of CEA's Hydro code\n");
    fprintf(stderr, "OPTIONS\n"
            "\t-h,--help\n\t    print this help message\n"
            "\t-v,--version\n\t    print configuration information\n"
            "\t--{no-}vtk-output\n\t    disable/enable vtk output (default: disabled)\n"
            "\t--output-interval <number>\n\t    set file output interval (default: infinity)\n"
            "\t--output <root>\n\t    set output root for timeseries\n"
            "\t-i,--input <file>\n\t    read parameters from <file>\n"
            "\t-q,--quiet\n\t    lower output to stdout\n"
            "\t--<key> <value>\n\t    if not captured by the above, try to add key/value to simulation parameters\n"
            );
}

int main(int argc, char *argv[])
{
    char   *input_file             = 0;
    char   *timeseries_output_root = 0;
    bool    vtk_output             = false;
    REAL_T  output_interval        = std::numeric_limits<REAL_T>::max();
    int     quiet                  = 0;

    const option opts[] =
    {
        {"help",            no_argument,       0, 'h'},
        {"version",         no_argument,       0, 'v'},
        {"vtk-output",      no_argument,       0, 'k'},
        {"no-vtk-output",   no_argument,       0, 'K'},
        {"output-interval", required_argument, 0, 'T'},
        {"output",          required_argument, 0, 'o'},
        {"input",           required_argument, 0, 'i'},
        {"quiet",           no_argument,       0, 'q'},
        {0, 0, 0, 0},
    };

    DECLARE_ARRAY_ALL(char *, candidate_kvs);
    INIT_ARRAY(candidate_kvs, 0);

    int opt;
    do
    {
        int in_ind = optind;
        opterr     = 0;
        opt        = getopt_long(argc, argv, "hvo:i:q", opts, 0);
        switch(opt)
        {
        case 0:
            break;
        case '?':
            if(optopt == 0)
            {
                EXTEND_ARRAY(candidate_kvs, 1);
                char *nondash = argv[in_ind];
                while(*nondash && *nondash == '-')
                    ++nondash;

                char *kvstr;
                int   klen = strlen(nondash);
                if(!strchr(nondash, '='))
                {
                    const int vlen = argv[in_ind+1] ? strlen(argv[in_ind+1]) : 0;
                    kvstr = (char*)malloc(klen + 1 + vlen + 1);
                    strcpy(kvstr, nondash);
                    kvstr[klen] = '=';
                    if(vlen)
                    {
                        strcpy(kvstr + klen + 1, argv[in_ind+1]);
                        ++optind;
                    }
                }
                else
                    kvstr = strdup(nondash);

                candidate_kvs[candidate_kvs_n++] = kvstr;
            }
            break;
        case 'k':
            vtk_output = true;
            break;
        case 'K':
            vtk_output = false;
            break;
        case 'T':
            output_interval = atof(optarg);
            break;
        case 'h':
            help(argv[0]);
            exit(0);
        case 'v':
            print_config(stderr);
            exit(0);
        case 'o':
            timeseries_output_root = strdup(optarg);
            break;
        case 'i':
            input_file = strdup(optarg);
            break;
        case 'q':
            ++quiet;
            break;
        default:
            usage(argv[0]);
        case -1:
            break;
        };
    }
    while(opt != -1);

    if(optind < argc)
       usage(argv[0]);

    fprintf(stderr, "==========================================\n");
    print_config(stderr);

    hydro H;
    load_hydro_params(&H, input_file, quiet);
    if(input_file)
        free(input_file);

    for(int i = 0; i < candidate_kvs_n; ++i)
    {
        if(quiet < 1)
            printf("[setup] Trying kv %s\n", candidate_kvs[i]);

        bool res = hydro_set_kv(&H, candidate_kvs[i]);
        if(!res && quiet < 2)
            printf("[setup] Hydro didn't accept option kv %s\n", candidate_kvs[i]);
        free(candidate_kvs[i]);
    }

    FREE_ARRAY_ALL(candidate_kvs);

    const int nthreads = omp_get_max_threads();

    init_hydro(&H);
    hydro_decomp HD;
    init_hydro_decomp(&HD, &H, nthreads, quiet);

    timeseries_writer  tw;
    if(timeseries_output_root)
    {
        tw.initialize(timeseries_output_root, 1024*1024*1024);

        tw.new_static("nx", sizeof(int));
        const int xw = H.global_n[0] + 2*2;
        tw.append(&xw, sizeof(int));
        tw.new_static("ny", sizeof(int));
        const int yw = H.global_n[1] + 2*2;
        tw.append(&yw, sizeof(int));

        gather(&H, &HD);
        write_hydro_ts(&tw, &H);
        if(quiet < 2)
            printf("[setup] output timeseries at %s (ts % 5d)\n", timeseries_output_root, -1);
    }

    REAL_T last_write  = 0.0;
    int    write_count = 0;

    REAL_T dt    = compute_timestep(&H) / (REAL_T) 2.0;
    REAL_T dt_dx = dt/H.dx;
    set_scheme(H.scheme, dt_dx);
    unsigned long long total_step_time = 0;

    REAL_T    *courantv     = (REAL_T*) _mm_malloc(64 * nthreads, 64);
    const int  cache_stride = 64/sizeof(REAL_T);

    print_params(stderr, &H);
    fprintf(stderr, "TIMESERIES_ROOT     = %s\n",  timeseries_output_root);
    fprintf(stderr, "------------------------------------------\n");

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();

        send_low_x_edge (HD.tiles + tid);
        send_high_x_edge(HD.tiles + tid);
        send_low_y_edge (HD.tiles + tid);
        send_high_y_edge(HD.tiles + tid);

        #pragma omp barrier
        // todo: handle > thread per tile
        while(H.step < H.nstepmax && H.t < H.tend)
        {
            unsigned long long step_start;
            #pragma omp barrier
            if(tid == 0)
            {
                H.t += dt;
                step_start = _rdtsc();
            }

            if(H.step % 2 == 0)
            {
                send_low_x_edge (HD.tiles + tid);
                send_high_x_edge(HD.tiles + tid);
                #pragma omp barrier
                x_step(&HD, &H, dt_dx, tid, false);
                #pragma omp barrier
                send_low_y_edge (HD.tiles + tid);
                send_high_y_edge(HD.tiles + tid);
                #pragma omp barrier
                y_step(&HD, &H, dt_dx, tid);
            }
            else
            {
                send_low_y_edge (HD.tiles + tid);
                send_high_y_edge(HD.tiles + tid);
                #pragma omp barrier
                y_step(&HD, &H, dt_dx, tid);
                #pragma omp barrier
                send_low_x_edge (HD.tiles + tid);
                send_high_x_edge(HD.tiles + tid);
                #pragma omp barrier
                courantv[tid*cache_stride] = x_step(&HD, &H, dt_dx, tid, true);
                #pragma omp barrier
                if(tid == 0)
                {
                    for(int i = 1; i < nthreads; ++i)
                        courantv[0] = std::max(courantv[i*cache_stride], courantv[0]);

                    dt            = H.courant_number * H.dx / courantv[0];
                    dt_dx         = dt/H.dx;
                    if(!set_scheme(H.scheme, dt_dx))
                        die("[error] Unknown limiter scheme! (%d)\n", H.scheme);
                }
            }

            #pragma omp barrier
            unsigned long long step_end;
            if(tid == 0)
            {
                step_end         = _rdtsc();
                total_step_time += step_end - step_start;

                ++H.step;
                switch(quiet)
                {
                case 1:
                    printf("\r");
                case 0:
                    printf("[run] step: %d, t: %le (dt: %le)", H.step, H.t, dt);
                    break;
                default:
                    break;
                }

                if(H.step % 2 == 0 && H.t - last_write > output_interval)
                {
                    int inc = 0;
                    if(vtk_output)
                    {
                        gather(&H, &HD);
                        vtkfile(write_count, H.q, H.global_n, 2, H.ystride, H.varstride, H.dx);
                        if(quiet < 2)
                            printf(" (vtk % 5d)", write_count);
                        inc        = 1;
                        last_write = H.t;
                    }
                    if(timeseries_output_root)
                    {
                        if(!inc)
                            gather(&H, &HD);
                        write_hydro_ts(&tw, &H);
                        if(quiet < 2)
                            printf(" (ts % 5d)", write_count);
                        inc        = 1;
                        last_write = H.t;
                    }
                    write_count += inc;
                }
                switch(quiet)
                {
                case 0:
                    printf("\n");
                    break;
                default:
                    break;
                }
            }

            #pragma omp barrier
        }

        #pragma omp barrier
        if(tid == 0)
        {
            if(quiet == 1)
                printf("\n");

            if(timeseries_output_root && H.t - last_write > 0.0)
            {
                gather(&H, &HD);
                write_hydro_ts(&tw, &H);
                if(quiet < 2)
                    printf("[setup] output timeseries at %s (ts % 5d)\n", timeseries_output_root, write_count);
            }
        }
    }

    fprintf(stderr, "%-15s = % llu\n", "TOTAL_CYCLES",     total_step_time);
    fprintf(stderr, "%-15s = % d\n",   "NSTEP",            H.step);
    fprintf(stderr, "%-15s = % le\n",  "CYCLES_PER_STEP", (double)total_step_time/(double)H.step);
    fprintf(stderr, "%-15s = % le\n",  "CYCLES_PER_CELL", (double)total_step_time/(double)(H.step*H.global_n[0]*H.global_n[1]));
    fprintf(stderr, "==========================================\n");

    if(timeseries_output_root)
        tw.finish();

    return 0;
}
