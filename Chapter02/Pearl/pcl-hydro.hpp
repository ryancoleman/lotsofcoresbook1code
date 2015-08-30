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
#ifndef __PCL_HYDRO_HPP__
#define __PCL_HYDRO_HPP__

#include "arch.hpp"
#include "timeseries.hpp"

static const REAL_T GAMMA     = 1.4;
static const REAL_T GAMMA6    = (GAMMA + 1) / (2.0 * GAMMA);
static const REAL_T SMALLC    = 1e-10;
static const REAL_T SMALLR    = 1e-10;
static const REAL_T SMALLP    = SMALLC*SMALLC / GAMMA;
static const REAL_T SMALLPP   = SMALLR * SMALLP;
static const REAL_T PRECISION = 1e-6;

static const int NITER_RIEMANN = 10;
extern REAL_T ZEROR;
extern REAL_T ZEROL;
extern REAL_T PROJECT;

struct hydro
{
    typedef enum { MUSCL = 1, PLMDE = 2, COLLELA = 3} hscheme;
    int global_n[2];

    int nxystep;

    int ystride;
    int varstride;

    int     testcase;
    hscheme scheme;

    int          step;
    unsigned int nstepmax;
    int          iorder;
    REAL_T       slope_type;
    REAL_T       inv_slope_type;

    REAL_T courant_number;
    REAL_T dx;
    REAL_T t;
    REAL_T tend;

    REAL_T *q;
};

// util functions
void init_hydro(hydro *h);
bool load_hydro_params(hydro *h, const char *file, int quiet);
bool hydro_set_kv(hydro *H, char *kvstr);

void write_hydro_ts(timeseries_writer *tw, const hydro *h);

REAL_T compute_timestep(const hydro *h);
bool   set_scheme(hydro::hscheme s, const REAL_T dt_dx);
void   vtkfile(int step, const REAL_T *q, const int n[2], const int padding, const int ystride, const int varstride, const double dx);

void set_boundaries(      REAL_T *restrict dest_base,
                    const REAL_T *restrict signs,
                    const int              width,
                    const int              stride,
                    const int              nv,
                    const int              vstride);

// serial core functions
void conservative_to_primitive(      REAL_T *restrict prim_rho, REAL_T *restrict inv_prim_rho,       REAL_T *restrict prim_u,          REAL_T *restrict prim_v,          REAL_T *restrict E_internal,
                                     const REAL_T           cons_rho,                                const REAL_T           cons_rhou, const REAL_T           cons_rhov, const REAL_T           cons_E);
REAL_T equation_of_state(const REAL_T rho,
                         const REAL_T E_internal);
REAL_T speed_of_sound(const REAL_T inv_rho,
                      const REAL_T p);
REAL_T slope(const REAL_T nbv_m, const REAL_T nbv_0, const REAL_T nbv_p,
             const REAL_T slope_type, const REAL_T inv_slope_type);
void flux(      REAL_T *restrict flux_rho,                              REAL_T *restrict flux_u,         REAL_T *restrict flux_v,         REAL_T *restrict flux_p,
          const REAL_T           rho,      const REAL_T inv_rho,  const REAL_T           u,        const REAL_T           v,        const REAL_T           p,
          const REAL_T           sp_m,     const REAL_T sp_0,     const REAL_T           sp_p,
          const REAL_T           alpha_m,  const REAL_T alpha_0r, const REAL_T           alpha_0v, const REAL_T           alpha_p,
          const REAL_T           c);
void trace(      REAL_T *restrict flux_rho_m,                             REAL_T *restrict flux_u_m,       REAL_T *restrict flux_v_m,       REAL_T *restrict flux_p_m,
                        REAL_T *restrict flux_rho_p,                             REAL_T *restrict flux_u_p,       REAL_T *restrict flux_v_p,       REAL_T *restrict flux_p_p,
           const REAL_T           rho,        const REAL_T inv_rho, const REAL_T           u,        const REAL_T           v,        const REAL_T           p,
           const REAL_T           drho,       const REAL_T du,      const REAL_T           dv,       const REAL_T           dp,
           const REAL_T           c,          const REAL_T inv_c,
           const REAL_T           dtdx);
void riemann(      REAL_T *restrict     gdnv_rho,       REAL_T *restrict gdnv_u,           REAL_T *restrict gdnv_v,           REAL_T *restrict gdnv_p,
             const REAL_T            in_left_rho, const REAL_T           in_left_u,  const REAL_T           in_left_v,  const REAL_T           in_left_p,
             const REAL_T           in_right_rho, const REAL_T           in_right_u, const REAL_T           in_right_v, const REAL_T           in_right_p);
void cmpflx(      REAL_T *restrict flux_rho,       REAL_T *restrict flux_rhou,       REAL_T *restrict flux_rhov,       REAL_T *restrict flux_E,
            const REAL_T           gdnv_rho, const REAL_T           gdnv_u,    const REAL_T           gdnv_v,    const REAL_T           gdnv_p);
REAL_T update(const REAL_T  in,
              const REAL_T  flux_left, const REAL_T  flux_right,
              const REAL_T  dtdx);
void courant(      REAL_T *restrict courantv,
             const REAL_T           u, const REAL_T  v,
             const REAL_T           c);

struct strip_work
{
    REAL_T flux      [4][2]; // flux at i-1/2, i+1/2
    REAL_T left_flux [4][2]; // left_flux at i, i+1
    REAL_T prim      [5][3]; // prim for i, i+1, i+2
};

void strip_prime(strip_work   *restrict sw,
                 const REAL_T *restrict rho,
                 const REAL_T *restrict rhou,
                 const REAL_T *restrict rhov,
                 const REAL_T *restrict E,
                 const hydro  *restrict h,
                 const int              stride,
                 const REAL_T           dtdx);

REAL_T strip_stable(const hydro *restrict h,
                    REAL_T      *restrict rho,
                    REAL_T      *restrict rhou,
                    REAL_T      *restrict rhov,
                    REAL_T      *restrict E,
                    strip_work  *restrict sw,
                    const int             i,
                    const int             stride,
                    const REAL_T          dtdx,
                    const bool            do_courant);

// vector core functions
void vconservative_to_primitive(      VREAL_T *restrict prim_rho, VREAL_T *restrict inv_prim_rho,       VREAL_T *restrict prim_u,          VREAL_T *restrict prim_v,          VREAL_T *restrict E_internal,
                                const VREAL_T           cons_rho,                                 const VREAL_T           cons_rhou, const VREAL_T           cons_rhov, const VREAL_T           cons_E);
VREAL_T vequation_of_state(const VREAL_T rho,
                           const VREAL_T E_internal);
VREAL_T vspeed_of_sound(const VREAL_T inv_rho,
                        const VREAL_T p);
VREAL_T vslope(const VREAL_T nbv_m,      const VREAL_T nbv_0, const VREAL_T nbv_p,
               const VREAL_T slope_type, const VREAL_T inv_slope_type);
void vflux(      VREAL_T *restrict flux_rho,                              VREAL_T *restrict flux_u,         VREAL_T *restrict flux_v,         VREAL_T *restrict flux_p,
           const VREAL_T           rho,      const VREAL_T inv_rho,  const VREAL_T           u,        const VREAL_T           v,        const VREAL_T           p,
           const VREAL_T           sp_m,     const VREAL_T sp_0,     const VREAL_T           sp_p,
           const VREAL_T           alpha_m,  const VREAL_T alpha_0r, const VREAL_T           alpha_0v, const VREAL_T           alpha_p,
           const VREAL_T           c);
void vtrace(      VREAL_T *restrict flux_rho_m,                              VREAL_T *restrict flux_u_m,       VREAL_T *restrict flux_v_m,       VREAL_T *restrict flux_p_m,
                  VREAL_T *restrict flux_rho_p,                              VREAL_T *restrict flux_u_p,       VREAL_T *restrict flux_v_p,       VREAL_T *restrict flux_p_p,
            const VREAL_T           rho,        const VREAL_T inv_rho, const VREAL_T           u,        const VREAL_T           v,        const VREAL_T           p,
            const VREAL_T           drho,       const VREAL_T du,      const VREAL_T           dv,       const VREAL_T           dp,
            const VREAL_T           c,          const VREAL_T inv_c,
            const VREAL_T           dtdx);
void vriemann(      VREAL_T *restrict     gdnv_rho,       VREAL_T *restrict gdnv_u,           VREAL_T *restrict gdnv_v,           VREAL_T *restrict gdnv_p,
              const VREAL_T            in_left_rho, const VREAL_T           in_left_u,  const VREAL_T           in_left_v,  const VREAL_T           in_left_p,
              const VREAL_T           in_right_rho, const VREAL_T           in_right_u, const VREAL_T           in_right_v, const VREAL_T           in_right_p);
void vcmpflx(      VREAL_T *restrict flux_rho,       VREAL_T *restrict flux_rhou,       VREAL_T *restrict flux_rhov,       VREAL_T *restrict flux_E,
             const VREAL_T           gdnv_rho, const VREAL_T           gdnv_u,    const VREAL_T           gdnv_v,    const VREAL_T           gdnv_p);
VREAL_T vupdate(const VREAL_T  in,
                const VREAL_T  flux_left, const VREAL_T  flux_right,
                const VREAL_T  dtdx);
void vcourant(      VREAL_T *restrict courantv,
              const VREAL_T           u, const VREAL_T  v,
              const VREAL_T           c,
              const VMASK_T           write_mask);

struct vstrip_work
{
    VREAL_T flux      [4][2]; // flux at i-1/2, i+1/2
    VREAL_T left_flux [4][2]; // left_flux at i, i+1
    VREAL_T prim      [5][3]; // prim for i, i+1, i+2
};

void vstrip_prime(      vstrip_work  *restrict sw,
                  const REAL_T       *restrict rho,
                  const REAL_T       *restrict rhou,
                  const REAL_T       *restrict rhov,
                  const REAL_T       *restrict E,
                  const hydro        *restrict h,
                  const int                    stride,
                  const VREAL_T                dtdx);

VREAL_T vstrip_stable(const hydro       *restrict h,
                            REAL_T      *restrict rho,
                            REAL_T      *restrict rhou,
                            REAL_T      *restrict rhov,
                            REAL_T      *restrict E,
                            vstrip_work *restrict sw,
                      const int                   i,
                      const int                   stride,
                      const VREAL_T               dtdx,
                      const VMASK_T               write_mask,
                      const bool                  do_courant);

VREAL_T hstrip_stable(const hydro       *restrict h,
                            REAL_T      *restrict rho,
                            REAL_T      *restrict rhou,
                            REAL_T      *restrict rhov,
                            REAL_T      *restrict E,
                            vstrip_work *restrict sw,
                      const int                   i,
                      const int                   stride,
                      const VREAL_T               dtdx,
                      const VMASK_T               write_mask,
                      const bool                  do_courant);

#endif /* __PCL_HYDRO_HPP__ */
