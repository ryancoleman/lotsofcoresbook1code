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

void vconservative_to_primitive(      VREAL_T *restrict prim_rho, VREAL_T *restrict inv_prim_rho,       VREAL_T *restrict prim_u,          VREAL_T *restrict prim_v,          VREAL_T *restrict E_internal,
                                const VREAL_T           cons_rho,                                 const VREAL_T           cons_rhou, const VREAL_T           cons_rhov, const VREAL_T           cons_E)
{
    *prim_rho     = std::max(cons_rho, VREAL_T(SMALLR));
    *inv_prim_rho = rcp(*prim_rho);

    *prim_u = cons_rhou * *inv_prim_rho;
    *prim_v = cons_rhov * *inv_prim_rho;

    const VREAL_T E_kinetic = VREAL_T((REAL_T)0.5) * (*prim_u * *prim_u + *prim_v * *prim_v);
    *E_internal             = cons_E * *inv_prim_rho - E_kinetic;

    // todo: handle passive terms here
};

VREAL_T vequation_of_state(const VREAL_T rho,
                           const VREAL_T E_internal)
{
    return std::max((VREAL_T(GAMMA - (REAL_T) 1.0)) * rho * E_internal, VREAL_T(SMALLP));
}

VREAL_T vspeed_of_sound(const VREAL_T inv_rho,
                        const VREAL_T p)
{
    return my_sqrt(VREAL_T(GAMMA) * p * inv_rho);
}

VREAL_T vslope(const VREAL_T nbv_m,      const VREAL_T nbv_0, const VREAL_T nbv_p,
               const VREAL_T slope_type, const VREAL_T inv_slope_type)

{
    const VREAL_T left    = slope_type * (nbv_0 - nbv_m);
    const VREAL_T right   = slope_type * (nbv_p - nbv_0);
    const VREAL_T center  = ((VREAL_T) 0.5) * (left + right) * inv_slope_type;
    const VREAL_T sign    = select_gt(center,         (VREAL_T) 0.0, (VREAL_T) 1.0, (VREAL_T) -1.0);
    const VREAL_T llftrgt = select_le((left * right), (VREAL_T) 0.0, (VREAL_T) 1.0, (VREAL_T) 0.0);
    const VREAL_T t1      = std::min(std::abs(left), std::abs(right));
    return sign * std::min((((VREAL_T) 1.0) - llftrgt) * t1, std::abs(center));
}

void vflux(      VREAL_T *restrict flux_rho,                              VREAL_T *restrict flux_u,         VREAL_T *restrict flux_v,         VREAL_T *restrict flux_p,
           const VREAL_T           rho,      const VREAL_T inv_rho,  const VREAL_T           u,        const VREAL_T           v,        const VREAL_T           p,
           const VREAL_T           sp_m,     const VREAL_T sp_0,     const VREAL_T           sp_p,
           const VREAL_T           alpha_m,  const VREAL_T alpha_0r, const VREAL_T           alpha_0v, const VREAL_T           alpha_p,
           const VREAL_T           c)
{
    const VREAL_T a_p  = ((VREAL_T)-0.5) * sp_p * alpha_p;
    const VREAL_T a_m  = ((VREAL_T)-0.5) * sp_m * alpha_m;
    const VREAL_T a_0r = ((VREAL_T)-0.5) * sp_0 * alpha_0r;
    const VREAL_T a_0v = ((VREAL_T)-0.5) * sp_0 * alpha_0v;

    *flux_rho = rho + (a_p + a_m + a_0r);
    *flux_u   = u   + (a_p - a_m) * c * inv_rho;
    *flux_v   = v   + a_0v;
    *flux_p   = p   + (a_p + a_m) * c * c;
}

void vtrace(      VREAL_T *restrict flux_rho_m,                              VREAL_T *restrict flux_u_m,       VREAL_T *restrict flux_v_m,       VREAL_T *restrict flux_p_m,
                  VREAL_T *restrict flux_rho_p,                              VREAL_T *restrict flux_u_p,       VREAL_T *restrict flux_v_p,       VREAL_T *restrict flux_p_p,
            const VREAL_T           rho,        const VREAL_T inv_rho, const VREAL_T           u,        const VREAL_T           v,        const VREAL_T           p,
            const VREAL_T           drho,       const VREAL_T du,      const VREAL_T           dv,       const VREAL_T           dp,
            const VREAL_T           c,          const VREAL_T inv_c,
            const VREAL_T           dtdx)
{
    const VREAL_T alpha_m  = ((VREAL_T) 0.5) * (dp * ( inv_rho * inv_c ) - du) * rho * inv_c;
    const VREAL_T alpha_p  = ((VREAL_T) 0.5) * (dp * ( inv_rho * inv_c ) + du) * rho * inv_c;
    const VREAL_T alpha_0r = drho - dp * (inv_c*inv_c);
    const VREAL_T alpha_0v = dv;

    const VREAL_T right_sp_m = select_gt(u - c, VREAL_T(ZEROR), VREAL_T(PROJECT), (u - c) * dtdx + (VREAL_T) 1.0);
    const VREAL_T right_sp_p = select_gt(u + c, VREAL_T(ZEROR), VREAL_T(PROJECT), (u + c) * dtdx + (VREAL_T) 1.0);
    const VREAL_T right_sp_0 = select_gt(u,     VREAL_T(ZEROR), VREAL_T(PROJECT),  u      * dtdx + (VREAL_T) 1.0);

    vflux(flux_rho_p,          flux_u_p,   flux_v_p,   flux_p_p,
          rho,        inv_rho, u,          v,          p,
          right_sp_m,          right_sp_0, right_sp_p,
          alpha_m,             alpha_0r,   alpha_0v,   alpha_p,
          c);

    const VREAL_T left_sp_m = select_le(u - c, VREAL_T(ZEROL), VREAL_T(-PROJECT), (u - c) * dtdx - (VREAL_T) 1.0);
    const VREAL_T left_sp_p = select_le(u + c, VREAL_T(ZEROL), VREAL_T(-PROJECT), (u + c) * dtdx - (VREAL_T) 1.0);
    const VREAL_T left_sp_0 = select_le(u,     VREAL_T(ZEROL), VREAL_T(-PROJECT),  u      * dtdx - (VREAL_T) 1.0);

    vflux(flux_rho_m,          flux_u_m,   flux_v_m,   flux_p_m,
          rho,        inv_rho, u,          v,          p,
          left_sp_m,           left_sp_0,  left_sp_p,
          alpha_m,             alpha_0r,   alpha_0v,   alpha_p,
          c);

    // todo: handle passive terms
}

void vriemann(      VREAL_T *restrict     gdnv_rho,       VREAL_T *restrict gdnv_u,           VREAL_T *restrict gdnv_v,           VREAL_T *restrict gdnv_p,
              const VREAL_T            in_left_rho, const VREAL_T           in_left_u,  const VREAL_T           in_left_v,  const VREAL_T           in_left_p,
              const VREAL_T           in_right_rho, const VREAL_T           in_right_u, const VREAL_T           in_right_v, const VREAL_T           in_right_p)

{
    const VREAL_T left_rho = std::max(in_left_rho, VREAL_T(SMALLR));
    const VREAL_T left_u   = in_left_u;
    const VREAL_T left_v   = in_left_v;
    const VREAL_T left_p   = std::max(in_left_p, left_rho * VREAL_T(SMALLP));
    const VREAL_T left_c   = VREAL_T(GAMMA) * left_p * left_rho;

    const VREAL_T right_rho = std::max(in_right_rho, VREAL_T(SMALLR));
    const VREAL_T right_u   = in_right_u;
    const VREAL_T right_v   = in_right_v;
    const VREAL_T right_p   = std::max(in_right_p, right_rho * VREAL_T(SMALLP));
    const VREAL_T right_c   = VREAL_T(GAMMA) * right_p * right_rho;

    VREAL_T p_star;
    VMASK_T goon;
    mask_true(&goon);
    {
        const VREAL_T left_w  = my_sqrt(left_c);
        const VREAL_T right_w = my_sqrt(right_c);
        p_star                = std::max( (right_w * left_p + left_w * right_p + left_w * right_w * (left_u - right_u)) * rcp(left_w + right_w), (VREAL_T) 0.0);
}

    for(int i = 0; i < NITER_RIEMANN && !all_zero(goon); ++i)
    {
        const VREAL_T left_ww2  = left_rho  * ((VREAL_T) 0.5) * (VREAL_T(GAMMA + (REAL_T) 1.0) * p_star + VREAL_T(GAMMA - (REAL_T) 1.0) * left_p);
        const VREAL_T left_ww   = my_sqrt(left_ww2);
        const VREAL_T right_ww2 = right_rho * ((VREAL_T) 0.5) * (VREAL_T(GAMMA + (REAL_T) 1.0) * p_star + VREAL_T(GAMMA - (REAL_T) 1.0) * right_p);
        const VREAL_T right_ww  = my_sqrt(right_ww2);
        const VREAL_T tmp_num   = ((VREAL_T)2.0) * left_ww2 * right_ww2 * (left_ww * right_ww * (left_u - right_u) - left_ww * (p_star - right_p) - right_ww * (p_star - left_p));
        const VREAL_T tmp_den   = right_ww2*right_ww * (left_ww2 + left_c) + left_ww2*left_ww * (right_ww2 + right_c);
        const VREAL_T tmp       = tmp_num * rcp(tmp_den);
        const VREAL_T deleft_p  = std::max(tmp, -p_star);

        p_star += select_true(goon, deleft_p, (VREAL_T) 0.0);

        const VREAL_T uo = std::abs(deleft_p * rcp(p_star + SMALLPP));
        goon             = mask_and(goon, mask_gt(uo, VREAL_T(PRECISION)));
    }

    const VREAL_T left_w2  = left_rho  * ((VREAL_T) 0.5) * (VREAL_T(GAMMA + (REAL_T) 1.0) * p_star + VREAL_T(GAMMA - (REAL_T) 1.0) * left_p);
    const VREAL_T left_w   = my_sqrt(left_w2);
    const VREAL_T right_w2 = right_rho * ((VREAL_T) 0.5) * (VREAL_T(GAMMA + (REAL_T) 1.0) * p_star + VREAL_T(GAMMA - (REAL_T) 1.0) * right_p);
    const VREAL_T right_w  = my_sqrt(right_w2);

    const VREAL_T u_star = (VREAL_T) 0.5 * (left_u + (left_p - p_star) * rcp(left_w) + right_u - (right_p - p_star) * rcp(right_w));

    const VREAL_T sgnm  = select_gt(u_star, (VREAL_T) 0.0, (VREAL_T) 1.0, (VREAL_T) -1.0);
    const VREAL_T rho_0 = select_gt(u_star, (VREAL_T) 0.0, left_rho,      right_rho);
    const VREAL_T u_0   = select_gt(u_star, (VREAL_T) 0.0, left_u,        right_u);
    const VREAL_T p_0   = select_gt(u_star, (VREAL_T) 0.0, left_p,        right_p);
    const VREAL_T w_0   = select_gt(u_star, (VREAL_T) 0.0, left_w,        right_w);
    const VREAL_T w2_0  = select_gt(u_star, (VREAL_T) 0.0, left_w2,       right_w2);

    const VREAL_T inv_rho_0 = rcp(rho_0);
    const VREAL_T c_0       = std::max(my_sqrt(std::abs(VREAL_T(GAMMA) * p_0 * inv_rho_0)),        VREAL_T(SMALLC));
    const VREAL_T rho_star  = std::max(w2_0*rho_0 * rcp(w2_0 + rho_0 * (p_0 - p_star)),            VREAL_T(SMALLR));
    const VREAL_T c_star    = std::max(my_sqrt(std::abs(VREAL_T(GAMMA) * p_star * rcp(rho_star))), VREAL_T(SMALLC));

    const VREAL_T ushock = w_0 * inv_rho_0 - sgnm * u_0;

    const VREAL_T spout  = select_ge(p_star, p_0, ushock, c_0    - sgnm * u_0);
    const VREAL_T spin   = select_ge(p_star, p_0, ushock, c_star - sgnm * u_star);

    const VREAL_T scr  = std::max(spout - spin, VREAL_T(SMALLC) + std::abs(spout + spin));
    VREAL_T       frac = std::max(std::min((((VREAL_T) 1.0) + (spout + spin) * rcp(scr)) * (VREAL_T) 0.5, (VREAL_T) 1.0), (VREAL_T) 0.0);

    const VMASK_T set_spout = mask_lt(spout, (VREAL_T) 0.0);
    const VMASK_T set_spin  = mask_gt(spin,  (VREAL_T) 0.0);
    frac                    = select_true(set_spout, (VREAL_T)0.0, frac);
    frac                    = select_true(mask_and(mask_not(set_spout), set_spin), (VREAL_T)1.0, frac);

    *gdnv_rho = (frac * rho_star + (((VREAL_T) 1.0) - frac) * rho_0);
    *gdnv_u   = (frac * u_star   + (((VREAL_T) 1.0) - frac) * u_0);
    *gdnv_v   = select_gt(u_star, (VREAL_T) 0.0, left_v, right_v);
    *gdnv_p   = (frac * p_star   + (((VREAL_T) 1.0) - frac) * p_0);

    // todo: handle passive
}

void vcmpflx(      VREAL_T *restrict flux_rho,       VREAL_T *restrict flux_rhou,       VREAL_T *restrict flux_rhov,       VREAL_T *restrict flux_E,
             const VREAL_T           gdnv_rho, const VREAL_T           gdnv_u,    const VREAL_T           gdnv_v,    const VREAL_T           gdnv_p)
{
    const VREAL_T mass_density = gdnv_rho * gdnv_u;
    *flux_rho                  = mass_density;
    *flux_rhou                 = mass_density * gdnv_u + gdnv_p;
    *flux_rhov                 = mass_density * gdnv_v;

    const VREAL_T E_kinetic = ((VREAL_T) 0.5) * gdnv_rho * (gdnv_u * gdnv_u + gdnv_v * gdnv_v);
    const VREAL_T E_total   = gdnv_p * VREAL_T((REAL_T) 1.0 / (REAL_T(GAMMA) - (REAL_T) 1.0)) + E_kinetic;
    *flux_E                 = gdnv_u * (E_total + gdnv_p);

    // todo: handle passive
}

VREAL_T vupdate(const VREAL_T  in,
                const VREAL_T  flux_left, const VREAL_T  flux_right,
                const VREAL_T  dtdx)
{
    return in  + (flux_left  - flux_right)  * dtdx;

    // todo: handle passive
}

void vcourant(      VREAL_T *restrict courantv,
              const VREAL_T           u, const VREAL_T  v,
              const VREAL_T           c,
              const VMASK_T           write_mask)
{
    const VREAL_T max_speed = std::max(c + std::abs(u), c + std::abs(v));
    const VMASK_T gt        = mask_gt(max_speed, *courantv);
    *courantv               = select_true(mask_and(gt, write_mask), max_speed, *courantv);

}

void vstrip_prime(      vstrip_work *restrict sw,
                  const REAL_T      *restrict rho,
                  const REAL_T      *restrict rhou,
                  const REAL_T      *restrict rhov,
                  const REAL_T      *restrict E,
                  const hydro       *restrict h,
                  const int                   stride,
                  const VREAL_T               dtdx)
{
    VREAL_T pre_prim[5][4]; // i-2, i-1, i, i+1
    for(int i = 0; i < 4; ++i)
    {
        const int src_offs = (- 2 + i)*stride;
        VREAL_T E_internal;
        vconservative_to_primitive(pre_prim[0] +i,       pre_prim[4] +i, pre_prim[1] +i,        pre_prim[2] +i, &E_internal,
                                   load(rho + src_offs),                 load(rhou + src_offs), load(rhov + src_offs), load(E + src_offs));
        pre_prim[3][i] = vequation_of_state(pre_prim[0][i], E_internal);
    }

    VREAL_T pre_dvar[4][2]; // i-1, i
    if(h->iorder != 1)
        for(int i = 0; i < 2; ++i)
            for(int v = 0; v < 4; ++ v)
                pre_dvar[v][i] = vslope(pre_prim[v][0+i], pre_prim[v][1+i], pre_prim[v][2+i], (VREAL_T)h->slope_type, (VREAL_T)h->inv_slope_type);

    VREAL_T pre_left_flux [4][2]; // i-1, i
    VREAL_T pre_right_flux[4][2]; // i-1, i
    for(int i = 0; i < 2; ++i)
    {
        const VREAL_T prim_c = vspeed_of_sound(pre_prim[4][i + 1], pre_prim[3][i+1]);
        vtrace(pre_left_flux [0] + i,                   pre_left_flux [1] + i,  pre_left_flux [2] + i, pre_left_flux [3] + i,
               pre_right_flux[0] + i,                   pre_right_flux[1] + i,  pre_right_flux[2] + i, pre_right_flux[3] + i,
               pre_prim[0]     [i+1], pre_prim[4][i+1], pre_prim[1]      [i+1], pre_prim[2]     [i+1], pre_prim[3]     [i+1],
               pre_dvar[0]       [i],                   pre_dvar[1]        [i], pre_dvar[2]       [i], pre_dvar[3]       [i],
               prim_c, rcp(prim_c),
               dtdx);
    }

    VREAL_T gdnv_rho, gdnv_u, gdnv_v, gdnv_p;
    vriemann(&gdnv_rho,        &gdnv_u,          &gdnv_v,          &gdnv_p,
             pre_left_flux [0][0], pre_left_flux [1][0], pre_left_flux [2][0], pre_left_flux [3][0],
             pre_right_flux[0][1], pre_right_flux[1][1], pre_right_flux[2][1], pre_right_flux[3][1]);

    vcmpflx(sw->flux[0] + 0, sw->flux[1] + 0, sw->flux[2] + 0, sw->flux[3] + 0,
           gdnv_rho,        gdnv_u,          gdnv_v,          gdnv_p);

    for(int v = 0; v < 4; ++ v)
        sw->left_flux[v][0] = pre_left_flux[v][1];

    for(int v = 0; v < 5; ++ v)
    {
        sw->prim[v][0] = pre_prim[v][2];
        sw->prim[v][1] = pre_prim[v][3];
    }
}

VREAL_T vstrip_stable(const hydro  *restrict h,
                      REAL_T       *restrict rho,
                      REAL_T       *restrict rhou,
                      REAL_T       *restrict rhov,
                      REAL_T       *restrict E,
                      vstrip_work  *restrict sw,
                      const int              i,
                      const int              stride,
                      const VREAL_T          dtdx,
                      const VMASK_T          write_mask,
                      const bool             do_courant)
{
    const int src_offs = (i + 2)*stride;

    VREAL_T E_internal;
    vconservative_to_primitive(sw->prim[0] + 2,      sw->prim[4] + 2, sw->prim[1] + 2,      sw->prim[2] + 2, &E_internal,
                               load(rho + src_offs),                  load(rhou+src_offs),  load(rhov + src_offs), load(E + src_offs));
    sw->prim[3][2] = vequation_of_state(sw->prim[0][2], E_internal);

    VREAL_T dvar[4];    // slope for i+1
    if(h->iorder != 1)
        for(int v = 0; v < 4; ++ v)
            dvar[v] = vslope(sw->prim[v][0], sw->prim[v][1], sw->prim[v][2], (VREAL_T) h->slope_type, (VREAL_T) h->inv_slope_type);

    VREAL_T right_flux[4];
    const VREAL_T prim_c = vspeed_of_sound(sw->prim[4][1], sw->prim[3][1]);
    vtrace(sw->left_flux [0] + 1,                 sw->left_flux [1] + 1, sw->left_flux [2] + 1, sw->left_flux [3] + 1,
           right_flux + 0,                        right_flux + 1,        right_flux + 2,        right_flux + 3,
           sw->prim[0]       [1], sw->prim[4][1], sw->prim[1]       [1], sw->prim[2]       [1], sw->prim[3]       [1],
           dvar[0],                                dvar[1],              dvar[2],               dvar[3],
           prim_c, rcp(prim_c),
           dtdx);

    VREAL_T gdnv_rho, gdnv_u, gdnv_v, gdnv_p;
    vriemann(&gdnv_rho,            &gdnv_u,              &gdnv_v,              &gdnv_p,
             sw->left_flux [0][0], sw->left_flux [1][0], sw->left_flux [2][0], sw->left_flux [3][0],
             right_flux[0],        right_flux[1],        right_flux[2],        right_flux[3]);

    vcmpflx(sw->flux[0] + 1, sw->flux[1] + 1, sw->flux[2] + 1, sw->flux[3] + 1,
            gdnv_rho,        gdnv_u,          gdnv_v,          gdnv_p);

    const VREAL_T new_rho  = vupdate(load(rho  + i*stride), sw->flux[0][0], sw->flux[0][1], VREAL_T(dtdx));
    const VREAL_T new_rhou = vupdate(load(rhou + i*stride), sw->flux[1][0], sw->flux[1][1], VREAL_T(dtdx));
    const VREAL_T new_rhov = vupdate(load(rhov + i*stride), sw->flux[2][0], sw->flux[2][1], VREAL_T(dtdx));
    const VREAL_T new_E    = vupdate(load(E    + i*stride), sw->flux[3][0], sw->flux[3][1], VREAL_T(dtdx));

    VREAL_T courantv = (VREAL_T) 0.0;
    if(do_courant)
    {
        VREAL_T prim_rho, prim_inv_rho, prim_u, prim_v, E_internal;
        vconservative_to_primitive(&prim_rho, &prim_inv_rho, &prim_u,  &prim_v,  &E_internal,
                                   new_rho,                  new_rhou, new_rhov, new_E);
        const VREAL_T prim_p = vequation_of_state(prim_rho,     E_internal);
        const VREAL_T prim_c = vspeed_of_sound   (prim_inv_rho, prim_p);
        vcourant(&courantv, prim_u, prim_v, prim_c, write_mask);
    }

    maskstore(rho  + i*stride, new_rho,  write_mask);
    maskstore(rhou + i*stride, new_rhou, write_mask);
    maskstore(rhov + i*stride, new_rhov, write_mask);
    maskstore(E    + i*stride, new_E,    write_mask);

    for(int v = 0; v < 4; ++ v)
    {
        sw->flux     [v][0] = sw->flux     [v][1];
        sw->left_flux[v][0] = sw->left_flux[v][1];
    }
    for(int v = 0; v < 5; ++ v)
    {
        sw->prim[v][0] = sw->prim[v][1];
        sw->prim[v][1] = sw->prim[v][2];
    }

    return courantv;
}

#ifdef HAVE_SIMD_TYPE
VREAL_T hstrip_stable(const hydro  *restrict h,
                      REAL_T       *restrict rho,
                      REAL_T       *restrict rhou,
                      REAL_T       *restrict rhov,
                      REAL_T       *restrict E,
                      vstrip_work  *restrict sw,
                      const int              i,
                      const int              stride,
                      const VREAL_T          dtdx,
                      const VMASK_T          write_mask,
                      const bool             do_courant)
{
    const int src_offs = (i + 2)*stride;

    VREAL_T E_internal;
    vconservative_to_primitive(sw->prim[0] + 2,      sw->prim[4] + 2, sw->prim[1] + 2,        sw->prim[2] + 2, &E_internal,
                               loadu(rho + src_offs),                 loadu(rhou + src_offs), loadu(rhov + src_offs), loadu(E + src_offs));
    sw->prim[3][2] = vequation_of_state(sw->prim[0][2], E_internal);

    for(int v = 0; v < 5; ++v)
    {
        rotate_left_wm2(sw->prim[v] + 0, sw->prim[v][2]);
        rotate_left_wm1(sw->prim[v] + 1, sw->prim[v][2]);
    }

    VREAL_T dvar[4];    // slope for i+1
    if(h->iorder != 1)
        for(int v = 0; v < 4; ++ v)
            dvar[v] = vslope(sw->prim[v][0], sw->prim[v][1], sw->prim[v][2], (VREAL_T) h->slope_type, (VREAL_T) h->inv_slope_type);

    VREAL_T right_flux[4];
    const VREAL_T prim_c = vspeed_of_sound(sw->prim[4][1], sw->prim[3][1]);
    vtrace(sw->left_flux [0] + 1,                 sw->left_flux [1] + 1, sw->left_flux [2] + 1, sw->left_flux [3] + 1,
           right_flux + 0,                        right_flux + 1,        right_flux + 2,        right_flux + 3,
           sw->prim[0]       [1], sw->prim[4][1], sw->prim[1]       [1], sw->prim[2]       [1], sw->prim[3]       [1],
           dvar[0],                               dvar[1],               dvar[2],               dvar[3],
           prim_c, rcp(prim_c),
           dtdx);

    for(int v = 0; v < 4; ++v)
        rotate_left_wm1(sw->left_flux[v] + 0, sw->left_flux[v][1]);

    VREAL_T gdnv_rho, gdnv_u, gdnv_v, gdnv_p;
    vriemann(&gdnv_rho,            &gdnv_u,              &gdnv_v,              &gdnv_p,
             sw->left_flux [0][0], sw->left_flux [1][0], sw->left_flux [2][0], sw->left_flux [3][0],
             right_flux[0],        right_flux[1],        right_flux[2],        right_flux[3]);

    vcmpflx(sw->flux[0] + 1, sw->flux[1] + 1, sw->flux[2] + 1, sw->flux[3] + 1,
            gdnv_rho,        gdnv_u,          gdnv_v,          gdnv_p);

    for(int v = 0; v < 4; ++v)
        rotate_left_wm1(sw->flux[v] + 0, sw->flux[v][1]);

    const VREAL_T new_rho  = vupdate(load(rho  + i*stride), sw->flux[0][0], sw->flux[0][1], VREAL_T(dtdx));
    const VREAL_T new_rhou = vupdate(load(rhou + i*stride), sw->flux[1][0], sw->flux[1][1], VREAL_T(dtdx));
    const VREAL_T new_rhov = vupdate(load(rhov + i*stride), sw->flux[2][0], sw->flux[2][1], VREAL_T(dtdx));
    const VREAL_T new_E    = vupdate(load(E    + i*stride), sw->flux[3][0], sw->flux[3][1], VREAL_T(dtdx));

    VREAL_T courantv = (VREAL_T) 0.0;
    if(do_courant)
    {
        VREAL_T prim_rho, prim_inv_rho, prim_u, prim_v, E_internal;
        vconservative_to_primitive(&prim_rho, &prim_inv_rho, &prim_u,  &prim_v,  &E_internal,
                                   new_rho,                  new_rhou, new_rhov, new_E);
        const VREAL_T prim_p = vequation_of_state(prim_rho,     E_internal);
        const VREAL_T prim_c = vspeed_of_sound   (prim_inv_rho, prim_p);
        vcourant(&courantv, prim_u, prim_v, prim_c, write_mask);
    }

    maskstore(rho  + i*stride, new_rho,  write_mask);
    maskstore(rhou + i*stride, new_rhou, write_mask);
    maskstore(rhov + i*stride, new_rhov, write_mask);
    maskstore(E    + i*stride, new_E,    write_mask);

    for(int v = 0; v < 4; ++ v)
    {
        sw->flux     [v][0] = sw->flux     [v][1];
        sw->left_flux[v][0] = sw->left_flux[v][1];
    }
    for(int v = 0; v < 5; ++ v)
    {
        sw->prim[v][0] = sw->prim[v][2];
        sw->prim[v][1] = sw->prim[v][2];
    }

    return courantv;
}
#endif
