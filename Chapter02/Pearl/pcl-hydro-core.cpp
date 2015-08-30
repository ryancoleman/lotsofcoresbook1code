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

REAL_T           ZEROR;
REAL_T           ZEROL;
REAL_T           PROJECT;

void conservative_to_primitive(      REAL_T *restrict prim_rho, REAL_T *restrict inv_prim_rho,       REAL_T *restrict prim_u,          REAL_T *restrict prim_v,          REAL_T *restrict E_internal,
                               const REAL_T           cons_rho,                                const REAL_T           cons_rhou, const REAL_T           cons_rhov, const REAL_T           cons_E)
{
    *prim_rho     = std::max(cons_rho, SMALLR);
    *inv_prim_rho = rcp(*prim_rho);

    *prim_u = cons_rhou * *inv_prim_rho;
    *prim_v = cons_rhov * *inv_prim_rho;

    const REAL_T E_kinetic = ((REAL_T) 0.5) * (*prim_u * *prim_u + *prim_v * *prim_v);
    *E_internal            = cons_E * *inv_prim_rho - E_kinetic;

    // todo: handle passive terms here
};

REAL_T equation_of_state(const REAL_T rho,
                         const REAL_T E_internal)
{
    return std::max((GAMMA - (REAL_T) 1.0) * rho * E_internal, SMALLP);
}

REAL_T speed_of_sound(const REAL_T inv_rho,
                      const REAL_T p)
{
    return my_sqrt(GAMMA * p * inv_rho);
}

REAL_T slope(const REAL_T nbv_m, const REAL_T nbv_0, const REAL_T nbv_p,
             const REAL_T slope_type, const REAL_T inv_slope_type)

{
    const REAL_T left    = slope_type * (nbv_0 - nbv_m);
    const REAL_T right   = slope_type * (nbv_p - nbv_0);
    const REAL_T center  = ((REAL_T) 0.5) * (left + right) * inv_slope_type;
    const REAL_T sign    = (center > (REAL_T) 0.0) ? (REAL_T) 1.0 : (REAL_T) -1.0;
    const REAL_T llftrgt = (left * right) <= (REAL_T) 0.0;
    const REAL_T t1      = std::min(std::abs(left), std::abs(right));
    return sign * std::min((((REAL_T) 1.0) - llftrgt) * t1, std::abs(center));
}

void flux(      REAL_T *restrict flux_rho,                              REAL_T *restrict flux_u,         REAL_T *restrict flux_v,         REAL_T *restrict flux_p,
          const REAL_T           rho,      const REAL_T inv_rho,  const REAL_T           u,        const REAL_T           v,        const REAL_T           p,
          const REAL_T           sp_m,     const REAL_T sp_0,     const REAL_T           sp_p,
          const REAL_T           alpha_m,  const REAL_T alpha_0r, const REAL_T           alpha_0v, const REAL_T           alpha_p,
          const REAL_T           c)
{
    const REAL_T a_p  = ((REAL_T)-0.5) * sp_p * alpha_p;
    const REAL_T a_m  = ((REAL_T)-0.5) * sp_m * alpha_m;
    const REAL_T a_0r = ((REAL_T)-0.5) * sp_0 * alpha_0r;
    const REAL_T a_0v = ((REAL_T)-0.5) * sp_0 * alpha_0v;

    *flux_rho = rho + (a_p + a_m + a_0r);
    *flux_u   = u   + (a_p - a_m) * c * inv_rho;
    *flux_v   = v   + a_0v;
    *flux_p   = p   + (a_p + a_m) * c * c;
}

void trace(      REAL_T *restrict flux_rho_m,                             REAL_T *restrict flux_u_m,       REAL_T *restrict flux_v_m,       REAL_T *restrict flux_p_m,
                 REAL_T *restrict flux_rho_p,                             REAL_T *restrict flux_u_p,       REAL_T *restrict flux_v_p,       REAL_T *restrict flux_p_p,
           const REAL_T           rho,        const REAL_T inv_rho, const REAL_T           u,        const REAL_T           v,        const REAL_T           p,
           const REAL_T           drho,       const REAL_T du,      const REAL_T           dv,       const REAL_T           dp,
           const REAL_T           c,          const REAL_T inv_c,
           const REAL_T           dtdx)
{
    const REAL_T alpha_m  = ((REAL_T) 0.5) * (dp * ( inv_rho * inv_c ) - du) * rho * inv_c;
    const REAL_T alpha_p  = ((REAL_T) 0.5) * (dp * ( inv_rho * inv_c ) + du) * rho * inv_c;
    const REAL_T alpha_0r = drho - dp * (inv_c*inv_c);
    const REAL_T alpha_0v = dv;

    const REAL_T right_sp_m = ((u - c) >= ZEROR) ? PROJECT : (u - c) * dtdx + ((REAL_T) 1.0);
    const REAL_T right_sp_p = ((u + c) >= ZEROR) ? PROJECT : (u + c) * dtdx + ((REAL_T) 1.0);
    const REAL_T right_sp_0 =  (u      >= ZEROR) ? PROJECT :  u      * dtdx + ((REAL_T) 1.0);

    flux(flux_rho_p,          flux_u_p,   flux_v_p,   flux_p_p,
         rho,        inv_rho, u,          v,          p,
         right_sp_m,          right_sp_0, right_sp_p,
         alpha_m,             alpha_0r,   alpha_0v,   alpha_p,
         c);

    const REAL_T left_sp_m = ((u - c) <= ZEROL) ? -PROJECT : (u - c) * dtdx - ((REAL_T) 1.0);
    const REAL_T left_sp_p = ((u + c) <= ZEROL) ? -PROJECT : (u + c) * dtdx - ((REAL_T) 1.0);
    const REAL_T left_sp_0 =  (u      <= ZEROL) ? -PROJECT :  u      * dtdx - ((REAL_T) 1.0);

    flux(flux_rho_m,          flux_u_m,   flux_v_m,   flux_p_m,
         rho,        inv_rho, u,          v,          p,
         left_sp_m,           left_sp_0,  left_sp_p,
         alpha_m,             alpha_0r,   alpha_0v,   alpha_p,
         c);

    // todo: handle passive terms
}

void riemann(      REAL_T *restrict     gdnv_rho,       REAL_T *restrict gdnv_u,           REAL_T *restrict gdnv_v,           REAL_T *restrict gdnv_p,
             const REAL_T            in_left_rho, const REAL_T           in_left_u,  const REAL_T           in_left_v,  const REAL_T           in_left_p,
             const REAL_T           in_right_rho, const REAL_T           in_right_u, const REAL_T           in_right_v, const REAL_T           in_right_p)

{
    const REAL_T left_rho = std::max(in_left_rho, SMALLR);
    const REAL_T left_u   = in_left_u;
    const REAL_T left_v   = in_left_v;
    const REAL_T left_p   = std::max(in_left_p, left_rho * SMALLP);
    const REAL_T left_c   = GAMMA * left_p * left_rho;

    const REAL_T right_rho = std::max(in_right_rho, SMALLR);
    const REAL_T right_u   = in_right_u;
    const REAL_T right_v   = in_right_v;
    const REAL_T right_p   = std::max(in_right_p, right_rho * SMALLP);
    const REAL_T right_c   = GAMMA * right_p * right_rho;

    REAL_T p_star;
    bool   goon = true;
    {
        const REAL_T left_w  = my_sqrt(left_c);
        const REAL_T right_w = my_sqrt(right_c);
        p_star               = std::max( (right_w * left_p + left_w * right_p + left_w * right_w * (left_u - right_u)) * rcp(left_w + right_w), (REAL_T) 0.0);
}

    for(int i = 0; i < NITER_RIEMANN; ++i)
    {
        if(goon)
        {
            const REAL_T left_ww2  = left_rho  * ((REAL_T) 0.5) * ((GAMMA + (REAL_T) 1.0) * p_star + (GAMMA - (REAL_T) 1.0) * left_p);
            const REAL_T left_ww   = my_sqrt(left_ww2);
            const REAL_T right_ww2 = right_rho * ((REAL_T) 0.5) * ((GAMMA + (REAL_T) 1.0) * p_star + (GAMMA - (REAL_T) 1.0) * right_p);
            const REAL_T right_ww  = my_sqrt(right_ww2);
            const REAL_T tmp_num   = ((REAL_T)2.0) * left_ww2 * right_ww2 * (left_ww * right_ww * (left_u - right_u) - left_ww * (p_star - right_p) - right_ww * (p_star - left_p));
            const REAL_T tmp_den   = right_ww2*right_ww * (left_ww2 + left_c) + left_ww2*left_ww * (right_ww2 + right_c);
            const REAL_T tmp       = tmp_num * rcp(tmp_den);
            const REAL_T deleft_p  = std::max(tmp, -p_star);

            p_star += deleft_p;

            const REAL_T uo = std::abs(deleft_p * rcp(p_star + SMALLPP));
            goon            = uo > PRECISION;
        }
    }

    const REAL_T left_w2  = left_rho  * ((REAL_T) 0.5) * ((GAMMA + (REAL_T) 1.0) * p_star + (GAMMA - (REAL_T) 1.0) * left_p);
    const REAL_T left_w   = my_sqrt(left_w2);
    const REAL_T right_w2 = right_rho * ((REAL_T) 0.5) * ((GAMMA + (REAL_T) 1.0) * p_star + (GAMMA - (REAL_T) 1.0) * right_p);
    const REAL_T right_w  = my_sqrt(right_w2);

    const REAL_T u_star = (REAL_T) 0.5 * (left_u + (left_p - p_star) * rcp(left_w) + right_u - (right_p - p_star) * rcp(right_w));

    const REAL_T sgnm  = (u_star > (REAL_T) 0.0) ? (REAL_T)1.0 : (REAL_T)-1.0;
    const REAL_T rho_0 = (u_star > (REAL_T) 0.0) ? left_rho    : right_rho;
    const REAL_T u_0   = (u_star > (REAL_T) 0.0) ? left_u      : right_u;
    const REAL_T p_0   = (u_star > (REAL_T) 0.0) ? left_p      : right_p;
    const REAL_T w_0   = (u_star > (REAL_T) 0.0) ? left_w      : right_w;
    const REAL_T w2_0  = (u_star > (REAL_T) 0.0) ? left_w2     : right_w2;

    const REAL_T inv_rho_0 = rcp(rho_0);
    const REAL_T c_0       = std::max(my_sqrt(std::abs(GAMMA * p_0 * inv_rho_0)),   SMALLC);
    const REAL_T rho_star  = std::max(w2_0*rho_0 * rcp(w2_0 + rho_0 * (p_0 - p_star)),   SMALLR);
    const REAL_T c_star    = std::max(my_sqrt(std::abs(GAMMA * p_star * rcp(rho_star))), SMALLC);

    const REAL_T ushock = w_0 * inv_rho_0 - sgnm * u_0;

    const REAL_T spout  = (p_star >= p_0) ? ushock : c_0    - sgnm * u_0;
    const REAL_T spin   = (p_star >= p_0) ? ushock : c_star - sgnm * u_star;

    REAL_T frac;
    if(spout < (REAL_T) 0.0)
        frac = (REAL_T)0.0;
    else if(spin > (REAL_T)0.0)
        frac = (REAL_T)1.0;
    else
    {
        const REAL_T scr = std::max(spout - spin, SMALLC + std::abs(spout + spin));
        frac             = std::max(std::min((((REAL_T) 1.0) + (spout + spin) * rcp(scr)) * ((REAL_T) 0.5), (REAL_T) 1.0), (REAL_T) 0.0);
    }

    *gdnv_rho = (frac * rho_star + (((REAL_T) 1.0) - frac) * rho_0);
    *gdnv_u   = (frac * u_star   + (((REAL_T) 1.0) - frac) * u_0);
    *gdnv_v   = (u_star > (REAL_T)0.0) ? left_v : right_v;
    *gdnv_p   = (frac * p_star   + (((REAL_T) 1.0) - frac) * p_0);

    // todo: handle passive
}

void cmpflx(      REAL_T *restrict flux_rho,       REAL_T *restrict flux_rhou,       REAL_T *restrict flux_rhov,       REAL_T *restrict flux_E,
            const REAL_T           gdnv_rho, const REAL_T           gdnv_u,    const REAL_T           gdnv_v,    const REAL_T           gdnv_p)
{
    const REAL_T mass_density = gdnv_rho * gdnv_u;
    *flux_rho                 = mass_density;
    *flux_rhou                = mass_density * gdnv_u + gdnv_p;
    *flux_rhov                = mass_density * gdnv_v;

    const REAL_T E_kinetic = ((REAL_T) 0.5) * gdnv_rho * (gdnv_u * gdnv_u + gdnv_v * gdnv_v);
    const REAL_T E_total   = gdnv_p * (((REAL_T) 1.0) / (GAMMA - ((REAL_T) 1.0))) + E_kinetic;
    *flux_E                = gdnv_u * (E_total + gdnv_p);

    // todo: handle passive
}

REAL_T update(const REAL_T  in,
              const REAL_T  flux_left, const REAL_T  flux_right,
              const REAL_T  dtdx)
{
    return in  + (flux_left  - flux_right)  * dtdx;

    // todo: handle passive
}

void courant(      REAL_T *restrict courantv,
             const REAL_T           u, const REAL_T  v,
             const REAL_T           c)
{
    *courantv = std::max(*courantv, std::max(c + std::abs(u), c + std::abs(v)));
}

void strip_prime(strip_work   *restrict sw,
                 const REAL_T *restrict rho,
                 const REAL_T *restrict rhou,
                 const REAL_T *restrict rhov,
                 const REAL_T *restrict E,
                 const hydro  *restrict h,
                 const int              stride,
                 const REAL_T           dtdx)
{
    REAL_T pre_prim[5][4]; // i-2, i-1, i, i+1
    for(int i = 0; i < 4; ++i)
    {
        const int src_offs = (- 2 + i)*stride;
        REAL_T E_internal;
        conservative_to_primitive(pre_prim[0] +i, pre_prim[4] +i, pre_prim[1] +i, pre_prim[2] +i, &E_internal,
                                   rho[src_offs],                 rhou[src_offs], rhov[src_offs], E[src_offs]);
        pre_prim[3][i] = equation_of_state(pre_prim[0][i], E_internal);
    }

    REAL_T pre_dvar[4][2]; // i-1, i
    if(h->iorder != 1)
        for(int i = 0; i < 2; ++i)
            for(int v = 0; v < 4; ++ v)
                pre_dvar[v][i] = slope(pre_prim[v][0+i], pre_prim[v][1+i], pre_prim[v][2+i], h->slope_type, h->inv_slope_type);

    REAL_T pre_left_flux [4][2]; // i-1, i
    REAL_T pre_right_flux[4][2]; // i-1, i
    for(int i = 0; i < 2; ++i)
    {
        const REAL_T prim_c = speed_of_sound(pre_prim[4][i + 1], pre_prim[3][i+1]);
        trace(pre_left_flux [0] + i,                   pre_left_flux [1] + i,  pre_left_flux [2] + i, pre_left_flux [3] + i,
              pre_right_flux[0] + i,                   pre_right_flux[1] + i,  pre_right_flux[2] + i, pre_right_flux[3] + i,
              pre_prim[0]     [i+1], pre_prim[4][i+1], pre_prim[1]      [i+1], pre_prim[2]     [i+1], pre_prim[3]     [i+1],
              pre_dvar[0]       [i],                   pre_dvar[1]        [i], pre_dvar[2]       [i], pre_dvar[3]       [i],
              prim_c, rcp(prim_c),
              dtdx);
    }

    REAL_T gdnv_rho, gdnv_u, gdnv_v, gdnv_p;
    riemann(&gdnv_rho,        &gdnv_u,          &gdnv_v,          &gdnv_p,
            pre_left_flux [0][0], pre_left_flux [1][0], pre_left_flux [2][0], pre_left_flux [3][0],
            pre_right_flux[0][1], pre_right_flux[1][1], pre_right_flux[2][1], pre_right_flux[3][1]);

    cmpflx(sw->flux[0] + 0, sw->flux[1] + 0, sw->flux[2] + 0, sw->flux[3] + 0,
           gdnv_rho,        gdnv_u,          gdnv_v,          gdnv_p);

    for(int v = 0; v < 4; ++ v)
        sw->left_flux[v][0] = pre_left_flux[v][1];

    for(int v = 0; v < 5; ++ v)
    {
        sw->prim[v][0] = pre_prim[v][2];
        sw->prim[v][1] = pre_prim[v][3];
    }
}

REAL_T strip_stable(const hydro *restrict h,
                    REAL_T      *restrict rho,
                    REAL_T      *restrict rhou,
                    REAL_T      *restrict rhov,
                    REAL_T      *restrict E,
                    strip_work  *restrict sw,
                    const int             i,
                    const int             stride,
                    const REAL_T          dtdx,
                    const bool            do_courant)
{
    const int src_offs = (i + 2)*stride;

    REAL_T E_internal;
    conservative_to_primitive(sw->prim[0] + 2, sw->prim[4] + 2, sw->prim[1] + 2, sw->prim[2] + 2, &E_internal,
                              rho[src_offs],                    rhou[src_offs],  rhov[src_offs], E[src_offs]);
    sw->prim[3][2] = equation_of_state(sw->prim[0][2], E_internal);

    REAL_T dvar[4];    // slope for i+1
    if(h->iorder != 1)
        for(int v = 0; v < 4; ++ v)
            dvar[v] = slope(sw->prim[v][0], sw->prim[v][1], sw->prim[v][2], h->slope_type, h->inv_slope_type);

    REAL_T right_flux[4];
    const REAL_T prim_c = speed_of_sound(sw->prim[4][1], sw->prim[3][1]);
    trace(sw->left_flux [0] + 1,                 sw->left_flux [1] + 1, sw->left_flux [2] + 1, sw->left_flux [3] + 1,
          right_flux + 0,                        right_flux + 1,        right_flux + 2,        right_flux + 3,
          sw->prim[0]       [1], sw->prim[4][1], sw->prim[1]       [1], sw->prim[2]       [1], sw->prim[3]       [1],
          dvar[0],                               dvar[1],              dvar[2],               dvar[3],
          prim_c, rcp(prim_c),
          dtdx);

    REAL_T gdnv_rho, gdnv_u, gdnv_v, gdnv_p;
    riemann(&gdnv_rho,            &gdnv_u,              &gdnv_v,              &gdnv_p,
            sw->left_flux [0][0], sw->left_flux [1][0], sw->left_flux [2][0], sw->left_flux [3][0],
            right_flux[0],        right_flux[1],        right_flux[2],        right_flux[3]);

    cmpflx(sw->flux[0] + 1, sw->flux[1] + 1, sw->flux[2] + 1, sw->flux[3] + 1,
           gdnv_rho,        gdnv_u,          gdnv_v,          gdnv_p);

    const REAL_T new_rho  = update(rho [i*stride], sw->flux[0][0], sw->flux[0][1], dtdx);
    const REAL_T new_rhou = update(rhou[i*stride], sw->flux[1][0], sw->flux[1][1], dtdx);
    const REAL_T new_rhov = update(rhov[i*stride], sw->flux[2][0], sw->flux[2][1], dtdx);
    const REAL_T new_E    = update(E   [i*stride], sw->flux[3][0], sw->flux[3][1], dtdx);

    REAL_T courantv = (REAL_T) 0.0;
    if(do_courant)
    {
        REAL_T prim_rho, prim_inv_rho, prim_u, prim_v, E_internal;
        conservative_to_primitive(&prim_rho, &prim_inv_rho, &prim_u,  &prim_v,  &E_internal,
                                  new_rho,                  new_rhou, new_rhov, new_E);
        const REAL_T prim_p = equation_of_state(prim_rho,     E_internal);
        const REAL_T prim_c = speed_of_sound   (prim_inv_rho, prim_p);
        courant(&courantv, prim_u, prim_v, prim_c);
    }

    rho [i*stride] = new_rho;
    rhou[i*stride] = new_rhou;
    rhov[i*stride] = new_rhov;
    E   [i*stride] = new_E;

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
