/*
  Helper files
  (C) Jason Sewall : Intel
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
#ifndef __ARCH_HPP__
#define __ARCH_HPP__

#undef HAVE_SIMD_TYPE

#ifdef KNC
#define HAVE_SIMD_TYPE
#define SIMD_STR "KNC"
#endif

#ifdef AVX3
#define HAVE_SIMD_TYPE
#define SIMD_STR "AVX3"
#endif

#ifdef AVX
#ifdef HAVE_SIMD_TYPE
#error "Multiple SIMD types given (AVX3, KNC, SSE, AVX)"
#endif
#define HAVE_SIMD_TYPE
#define SIMD_STR "AVX"
#endif

#ifdef SSE
#ifdef HAVE_SIMD_TYPE
#error "Multiple SIMD types given (AVX3, KNC, SSE, AVX)"
#endif
#define HAVE_SIMD_TYPE
#define SIMD_STR "SSE"
#endif

#ifndef HAVE_SIMD_TYPE
#warning "No SIMD type given (AVX3, KNC, SSE, AVX)"
#define SIMD_STR "NONE"
#endif

#undef HAVE_PRECISION

#ifdef SINGLE
#define HAVE_PRECISION
#endif

#ifdef DOUBLE
#ifdef HAVE_PRECISION
#error "Multiple precisions given (SINGLE, DOUBLE)"
#endif
#define HAVE_PRECISION
#endif

#ifndef HAVE_PRECISION
#error "No precision given (SINGLE, DOUBLE)"
#endif

#include <cmath>
#include <algorithm>

#ifdef SINGLE
typedef float REAL_T;
#define REAL_FMT "%f"

inline float rcp(const float x)
{
    float res;
    _mm_store_ss(&res, _mm_rcp_ss(_mm_load_ss(&x)));
    return res;
}

inline float my_sqrt(const float z)
{
    union
    {
        int tmp;
        float f;
    } u;

    u.f     = z;
    u.tmp  -= 1 << 23;          /* Subtract 2^m. */
    u.tmp >>= 1;                /* Divide by 2. */
    u.tmp  += 1 << 29;          /* Add ((b + 1) / 2) * 2^m. */

    return u.f;
}

#ifndef HAVE_SIMD_TYPE
typedef float VREAL_T;
typedef bool  VMASK_T;
#define SIMD_WIDTH 1
#include "xmmintrin.h"

#endif // HAVE_SIMD_TYPE

#ifdef AVX
typedef F32vec8 VREAL_T;
typedef F32vec8 VMASK_T;
#define SIMD_WIDTH 8
#endif // AVX

#ifdef KNC
typedef F32vec16 VREAL_T;
typedef F32vec16 VMASK_T;
#define SIMD_WIDTH 16
#endif // KNC
#endif // SINGLE

#ifdef DOUBLE
typedef double REAL_T;
#define REAL_FMT "%lf"

inline double rcp(const double x)
{
    return 1.0/x;
}

inline double my_sqrt(const double z)
{
    return std::sqrt(z);
}

#ifndef HAVE_SIMD_TYPE
typedef double VREAL_T;
typedef bool   VMASK_T;
#define SIMD_WIDTH 1


#endif // HAVE_SIMD_TYPE

#ifndef HAVE_SIMD_TYPE
typedef int VINT_T;

inline VREAL_T load(const REAL_T *r)
{
    return *r;
}

inline VREAL_T loadu(const REAL_T *r)
{
    return *r;
}

inline void maskstore(REAL_T *d, const VREAL_T s)
{
    *d = s;
}

inline void maskstore(REAL_T *d, const VREAL_T s, const VMASK_T m)
{
    if(m)
        *d = s;
}

inline void mask_true(bool *mask)
{
    *mask = true;
}

inline bool mask_not(bool mask)
{
    return !mask;
}

inline void linear_offset(VINT_T *linear)
{
    *linear = 0;
}

inline REAL_T select_true(const bool &mask, const REAL_T &iftrue, const REAL_T &iffalse)
{
    return mask ? iftrue : iffalse;
}

inline bool mask_gt(const REAL_T &l, const REAL_T &r)
{
    return l > r;
}

inline bool mask_lt(const REAL_T &l, const REAL_T &r)
{
    return l < r;
}

inline bool mask_and(const bool &l, const bool &r)
{
    return l && r;
}

inline bool mask_or(const bool &l, const bool &r)
{
    return l || r;
}

inline bool all_zero(const bool &l)
{
    return l == 0;
}

inline REAL_T select_lt(const REAL_T &a, const REAL_T &b, const REAL_T &c, const REAL_T &d)
{ return (a < b) ? c : d; }

inline REAL_T select_le(const REAL_T &a, const REAL_T &b, const REAL_T &c, const REAL_T &d)
{ return (a <= b) ? c : d; }

inline REAL_T select_gt(const REAL_T &a, const REAL_T &b, const REAL_T &c, const REAL_T &d)
{ return (a > b) ? c : d; }

inline REAL_T select_ge(const REAL_T &a, const REAL_T &b, const REAL_T &c, const REAL_T &d)
{ return (a >= b) ? c : d; }

#endif // HAVE_SIMD_TYPE

#ifdef SSE
#include <dvec.h>
#include "immintrin.h"

typedef F64vec2 VREAL_T;
typedef F64vec2 VMASK_T;
typedef F64vec2 VINT_T;
#define SIMD_WIDTH 2

namespace std
{
    inline F64vec2 sqrt(const F64vec2 &v)
    {
        return _mm_sqrt_pd(v);
    }

    inline F64vec2 abs(const F64vec2 &a)
    {
        static const union
        {
            int i[4];
            __m128d m;
        } __f64vec2_abs_mask = { 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff };
        return _mm_and_pd(a, __f64vec2_abs_mask.m);
    }

    inline F64vec2 max(const F64vec2 &l, const F64vec2 &r)
    {
        return _mm_max_pd(l, r);
    }

    inline F64vec2 min(const F64vec2 &l, const F64vec2 &r)
    {
        return _mm_min_pd(l, r);
    }
}

inline F64vec2 rcp(const F64vec2 &v)
{
    return F64vec2(1.0)/v;
}

inline F64vec2 my_sqrt(const F64vec2 &v)
{
    return std::sqrt(v);
}

inline void mask_true(F64vec2 *mask)
{
    static const union
    {
        int i[4];
        __m128d m;
    } __f64vec2_true = { 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff };
    *mask = __f64vec2_true.m;
}

inline void linear_offset(F64vec2 *linear)
{
    *linear = F64vec2(1, 0);
}

inline bool all_zero(const F64vec2 &l)
{
    return l.is_zero();
}

inline F64vec2 select_true(const F64vec2 &mask, const F64vec2 &iftrue, const F64vec2 &iffalse)
{
    return _mm_blendv_pd(iffalse, iftrue, mask);
}

inline F64vec2 mask_not(const F64vec2 &l)
{
    static const union
    {
        int i[4];
        __m128d m;
    } __f64vec2_true = { 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff };

    return _mm_andnot_pd(l, __f64vec2_true.m);
}

inline F64vec2 mask_gt(const F64vec2 &l, const F64vec2 &r)
{
    return _mm_cmp_pd(l, r, _CMP_GT_OS);
}

inline F64vec2 mask_lt(const F64vec2 &l, const F64vec2 &r)
{
    return _mm_cmp_pd(l, r, _CMP_LT_OS);
}

inline F64vec2 mask_and(const F64vec2 &l, const F64vec2 &r)
{
    return _mm_and_pd(l, r);
}

inline F64vec2 mask_or(const F64vec2 &l, const F64vec2 &r)
{
    return _mm_or_pd(l, r);
}

inline F64vec2 loadu(const double *r)
{
    F64vec2 res;
    loadu(res, const_cast<REAL_T*>(r));
    return res;
}

inline F64vec2 load(const double *r)
{
    return _mm_load_pd(r);
}

inline void store(double *r, const F64vec2 v)
{
    _mm_store_pd(r, v);
}

inline void rotate_left_wm2(F64vec2 *v0, const F64vec2 v1)
{
    //v0 {1.0,  2.0};
    //v1 {3.0,  4.0};

    //v0 {1.0, 2.0};
}

inline void rotate_left_wm1(F64vec2 *v0, const F64vec2 v1)
{
    //v0 {1.0,  2.0};
    //v1 {3.0,  4.0};

    //v0 {2.0, 3.0, 4.0};
    *v0 = _mm_castsi128_pd(_mm_alignr_epi8(_mm_castpd_si128(v1), _mm_castpd_si128(*v0), 8));
}

#endif // SSE

#ifdef AVX
#include <dvec.h>
#include "immintrin.h"

typedef F64vec4 VREAL_T;
typedef F64vec4 VMASK_T;
typedef F64vec4 VINT_T;
#define SIMD_WIDTH 4

namespace std
{
    inline F64vec4 sqrt(const F64vec4 &v)
    {
        return _mm256_sqrt_pd(v);
    }

    inline F64vec4 abs(const F64vec4 &a)
    {
        static const union
        {
            int i[8];
            __m256d m;
        } __f64vec4_abs_mask = { 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff};
        return _mm256_and_pd(a, __f64vec4_abs_mask.m);
    }

    inline F64vec4 max(const F64vec4 &l, const F64vec4 &r)
    {
        return _mm256_max_pd(l, r);
    }

    inline F64vec4 min(const F64vec4 &l, const F64vec4 &r)
    {
        return _mm256_min_pd(l, r);
    }
}

inline F64vec4 rcp(const F64vec4 &v)
{
    return F64vec4(1.0)/v;
}

inline F64vec4 my_sqrt(const F64vec4 &v)
{
    return std::sqrt(v);
}

inline void mask_true(F64vec4 *mask)
{
    static const union
    {
        int i[8];
        __m256d m;
    } __f64vec4_true = { 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff};

    *mask = __f64vec4_true.m;
}

inline void linear_offset(F64vec4 *linear)
{
    *linear = F64vec4(3, 2, 1, 0);
}

inline bool all_zero(const F64vec4 &l)
{
    return l.is_zero();
}

inline F64vec4 select_true(const F64vec4 &mask, const F64vec4 &iftrue, const F64vec4 &iffalse)
{
    return _mm256_blendv_pd(iffalse, iftrue, mask);
}

inline F64vec4 mask_not(const F64vec4 &l)
{
    static const union
    {
        int i[8];
        __m256d m;
    } __f64vec4_true = { 0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff,
                         0xffffffff, 0xffffffff, 0xffffffff, 0xffffffff};

    return _mm256_andnot_pd(l, __f64vec4_true.m);
}

inline F64vec4 mask_gt(const F64vec4 &l, const F64vec4 &r)
{
    return _mm256_cmp_pd(l, r, _CMP_GT_OS);
}

inline F64vec4 mask_lt(const F64vec4 &l, const F64vec4 &r)
{
    return _mm256_cmp_pd(l, r, _CMP_LT_OS);
}

inline F64vec4 mask_and(const F64vec4 &l, const F64vec4 &r)
{
    return _mm256_and_pd(l, r);
}

inline F64vec4 mask_or(const F64vec4 &l, const F64vec4 &r)
{
    return _mm256_or_pd(l, r);
}

inline F64vec4 loadu(const double *r)
{
    F64vec4 res;
    loadu(res, const_cast<REAL_T*>(r));
    return res;
}

inline F64vec4 load(const double *r)
{
    return _mm256_load_pd(r);
}

inline void store(double *r, const F64vec4 v)
{
    _mm256_store_pd(r, v);
}

inline void rotate_left_wm2(F64vec4 *v0, const F64vec4 v1)
{
    *v0 = _mm256_castpd128_pd256(_mm256_extractf128_pd(*v0, 1));
    *v0 = _mm256_insertf128_pd(*v0, _mm256_castpd256_pd128(v1), 1);
}

inline void rotate_left_wm1(F64vec4 *v0, const F64vec4 v1)
{
    // {1.0, 2.0, 3.0, 4.0};
    // {5.0, 6.0, 7.0, 8.0};

    const __m128d hiv0      = _mm256_extractf128_pd(*v0, 1); // {3.0, 4.0}
    const __m128d phiv0     = _mm_permute_pd(hiv0, 0x1); // {4.0, 3.0}
    const __m256d shufv1    = _mm256_permute_pd(v1, 0x1); // {6.0, 5.0, 8.0, 7.0};
    const __m128d shufv1_lo = _mm256_extractf128_pd(shufv1, 0); // {6.0, 5.0}
    const __m128d shufv1_hi = _mm256_extractf128_pd(shufv1, 1); // {8.0, 7.0}
    const __m128d v1_blend  = _mm_blend_pd(shufv1_lo, shufv1_hi, 0x2); // blend   {6.0, 7.0};
    const __m256d inserted  = _mm256_insertf128_pd(shufv1, v1_blend, 1); // insert  {6.0, 5.0, 6.0, 7.0};
    const __m256d blended   = _mm256_blend_pd(_mm256_castpd128_pd256(phiv0), inserted, 0xE);
    *v0                     = blended;
}

#endif // AVX

#ifdef KNC
#include <micvec.h>
#include "zmmintrin.h"

typedef F64vec8 VREAL_T;
typedef __mmask8 VMASK_T;
typedef F64vec8 VINT_T;
#define SIMD_WIDTH 8

namespace std
{
    inline F64vec8 sqrt(const F64vec8 &v)
    {
        return _mm512_sqrt_pd(v);
    }

    inline F64vec8 abs(const F64vec8 &a)
    {
        static const union
        {
            int i[16];
            __m512 m;
        } __i64vec8_abs_mask = { 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff};
        return _mm512_and_epi32(_mm512_castpd_si512(a), __i64vec8_abs_mask.m);
    }

    inline F64vec8 max(const F64vec8 &l, const F64vec8 &r)
    {
        return _mm512_max_pd(l, r);
    }

    inline F64vec8 min(const F64vec8 &l, const F64vec8 &r)
    {
        return _mm512_min_pd(l, r);
    }
}

inline F64vec8 rcp(const F64vec8 &v)
{
    return F64vec8(1.0)/v;
}

inline F64vec8 my_sqrt(const F64vec8 &v)
{
    return std::sqrt(v);
}

inline void mask_true(__mmask8 *mask)
{
    *mask = 0xFF;
}

inline void maskstore(double *d, const F64vec8 s, const __mmask8 m)
{
    _mm512_mask_store_pd(d, m, s);
}

inline void linear_offset(F64vec8 *linear)
{
    *linear = F64vec8(7, 6, 5, 4, 3, 2, 1, 0);
}

inline bool all_zero(const __mmask8 &l)
{
    return l == 0;
}

inline F64vec8 select_true(const __mmask8 &mask, const F64vec8 &iftrue, const F64vec8 &iffalse)
{
    return _mm512_castsi512_pd(_mm512_mask_or_epi64(_mm512_castpd_si512(iffalse), mask, _mm512_castpd_si512(iftrue), _mm512_castpd_si512(iftrue)));
}

inline F64vec8 select_gt(const F64vec8 &a, const F64vec8 &b, const F64vec8 &c, const F64vec8 &d)
{
    return select_lt(b, a, c, d);
}

inline F64vec8 select_ge(const F64vec8 &a, const F64vec8 &b, const F64vec8 &c, const F64vec8 &d)
{
    return select_le(b, a, c, d);
}

inline __mmask8 mask_not(const __mmask8 &l)
{
    return ~l;
}

inline __mmask8 mask_gt(const F64vec8 &l, const F64vec8 &r)
{
    return _mm512_cmp_pd_mask(l, r, _CMP_GT_OS);
}

inline __mmask8 mask_lt(const F64vec8 &l, const F64vec8 &r)
{
    return _mm512_cmp_pd_mask(l, r, _CMP_LT_OS);
}

inline __mmask8 mask_and(const __mmask8 &l, const __mmask8 &r)
{
    return l & r;
}

inline __mmask8 mask_or(const __mmask8 &l, const __mmask8 &r)
{
    return l | r;
}

inline F64vec8 load(const double *r)
{
    return _mm512_load_pd(r);
}

#pragma warning disable 592
inline F64vec8 loadu(const double *d)
{
    __m512d vdst;
    vdst = _mm512_loadunpacklo_pd(vdst, d);
    vdst = _mm512_loadunpackhi_pd(vdst, d+8);
    return vdst;
}

inline void store(double *r, const F64vec8 v)
{
    _mm512_store_pd(r, v);
}

inline void rotate_left_wm2(F64vec8 *v0, const F64vec8 v1)
{
    //v0 {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0};
    //v1 {9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    //v0 {7.0,  8.0,  9.0,  10.0,  11.0,  12.0,  13.0, 14.0};

    static const I32vec16 shift2(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12);
    *v0 = _mm512_permutevar_epi32     (                                   shift2, _mm512_castpd_si512(*v0));
    *v0 = _mm512_mask_permutevar_epi32(_mm512_castpd_si512(*v0), 0xFFF0U, shift2, _mm512_castpd_si512(v1));
}

inline void rotate_left_wm1(F64vec8 *v0, const F64vec8 v1)
{
    //v0 {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0};
    //v1 {9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    //v0 {8, 9, 10, 11, 12, 13, 14, 15};

    static const I32vec16 shift1(13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14);
    *v0 = _mm512_permutevar_epi32     (                                   shift1, _mm512_castpd_si512(*v0));
    *v0 = _mm512_mask_permutevar_epi32(_mm512_castpd_si512(*v0), 0xFFFCU, shift1, _mm512_castpd_si512(v1));
}
#endif // KNC

#ifdef AVX3
#include "micvec.h"
#include "ia32intrin.h"

typedef F64vec8 VREAL_T;
typedef __mmask8 VMASK_T;
typedef F64vec8 VINT_T;
#define SIMD_WIDTH 8

namespace std
{
    inline F64vec8 sqrt(const F64vec8 &v)
    {
        return _mm512_sqrt_pd(v);
    }

    inline F64vec8 abs(const F64vec8 &a)
    {
        static const union
        {
            int i[16];
            __m512 m;
        } __i64vec8_abs_mask = { 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff,
                                 0xffffffff, 0x7fffffff, 0xffffffff, 0x7fffffff};
        return _mm512_and_epi32(_mm512_castpd_si512(a), __i64vec8_abs_mask.m);
    }

    inline F64vec8 max(const F64vec8 &l, const F64vec8 &r)
    {
        return _mm512_max_pd(l, r);
    }

    inline F64vec8 min(const F64vec8 &l, const F64vec8 &r)
    {
        return _mm512_min_pd(l, r);
    }
}

inline F64vec8 rcp(const F64vec8 &v)
{
    return F64vec8(1.0)/v;
}

inline F64vec8 my_sqrt(const F64vec8 &v)
{
    return std::sqrt(v);
}

inline void mask_true(__mmask8 *mask)
{
    *mask = 0xFF;
}

inline void maskstore(double *d, const F64vec8 s, const __mmask8 m)
{
    _mm512_mask_store_pd(d, m, s);
}

inline void linear_offset(F64vec8 *linear)
{
    *linear = F64vec8(7, 6, 5, 4, 3, 2, 1, 0);
}

inline bool all_zero(const __mmask8 &l)
{
    return l == 0;
}

inline F64vec8 select_true(const __mmask8 &mask, const F64vec8 &iftrue, const F64vec8 &iffalse)
{
    return _mm512_castsi512_pd(_mm512_mask_or_epi64(_mm512_castpd_si512(iffalse), mask, _mm512_castpd_si512(iftrue), _mm512_castpd_si512(iftrue)));
}

inline F64vec8 select_gt(const F64vec8 &a, const F64vec8 &b, const F64vec8 &c, const F64vec8 &d)
{
    return select_lt(b, a, c, d);
}

inline F64vec8 select_ge(const F64vec8 &a, const F64vec8 &b, const F64vec8 &c, const F64vec8 &d)
{
    return select_le(b, a, c, d);
}

inline __mmask8 mask_not(const __mmask8 &l)
{
    return ~l;
}

inline __mmask8 mask_gt(const F64vec8 &l, const F64vec8 &r)
{
    return _mm512_cmp_pd_mask(l, r, _CMP_GT_OS);
}

inline __mmask8 mask_lt(const F64vec8 &l, const F64vec8 &r)
{
    return _mm512_cmp_pd_mask(l, r, _CMP_LT_OS);
}

inline __mmask8 mask_and(const __mmask8 &l, const __mmask8 &r)
{
    return l & r;
}

inline __mmask8 mask_or(const __mmask8 &l, const __mmask8 &r)
{
    return l | r;
}

inline F64vec8 load(const double *r)
{
    return _mm512_load_pd(r);
}

inline F64vec8 loadu(const double *d)
{
   return _mm512_loadu_pd(d);
}

inline void store(double *r, const F64vec8 v)
{
    _mm512_store_pd(r, v);
}

inline void rotate_left_wm2(F64vec8 *v0, const F64vec8 v1)
{
    //v0 {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0};
    //v1 {9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    //v0 {7.0,  8.0,  9.0,  10.0,  11.0,  12.0,  13.0, 14.0};

    static const I32vec16 shift2(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12);
    *v0 = _mm512_permutevar_epi32     (                                   shift2, _mm512_castpd_si512(*v0));
    *v0 = _mm512_mask_permutevar_epi32(_mm512_castpd_si512(*v0), 0xFFF0U, shift2, _mm512_castpd_si512(v1));
}

inline void rotate_left_wm1(F64vec8 *v0, const F64vec8 v1)
{
    //v0 {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0};
    //v1 {9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    //v0 {8, 9, 10, 11, 12, 13, 14, 15};

    static const I32vec16 shift1(13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 15, 14);
    *v0 = _mm512_permutevar_epi32     (                                   shift1, _mm512_castpd_si512(*v0));
    *v0 = _mm512_mask_permutevar_epi32(_mm512_castpd_si512(*v0), 0xFFFCU, shift1, _mm512_castpd_si512(v1));
}
#endif // AVX3
#endif // DOUBLE

#endif /* __ARCH_HPP__ */
