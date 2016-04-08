!Copyright (c) 2014, Per Berg and Jacob Weismann Poulsen, DMI
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met: 
! 
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer. 
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution. 
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
!The views and conclusions contained in the software and documentation are those
!of the authors and should not be interpreted as representing official policies, 
!either expressed or implied, of the FreeBSD Project.

module tflow_simd_srf_col
  use constants, only : qrt, half, zero, one, two,                             &
                        onethird, twothird, fourthird
  implicit none
  private

  !- Private work arrays -------------------------------------------------------
  !
  !  Now 9 work arrays, previous it was 15.
  !  Substitution as compared to ealier versions of tflow:
  !    tu/v/w  --> t1/2/3
  !    dx/y/zt --> t4/5/6
  !    dtu/v/w --> t1/2/3
  !    rin/out --> t7/8
  !    cx/y/z  --> t4/5/6
  !
  !  Computational flowchart (t,u,v,w ao args are implicit):
  !  A: c_tu/v/w --> t1/2/3 --> halo(t1/2) --> c_tt --> tt (release t1/2/3)
  !  B: c_delta --> t4/5/6 --> halo(t4/5/6) --> c_dt --> t1/2/3 (release t4/5/6)
  !  C: tt, t1/2/3 --> c_rin_rout --> t7/8 
  !  D: t1/2/3, t7/8 --> c_cx_cy_cz --> t1/2/3 (release t7/8)
  !  E: tt, t1/2/3, t4/5/6 --> t (release tt, t1/2/3/4/5/6/7)
  !
  !-----------------------------------------------------------------------------
  real(8), allocatable, private      :: t101(:), t102(:), t201(:), t202(:)
  real(8), allocatable, public, save :: t301(:), t302(:)
  real(8), allocatable, private      :: t401(:), t402(:), t501(:), t502(:)
  real(8), allocatable, private      :: t601(:), t602(:), t701(:), t702(:)
  real(8), allocatable, private      :: t801(:), t802(:), tt01(:), tt02(:)
#ifdef _OPENACC
  integer(4), private :: s_idx(1:2,1:16)=-1 ! FIXME hardcoding, 16 streams
#endif
  
  !- Other private module data -------------------------------------------------
  real(8),              private :: dx, dy, dt, dxdy, ddxdy, dtdx, dtdy, dtdxdy
  integer(4), save,     private :: n3dmax = 0

  public  :: tflow_int_simd, tflow_alloc_simd, tflow_ft_simd, tflow_up_ext
  private :: c_rin_rout, c_tu, c_tv, c_tw, c_tt, c_delta, c_dtx, c_cx_cy_cz,   &
             advection, c_dty, c_dtz, diffusion_vi, diffusion_hx, copy_t2tt,   &
             tflow_assign

contains
    
  ! ----------------------------------------------------------------------------

  subroutine tflow_assign ( n2d, n2dhalo, kh, mcol, idx, const, arr )
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: n2d, n2dhalo, kh(0:), mcol(0:)
    integer(4), intent(inout) :: idx(1:,1:)
    real(8),    intent(inout) :: arr(0:)
    real(8),    intent(in)    :: const

    integer(4) :: n, n2dl, n2du, ml, mu, kb

    include 'tflow_assign.inc'

    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    arr(n2dl:n2du) = zero

    do n=n2dl,n2du
      kb = kh(n)
      if (kb <= 1) cycle
      ml = mcol(n)
      mu = ml + kb - 2
      arr(ml:mu) = const
    enddo

!$OMP MASTER
    ! handle the halo:
    if (n2dhalo > 0) then
      n2dl = n2d+1
      n2du = n2d+n2dhalo
      arr(n2dl:n2du) = const
      do n=n2dl,n2du
        kb = kh(n)
        if (kb > 1) then
          ml = mcol(n)
          mu = ml + kb - 2
          arr(ml:mu) = const
        endif
      enddo
    endif
!$OMP END MASTER
  end subroutine tflow_assign

  ! ----------------------------------------------------------------------------

  subroutine tflow_alloc_simd (n3d,ierr) 

    implicit none

    integer(4), intent(in)    :: n3d
    integer(4), intent(inout) :: ierr

    if (n3d > n3dmax) then
      n3dmax = n3d  
    endif

    allocate (tt01(0:n3d), tt02(0:n3d),                                        &
              t101(0:n3d), t102(0:n3d),                                        &
              t201(0:n3d), t202(0:n3d),                                        &
              t301(0:n3d), t302(0:n3d),                                        &
              t401(0:n3d), t402(0:n3d),                                        &
              t501(0:n3d), t502(0:n3d),                                        &
              t601(0:n3d), t602(0:n3d),                                        &
              t701(0:n3d), t702(0:n3d),                                        &
              t801(0:n3d), t802(0:n3d),                                        &
              stat=ierr                      )
  end subroutine tflow_alloc_simd


  subroutine tflow_ft_simd (mcol, n2d, kh, idx, n3d)
    ! Let the domain with the largest Ir determine the numa first-touch
    ! of tflow tmp-arrays in the wet points of that domain on the present task.
    ! The remaining points (i.e. land point, halo points, points in the largest
    ! domain) of the arrays initialised with zeroes from the MASTER thread.

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, n3d, mcol(0:), kh(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_ft_simd.inc'

    !- local vars --------------------------------------------------------------
    integer(4) :: n, n2dl, n2du, kb, ml, mu

    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    ! NUMA first-touch in surface:
    t101 (n2dl:n2du) = zero
    t102 (n2dl:n2du) = zero
    t201 (n2dl:n2du) = zero
    t202 (n2dl:n2du) = zero
    t301 (n2dl:n2du) = zero
    t302 (n2dl:n2du) = zero
    t401 (n2dl:n2du) = zero
    t402 (n2dl:n2du) = zero
    t501 (n2dl:n2du) = zero
    t502 (n2dl:n2du) = zero
    t601 (n2dl:n2du) = zero
    t602 (n2dl:n2du) = zero
    t701 (n2dl:n2du) = zero
    t702 (n2dl:n2du) = zero
    t801 (n2dl:n2du) = zero
    t802 (n2dl:n2du) = zero
    tt01 (n2dl:n2du) = zero
    tt02 (n2dl:n2du) = zero

    ! NUMA first-touch in sub-surface:
    do n=n2dl,n2du
      kb = kh(n)
      if (kb <= 1) cycle
  
      ml = mcol(n)
      mu = ml + kb -2

      t101 (ml:mu) = zero
      t102 (ml:mu) = zero
      t201 (ml:mu) = zero
      t202 (ml:mu) = zero
      t301 (ml:mu) = zero
      t302 (ml:mu) = zero
      t401 (ml:mu) = zero
      t402 (ml:mu) = zero
      t501 (ml:mu) = zero
      t502 (ml:mu) = zero
      t601 (ml:mu) = zero
      t602 (ml:mu) = zero
      t701 (ml:mu) = zero
      t702 (ml:mu) = zero
      t801 (ml:mu) = zero
      t802 (ml:mu) = zero
      tt01 (ml:mu) = zero
      tt02 (ml:mu) = zero
    enddo

    ! make sure to have nice initial values in land point, 
    ! and possibly the halo:
!$OMP MASTER
    t101 (0) = zero
    t102 (0) = zero
    t201 (0) = zero
    t202 (0) = zero
    t301 (0) = zero
    t302 (0) = zero
    t401 (0) = zero
    t402 (0) = zero
    t501 (0) = zero
    t502 (0) = zero
    t601 (0) = zero
    t602 (0) = zero
    t701 (0) = zero
    t702 (0) = zero
    t801 (0) = zero
    t802 (0) = zero
    tt01 (0) = zero
    tt02 (0) = zero

    if (n3d < n3dmax) then
      ml = n3d + 1
      mu = n3dmax

      t101 (ml:mu) = zero
      t102 (ml:mu) = zero
      t201 (ml:mu) = zero
      t202 (ml:mu) = zero
      t301 (ml:mu) = zero
      t302 (ml:mu) = zero
      t401 (ml:mu) = zero
      t402 (ml:mu) = zero
      t501 (ml:mu) = zero
      t502 (ml:mu) = zero
      t601 (ml:mu) = zero
      t602 (ml:mu) = zero
      t701 (ml:mu) = zero
      t702 (ml:mu) = zero
      t801 (ml:mu) = zero
      t802 (ml:mu) = zero
      tt01 (ml:mu) = zero
      tt02 (ml:mu) = zero
    endif
!$OMP END MASTER

  end subroutine tflow_ft_simd

  !-----------------------------------------------------------------------------
    
  subroutine diffusion_hx(msrf,mcol,ind,n2d,kh,hn,hx,hy,eddyh,cosphi,idx,uvdam,&
                          t1d01,t1d02)
    ! NOTE: it uses the global module variables: t1d and tt
    ! NOTE: assumes eddyh and tt are zero on land (these were initialized to 0).

    !- modules -----------------------------------------------------------------
    use dmi_omp,   only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d
    integer(4), intent(in)    :: msrf(0:,0:), mcol(0:), ind(:,:), kh(0:)
    real(8),    intent(in)    :: hn(0:), hx(0:), hy(0:)
    real(8),    intent(in)    :: cosphi(:,0:), eddyh(0:), uvdam(:,0:)
    real(8),    intent(inout) :: t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'diffusion_hx.inc'

    !  simple vars:
    integer(4) :: n, k, kb, i, j, mi, me, mn, mw, ms
    integer(4) :: n2dl, n2du, mi0, me0, mw0, ms0, mn0, kmx
    real(8)    :: ahx, ahy, dh, dp, fx0, fy0, fx, fy, dahe, dahw, dahn, dahs
    real(8)    :: dtl2, ede, edw, eds, edn, lf, ude, udw, vdn, vds

    real(8), parameter :: p2 = 0.2_8, p2s = p2*p2

    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    !- Horizontal diffusion ----------------------------------------------------
    !   The effective horizontal eddy diffusivity, D_effective, is obtained from
    !   the horizontal eddy viscosity, Eh, by a relaxation according to:
    !       D_effective = 0.25*Eh/(0.1 + Eh*dt/L**2)
    !   where L is the length scale applied in the Smagorinsky model,
    !       L**2 = Csmag**2 * dx * dy,    Csmag = 0.2
    !   This yields a dimensionless diffusivity with maximum value approx 
    !       0.25 * Csmag**2 = 0.01
    !
    !   FIXME: Consider to make the magic numbers 0.1, 0.25 and 0.2 configurable
    !          i.e. through the argument list. At least, 0.2 should be exactly
    !          the same as the Csmag value use elsewhere in model.
    ! 

    !- some constants:
    dtl2 = dt/(p2s*dx*dy)
    ahx  = qrt*dt/(dx*dx)
    ahy  = qrt*dt/(dy*dy)

    ! k=1 un-rolled:
    do n=n2dl,n2du
      kb = kh(n)
      if (kb < 1) cycle
      i = ind(1,n)
      j = ind(2,n)

      dp  = one/cosphi(1,i)
      fx0 = ahx*dp*dp
      fy0 = ahy*dp
      lf  = dtl2*dp

      mi = n
      mn = msrf(i-1,j  )
      mw = msrf(i,  j-1)
      me = msrf(i,  j+1)
      ms = msrf(i+1,j  )

      ude = uvdam(1,mi)
      udw = uvdam(1,mw)
      vdn = uvdam(2,mn)
      vds = uvdam(2,mi)

      dh   = one/hn(mi)

      fx   = fx0*dh
      ede  = (eddyh(me)+eddyh(mi))*ude
      dahe = fx*hx(mi)*ede/(p2 + lf*ede)
      edw  = (eddyh(mw)+eddyh(mi))*udw
      dahw = fx*hx(mw)*edw/(p2 + lf*edw)

      fy   = fy0*dh
      eds  = (eddyh(ms)+eddyh(mi))*vds
      dahs = fy*hy(mi)*cosphi(2,i  )*eds/(p2 + lf*eds)
      edn  = (eddyh(mn)+eddyh(mi))*vdn
      dahn = fy*hy(mn)*cosphi(2,i-1)*edn/(p2 + lf*edn)

      t1d01(mi) = t1d01(mi)                                                    &
                + dahe*(tt01(me)-tt01(mi)) - dahw*(tt01(mi)-tt01(mw))          &
                + dahn*(tt01(mn)-tt01(mi)) - dahs*(tt01(mi)-tt01(ms)) 
      t1d02(mi) = t1d02(mi)                                                    &
                + dahe*(tt02(me)-tt02(mi)) - dahw*(tt02(mi)-tt02(mw))          &
                + dahn*(tt02(mn)-tt02(mi)) - dahs*(tt02(mi)-tt02(ms)) 
    enddo

    do n=n2dl,n2du
      kb = kh(n)
      if (kb <= 1) cycle
      i = ind(1,n)
      j = ind(2,n)

      dp  = one/cosphi(1,i)
      fx0 = ahx*dp*dp
      fy0 = ahy*dp
      lf  = dtl2*dp

      ude = uvdam(1,n)
      udw = uvdam(1,msrf(i,j-1))
      vdn = uvdam(2,msrf(i-1,j))
      vds = uvdam(2,n)

      mi0 = mcol(n) - 2
      me0 = mcol(msrf(i,  j+1)) - 2
      mw0 = mcol(msrf(i,  j-1)) - 2
      mn0 = mcol(msrf(i-1,j  )) - 2
      ms0 = mcol(msrf(i+1,j  )) - 2

      fy = fy0*cosphi(2,i)

      kmx = min(kb, kh(msrf(i,j+1)), kh(msrf(i,j-1)), kh(msrf(i-1,j)),         &
                kh(msrf(i+1,j)))

      ! mainloop: ld=9+2*5+1=20, ld+st=2, st=0, op=20+2*15=50, op/ldst=50/(20+4)
      do k=2,kmx
        mi = mi0 + k
        me = me0 + k
        mw = mw0 + k
        mn = mn0 + k
        ms = ms0 + k

        ede  = eddyh(me)+eddyh(mi)
        dahe = ude*hx(mi)*ede/(p2 + lf*ede)

        edw  = eddyh(mw)+eddyh(mi)
        dahw = udw*hx(mw)*edw/(p2 + lf*edw)

        edn  = eddyh(mn)+eddyh(mi)
        dahn = vdn*hy(mn)*edn/(p2 + lf*edn)

        eds  = eddyh(ms)+eddyh(mi)
        dahs = vds*hy(mi)*eds/(p2 + lf*eds)

        t1d01(mi) = t1d01(mi) + ( fx0*( dahe*(tt01(me)-tt01(mi))               &
                                       -dahw*(tt01(mi)-tt01(mw)))              &
                                 + fy*( dahn*(tt01(mn)-tt01(mi))               &
                                       -dahs*(tt01(mi)-tt01(ms))) )/hn(mi)
        t1d02(mi) = t1d02(mi) + ( fx0*( dahe*(tt02(me)-tt02(mi))               &
                                       -dahw*(tt02(mi)-tt02(mw)))              &
                                 + fy*( dahn*(tt02(mn)-tt02(mi))               &
                                       -dahs*(tt02(mi)-tt02(ms))) )/hn(mi)
      enddo

      ! remainder loops:
      if (ude > zero) then
        do k=max(2,kmx+1),min(kb,kh(msrf(i,j+1)))
          mi = mi0 + k
          me = me0 + k
          ede  = eddyh(me)+eddyh(mi)
          dahe = fx0*hx(mi)*ede/(p2 + lf*ede)/hn(mi)
          t1d01(mi) = t1d01(mi) + dahe*(tt01(me)-tt01(mi))
          t1d02(mi) = t1d02(mi) + dahe*(tt02(me)-tt02(mi))
        enddo
      endif
      if (udw > zero) then
        do k=max(2,kmx+1),min(kb,kh(msrf(i,j-1)))
          mi = mi0 + k
          mw = mw0 + k
          edw  = eddyh(mw)+eddyh(mi)
          dahw = fx0*hx(mw)*edw/(p2 + lf*edw)/hn(mi)
          t1d01(mi) = t1d01(mi) - dahw*(tt01(mi)-tt01(mw))
          t1d02(mi) = t1d02(mi) - dahw*(tt02(mi)-tt02(mw))
        enddo
      endif
      if (vdn > zero) then
        do k=max(2,kmx+1),min(kb,kh(msrf(i-1,j)))
          mi = mi0 + k
          mn = mn0 + k
          edn  = eddyh(mn)+eddyh(mi)
          dahn = fy*hy(mn)*edn/(p2 + lf*edn)/hn(mi)
          t1d01(mi) = t1d01(mi) + dahn*(tt01(mn)-tt01(mi))
          t1d02(mi) = t1d02(mi) + dahn*(tt02(mn)-tt02(mi))
        enddo
      endif
      if (vds > zero) then
        do k=max(2,kmx+1),min(kb,kh(msrf(i+1,j)))
          mi = mi0 + k
          ms = ms0 + k
          eds  = eddyh(ms)+eddyh(mi)
          dahs = fy*hy(mi)*eds/(p2 + lf*eds)/hn(mi)
          t1d01(mi) = t1d01(mi) - dahs*(tt01(mi)-tt01(ms))
          t1d02(mi) = t1d02(mi) - dahs*(tt02(mi)-tt02(ms))
        enddo
      endif
    enddo
    !- Horizontal diffusion end ------------------------------------------------

  end subroutine diffusion_hx

  !-----------------------------------------------------------------------------

  subroutine diffusion_vi(kmax,mcol,n2d,kh,h,hn,avt,avs,idx,t1d01,t1d02,q)
    !  inlined version over 2 components

    !- modules -----------------------------------------------------------------
    use constants, only : cpw, rhow, hnull
    use dmi_omp,   only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)           :: kmax, n2d, mcol(0:), kh(0:)
    real(8),    intent(in)           :: h(0:), avt(0:), avs(0:), hn(0:)
    real(8),    intent(inout)        :: t1d01(0:), t1d02(0:)
    real(8),    intent(in), optional :: q(0:)
    integer(4), intent(inout)        :: idx(1:,1:)

    include 'diffusion_vi.inc'

    !  automatic arrays:
    real(8), dimension(kmax+1) :: qa01, qa02, qb01, qb02, qc01, qc02,          &
                                  dtiav01, dtiav02, hm,                        &
                                  pgka01, pgka02, pgkb01, pgkb02
    
    !  simple vars:
    integer(4) :: n, k, kb, no, nu, kbp1, kp1, kbm2, ki, n2dl, n2du
    real(8)    :: dti2, fac, qfac, cn01, cn02

    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    !- some preparations:
    dti2 = two*dt
    if (present(q)) then
      qfac = dt/(cpw*rhow)
      t1d01(n2dl:n2du) = t1d01(n2dl:n2du) + qfac*q(n2dl:n2du)/h(n2dl:n2du)
    endif

    !- Vertical diffusion ------------------------------------------------------
    do n=n2dl,n2du
      if (h(n) < hnull) cycle

      kb=kh(n)
      if (kb <= 1) cycle
    
      !- Implicit part:
      hm(1)      = hn(n)
      dtiav01(1) = zero
      dtiav02(1) = zero

      no = mcol(n) - 2

      do k=2,kb
        nu = no + k
        hm(k) = hn(nu)
        fac   = dti2/(hm(k-1)+hm(k))
        dtiav01(k) = fac*avt(nu)
        dtiav02(k) = fac*avs(nu)
      enddo
      kbp1 = kb+1
      dtiav01(kbp1) = zero
      dtiav02(kbp1) = zero
      do k=1,kb
        kp1=k+1
        qa01(k) = dtiav01(k  )/hm(k)
        qc01(k) = dtiav01(kp1)/hm(k)
        qb01(k) = one + qa01(k) + qc01(k)
        qa02(k) = dtiav02(k  )/hm(k)
        qc02(k) = dtiav02(kp1)/hm(k)
        qb02(k) = one + qa02(k) + qc02(k)
      enddo
      pgka01(1)=zero
      pgkb01(1)=zero
      pgka02(1)=zero
      pgkb02(1)=zero
      k   = 1
      kp1 = k+1
      cn01        = one/(qb01(k) - qa01(k)*pgka01(k))
      pgka01(kp1) = qc01(k)*cn01
      pgkb01(kp1) = (qa01(k)*pgkb01(k) + t1d01(n))*cn01
      cn02        = one/(qb02(k) - qa02(k)*pgka02(k))
      pgka02(kp1) = qc02(k)*cn02
      pgkb02(kp1) = (qa02(k)*pgkb02(k) + t1d02(n))*cn02

      do k=2,kb
        kp1  = k+1
        cn01        = one/(qb01(k) - qa01(k)*pgka01(k))
        pgka01(kp1) = qc01(k)*cn01
        pgkb01(kp1) = (qa01(k)*pgkb01(k) + t1d01(no+k))*cn01
        cn02        = one/(qb02(k) - qa02(k)*pgka02(k))
        pgka02(kp1) = qc02(k)*cn02
        pgkb02(kp1) = (qa02(k)*pgkb02(k) + t1d02(no+k))*cn02
      enddo

      t1d01(no+kb) = pgkb01(kbp1)
      t1d02(no+kb) = pgkb02(kbp1)
      kbm2 = kb-2
      do ki=1,kbm2
        k   = kb-ki
        kp1 = k+1
        t1d01(no+k) = pgka01(kp1)*t1d01(no+kp1) + pgkb01(kp1)
        t1d02(no+k) = pgka02(kp1)*t1d02(no+kp1) + pgkb02(kp1)
      enddo
      k   = 1
      kp1 = k+1
      t1d01(n) = pgka01(kp1)*t1d01(no+kp1) + pgkb01(kp1)
      t1d02(n) = pgka02(kp1)*t1d02(no+kp1) + pgkb02(kp1)
    enddo
    
    !- Vertical diffusion ended ------------------------------------------------

  end subroutine diffusion_vi

  !-----------------------------------------------------------------------------
    
  subroutine copy_t2tt (mcol,n2d,n2dhalo,kh,idx,t1d01,t1d02)
    ! uses t1d, modifies tt

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, n2dhalo, mcol(0:), kh(0:)
    integer(4), intent(inout) :: idx(1:,1:)
    real(8),    intent(in)    :: t1d01(0:), t1d02(0:)

    include 'copy_t2tt.inc'

    !- local vars --------------------------------------------------------------
    integer(4) :: n, n2dl, n2du, kb, ml, mu

    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    tt01(n2dl:n2du) = t1d01(n2dl:n2du)
    tt02(n2dl:n2du) = t1d02(n2dl:n2du)
    do n=n2dl,n2du
      kb = kh(n)
      if (kb > 1) then
        ml = mcol(n)
        mu = ml + kb - 2
        tt01(ml:mu) = t1d01(ml:mu)
        tt02(ml:mu) = t1d02(ml:mu)
      endif
    enddo

!$OMP MASTER
    ! handle the halo:
    if (n2dhalo > 0) then
      n2dl = n2d+1
      n2du = n2d+n2dhalo
      tt01(n2dl:n2du) = t1d01(n2dl:n2du)
      tt02(n2dl:n2du) = t1d02(n2dl:n2du)
      do n=n2dl,n2du
        kb = kh(n)
        if (kb > 1) then
          ml = mcol(n)
          mu = ml + kb - 2
          tt01(ml:mu) = t1d01(ml:mu)
          tt02(ml:mu) = t1d02(ml:mu)
        endif
      enddo
    endif
!$OMP END MASTER
  end subroutine copy_t2tt

  !-----------------------------------------------------------------------------

  subroutine c_delta (msrf,mcol,ind,n2d,kh,khu,khv,h,                          &
                      nbpz,krz,nbpu,kru,nbpv,krv,idx,t1d01,t1d02)

    ! NOTE: it modifies the global module variables: t4, t5 ,t6 

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d,nbpz,nbpu,nbpv
    integer(4), intent(in)    :: msrf(0:,0:), mcol(0:), ind(:,:)
    real(8),    intent(in)    :: h(0:), t1d01(0:), t1d02(0:)
    integer(4), intent(in)    :: krz(:,0:),kru(:,0:),krv(:,0:)
    integer(4), intent(in)    :: khu(0:),khv(0:),kh(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_delta.inc'

    !  simple vars:
    real(8)    :: fac, qqa, qqb
    integer(4) :: n,ntyp,k,kb,kbuw,kbue,kbvn,kbvs,i,j,mi,mm,mp,ss,ee,ww
    integer(4) :: mi0, mm0, mp0, nn, mw0, me0, ms0, mn0, me2, ms2, e2, s2
    integer(4) :: n2dl, n2du,nbpzl, nbpzu,nbpul, nbpuu,nbpvl, nbpvu
#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,kru,msrf,mcol,khu,krv,khv,h,kh,krz)                        &
!$ACC   pcreate(t1d01,t1d02,t401,t402,t501,t502,t601,t602) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,mi,n,i,j,mm,mp,kb,mm0,mp0,mi0) 
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif
    ! compute t4  --------------------------------------------------------------
    do n = n2dl, n2du
      kb=kh(n)
      i=ind(1,n)
      j=ind(2,n)

      ! unroll k=1
      mi = n
      mm = msrf(i,j-1)
      mp = msrf(i,j+1)
      if (mm > 0 .and. mp > 0) then
        t401(mi) = half*(t1d01(mp)-t1d01(mm))
        t402(mi) = half*(t1d02(mp)-t1d02(mm))
      elseif (mm > 0) then
        t401(mi) = t1d01(mi)-t1d01(mm)
        t402(mi) = t1d02(mi)-t1d02(mm)
      elseif (mp > 0) then
        t401(mi) = t1d01(mp)-t1d01(mi)
        t402(mi) = t1d02(mp)-t1d02(mi)
      else
        t401(mi) = zero
        t402(mi) = zero
      endif

      if (kb < 2) cycle

      mi0 = mcol(n) - 2
      mm0 = mcol(msrf(i,j-1)) - 2
      mp0 = mcol(msrf(i,j+1)) - 2
      do k=2,min(kb,kh(msrf(i,j-1)),kh(msrf(i,j+1)))
        t401(mi0+k) = half*(t1d01(mp0+k)-t1d01(mm0+k))
        t402(mi0+k) = half*(t1d02(mp0+k)-t1d02(mm0+k))
      enddo
      if (kh(msrf(i,j-1)) < kh(msrf(i,j+1))) then
        do k=max(1,kh(msrf(i,j-1)))+1,min(kb,kh(msrf(i,j+1)))
          t401(mi0+k) = t1d01(mp0+k)-t1d01(mi0+k)
          t402(mi0+k) = t1d02(mp0+k)-t1d02(mi0+k)
        enddo
      else
        do k=max(1,kh(msrf(i,j+1)))+1,min(kb,kh(msrf(i,j-1)))
          t401(mi0+k) = t1d01(mi0+k)-t1d01(mm0+k)
          t402(mi0+k) = t1d02(mi0+k)-t1d02(mm0+k)
        enddo
      endif
      do k=max(1,kh(msrf(i,j-1)),kh(msrf(i,j+1)))+1,kb
        t401(mi0+k) = zero
        t402(mi0+k) = zero
      enddo
    enddo
    
    ! compute t5  --------------------------------------------------------------
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (i,j,k,kb,mi,mi0,mm,mm0,mp,mp0,n)                       
#endif
    do n=n2dl, n2du
      kb=kh(n)
      i=ind(1,n)
      j=ind(2,n)

      ! unroll k=1
      mi = n
      mm = msrf(i-1,j)
      mp = msrf(i+1,j)
      if (mm > 0 .and. mp > 0) then
        t501(mi) = half*(t1d01(mm)-t1d01(mp))
        t502(mi) = half*(t1d02(mm)-t1d02(mp))
      elseif (mm > 0) then
        t501(mi) = t1d01(mm)-t1d01(mi)
        t502(mi) = t1d02(mm)-t1d02(mi)
      elseif (mp > 0) then
        t501(mi) = t1d01(mi)-t1d01(mp)
        t502(mi) = t1d02(mi)-t1d02(mp)
      else
        t501(mi) = zero
        t502(mi) = zero
      endif

      if (kb < 2) cycle

      mi0 = mcol(n) - 2
      mm0 = mcol(msrf(i-1,j)) - 2
      mp0 = mcol(msrf(i+1,j)) - 2
      do k=2,min(kb,kh(msrf(i-1,j)),kh(msrf(i+1,j)))
        t501(mi0+k) = half*(t1d01(mm0+k)-t1d01(mp0+k))
        t502(mi0+k) = half*(t1d02(mm0+k)-t1d02(mp0+k))
      enddo
      if (kh(msrf(i-1,j)) < kh(msrf(i+1,j))) then
        do k=max(1,kh(msrf(i-1,j)))+1,min(kb,kh(msrf(i+1,j)))
          t501(mi0+k) = t1d01(mi0+k)-t1d01(mp0+k)
          t502(mi0+k) = t1d02(mi0+k)-t1d02(mp0+k)
        enddo
      else
        do k=max(1,kh(msrf(i+1,j)))+1,min(kb,kh(msrf(i-1,j)))
          t501(mi0+k) = t1d01(mm0+k)-t1d01(mi0+k)
          t502(mi0+k) = t1d02(mm0+k)-t1d02(mi0+k)
        enddo
      endif
      do k=max(1,kh(msrf(i-1,j)),kh(msrf(i+1,j)))+1,kb
        t501(mi0+k) = zero
        t502(mi0+k) = zero
      enddo
    enddo
    

    ! compute t6  --------------------------------------------------------------
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private(fac,i,j,k,kb,mi,mi0,mm,mp,n,qqa,qqb)                   
#endif
    do n=n2dl, n2du
      kb = kh(n)

      ! un-roll k=1:
      if (kb == 1) then
        t601(n) = zero
        t602(n) = zero
        cycle
      endif
      mp = mcol(n)
      t601(n) = t1d01(n)-t1d01(mp)
      t602(n) = t1d02(n)-t1d02(mp)

      if (kb > 2) then
        ! un-roll k=2:
        mm = n
        mi = mcol(n)
        mp = mp + 1

        fac = one/(h(mm)+two*h(mi)+h(mp))
        qqa = (h(mm)+h(mi))*fac
        qqb = (h(mi)+h(mp))*fac

        t601(mi) = (one-qqa)*(t1d01(mm)-t1d01(mi))                             &
                 + (one-qqb)*(t1d01(mi)-t1d01(mp))
        t602(mi) = (one-qqa)*(t1d02(mm)-t1d02(mi))                             &
                 + (one-qqb)*(t1d02(mi)-t1d02(mp))
      endif

      mi0 = mcol(n) - 2

      ! un-roll k=kb:
      mi = mi0 + kb
      if (kb == 2) then
        mm = n
      else
        mm = mi - 1
      endif
      t601(mi) = t1d01(mm)-t1d01(mi)
      t602(mi) = t1d02(mm)-t1d02(mi)

      ! loop except k=1, k=2 and k=kb:
      do k=3,kb-1
        mi = mi0 + k
        mp = mi + 1
        mm = mi - 1

        fac = one/(h(mm)+two*h(mi)+h(mp))
        qqa = (h(mm)+h(mi))*fac
        qqb = (h(mi)+h(mp))*fac

        t601(mi) = (one-qqa)*(t1d01(mm)-t1d01(mi))                             &
                 + (one-qqb)*(t1d01(mi)-t1d01(mp))
        t602(mi) = (one-qqa)*(t1d02(mm)-t1d02(mi))                             &
                 + (one-qqb)*(t1d02(mi)-t1d02(mp))
      enddo
    enddo
#ifdef _OPENACC
    enddo ! streams
!$ACC WAIT
#endif

    
    ! compute boundary contributions to t4  and t5   ---------------------------
    ! - do this as the last actions in this routine to gather the barriers
    !   and push these as far as posible.
    if (nbpz > 0 .or. nbpu > 0 .or. nbpv > 0) then
      ! due to barrier(s)

      if (nbpz > 0) then
!$OMP BARRIER
! barrier since n-loop has changed
! MPI: NO need for halo swap
        call domp_get_domain(1, nbpz, nbpzl, nbpzu)
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (ee,i,j,k,kbue,kbuw,kbvn,kbvs,me0,mi,mi0,mn0,ms0,mw0,n,        &
!$ACC            nn,ntyp,ss,ww)                                         
#endif
        do n=nbpzl,nbpzu
          i = krz(1,n)
          j = krz(2,n)
          mi  = msrf(i,j)
          mi0 = mcol(mi) - 2
          ntyp = abs(krz(3,n))
          kbuw = khu(msrf(i,  j-1))
          kbue = khu(msrf(i,  j  ))
          kbvn = khv(msrf(i-1,j  ))
          kbvs = khv(msrf(i,  j  ))
          if     (ntyp == 3 .and. kbuw > 0) then
            ww = msrf(i,j-1)
            t401(mi) = t1d01(mi)-t1d01(ww)
            t402(mi) = t1d02(mi)-t1d02(ww)
            mw0 = mcol(ww) - 2
            do k=2,kbuw
              mi = mi0 + k
              ww = mw0 + k
              t401(mi) = t1d01(mi)-t1d01(ww)
              t402(mi) = t1d02(mi)-t1d02(ww)
            enddo
          elseif (ntyp == 1 .and. kbue > 0) then
            ee = msrf(i,j+1)
            t401(mi) = t1d01(ee)-t1d01(mi)
            t402(mi) = t1d02(ee)-t1d02(mi)
            me0 = mcol(ee) - 2
            do k=2,kbue
              mi = mi0 + k
              ee = me0 + k
              t401(mi) = t1d01(ee)-t1d01(mi)
              t402(mi) = t1d02(ee)-t1d02(mi)
            enddo
          elseif (ntyp == 4 .and. kbvn > 0) then
            nn = msrf(i-1,j)
            t501(mi) = t1d01(nn)-t1d01(mi)
            t502(mi) = t1d02(nn)-t1d02(mi)
            mn0 = mcol(nn) - 2
            do k=2,kbvn
              mi = mi0 + k
              nn = mn0 + k
              t501(mi) = t1d01(nn)-t1d01(mi)
              t502(mi) = t1d02(nn)-t1d02(mi)
            enddo
          elseif (ntyp == 2  .and. kbvs > 0) then
            ss = msrf(i+1,j)
            t501(mi) = t1d01(mi)-t1d01(ss)
            t502(mi) = t1d02(mi)-t1d02(ss)
            ms0 = mcol(ss) - 2
            do k=2,kbvs
              mi = mi0 + k
              ss = ms0 + k
              t501(mi) = t1d01(mi)-t1d01(ss)
              t502(mi) = t1d02(mi)-t1d02(ss)
            enddo
          endif
        enddo
      endif
    
      if (nbpu > 0) then
!$OMP BARRIER
! barrier since n-loop has changed
! MPI: NO need for halo swap
        call domp_get_domain(1, nbpu, nbpul, nbpuu)
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (e2,ee,i,j,k,kbue,kbuw,me0,me2,mi,mi0,mw0,n,ntyp,ww)    
#endif
        do n=nbpul,nbpuu
          i = kru(1,n)
          j = kru(2,n)
          mi  = msrf(i,j)
          mi0 = mcol(mi) - 2
          ntyp = kru(3,n)
          kbuw = khu(msrf(i,j-1))
          kbue = khu(msrf(i,j+1))
          if     (ntyp == 3 .and. kbuw > 0) then
            ww = msrf(i,j-1)
            t401(mi) = t1d01(mi)-t1d01(ww)
            t402(mi) = t1d02(mi)-t1d02(ww)
            mw0 = mcol(ww) - 2
            do k=2,kbuw
              mi = mi0 + k
              ww = mw0 + k
              t401(mi) = t1d01(mi)-t1d01(ww)
              t402(mi) = t1d02(mi)-t1d02(ww)
            enddo
          elseif (ntyp == 1 .and. kbue > 0) then
            ee = msrf(i,j+1)
            e2 = msrf(i,j+2)
            if (e2 > 0) then
              t401(ee) = t1d01(e2)-t1d01(ee)
              t402(ee) = t1d02(e2)-t1d02(ee)
            else
              t401(ee) = -t1d01(ee)
              t402(ee) = -t1d02(ee)
            endif
            me0 = mcol(ee) - 2
            me2 = mcol(e2) - 2
            do k=2,min(kbue,kh(msrf(i,j+2)))
              ee = me0 + k
              e2 = me2 + k
              t401(ee) = t1d01(e2)-t1d01(ee)
              t402(ee) = t1d02(e2)-t1d02(ee)
            enddo
            do k=min(kbue,max(kh(msrf(i,j+2)),1))+1,kbue
              ee = me0 + k
              t401(ee) = -t1d01(ee)
              t402(ee) = -t1d02(ee)
            enddo
          endif
        enddo
      endif

      if (nbpv > 0) then
!$OMP BARRIER
! barrier since n-loop has changed
! MPI: NO need for halo swap
        call domp_get_domain(1, nbpv, nbpvl, nbpvu)
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (nn,s2,ss,i,j,k,kbvn,kbvs,ms0,ms2,mi,mi0,mw0,n,ntyp,ww)    
#endif
        do n=nbpvl,nbpvu
          i = krv(1,n)
          j = krv(2,n)
          mi  = msrf(i,j)
          mi0 = mcol(mi) - 2
          ntyp = krv(3,n)
          kbvn = khv(msrf(i-1,j))
          kbvs = khv(msrf(i+1,j))
          if (ntyp == 4 .and. kbvn > 0) then
            nn = msrf(i-1,j)
            t501(mi) = t1d01(nn)-t1d01(mi)
            t502(mi) = t1d02(nn)-t1d02(mi)
            mn0 = mcol(nn) - 2
            do k=2,kbvn
              mi = mi0 + k
              nn = mn0 + k
              t501(mi) = t1d01(nn)-t1d01(mi)
              t502(mi) = t1d02(nn)-t1d02(mi)
            enddo
          elseif (ntyp == 2 .and. kbvs > 0) then
            ss = msrf(i+1,j) - 2
            s2 = msrf(i+2,j) - 2
            if (s2 > 0) then
              t501(ss) = t1d01(ss)-t1d01(s2)
              t502(ss) = t1d02(ss)-t1d02(s2)
            else
              t501(ss) = t1d01(ss)
              t502(ss) = t1d02(ss)
            endif
            ms0 = mcol(ss) - 2
            ms2 = mcol(s2) - 2
            do k=2,min(kbvs,kh(msrf(i+2,j)))
              ss = ms0 + k
              s2 = ms2 + k
              t501(ss) = t1d01(ss)-t1d01(s2)
              t502(ss) = t1d02(ss)-t1d02(s2)
            enddo
            do k=min(kbvs,max(kh(msrf(i+2,j)),1))+1,kbvs
              ss = ms0 + k
              t501(ss) = t1d01(ss)
              t502(ss) = t1d02(ss)
            enddo
          endif
        enddo
      endif
    endif
#ifdef _OPENACC
!$ACC END DATA
#endif

  end subroutine c_delta
    
  subroutine c_tt (kmax,msrf,mcol,ind,n2d,kh,khu,khv,hx,hy,h_old,h_new,u,v,w,  &
                   nbpz,krz,nbpq,krq,rwqk,rwzkt,rwqkt,cosphi,idx,t1d01,t1d02)
    
    ! NOTE: it modifies the global module variables: tt and uses t1, t2, t3

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain
    
    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: kmax, n2d
    integer(4), intent(in)    :: msrf(0:,0:), mcol(0:), ind(:,:)
    integer(4), intent(in)    :: kh(0:), khu(0:), khv(0:)
    integer(4), intent(in)    :: nbpz, nbpq
    integer(4), intent(in)    :: krq(:,0:), krz(:,0:)
    real(8),    intent(in)    :: hx(0:), hy(0:), h_old(0:), h_new(0:)
    real(8),    intent(in)    :: u(0:), v(0:), w(0:)
    real(8),    intent(in)    :: rwqk(0:,:), rwzkt(:,:,0:)
    real(8),    intent(in)    :: rwqkt(:,:,0:), cosphi(:,0:)
    real(8),    intent(in)    :: t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)
    
    include 'tflow_tt.inc'

    !- local vars --------------------------------------------------------------
    !  automatic:
    real(8) :: rwzkin(2:max(2,kmax)), rwzkout(2:max(2,kmax)), rwzkin1, rwzkout1
    
    !  simple vars:
    integer(4) :: nuw, nue, nvn, nvs, kbuw, kbue, kbvn, kbvs, kbm
    integer(4) :: n, i, j, k, kb
    integer(4) :: mi, mn, mw, me, ms, md, mu, mi0, mn0, mw0, me0, ms0
    integer(4) :: n2dl, n2du, nbpql, nbpqu, nbpzl, nbpzu
    real(8)    :: fac0, fac, facx, facy1, facy2, factxy, uuu, vvv, wh
    real(8)    :: ws, ddx, ddy, ddt, dtf0, dtdx1, dtdx2
    real(8)    :: u1, u2, u3, u4, v1, v2, v3, v4, w1, w2, w3, w4
    real(8)    :: humi, humw, hvmi, hvmn, tmp1, tmp2
#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(kh,ind,cosphi,msrf,h_old,u,hx,v,hy,w,khv,khu,krz,h_new,rwzkt,  &
!$ACC           krq,rwqk,rwqkt,mcol)                                           &
!$ACC   pcreate(t101,t102,t1d01,t1d02,t201,t202,                               &
!$ACC           t301,t302,tt01,tt02) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,n,tmp1,tmp2,u1,u2,u3,u4,v1,v2,v3,v4,w2,kb,i,fac0,j,         &
!$ACC            fac,humi,humw,dtdx2,hvmn,dtdx1,ms,hvmi,w3,dtf0,kbm,           &
!$ACC            mw0,me0,mn0,ms0,w1,MI0,md,me,mi,mn,mu,mw,w4)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif
    do n=n2dl,n2du
      kb=kh(n)
      if (kb < 1) cycle
      i=ind(1,n)
      j=ind(2,n)
      fac0  = dxdy*cosphi(1,i)
      dtf0  = dt*fac0
      dtdx1 = dtdx*cosphi(2,i)
      dtdx2 = dtdx*cosphi(2,i-1)

      ! k=1 unrolled:
      mi = n
      mw = msrf(i,j-1)
      mn = msrf(i-1,j)
      me = msrf(i,j+1)
      ms = msrf(i+1,j)

      fac = fac0*h_old(mi)
      tmp1 = t1d01(mi)*fac
      tmp2 = t1d02(mi)*fac
      u1 = u(mi)
      if (u1 > zero) then
        humi = dtdy*hx(mi)*u1
        tmp1 = tmp1 - humi*t101(mi)
        tmp2 = tmp2 - humi*t102(mi)
      elseif (u1 < zero) then
        humi = dtdy*hx(mi)*u1
        tmp1 = tmp1 - humi*t101(me)
        tmp2 = tmp2 - humi*t102(me)
      endif
      u3 = u(mw)
      if (u3 > zero) then
        humw = dtdy*hx(mw)*u3
        tmp1 = tmp1 + humw*t101(mw)
        tmp2 = tmp2 + humw*t102(mw)
      elseif (u3 < zero) then
        humw = dtdy*hx(mw)*u3
        tmp1 = tmp1 + humw*t101(mi)
        tmp2 = tmp2 + humw*t102(mi)
      endif
      v1 = v(mn)
      if (v1 > zero) then
        hvmn = dtdx2*hy(mn)*v1
        tmp1 = tmp1 - hvmn*t201(mi)
        tmp2 = tmp2 - hvmn*t202(mi)
      elseif (v1 < zero) then
        hvmn = dtdx2*hy(mn)*v1
        tmp1 = tmp1 - hvmn*t201(mn)
        tmp2 = tmp2 - hvmn*t202(mn)
      endif
      v3 = v(mi)
      if (v3 > zero) then
        hvmi = dtdx1*hy(mi)*v3
        tmp1 = tmp1 + hvmi*t201(ms)
        tmp2 = tmp2 + hvmi*t202(ms)
      elseif (v3 < zero) then
        hvmi = dtdx1*hy(mi)*v3
        tmp1 = tmp1 + hvmi*t201(mi)
        tmp2 = tmp2 + hvmi*t202(mi)
      endif

      tt01(mi) = tmp1
      tt02(mi) = tmp2

      if (kb > 1) then
        md = mcol(n)
        w3 = w(md)
        if (w3 > zero) then
          w4 = dtf0*w3
          tt01(mi) = tt01(mi) + w4*t301(md)
          tt02(mi) = tt02(mi) + w4*t302(md)
        elseif (w3 < zero) then
          w4 = dtf0*w3
          tt01(mi) = tt01(mi) + w4*t301(mi)
          tt02(mi) = tt02(mi) + w4*t302(mi)
        endif

        mi0 = mcol(n) - 2
        mn0 = mcol(msrf(i-1,j  )) - 2
        mw0 = mcol(msrf(i,  j-1)) - 2
        ms0 = mcol(msrf(i+1,j  )) - 2
        me0 = mcol(msrf(i,  j+1)) - 2
      else
        cycle
      endif

      kbm = min(kb-1,khu(msrf(i,j-1)),khv(msrf(i-1,j)),khu(n),khv(n))

      ! unroll k=2:
      k = 2
      if (k < kb) then
        mi = mi0 + k
        md = mi + 1
        mu = n

        fac = fac0*h_old(mi)
        tmp1 = t1d01(mi)*fac
        tmp2 = t1d02(mi)*fac

        if (khu(n) >= k) then 
          u1 = u(mi)
          if (u1 > zero) then
            humi = dtdy*hx(mi)*u1
            tmp1 = tmp1 - humi*t101(mi)
            tmp2 = tmp2 - humi*t102(mi)
          elseif (u1 < zero) then
            me = me0 + k
            humi = dtdy*hx(mi)*u1
            tmp1 = tmp1 - humi*t101(me)
            tmp2 = tmp2 - humi*t102(me)
          endif
        endif
        if (khu(msrf(i,j-1)) >= k) then 
          mw = mw0 + k
          u3 = u(mw)
          if (u3 > zero) then
            humw = dtdy*hx(mw)*u3
            tmp1 = tmp1 + humw*t101(mw)
            tmp2 = tmp2 + humw*t102(mw)
          elseif (u3 < zero) then
            humw = dtdy*hx(mw)*u3
            tmp1 = tmp1 + humw*t101(mi)
            tmp2 = tmp2 + humw*t102(mi)
          endif
        endif
        if (khv(msrf(i-1,j)) >= k) then 
          mn = mn0 + k
          v1 = v(mn)
          if (v1 > zero) then
            hvmn = dtdx2*hy(mn)*v1
            tmp1 = tmp1 - hvmn*t201(mi)
            tmp2 = tmp2 - hvmn*t202(mi)
          elseif (v1 < zero) then
            hvmn = dtdx2*hy(mn)*v1
            tmp1 = tmp1 - hvmn*t201(mn)
            tmp2 = tmp2 - hvmn*t202(mn)
          endif
        endif
        if (khv(n) >= k) then 
          v3 = v(mi)
          if (v3 > zero) then
            ms = ms0 + k
            hvmi = dtdx1*hy(mi)*v3
            tmp1 = tmp1 + hvmi*t201(ms)
            tmp2 = tmp2 + hvmi*t202(ms)
          elseif (v3 < zero) then
            hvmi = dtdx1*hy(mi)*v3
            tmp1 = tmp1 + hvmi*t201(mi)
            tmp2 = tmp2 + hvmi*t202(mi)
          endif
        endif
        w1 = w(mi)
        if (w1 > zero) then
          w2   = dtf0*w1
          tmp1 = tmp1 - w2*t301(mi)
          tmp2 = tmp2 - w2*t302(mi)
        elseif (w1 < zero) then
          w2   = dtf0*w1
          tmp1 = tmp1 - w2*t301(mu)
          tmp2 = tmp2 - w2*t302(mu)
        endif
        w3 = w(md)
        if (w3 > zero) then
          w4   = dtf0*w3
          tmp1 = tmp1 + w4*t301(md)
          tmp2 = tmp2 + w4*t302(md)
        elseif (w3 < zero) then
          w4   = dtf0*w3
          tmp1 = tmp1 + w4*t301(mi)
          tmp2 = tmp2 + w4*t302(mi)
        endif
        tt01(mi) = tmp1
        tt02(mi) = tmp2
      endif

      ! ld/st=11+2*(2+10)=36, flop=21+2*23=67, I=67/36 ~ 1.9
      do k=3,kbm
        mi = mi0 + k
        mw = mw0 + k
        me = me0 + k
        mn = mn0 + k
        ms = ms0 + k
        md = mi + 1
        mu = mi - 1

        fac  = fac0*h_old(mi)
        humi = hx(mi)*u(mi)*dtdy
        humw = hx(mw)*u(mw)*dtdy
        hvmn = hy(mn)*v(mn)*dtdx2
        hvmi = hy(mi)*v(mi)*dtdx1
        u1 = max(humi,zero)
        u2 = min(humi,zero)
        u3 = max(humw,zero)
        u4 = min(humw,zero)
        v1 = max(hvmn,zero)
        v2 = min(hvmn,zero)
        v3 = max(hvmi,zero)
        v4 = min(hvmi,zero)
        w1 = max(w(mi),zero)
        w2 = min(w(mi),zero)
        w3 = max(w(md),zero)
        w4 = min(w(md),zero)

        tt01(mi) = t1d01(mi)*fac                                               &
                 -        u2*t101(me) - (u1-u4)*t101(mi) + u3*t101(mw)         &
                 -        v2*t201(mn) - (v1-v4)*t201(mi) + v3*t201(ms)         &
                 - dtf0 *(w2*t301(mu) + (w1-w4)*t301(mi) - w3*t301(md))
        tt02(mi) = t1d02(mi)*fac                                               &
                 -        u2*t102(me) - (u1-u4)*t102(mi) + u3*t102(mw)         &
                 -        v2*t202(mn) - (v1-v4)*t202(mi) + v3*t202(ms)         &
                 - dtf0 *(w2*t302(mu) + (w1-w4)*t302(mi) - w3*t302(md))
      enddo

      ! remainder loops:
      do k=max(3,max(kbm,1)+1),kb-1
        mi = mi0 + k
        md = mi + 1
        mu = mi - 1

        fac  = fac0*h_old(mi)
        w1 = max(w(mi),zero)
        w2 = min(w(mi),zero)
        w3 = max(w(md),zero)
        w4 = min(w(md),zero)

        tt01(mi) = t1d01(mi)*fac                                               &
                 - dtf0 *(w2*t301(mu) + (w1-w4)*t301(mi) - w3*t301(md))
        tt02(mi) = t1d02(mi)*fac                                               &
                 - dtf0 *(w2*t302(mu) + (w1-w4)*t302(mi) - w3*t302(md))
      enddo
      do k=max(3,max(kbm,1)+1),min(kb-1,khu(n))
        mi = mi0 + k
        me = me0 + k

        humi = hx(mi)*u(mi)*dtdy
        u1 = max(humi,zero)
        u2 = min(humi,zero)

        tt01(mi) = tt01(mi) - u2*t101(me) - u1*t101(mi)
        tt02(mi) = tt02(mi) - u2*t102(me) - u1*t102(mi)
      enddo
      do k=max(3,max(kbm,1)+1),min(kb-1,khu(msrf(i,j-1)))
        mi = mi0 + k
        mw = mw0 + k

        humw = hx(mw)*u(mw)*dtdy
        u3 = max(humw,zero)
        u4 = min(humw,zero)

        tt01(mi) = tt01(mi) + u4*t101(mi) + u3*t101(mw)
        tt02(mi) = tt02(mi) + u4*t102(mi) + u3*t102(mw)
      enddo
      do k=max(3,max(kbm,1)+1),min(kb-1,khv(n))
        mi = mi0 + k
        ms = ms0 + k

        hvmi = hy(mi)*v(mi)*dtdx1
        v3 = max(hvmi,zero)
        v4 = min(hvmi,zero)

        tt01(mi) = tt01(mi) + v4*t201(mi) + v3*t201(ms)
        tt02(mi) = tt02(mi) + v4*t202(mi) + v3*t202(ms)
      enddo
      do k=max(3,max(kbm,1)+1),min(kb-1,khv(msrf(i-1,j)))
        mi = mi0 + k
        mn = mn0 + k

        hvmn = hy(mn)*v(mn)*dtdx2
        v1 = max(hvmn,zero)
        v2 = min(hvmn,zero)

        tt01(mi) = tt01(mi) - v2*t201(mn) - v1*t201(mi)
        tt02(mi) = tt02(mi) - v2*t202(mn) - v1*t202(mi)
      enddo

      ! k=kb unrolled
      mi = mi0 + kb
      if (kb > 2) then
        mu = mi - 1
      else
        mu = n
      endif

      fac  = fac0*h_old(mi)
      tmp1 = t1d01(mi)*fac
      tmp2 = t1d02(mi)*fac
      if (khu(n) >= kb) then 
        u1 = u(mi)
        if (u1 > zero) then
          humi = dtdy*hx(mi)*u1
          tmp1 = tmp1 - humi*t101(mi)
          tmp2 = tmp2 - humi*t102(mi)
        elseif (u1 < zero) then
          me = me0 + kb
          humi = dtdy*hx(mi)*u1
          tmp1 = tmp1 - humi*t101(me)
          tmp2 = tmp2 - humi*t102(me)
        endif
      endif
      if (khu(msrf(i,j-1)) >= kb) then 
        mw = mw0 + kb
        u3 = u(mw)
        if (u3 > zero) then
          humw = dtdy*hx(mw)*u3
          tmp1 = tmp1 + humw*t101(mw)
          tmp2 = tmp2 + humw*t102(mw)
        elseif (u3 < zero) then
          humw = dtdy*hx(mw)*u3
          tmp1 = tmp1 + humw*t101(mi)
          tmp2 = tmp2 + humw*t102(mi)
        endif
      endif
      if (khv(msrf(i-1,j)) >= kb) then 
        mn = mn0 + kb
        v1 = v(mn)
        if (v1 > zero) then
          hvmn = dtdx2*hy(mn)*v1
          tmp1 = tmp1 - hvmn*t201(mi)
          tmp2 = tmp2 - hvmn*t202(mi)
        elseif (v1 < zero) then
          hvmn = dtdx2*hy(mn)*v1
          tmp1 = tmp1 - hvmn*t201(mn)
          tmp2 = tmp2 - hvmn*t202(mn)
        endif
      endif
      if (khv(n) >= kb) then 
        v3 = v(mi)
        if (v3 > zero) then
          ms = ms0 + kb
          hvmi = dtdx1*hy(mi)*v3
          tmp1 = tmp1 + hvmi*t201(ms)
          tmp2 = tmp2 + hvmi*t202(ms)
        elseif (v3 < zero) then
          hvmi = dtdx1*hy(mi)*v3
          tmp1 = tmp1 + hvmi*t201(mi)
          tmp2 = tmp2 + hvmi*t202(mi)
        endif
      endif
      w1 = w(mi)
      if (w1 > zero) then
        w2   = dtf0*w1
        tmp1 = tmp1 - w2*t301(mi)
        tmp2 = tmp2 - w2*t302(mi)
      elseif (w1 < zero) then
        w2   = dtf0*w1
        tmp1 = tmp1 - w2*t301(mu)
        tmp2 = tmp2 - w2*t302(mu)
      endif
      tt01(mi) = tmp1
      tt02(mi) = tmp2
    enddo
#ifdef _OPENACC
    enddo
!$ACC WAIT
#endif
    

! fixme - consider fusion of the two decompos to get better balancing:
    if (nbpq > 0 .or. nbpz > 0) then
      ! due to barrier

      if (nbpq > 0) then
!$OMP BARRIER
        call domp_get_domain(1, nbpq, nbpql, nbpqu)
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (n,i,fac,mi)                                            
#endif
        do n=nbpql,nbpqu
          i   = krq(1,n)
          mi  = msrf(i,krq(2,n))
          fac = rwqk(n,1)*dtdxdy*cosphi(1,i)
          tt01(mi) = tt01(mi) + fac*rwqkt(1,1,n)
          tt02(mi) = tt02(mi) + fac*rwqkt(2,1,n)
        enddo
#ifdef _OPENACC
!$ACC  parallel loop vector_length(32)                                         &
!$ACC   private (n,i,j,fac0,mi,k,fac)
#endif
        do n=nbpql,nbpqu
          i = krq(1,n)
          j = krq(2,n)
          fac0 = dtdxdy*cosphi(1,i)
          do k=2,kh(msrf(i,j))
            fac = rwqk(n,k)*fac0
            mi  = mcol(msrf(i,j)) + k - 2
            tt01(mi) = tt01(mi) + fac*rwqkt(1,k,n)
            tt02(mi) = tt02(mi) + fac*rwqkt(2,k,n)
          enddo
        enddo
      endif
    
      if (nbpz > 0) then
!$OMP BARRIER
        call domp_get_domain(1, nbpz, nbpzl, nbpzu)
        ddx = one/dx
        ddy = one/dy
        ddt = one/dt
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (k,n,i,j,kb,nue,kbuw,kbue,kbvn,kbvs,facx,fac,nuw,facy1,facy2,  &
!$ACC           nvn,nvs,uuu,vvv,mi,wh,ws,factxy,rwzkin,rwzkout,rwzkin1,rwzkout1)
#endif
        do n=nbpzl,nbpzu
          i = krz(1,n)
          j = krz(2,n)
          if (n>1 .and. i==krz(1,n-1) .and. j==krz(2,n-1)) cycle
          kb = kh(msrf(i,j))
          if (kb < 1) cycle
          kbuw = khu(msrf(i,j-1))
          kbue = khu(msrf(i, j ))
          kbvn = khv(msrf(i-1,j))
          kbvs = khv(msrf(i, j ))
          rwzkin1      =zero
          rwzkout1     =zero
          rwzkin (2:kb)=zero
          rwzkout(2:kb)=zero

          facx   = ddx/cosphi(1,i)
          facy1  = ddy*cosphi(2,i-1)/cosphi(1,i)
          facy2  = ddy*cosphi(2, i )/cosphi(1,i)
          factxy = dtdxdy*cosphi(1,i)

          ! k=1 unrolled:
          nuw = msrf(i,j-1)
          nue = msrf(i, j )
          nvn = msrf(i-1,j)
          nvs = msrf(i, j )
          if (1 > kbuw) nuw=0
          if (1 > kbue) nue=0
          if (1 > kbvn) nvn=0
          if (1 > kbvs) nvs=0
          fac = hx(nue)*facx
          if (abs(krz(3,n)) == 3 .and. u(nue) > zero) then
            !rwzkout(1,n)=rwzkout(1,n)
            rwzkin1 = rwzkin1 + fac*u(nue)
          else
            rwzkin1  = rwzkin1  + fac*min(u(nue),zero)
            rwzkout1 = rwzkout1 + fac*max(u(nue),zero)
          endif
          fac = hx(nuw)*facx
          if (abs(krz(3,n)) == 1 .and. u(nuw) < zero) then
            !rwzkout(1)=rwzkout(1)
            rwzkin1 = rwzkin1 - fac*u(nuw)
          else
            rwzkin1  = rwzkin1  - fac*max(u(nuw),zero)
            rwzkout1 = rwzkout1 - fac*min(u(nuw),zero)
          endif
          fac = hy(nvn)*facy1
          if (abs(krz(3,n)) == 2 .and. v(nvn) > zero) then
            !rwzkout(1)=rwzkout(1)
            rwzkin1 = rwzkin1 + fac*v(nvn)
          else
            rwzkin1  = rwzkin1  + fac*min(v(nvn),zero)
            rwzkout1 = rwzkout1 + fac*max(v(nvn),zero)
          endif
          fac = hy(nvs)*facy2
          if (abs(krz(3,n)) == 4 .and. v(nvs) < zero) then
            !rwzkout(1)=rwzkout(1)
            rwzkin1 = rwzkin1 - fac*v(nvs)
          else
            rwzkin1  = rwzkin1  - fac*max(v(nvs),zero)
            rwzkout1 = rwzkout1 - fac*min(v(nvs),zero)
          endif

          ! fixme: consider moving IF(abs(krz...) branches outside loop
          !        and merge with above IF branches
          do k=2,min(kb,kh(msrf(i,j+1)))
            nue = mcol(msrf(i,j)) + k - 2
            fac = hx(nue)*facx
            uuu = u(nue)
            if (abs(krz(3,n)) == 3 .and. uuu > zero) then
              !rwzkout(k,n)=rwzkout(k,n)
              rwzkin (k)=rwzkin (k) + fac*uuu
            else
              rwzkin (k)=rwzkin (k) + fac*min(uuu,zero)
              rwzkout(k)=rwzkout(k) + fac*max(uuu,zero)
            endif
          enddo
          do k=2,min(kb,kh(msrf(i,j-1)))
            nuw = mcol(msrf(i,j-1)) + k - 2
            fac = hx(nuw)*facx
            uuu = u(nuw)
            if (abs(krz(3,n)) == 1 .and. uuu < zero) then
              !rwzkout(k)=rwzkout(k)
              rwzkin (k)=rwzkin (k) - fac*uuu
            else
              rwzkin (k)=rwzkin (k) - fac*max(uuu,zero)
              rwzkout(k)=rwzkout(k) - fac*min(uuu,zero)
            endif
          enddo
          do k=2,min(kb,kh(msrf(i-1,j)))
            nvn = mcol(msrf(i-1,j)) + k - 2
            fac = hy(nvn)*facy1
            vvv = v(nvn)
            if (abs(krz(3,n)) == 2 .and. vvv > zero) then
              !rwzkout(k)=rwzkout(k)
              rwzkin (k)=rwzkin (k) + fac*vvv
            else
              rwzkin (k)=rwzkin (k) + fac*min(vvv,zero)
              rwzkout(k)=rwzkout(k) + fac*max(vvv,zero)
            endif
          enddo
          do k=2,min(kb,kh(msrf(i+1,j)))
            nvs = mcol(msrf(i,j)) + k - 2
            fac = hy(nvs)*facy2
            vvv = v(nvs)
            if (abs(krz(3,n)) == 4 .and. vvv < zero) then
              !rwzkout(k)=rwzkout(k)
              rwzkin (k)=rwzkin (k) - fac*vvv
            else
              rwzkin (k)=rwzkin (k) - fac*max(vvv,zero)
              rwzkout(k)=rwzkout(k) - fac*min(vvv,zero)
            endif
          enddo

          mi = msrf(i,j)
          wh = ddt*(h_new(mi)-h_old(mi))
          rwzkin1  = rwzkin1  + min(wh,zero)
          rwzkout1 = rwzkout1 + max(wh,zero)
          do k=2,kb
            mi = mcol(msrf(i,j)) + k - 2
            wh = ddt*(h_new(mi)-h_old(mi))
            ws = w(mi)
            rwzkout(k) = rwzkout(k) + max(ws,zero) + max(wh,zero)
            rwzkin (k) = rwzkin (k) + min(ws,zero) + min(wh,zero)
          enddo
          if (kb > 1) then
            ws = w(mcol(msrf(i,j)))
            rwzkin1  = rwzkin1  - max(ws,zero)
            rwzkout1 = rwzkout1 - min(ws,zero)
          endif
          do k=2,kb-1
            ws = w(mcol(msrf(i,j)) + k - 1)
            rwzkin (k) = rwzkin (k) - max(ws,zero)
            rwzkout(k) = rwzkout(k) - min(ws,zero)
          enddo

          mi = msrf(i,j)
          tt01(mi) = tt01(mi)                                                  &
                   + factxy*(rwzkin1*t1d01(mi) + rwzkout1*rwzkt(1,1,n))
          tt02(mi) = tt02(mi)                                                  &
                   + factxy*(rwzkin1*t1d02(mi) + rwzkout1*rwzkt(2,1,n))
          do k=2,kb
            mi = mcol(msrf(i,j)) + k - 2
            tt01(mi) = tt01(mi)                                                &
                     + factxy*(rwzkin(k)*t1d01(mi) + rwzkout(k)*rwzkt(1,k,n))
            tt02(mi) = tt02(mi)                                                &
                     + factxy*(rwzkin(k)*t1d02(mi) + rwzkout(k)*rwzkt(2,k,n))
          enddo
        enddo
      endif

!$OMP BARRIER
    endif
    
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (n,fac)                                                 
#endif
    do n=n2dl,n2du
      fac = ddxdy/cosphi(1,ind(1,n))/h_new(n)
      tt01(n) = tt01(n)*fac
      tt02(n) = tt02(n)*fac
    enddo
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   private (mi,n,kb,i,j,fac0,fac)                                  
#endif
    do n=n2dl,n2du
      kb=kh(n)
      if (kb < 2) cycle
      fac0 = ddxdy/cosphi(1,ind(1,n))
      mi0 = mcol(n) - 2
      do k=2,kb
        mi = mi0 + k
        fac = fac0/h_new(mi)
        tt01(mi) = tt01(mi)*fac
        tt02(mi) = tt02(mi)*fac
      enddo
    enddo
#ifdef _OPENACC
!$ACC END DATA
#endif
    
  end subroutine c_tt
    
  ! ----------------------------------------------------------------------------

  subroutine c_tw (msrf,mcol,ind,n2d,kh,u,v,cosphi,idx,t1d01,t1d02)
      
    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, msrf(0:,0:), mcol(0:), ind(:,:), kh(0:)
    real(8),    intent(in)    :: u(0:), v(0:), cosphi(:,0:)
    real(8),    intent(in)    :: t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow.c_tw.inc'
    
    !- local vars --------------------------------------------------------------
    integer(4) :: n, i, j, k, kb, mi, nn, ss, ee, ww, nw, ne, sw, se, kbm, klow
    integer(4) :: mi0, nn0, ss0, ee0, ww0, nw0, ne0, sw0, se0, n2dl, n2du
    real(8)    :: dtddy, dtddx, facy, facx, mhdtddy, mhdtddx, dtddx0
    real(8)    :: u1, u2, u3, u4, u5, u6, v1, v2, v3, v4, v5, v6
    
    ! --------------------------------------------------------------------------

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,cosphi,kh,msrf,mcol,u,v)                                   &
!$ACC   pcreate(t1d01,t1d02,t301,t302) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif
    dtddy   = dt/dy
    facy    = -onethird*dtddy
    mhdtddy = -half*dtddy
    dtddx0  = dt/dx

#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,klow,mi,n,i,dtddx,j,v5,v6,v2,u1,v4,v3,u2,mhdtddx,u6,        &
!$ACC            u5,facx,u3,u4,kb,kbm,mi0,ww0,nn0,nw0,ss0,sw0,ee0,ne0,         &
!$ACC            se0,ee,ne,nn,nw,se,ss,sw,v1,ww)
#endif
    do n=n2dl,n2du
      kb = kh(n)

      i  = ind(1,n)
      j  = ind(2,n)

      dtddx   = dtddx0/cosphi(1,i)
      facx    = -onethird*dtddx
      mhdtddx = -half*dtddx
      
      ! unroll k=1:
      mi = n
      nn = msrf(i-1,j  )
      ss = msrf(i+1,j  )
      ee = msrf(i,  j+1)
      ww = msrf(i,  j-1)
      nw = msrf(i-1,j-1)
      sw = msrf(i+1,j-1)
      se = msrf(i+1,j+1)
      ne = msrf(i-1,j+1)

      u1 = min(u(mi),zero)
      u2 = max(u(ww),zero)
      u3 = min(u(nn),zero)
      u4 = max(u(nw),zero)
      u5 = min(u(ss),zero)
      u6 = max(u(sw),zero)

      v1 = max(v(mi),zero)
      v2 = min(v(nn),zero)
      v3 = max(v(ww),zero)
      v4 = min(v(nw),zero)
      v5 = max(v(ee),zero)
      v6 = min(v(ne),zero)
    
      t301(mi) = t1d01(mi)                                                     &
               + mhdtddx*(  u2*(  t1d01(mi) - t1d01(ww)                        &
                                + facy*(  v2*(t1d01(nn)-t1d01(mi))             &
                                        + v1*(t1d01(mi)-t1d01(ss))             &
                                        - v4*(t1d01(nw)-t1d01(ww))             &
                                        - v3*(t1d01(ww)-t1d01(sw)) ))          &
                          + u1*(  t1d01(ee) - t1d01(mi)                        &
                                + facy*(  v6*(t1d01(ne)-t1d01(ee))             &
                                        + v5*(t1d01(ee)-t1d01(se))             &
                                        - v2*(t1d01(nn)-t1d01(mi))             &
                                        - v1*(t1d01(mi)-t1d01(ss)) )))         &
               + mhdtddy*(  v2*( t1d01(nn) - t1d01(mi)                         &
                                + facx*(  u4*(t1d01(nn)-t1d01(nw))             &
                                        + u3*(t1d01(ne)-t1d01(nn))             &
                                        - u2*(t1d01(mi)-t1d01(ww))             &
                                        - u1*(t1d01(ee)-t1d01(mi)) ))          &
                          + v1*( t1d01(mi) - t1d01(ss)                         &
                                + facx*(  u2*(t1d01(mi)-t1d01(ww))             &
                                        + u1*(t1d01(ee)-t1d01(mi))             &
                                        - u6*(t1d01(ss)-t1d01(sw))             &
                                        - u5*(t1d01(se)-t1d01(ss)) )))
      t302(mi) = t1d02(mi)                                                     &
               + mhdtddx*(  u2*(  t1d02(mi) - t1d02(ww)                        &
                                + facy*(  v2*(t1d02(nn)-t1d02(mi))             &
                                        + v1*(t1d02(mi)-t1d02(ss))             &
                                        - v4*(t1d02(nw)-t1d02(ww))             &
                                        - v3*(t1d02(ww)-t1d02(sw)) ))          &
                          + u1*(  t1d02(ee) - t1d02(mi)                        &
                                + facy*(  v6*(t1d02(ne)-t1d02(ee))             &
                                        + v5*(t1d02(ee)-t1d02(se))             &
                                        - v2*(t1d02(nn)-t1d02(mi))             &
                                        - v1*(t1d02(mi)-t1d02(ss)) )))         &
               + mhdtddy*(  v2*( t1d02(nn) - t1d02(mi)                         &
                                + facx*(  u4*(t1d02(nn)-t1d02(nw))             &
                                        + u3*(t1d02(ne)-t1d02(nn))             &
                                        - u2*(t1d02(mi)-t1d02(ww))             &
                                        - u1*(t1d02(ee)-t1d02(mi)) ))          &
                          + v1*( t1d02(mi) - t1d02(ss)                         &
                                + facx*(  u2*(t1d02(mi)-t1d02(ww))             &
                                        + u1*(t1d02(ee)-t1d02(mi))             &
                                        - u6*(t1d02(ss)-t1d02(sw))             &
                                        - u5*(t1d02(se)-t1d02(ss)) )))

      if (kb < 2) cycle

      ! do the water column, k=2:kb ...
      ! the nine index offsets:
      mi0 = mcol(n) - 2
      ww0 = mcol(msrf(i,  j-1)) - 2
      nw0 = mcol(msrf(i-1,j-1)) - 2
      sw0 = mcol(msrf(i+1,j-1)) - 2
      nn0 = mcol(msrf(i-1,j  )) - 2
      ss0 = mcol(msrf(i+1,j  )) - 2
      ee0 = mcol(msrf(i,  j+1)) - 2
      ne0 = mcol(msrf(i-1,j+1)) - 2
      se0 = mcol(msrf(i+1,j+1)) - 2

      kbm = min(kb,kh(msrf(i,j-1)),kh(msrf(i-1,j-1)),kh(msrf(i+1,j-1)),        &
                                   kh(msrf(i-1,j)),  kh(msrf(i+1,j)),          &
                   kh(msrf(i,j+1)),kh(msrf(i-1,j+1)),kh(msrf(i+1,j+1)))

      ! ld/st = 12+2*(2+9) = 34; flop = 12+2*66 = 134; I = 134/34 ~ 3.94
      do k=2,kbm
        mi = mi0 + k
        ww = ww0 + k
        nn = nn0 + k
        nw = nw0 + k
        ss = ss0 + k
        sw = sw0 + k
        ee = ee0 + k
        ne = ne0 + k
        se = se0 + k

        u1 = min(u(mi),zero)
        u2 = max(u(ww),zero)
        u3 = min(u(nn),zero)
        u4 = max(u(nw),zero)
        u5 = min(u(ss),zero)
        u6 = max(u(sw),zero)

        v1 = max(v(mi),zero)
        v2 = min(v(nn),zero)
        v3 = max(v(ww),zero)
        v4 = min(v(nw),zero)
        v5 = max(v(ee),zero)
        v6 = min(v(ne),zero)
    
        t301(mi) = t1d01(mi)                                                   &
                 + mhdtddx*(  u2*(  t1d01(mi) - t1d01(ww)                      &
                                  + facy*(  v2*(t1d01(nn)-t1d01(mi))           &
                                          + v1*(t1d01(mi)-t1d01(ss))           &
                                          - v4*(t1d01(nw)-t1d01(ww))           &
                                          - v3*(t1d01(ww)-t1d01(sw)) ))        &
                            + u1*(  t1d01(ee) - t1d01(mi)                      &
                                  + facy*(  v6*(t1d01(ne)-t1d01(ee))           &
                                          + v5*(t1d01(ee)-t1d01(se))           &
                                          - v2*(t1d01(nn)-t1d01(mi))           &
                                          - v1*(t1d01(mi)-t1d01(ss)) )))       &
                 + mhdtddy*(  v2*( t1d01(nn) - t1d01(mi)                       &
                                  + facx*(  u4*(t1d01(nn)-t1d01(nw))           &
                                          + u3*(t1d01(ne)-t1d01(nn))           &
                                          - u2*(t1d01(mi)-t1d01(ww))           &
                                          - u1*(t1d01(ee)-t1d01(mi)) ))        &
                            + v1*( t1d01(mi) - t1d01(ss)                       &
                                  + facx*(  u2*(t1d01(mi)-t1d01(ww))           &
                                          + u1*(t1d01(ee)-t1d01(mi))           &
                                          - u6*(t1d01(ss)-t1d01(sw))           &
                                          - u5*(t1d01(se)-t1d01(ss)) )))
        t302(mi) = t1d02(mi)                                                   &
                 + mhdtddx*(  u2*(  t1d02(mi) - t1d02(ww)                      &
                                  + facy*(  v2*(t1d02(nn)-t1d02(mi))           &
                                          + v1*(t1d02(mi)-t1d02(ss))           &
                                          - v4*(t1d02(nw)-t1d02(ww))           &
                                          - v3*(t1d02(ww)-t1d02(sw)) ))        &
                            + u1*(  t1d02(ee) - t1d02(mi)                      &
                                  + facy*(  v6*(t1d02(ne)-t1d02(ee))           &
                                          + v5*(t1d02(ee)-t1d02(se))           &
                                          - v2*(t1d02(nn)-t1d02(mi))           &
                                          - v1*(t1d02(mi)-t1d02(ss)) )))       &
                 + mhdtddy*(  v2*( t1d02(nn) - t1d02(mi)                       &
                                  + facx*(  u4*(t1d02(nn)-t1d02(nw))           &
                                          + u3*(t1d02(ne)-t1d02(nn))           &
                                          - u2*(t1d02(mi)-t1d02(ww))           &
                                          - u1*(t1d02(ee)-t1d02(mi)) ))        &
                            + v1*( t1d02(mi) - t1d02(ss)                       &
                                  + facx*(  u2*(t1d02(mi)-t1d02(ww))           &
                                          + u1*(t1d02(ee)-t1d02(mi))           &
                                          - u6*(t1d02(ss)-t1d02(sw))           &
                                          - u5*(t1d02(se)-t1d02(ss)) )))
      enddo

      ! remainder loops:
      klow = max(kbm,1) + 1

      do k=klow,kb
        mi = mi0 + k
        t301(mi) = t1d01(mi)
        t302(mi) = t1d02(mi)
      enddo

      do k=klow,min(kb,kh(msrf(i,j-1)))
        mi = mi0 + k
        ww = ww0 + k
        u2 = mhdtddx*max(u(ww),zero)
        t301(mi) = t301(mi) + u2*(t1d01(mi) - t1d01(ww))
        t302(mi) = t302(mi) + u2*(t1d02(mi) - t1d02(ww))
      enddo
      do k=klow,min(kb,kh(msrf(i,j-1)),kh(msrf(i-1,j-1)))
        mi = mi0 + k
        ww = ww0 + k
        nw = nw0 + k
        u2 = mhdtddx*max(u(ww),zero)*facy*min(v(nw),zero)
        t301(mi) = t301(mi) - u2*(t1d01(nw)-t1d01(ww)) 
        t302(mi) = t302(mi) - u2*(t1d02(nw)-t1d02(ww)) 
      enddo
      do k=klow,min(kb,kh(msrf(i,j-1)),kh(msrf(i+1,j-1)))
        mi = mi0 + k
        ww = ww0 + k
        sw = sw0 + k
        u2 = mhdtddx*max(u(ww),zero)*facy*max(v(ww),zero)
        t301(mi) = t301(mi) - u2*(t1d01(ww)-t1d01(sw)) 
        t302(mi) = t302(mi) - u2*(t1d02(ww)-t1d02(sw)) 
      enddo
      do k=klow,min(kb,kh(msrf(i,j-1)),kh(msrf(i-1,j)))
        mi = mi0 + k
        ww = ww0 + k
        nn = nn0 + k
        u2 = mhdtddx*max(u(ww),zero)*facy*min(v(nn),zero)
        t301(mi) = t301(mi) + u2*(t1d01(nn)-t1d01(mi)) 
        t302(mi) = t302(mi) + u2*(t1d02(nn)-t1d02(mi)) 
      enddo
      do k=klow,min(kb,kh(msrf(i,j-1)),kh(msrf(i+1,j)))
        mi = mi0 + k
        ww = ww0 + k
        ss = ss0 + k
        u2 = mhdtddx*max(u(ww),zero)*facy*max(v(mi),zero)
        t301(mi) = t301(mi) + u2*(t1d01(mi)-t1d01(ss)) 
        t302(mi) = t302(mi) + u2*(t1d02(mi)-t1d02(ss)) 
      enddo

      do k=klow,min(kb,kh(msrf(i,j+1)))
        mi = mi0 + k
        ee = ee0 + k
        u1 = mhdtddx*min(u(mi),zero)
        t301(mi) = t301(mi) + u1*(t1d01(ee)-t1d01(mi))
        t302(mi) = t302(mi) + u1*(t1d02(ee)-t1d02(mi))
      enddo
      do k=klow,min(kb,kh(msrf(i,j+1)),kh(msrf(i-1,j)))
        mi = mi0 + k
        nn = nn0 + k
        u1 = mhdtddx*min(u(mi),zero)*facy*min(v(nn),zero)
        t301(mi) = t301(mi) - u1*(t1d01(nn)-t1d01(mi))
        t302(mi) = t302(mi) - u1*(t1d02(nn)-t1d02(mi))
      enddo
      do k=klow,min(kb,kh(msrf(i,j+1)),kh(msrf(i+1,j)))
        mi = mi0 + k
        ss = ss0 + k
        u1 = mhdtddx*min(u(mi),zero)*facy*max(v(mi),zero)
        t301(mi) = t301(mi) - u1*(t1d01(mi)-t1d01(ss))
        t302(mi) = t302(mi) - u1*(t1d02(mi)-t1d02(ss))
      enddo
      do k=klow,min(kb,kh(msrf(i,j+1)),kh(msrf(i-1,j+1)))
        mi = mi0 + k
        ne = ne0 + k
        ee = ee0 + k
        u1 = mhdtddx*min(u(mi),zero)*facy*min(v(ne),zero)
        t301(mi) = t301(mi) + u1*(t1d01(ne)-t1d01(ee))
        t302(mi) = t302(mi) + u1*(t1d02(ne)-t1d02(ee))
      enddo
      do k=klow,min(kb,kh(msrf(i,j+1)),kh(msrf(i+1,j+1)))
        mi = mi0 + k
        se = se0 + k
        ee = ee0 + k
        u1 = mhdtddx*min(u(mi),zero)*facy*max(v(ee),zero)
        t301(mi) = t301(mi) + u1*(t1d01(ee)-t1d01(se))
        t302(mi) = t302(mi) + u1*(t1d02(ee)-t1d02(se))
      enddo

      do k=klow,min(kb,kh(msrf(i-1,j)))
        mi = mi0 + k
        nn = nn0 + k
        v2 = mhdtddy*min(v(nn),zero)
        t301(mi) = t301(mi) + v2*(t1d01(nn)-t1d01(mi))
        t302(mi) = t302(mi) + v2*(t1d02(nn)-t1d02(mi))
      enddo
      do k=klow,min(kb,kh(msrf(i-1,j)),kh(msrf(i-1,j-1)))
        mi = mi0 + k
        nn = nn0 + k
        nw = nw0 + k
        v2 = mhdtddy*min(v(nn),zero)*facx*max(u(nw),zero)
        t301(mi) = t301(mi) + v2*(t1d01(nn)-t1d01(nw))
        t302(mi) = t302(mi) + v2*(t1d02(nn)-t1d02(nw))
      enddo
      do k=klow,min(kb,kh(msrf(i-1,j)),kh(msrf(i-1,j+1)))
        mi = mi0 + k
        nn = nn0 + k
        ne = ne0 + k
        v2 = mhdtddy*min(v(nn),zero)*facx*min(u(nn),zero)
        t301(mi) = t301(mi) + v2*(t1d01(ne)-t1d01(nn))
        t302(mi) = t302(mi) + v2*(t1d02(ne)-t1d02(nn))
      enddo
      do k=klow,min(kb,kh(msrf(i-1,j)),kh(msrf(i,j-1)))
        mi = mi0 + k
        nn = nn0 + k
        ww = ww0 + k
        v2 = mhdtddy*min(v(nn),zero)*facx*max(u(ww),zero)
        t301(mi) = t301(mi) - v2*(t1d01(mi)-t1d01(ww))
        t302(mi) = t302(mi) - v2*(t1d02(mi)-t1d02(ww))
      enddo
      do k=klow,min(kb,kh(msrf(i-1,j)),kh(msrf(i,j+1)))
        mi = mi0 + k
        nn = nn0 + k
        ee = ee0 + k
        v2 = mhdtddy*min(v(nn),zero)*facx*min(u(mi),zero)
        t301(mi) = t301(mi) - v2*(t1d01(ee)-t1d01(mi))
        t302(mi) = t302(mi) - v2*(t1d02(ee)-t1d02(mi))
      enddo

      do k=klow,min(kb,kh(msrf(i+1,j)))
        mi = mi0 + k
        ss = ss0 + k
        v1 =  mhdtddy*max(v(mi),zero)
        t301(mi) = t301(mi) + v1*(t1d01(mi)-t1d01(ss))
        t302(mi) = t302(mi) + v1*(t1d02(mi)-t1d02(ss))
      enddo
      do k=klow,min(kb,kh(msrf(i+1,j)),kh(msrf(i,j-1)))
        mi = mi0 + k
        ww = ww0 + k
        v1 =  mhdtddy*max(v(mi),zero)*facx*max(u(ww),zero)
        t301(mi) = t301(mi) + v1*(t1d01(mi)-t1d01(ww))
        t302(mi) = t302(mi) + v1*(t1d02(mi)-t1d02(ww))
      enddo
      do k=klow,min(kb,kh(msrf(i+1,j)),kh(msrf(i,j+1)))
        mi = mi0 + k
        ee = ee0 + k
        v1 =  mhdtddy*max(v(mi),zero)*facx*min(u(mi),zero)
        t301(mi) = t301(mi) + v1*(t1d01(ee)-t1d01(mi))
        t302(mi) = t302(mi) + v1*(t1d02(ee)-t1d02(mi))
      enddo
      do k=klow,min(kb,kh(msrf(i+1,j)),kh(msrf(i+1,j-1)))
        mi = mi0 + k
        sw = sw0 + k
        ss = ss0 + k
        v1 =  mhdtddy*max(v(mi),zero)*facx*max(u(sw),zero)
        t301(mi) = t301(mi) - v1*(t1d01(ss)-t1d01(sw))
        t302(mi) = t302(mi) - v1*(t1d02(ss)-t1d02(sw))
      enddo
      do k=klow,min(kb,kh(msrf(i+1,j)),kh(msrf(i+1,j+1)))
        mi = mi0 + k
        se = se0 + k
        ss = ss0 + k
        v1 =  mhdtddy*max(v(mi),zero)*facx*min(u(ss),zero)
        t301(mi) = t301(mi) - v1*(t1d01(se)-t1d01(ss))
        t302(mi) = t302(mi) - v1*(t1d02(se)-t1d02(ss))
      enddo

    enddo
#ifdef _OPENACC
    enddo ! streams
!$ACC WAIT
!$ACC END DATA
#endif

  end subroutine c_tw
    
  ! ----------------------------------------------------------------------------

  subroutine c_tu (kmax,msrf,mcol,ind,n2d,advecstab,h_new,khu,v,w,kh,idx,      &
                   t1d01,t1d02)
    
    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in) :: kmax,n2d
    integer(4), intent(in) :: msrf(0:,0:),mcol(0:),ind(:,:),khu(0:),kh(0:)
    real(8), intent(in)    :: advecstab, h_new(0:)
    real(8), intent(in)    :: v(0:),w(0:), t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)
    
    include 'tflow.c_tu.inc'

    !- local vars --------------------------------------------------------------
    integer(4) :: n, i, j, k, mi, kb, kbn, kbs, kbb, kb2, n2dl, n2du
    integer(4) :: nn, ss, nd, md, sd, nu, mu, su, mi0, mn0, ms0
    real(8)    :: fac, dtddz, dtddy, facz, facy, mhdtddz, mhdtddy
    real(8)    :: f1, f2, f3, f4, f5, hdt, v1, v2, v3, v4, v5, v6, wd1, wd2
    real(8)    :: w1(3:max(kmax,3)), w2(3:max(kmax,3)), w1srf, w2srf

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,khu,h_new,msrf,mcol,v,w,kh)                                &
!$ACC   pcreate(t101,t102,t1d01,t1d02,w1,w2) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif

    dtddy   = dt/dy
    facy    = -onethird*dtddy
    mhdtddy = -half*dtddy
    hdt     = half*dt

#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (f1,f2,k,kbb,mi,mu,n,su,kb,i,j,dtddz,facz,                     &
!$ACC            fac,mhdtddz,mi0,mn0,ms0,v2,v1,v3,v4,kbs,kbn,f3,wd1,f4,        &
!$ACC            v6,v5,f5,wd2,md,nd,nn,nu,sd,ss)
#endif
    do n=n2dl, n2du
      kb  = kh(n)
      i   = ind(1,n)
      j   = ind(2,n)
      kbb = max(khu(n),khu(msrf(i,j-1)))

      ! make sure to have nice values throughout
      if (kb > 1 .and. kb > kbb) then
        mi0 = mcol(n) - 2
        t101(mi0+max(2,kbb+1):mi0+kb) = zero
        t102(mi0+max(2,kbb+1):mi0+kb) = zero
      endif
      if (kbb == 0) then
        t101(n) = zero
        t102(n) = zero
        cycle
      endif

      ! First do k=1 -----------------------------------------------------------
      k = 1
      mi = n
      nn = msrf(i-1,j)
      ss = msrf(i+1,j)

      t101(mi) = t1d01(mi)
      t102(mi) = t1d02(mi)
    
      dtddz = dt/(h_new(mi) - advecstab)
      facz  = -onethird*dtddz
      mhdtddz = -half*dtddz

      if (k < kmax) then  ! k=1<kmax
        md = mcol(n)
        if (k < kh(msrf(i+1,j))) then
          sd = mcol(msrf(i+1,j))
        else
          sd = 0
        endif
        if (k < kh(msrf(i-1,j))) then
          nd = mcol(msrf(i-1,j))
        else
          nd = 0
        endif

        if (v(nn) < zero) then
          f1 = t1d01(nn) - t1d01(mi)
          f2 = t1d02(nn) - t1d02(mi)
          if (w(nd) > zero) then
            fac = facz*w(nd)
            f1 = f1 + fac*(t1d01(nn)-t1d01(nd))
            f2 = f2 + fac*(t1d02(nn)-t1d02(nd))
          endif
          if (w(md) > zero) then
            fac = facz*w(md)
            f1 = f1 - fac*(t1d01(mi)-t1d01(md))
            f2 = f2 - fac*(t1d02(mi)-t1d02(md))
          endif
          fac = mhdtddy*v(nn)
          t101(mi) = t101(mi) + fac*f1
          t102(mi) = t102(mi) + fac*f2
        endif
    
        if (v(mi) > zero) then
          f1 = t1d01(mi) - t1d01(ss)
          f2 = t1d02(mi) - t1d02(ss)
          if (w(md) > zero) then
            fac = facz*w(md)
            f1 = f1 + fac*(t1d01(mi)-t1d01(md))
            f2 = f2 + fac*(t1d02(mi)-t1d02(md))
          endif
          if (w(sd) > zero) then
            fac = facz*w(sd)
            f1 = f1 - fac*(t1d01(ss)-t1d01(sd))
            f2 = f2 - fac*(t1d02(ss)-t1d02(sd))
          endif
          fac = mhdtddy*v(mi)
          t101(mi) = t101(mi) + fac*f1
          t102(mi) = t102(mi) + fac*f2
        endif
    
        if (w(md) > zero) then
          f1 = t1d01(mi) - t1d01(md)
          f2 = t1d02(mi) - t1d02(md)
          if (v(nn) < zero) then
            fac = facy*v(nn)
            f1 = f1 + fac*(t1d01(nn)-t1d01(mi))
            f2 = f2 + fac*(t1d02(nn)-t1d02(mi))
          endif
          if (v(mi) > zero) then
            fac = facy*v(mi)
            f1 = f1 + fac*(t1d01(mi)-t1d01(ss))
            f2 = f2 + fac*(t1d02(mi)-t1d02(ss))
          endif
          if (v(nd) < zero) then
            fac = facy*v(nd)
            f1 = f1 - fac*(t1d01(nd)-t1d01(md))
            f2 = f2 - fac*(t1d02(nd)-t1d02(md))
          endif
          if (v(md) > zero) then
            fac = facy*v(md)
            f1 = f1 - fac*(t1d01(md)-t1d01(sd))
            f2 = f2 - fac*(t1d02(md)-t1d02(sd))
          endif
          fac = mhdtddz*w(md)
          t101(mi) = t101(mi) + fac*f1
          t102(mi) = t102(mi) + fac*f2
        endif

      else ! k=1=kmax
        if (v(nn) < zero) then
          fac = mhdtddy*v(nn)
          t101(mi) = t101(mi) + fac*(t1d01(nn)-t1d01(mi))
          t102(mi) = t102(mi) + fac*(t1d02(nn)-t1d02(mi))
        endif
        if (v(mi) > zero) then
          fac = mhdtddy*v(mi)
          t101(mi) = t101(mi) + fac*(t1d01(mi)-t1d01(ss))
          t102(mi) = t102(mi) + fac*(t1d02(mi)-t1d02(ss))
        endif

      endif
      ! This ends treatment of k=1 ---------------------------------------------
      if (kbb < 2) cycle

      ! Do all k in-between k=1 and kmax ---------------------------------------

      mi0 = mcol(n) - 2
      ms0 = mcol(msrf(i+1,j)) - 2
      mn0 = mcol(msrf(i-1,j)) - 2
      kbn = kh(msrf(i-1,j))
      kbs = kh(msrf(i+1,j)) 
      kb2 = min(kbb+1,kb)

      !## PART II:
      k = 2
      mi = mi0 + k
      nn = mn0 + k
      ss = ms0 + k
      mu = n
      nu = msrf(i-1,j)
      su = msrf(i+1,j)
      v1 = min(v(nu),zero)
      v2 = max(v(mu),zero)
      v3 = min(v(nn),zero)
      v4 = max(v(mi),zero)
      w1srf = t1d01(mu)-t1d01(mi) + facy*( v1*(t1d01(nu)-t1d01(mu))            &
                                          +v2*(t1d01(mu)-t1d01(su))            &
                                          -v3*(t1d01(nn)-t1d01(mi))            &
                                          -v4*(t1d01(mi)-t1d01(ss)))
      w2srf = t1d02(mu)-t1d02(mi) + facy*( v1*(t1d02(nu)-t1d02(mu))            &
                                          +v2*(t1d02(mu)-t1d02(su))            &
                                          -v3*(t1d02(nn)-t1d02(mi))            &
                                          -v4*(t1d02(mi)-t1d02(ss)))
      do k=3,min(kb2,kbn,kbs)
        ! Flops  (4*min/max, 28 *Mult/Add/Sub); 32
        ! Memory (4 + 7 + 7 +2s): 20 / 22
        ! Intensity: 32/20 or 32/22
        mi = mi0 + k
        nn = mn0 + k
        ss = ms0 + k
        mu = mi - 1
        nu = nn - 1
        su = ss - 1
        v1 = min(v(nu),zero)
        v2 = max(v(mu),zero)
        v3 = min(v(nn),zero)
        v4 = max(v(mi),zero)
        w1(k) = t1d01(mu)-t1d01(mi) + facy*( v1*(t1d01(nu)-t1d01(mu))          &
                                            +v2*(t1d01(mu)-t1d01(su))          &
                                            -v3*(t1d01(nn)-t1d01(mi))          &
                                            -v4*(t1d01(mi)-t1d01(ss)))
        w2(k) = t1d02(mu)-t1d02(mi) + facy*( v1*(t1d02(nu)-t1d02(mu))          &
                                            +v2*(t1d02(mu)-t1d02(su))          &
                                            -v3*(t1d02(nn)-t1d02(mi))          &
                                            -v4*(t1d02(mi)-t1d02(ss)))
      enddo
      if (kbs <= min(kb2,kbn)) then
        do k=max(3,kbs+1),min(kb2,kbn)
          ! Flops  (2*min/max, 2*8 Mult/Add/Sub); 18
          ! Memory (2 + 4 + 4 +2s): 12 / 14
          ! Intensity: 18/12 or 18/14
          mi = mi0 + k
          nn = mn0 + k
          mu = mi - 1
          nu = nn - 1
          v1 = min(v(nu),zero)
          v3 = min(v(nn),zero)
          w1(k) = t1d01(mu)-t1d01(mi) + facy*( v1*(t1d01(nu)-t1d01(mu))        &
                                              -v3*(t1d01(nn)-t1d01(mi)))
          w2(k) = t1d02(mu)-t1d02(mi) + facy*( v1*(t1d02(nu)-t1d02(mu))        &
                                              -v3*(t1d02(nn)-t1d02(mi)))
        enddo
        do k=max(3,kbn+1),kb2
          ! Flops  (2 Mult/Add/Sub); 2
          ! Memory (2 + 2 +2s): 6 / 8
          ! Intensity: 2/6 or 2/8
          mi = mi0 + k
          mu = mi - 1
          w1(k) = t1d01(mu)-t1d01(mi)
          w2(k) = t1d02(mu)-t1d02(mi)
        enddo
      elseif (kbn <= min(kb2,kbs)) then
        do k=max(3,kbn+1),min(kb2,kbs)
          ! Flops  (2*min/max, 2*8 Mult/Add/Sub); 18
          ! Memory (2 + 4 + 4 +2s): 12 / 14
          ! Intensity: 18/12 or 18/14
          mi = mi0 + k
          ss = ms0 + k
          mu = mi -1
          su = ss - 1
          v2 = max(v(mu),zero)
          v4 = max(v(mi),zero)
          w1(k) = t1d01(mu)-t1d01(mi) + facy*( v2*(t1d01(mu)-t1d01(su))        &
                                              -v4*(t1d01(mi)-t1d01(ss)))
          w2(k) = t1d02(mu)-t1d02(mi) + facy*( v2*(t1d02(mu)-t1d02(su))        &
                                              -v4*(t1d02(mi)-t1d02(ss)))
        enddo
        do k=max(3,kbs+1),kb2
          ! Flops  (2 Mult/Add/Sub); 2
          ! Memory (2 + 2 +2s): 6 / 8
          ! Intensity: 2/6 or 2/8
          mi = mi0 + k
          mu = mi - 1
          w1(k) = t1d01(mu)-t1d01(mi)
          w2(k) = t1d02(mu)-t1d02(mi)
        enddo
      endif

      !## PART I:

      ! k=2 unrolled:
      k = 2
      mi = mi0 + k
      nn = mn0 + k
      ss = ms0 + k
      mu = n
      nu = msrf(i-1,j)
      su = msrf(i+1,j)
      md = 0
      nd = 0
      sd = 0
      if (kb  > 2) md = mi+1
      if (kbn > 2) nd = nn+1
      if (kbs > 2) sd = ss+1
      if (kb2 > k) then
        wd1 = w1(k+1)
        wd2 = w2(k+1)
      else
        wd1 = zero
        wd2 = zero
      endif
      f1 = mhdtddy*min(v(nn),zero)
      f2 = mhdtddy*max(v(mi),zero)
      f3 = -onethird*dt/h_new(mi)
      f4 =-hdt/h_new(mi)
      f5 = (f2-f1)*f3
      v1 = min(w(nn),zero)
      v2 = max(w(nd),zero)
      v3 = min(w(mi),zero)
      v4 = max(w(md),zero)
      v5 = min(w(ss),zero)
      v6 = max(w(sd),zero)
      t101(mi) = t1d01(mi) + f4*(v3*w1srf + v4*wd1)                            &
               + f1*(t1d01(nn)-t1d01(mi) + f3*(  v1*(t1d01(nu)-t1d01(nn))      &
                                               + v2*(t1d01(nn)-t1d01(nd)) ))   &
               + f2*(t1d01(mi)-t1d01(ss) + f3*(- v5*(t1d01(su)-t1d01(ss))      &
                                               - v6*(t1d01(ss)-t1d01(sd)) ))   &
               + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
      t102(mi) = t1d02(mi) + f4*(v3*w2srf + v4*wd2)                            &
               + f1*(t1d02(nn)-t1d02(mi) + f3*(  v1*(t1d02(nu)-t1d02(nn))      &
                                               + v2*(t1d02(nn)-t1d02(nd)) ))   &
               + f2*(t1d02(mi)-t1d02(ss) + f3*(- v5*(t1d02(su)-t1d02(ss))      &
                                               - v6*(t1d02(ss)-t1d02(sd)) ))   &
               + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))

      ! main loop:
      do k=3,min(kbb,kbn,kbs)-1
        ! Flops  (8*min/max, 7+2*33+ *Div/Mult/Add/Sub); 81
        ! Memory (10 + 2*13 + 2s): 38 / 40
        ! Intensity: 81/38 or 81/40
        mi = mi0 + k
        nn = mn0 + k
        ss = ms0 + k
        mu = mi - 1
        md = mi + 1
        nu = nn - 1
        nd = nn + 1
        su = ss - 1
        sd = ss + 1
        f1 = mhdtddy*min(v(nn),zero)
        f2 = mhdtddy*max(v(mi),zero)
        f3 = -onethird*dt/h_new(mi)
        f4 =-hdt/h_new(mi)
        f5 = (f2-f1)*f3
        v1 = min(w(nn),zero)
        v2 = max(w(nd),zero)
        v3 = min(w(mi),zero)
        v4 = max(w(md),zero)
        v5 = min(w(ss),zero)
        v6 = max(w(sd),zero)
        t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                      &
                 + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn))     &
                                                 +v2*(t1d01(nn)-t1d01(nd))))   &
                 + f2*(t1d01(mi)-t1d01(ss) + f3*(-v5*(t1d01(su)-t1d01(ss))     &
                                                 -v6*(t1d01(ss)-t1d01(sd))))   &
                 + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
        t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                      &
                 + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn))     &
                                                 +v2*(t1d02(nn)-t1d02(nd))))   &
                 + f2*(t1d02(mi)-t1d02(ss) + f3*(-v5*(t1d02(su)-t1d02(ss))     &
                                                 -v6*(t1d02(ss)-t1d02(sd))))   &
                 + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
      enddo

      ! remainder stuff:
      if (kbb <= min(kbn,kbs)) then
        if (kbb >= 3) then
          k = kbb
          mi = mi0 + k
          nn = mn0 + k
          ss = ms0 + k
          mu = mi - 1
          nu = nn - 1
          su = ss - 1
          md = 0
          nd = 0
          sd = 0
          if (kb  > kbb) md = mi + 1
          if (kbn > kbb) nd = nn + 1
          if (kbs > kbb) sd = ss + 1
          if (kb2 > k) then
            wd1 = w1(k+1)
            wd2 = w2(k+1)
          else
            wd1 = zero
            wd2 = zero
          endif
          f1 = mhdtddy*min(v(nn),zero)
          f2 = mhdtddy*max(v(mi),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          f5 = (f2-f1)*f3
          v1 = min(w(nn),zero)
          v2 = max(w(nd),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ss),zero)
          v6 = max(w(sd),zero)
          t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                        &
                   + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn))   &
                                                   +v2*(t1d01(nn)-t1d01(nd)))) &
                   + f2*(t1d01(mi)-t1d01(ss) + f3*(-v5*(t1d01(su)-t1d01(ss))   &
                                                   -v6*(t1d01(ss)-t1d01(sd)))) &
                   + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
          t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                        &
                   + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn))   &
                                                   +v2*(t1d02(nn)-t1d02(nd)))) &
                   + f2*(t1d02(mi)-t1d02(ss) + f3*(-v5*(t1d02(su)-t1d02(ss))   &
                                                   -v6*(t1d02(ss)-t1d02(sd)))) &
                   + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
        endif
      elseif (kbn <= min(kbb,kbs)) then
        if (kbn >= 3) then
          k = kbn
          mi = mi0 + k
          nn = mn0 + k
          ss = ms0 + k
          mu = mi - 1
          nu = nn - 1
          su = ss - 1
          md = 0
          sd = 0
          if (kb  > kbn) md = mi + 1
          if (kbs > kbn) sd = ss + 1
          if (kb2 > k) then
            wd1 = w1(k+1)
            wd2 = w2(k+1)
          else
            wd1 = zero
            wd2 = zero
          endif
          f1 = mhdtddy*min(v(nn),zero)
          f2 = mhdtddy*max(v(mi),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          f5 = (f2-f1)*f3
          v1 = min(w(nn),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ss),zero)
          v6 = max(w(sd),zero)
          t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                        &
                   + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn)))) &
                   + f2*(t1d01(mi)-t1d01(ss) + f3*(-v5*(t1d01(su)-t1d01(ss))   &
                                                   -v6*(t1d01(ss)-t1d01(sd)))) &
                   + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
          t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                        &
                   + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn)))) &
                   + f2*(t1d02(mi)-t1d02(ss) + f3*(-v5*(t1d02(su)-t1d02(ss))   &
                                                   -v6*(t1d02(ss)-t1d02(sd)))) &
                   + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
        endif
        do k=max(3,kbn+1),min(kbb,kbs)-1
          mi = mi0 + k
          ss = ms0 + k
          mu = mi - 1
          md = mi + 1
          su = ss - 1
          sd = ss + 1
          f2 = mhdtddy*max(v(mi),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ss),zero)
          v6 = max(w(sd),zero)
          t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                    &
                   + f2*(t1d01(mi)-t1d01(ss) + f3*( v3*(t1d01(mu)-t1d01(mi))   &
                                                   +v4*(t1d01(mi)-t1d01(md))   &
                                                   -v5*(t1d01(su)-t1d01(ss))   &
                                                   -v6*(t1d01(ss)-t1d01(sd))))
          t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                    &
                   + f2*(t1d02(mi)-t1d02(ss) + f3*( v3*(t1d02(mu)-t1d02(mi))   &
                                                   +v4*(t1d02(mi)-t1d02(md))   &
                                                   -v5*(t1d02(su)-t1d02(ss))   &
                                                   -v6*(t1d02(ss)-t1d02(sd))))
        enddo
        if (kbb <= kbs) then
          if (kbb >= 3 .and. kbn < kbb) then
            k = kbb
            mi = mi0 + k
            ss = ms0 + k
            mu = mi - 1
            su = ss - 1
            md = 0
            sd = 0
            if (kb  > kbb) md = mi + 1
            if (kbs > kbb) sd = ss + 1
            if (kb2 > k) then
              wd1 = w1(k+1)
              wd2 = w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f2 = mhdtddy*max(v(mi),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            v5 = min(w(ss),zero)
            v6 = max(w(sd),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                      &
                     + f2*(t1d01(mi)-t1d01(ss) + f3*( v3*(t1d01(mu)-t1d01(mi)) &
                                                     +v4*(t1d01(mi)-t1d01(md)) &
                                                     -v5*(t1d01(su)-t1d01(ss)) &
                                                     -v6*(t1d01(ss)-t1d01(sd))))
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                      &
                     + f2*(t1d02(mi)-t1d02(ss) + f3*( v3*(t1d02(mu)-t1d02(mi)) &
                                                     +v4*(t1d02(mi)-t1d02(md)) &
                                                     -v5*(t1d02(su)-t1d02(ss)) &
                                                     -v6*(t1d02(ss)-t1d02(sd))))
          endif
        else ! kbs < kbb
          if (kbs >= 3 .and. kbn < kbs) then
            k = kbs
            mi = mi0 + k
            ss = ms0 + k
            mu = mi - 1
            su = ss - 1
            md = 0
            sd = 0
            if (kb > kbs) md = mi + 1
            f2 = mhdtddy*max(v(mi),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            v5 = min(w(ss),zero)
            v6 = max(w(sd),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                  &
                     + f2*(t1d01(mi)-t1d01(ss) + f3*( v3*(t1d01(mu)-t1d01(mi)) &
                                                     +v4*(t1d01(mi)-t1d01(md)) &
                                                     -v5*(t1d01(su)-t1d01(ss)) &
                                                     -v6*(t1d01(ss)-t1d01(sd))))
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                  &
                     + f2*(t1d02(mi)-t1d02(ss) + f3*( v3*(t1d02(mu)-t1d02(mi)) &
                                                     +v4*(t1d02(mi)-t1d02(md)) &
                                                     -v5*(t1d02(su)-t1d02(ss)) &
                                                     -v6*(t1d02(ss)-t1d02(sd))))
          endif
          do k=max(3,kbs+1),kbb-1
            mi = mi0 + k
            md = mi + 1
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))
          enddo
          if (kbb >= 3) then
            k = kbb
            mi = mi0 + k
            if (kb2 > k) then
              md = mi + 1
              wd1 = max(w(md),zero)*w1(k+1)
              wd2 = max(w(md),zero)*w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + wd1)
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + wd2)
          endif
        endif

      elseif (kbs <= min(kbb,kbn)) then
        if (kbs >= 3) then
          k = kbs
          mi = mi0 + k
          nn = mn0 + k
          ss = ms0 + k
          mu = mi - 1
          nu = nn - 1
          su = ss - 1
          md = 0
          nd = 0
          if (kb  > kbs) md = mi + 1
          if (kbn > kbs) nd = nn + 1
          if (kb2 > k) then
            wd1 = w1(k+1)
            wd2 = w2(k+1)
          else
            wd1 = zero
            wd2 = zero
          endif
          f1 = mhdtddy*min(v(nn),zero)
          f2 = mhdtddy*max(v(mi),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          f5 = (f2-f1)*f3
          v1 = min(w(nn),zero)
          v2 = min(w(nd),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ss),zero)
          t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                        &
                   + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn))   &
                                                   +v2*(t1d01(nn)-t1d01(nd)))) &
                   + f2*(t1d01(mi)-t1d01(ss) + f3*(-v5*(t1d01(su)-t1d01(ss)))) &
                   + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
          t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                        &
                   + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn))   &
                                                   +v2*(t1d02(nn)-t1d02(nd)))) &
                   + f2*(t1d02(mi)-t1d02(ss) + f3*(-v5*(t1d02(su)-t1d02(ss)))) &
                   + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
        endif
        do k=max(3,kbs+1),min(kbb,kbn)-1
          mi = mi0 + k
          nn = mn0 + k
          mu = mi - 1
          md = mi + 1
          nu = nn - 1
          nd = nn + 1
          f1 = mhdtddy*min(v(nn),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          v1 = min(w(nn),zero)
          v2 = min(w(nd),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                    &
                   + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn))   &
                                                   +v2*(t1d01(nn)-t1d01(nd))   &
                                                   -v3*(t1d01(mu)-t1d01(mi))   &
                                                   -v4*(t1d01(mi)-t1d01(md))))
          t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                    &
                   + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn))   &
                                                   +v2*(t1d02(nn)-t1d02(nd))   &
                                                   -v3*(t1d02(mu)-t1d02(mi))   &
                                                   -v4*(t1d02(mi)-t1d02(md))))
        enddo
        if (kbb <= kbn) then
          if (kbb >= 3 .and. kbs < kbb) then
            k = kbb
            mi = mi0 + k
            nn = mn0 + k
            mu = mi - 1
            nu = nn - 1
            md = 0
            nd = 0
            if (kb  > kbb) md = mi + 1
            if (kbn > kbb) nd = nn + 1
            if (kb2 > k) then
              wd1 = w1(k+1)
              wd2 = w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f1 = mhdtddy*max(v(nn),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v1 = min(w(nn),zero)
            v2 = min(w(nd),zero)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                      &
                     + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn)) &
                                                     +v2*(t1d01(nn)-t1d01(nd)) &
                                                     -v3*(t1d01(mu)-t1d01(mi)) &
                                                     -v4*(t1d01(mi)-t1d01(md))))
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                      &
                     + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn)) &
                                                     +v2*(t1d02(nn)-t1d02(nd)) &
                                                     -v3*(t1d02(mu)-t1d02(mi)) &
                                                     -v4*(t1d02(mi)-t1d02(md))))
          endif
        else ! kbn < kbb
          if (kbn >= 3 .and. kbs < kbn) then
            k = kbn
            mi = mi0 + k
            nn = mn0 + k
            mu = mi - 1
            nu = nn - 1
            md = 0
            nd = 0
            if (kb  > kbn) md = mi + 1
            f1 = mhdtddy*max(v(nn),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v1 = min(w(nn),zero)
            v2 = min(w(nd),zero)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                  &
                     + f1*(t1d01(nn)-t1d01(mi) + f3*( v1*(t1d01(nu)-t1d01(nn)) &
                                                     +v2*(t1d01(nn)-t1d01(nd)) &
                                                     -v3*(t1d01(mu)-t1d01(mi)) &
                                                     -v4*(t1d01(mi)-t1d01(md))))
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                  &
                     + f1*(t1d02(nn)-t1d02(mi) + f3*( v1*(t1d02(nu)-t1d02(nn)) &
                                                     +v2*(t1d02(nn)-t1d02(nd)) &
                                                     -v3*(t1d02(mu)-t1d02(mi)) &
                                                     -v4*(t1d02(mi)-t1d02(md))))
          endif
          do k=max(3,kbn+1),kbb-1
            mi = mi0 + k
            md = mi + 1
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))
          enddo
          if (kbb >= 3) then
            k = kbb
            mi = mi0 + k
            if (kb2 > k) then
              md = mi + 1
              wd1 = max(w(md),zero)*w1(k+1)
              wd2 = max(w(md),zero)*w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            t101(mi) = t1d01(mi) + f4*(v3*w1(k) + wd1)
            t102(mi) = t1d02(mi) + f4*(v3*w2(k) + wd2)
          endif
        endif

      endif

    enddo ! n loop

#ifdef _OPENACC
    enddo ! stream
!$ACC WAIT
!$ACC END DATA
#endif

  end subroutine c_tu

  ! ----------------------------------------------------------------------------

  subroutine c_tv (kmax,msrf,mcol,ind,n2d,advecstab,h_new,khv,u,w,cosphi,kh,   &
                   idx,t1d01,t1d02)
    
    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: kmax,n2d
    integer(4), intent(in)    :: msrf(0:,0:),mcol(0:),ind(:,:),khv(0:),kh(0:)
    real(8),    intent(in)    :: advecstab, h_new(0:), cosphi(:,0:)
    real(8),     intent(in)   :: u(0:), w(0:), t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)
    
    include 'tflow.c_tv.inc'

    !- local vars --------------------------------------------------------------
    integer(4) :: n, i, j, k, mi, kb, kbe, kbw, kbb, kb2, n2dl, n2du
    integer(4) :: ee, ww, ed, md, wd, eu, mu, wu, mi0, me0, mw0
    real(8)    :: fac, dtddz, dtddx, dtddx0, facz, facx, mhdtddz, mhdtddx
    real(8)    :: f1, f2, f3, f4, f5, hdt, v1, v2, v3, v4, v5, v6, wd1, wd2
    real(8)    :: w1(3:max(kmax,3)), w2(3:max(kmax,3)), w1srf, w2srf
    
    ! --------------------------------------------------------------------------

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,khv,h_new,khv,msrf,mcol,u,w,kh,cosphi)                     &
!$ACC   pcreate(t1d01,t1d02,t201,t202,w1,w2) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif

    dtddx0 = dt/dx
    hdt    = half*dt

#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (f1,f2,k,kbb,mi,mu,n,wu,kb,i,j,dtddx,dtddz,                    &
!$ACC            facz,fac,mhdtddx,facx,mhdtddz,mi0,me0,mw0,v2,v1,v3,v4,        &
!$ACC            kbw,kbe,f3,wd1,f4,v6,v5,f5,wd2,ed,ee,eu,md,wd,ww)
#endif
    do n=n2dl, n2du
      kb  = kh(n)
      i   = ind(1,n)
      j   = ind(2,n)
      kbb = max(khv(n),khv(msrf(i-1,j)))

      ! make sure to have nice values throughout
      if (kb > 1 .and. kb > kbb) then
        mi0 = mcol(n) - 2
        t201(mi0+max(2,kbb+1):mi0+kb) = zero
        t202(mi0+max(2,kbb+1):mi0+kb) = zero
      endif
      if (kbb == 0) then
        t201(n) = zero
        t202(n) = zero
        cycle
      endif

      dtddx   = dtddx0/cosphi(1,i)
      facx    = -onethird*dtddx
      mhdtddx = -half*dtddx

      ! First do k=1 -----------------------------------------------------------
      k = 1
      mi = n
      ee = msrf(i,j+1)
      ww = msrf(i,j-1)

      t201(mi) = t1d01(mi)
      t202(mi) = t1d02(mi)
    
      dtddz = dt/(h_new(mi) - advecstab)
      facz  = -onethird*dtddz
      mhdtddz = -half*dtddz

      if (k < kmax) then  ! k=1<kmax
        md = mcol(n)
        if (k < kh(msrf(i,j+1))) then
          ed = mcol(msrf(i,j+1))
        else
          ed = 0
        endif
        if (k < kh(msrf(i,j-1))) then
          wd = mcol(msrf(i,j-1))
        else
          wd = 0
        endif

        if (u(mi) < zero) then
          f1 = t1d01(ee) - t1d01(mi)
          f2 = t1d02(ee) - t1d02(mi)
          if (w(ed) > zero) then
            fac = facz*w(ed)
            f1 = f1 + fac*(t1d01(ee)-t1d01(ed))
            f2 = f2 + fac*(t1d02(ee)-t1d02(ed))
          endif
          if (w(md) > zero) then
            fac = facz*w(md)
            f1 = f1 - fac*(t1d01(mi)-t1d01(md))
            f2 = f2 - fac*(t1d02(mi)-t1d02(md))
          endif
          fac = mhdtddx*u(mi)
          t201(mi) = t201(mi) + fac*f1
          t202(mi) = t202(mi) + fac*f2
        endif
    
        if (u(ww) > zero) then
          f1 = t1d01(mi) - t1d01(ww)
          f2 = t1d02(mi) - t1d02(ww)
          if (w(md) > zero) then
            fac = facz*w(md)
            f1 = f1 + fac*(t1d01(mi)-t1d01(md))
            f2 = f2 + fac*(t1d02(mi)-t1d02(md))
          endif
          if (w(wd) > zero) then
            fac = facz*w(wd)
            f1 = f1 - fac*(t1d01(ww)-t1d01(wd))
            f2 = f2 - fac*(t1d02(ww)-t1d02(wd))
          endif
          fac = mhdtddx*u(ww)
          t201(mi) = t201(mi) + fac*f1
          t202(mi) = t202(mi) + fac*f2
        endif
    
        if (w(md) > zero) then
          f1 = t1d01(mi) - t1d01(md)
          f2 = t1d02(mi) - t1d02(md)
          if (u(mi) < zero) then
            fac = facx*u(mi)
            f1 = f1 + fac*(t1d01(ee)-t1d01(mi))
            f2 = f2 + fac*(t1d02(ee)-t1d02(mi))
          endif
          if (u(ww) > zero) then
            fac = facx*u(ww)
            f1 = f1 + fac*(t1d01(mi)-t1d01(ww))
            f2 = f2 + fac*(t1d02(mi)-t1d02(ww))
          endif
          if (u(md) < zero) then
            fac = facx*u(md)
            f1 = f1 - fac*(t1d01(ed)-t1d01(md))
            f2 = f2 - fac*(t1d02(ed)-t1d02(md))
          endif
          if (u(wd) > zero) then
            fac = facx*u(wd)
            f1 = f1 - fac*(t1d01(md)-t1d01(wd))
            f2 = f2 - fac*(t1d02(md)-t1d02(wd))
          endif
          fac = mhdtddz*w(md)
          t201(mi) = t201(mi) + fac*f1
          t202(mi) = t202(mi) + fac*f2
        endif

      else ! k=1=kmax
        if (u(mi) < zero) then
          fac = mhdtddx*u(mi)
          t201(mi) = t201(mi) + fac*(t1d01(ee)-t1d01(mi))
          t202(mi) = t202(mi) + fac*(t1d02(ee)-t1d02(mi))
        endif
        if (u(ww) > zero) then
          fac = mhdtddx*u(ww)
          t201(mi) = t201(mi) + fac*(t1d01(mi)-t1d01(ww))
          t202(mi) = t202(mi) + fac*(t1d02(mi)-t1d02(ww))
        endif

      endif
      ! This ends treatment of k=1 ---------------------------------------------
      if (kbb < 2) cycle

      ! Do all k in-between k=1 and kmax ---------------------------------------

      mi0 = mcol(n) - 2
      mw0 = mcol(msrf(i,j-1)) - 2
      me0 = mcol(msrf(i,j+1)) - 2
      kbe = kh(msrf(i,j+1))
      kbw = kh(msrf(i,j-1)) 
      kb2 = min(kbb+1,kb)

      !## PART II:
      k = 2
      mi = mi0 + k
      ee = me0 + k
      ww = mw0 + k
      mu = n
      eu = msrf(i,j+1)
      wu = msrf(i,j-1)
      v1 = min(u(mu),zero)
      v2 = max(u(wu),zero)
      v3 = min(u(mi),zero)
      v4 = max(u(ww),zero)
      w1srf = t1d01(mu)-t1d01(mi) + facx*( v1*(t1d01(eu)-t1d01(mu))            &
                                          +v2*(t1d01(mu)-t1d01(wu))            &
                                          -v3*(t1d01(ee)-t1d01(mi))            &
                                          -v4*(t1d01(mi)-t1d01(ww)))
      w2srf = t1d02(mu)-t1d02(mi) + facx*( v1*(t1d02(eu)-t1d02(mu))            &
                                          +v2*(t1d02(mu)-t1d02(wu))            &
                                          -v3*(t1d02(ee)-t1d02(mi))            &
                                          -v4*(t1d02(mi)-t1d02(ww)))
      do k=3,min(kb2,kbe,kbw)
        ! Flops  (4*min/max, 28 *Mult/Add/Sub); 32
        ! Memory (4 + 7 + 7 +2s): 20 / 22
        ! Intensity: 32/20 or 32/22
        mi = mi0 + k
        ee = me0 + k
        ww = mw0 + k
        mu = mi - 1
        eu = ee - 1
        wu = ww - 1
        v1 = min(u(mu),zero)
        v2 = max(u(wu),zero)
        v3 = min(u(mi),zero)
        v4 = max(u(ww),zero)
        w1(k) = t1d01(mu)-t1d01(mi) + facx*( v1*(t1d01(eu)-t1d01(mu))          &
                                            +v2*(t1d01(mu)-t1d01(wu))          &
                                            -v3*(t1d01(ee)-t1d01(mi))          &
                                            -v4*(t1d01(mi)-t1d01(ww)))
        w2(k) = t1d02(mu)-t1d02(mi) + facx*( v1*(t1d02(eu)-t1d02(mu))          &
                                            +v2*(t1d02(mu)-t1d02(wu))          &
                                            -v3*(t1d02(ee)-t1d02(mi))          &
                                            -v4*(t1d02(mi)-t1d02(ww)))
      enddo
      if (kbw <= min(kb2,kbe)) then
        do k=max(3,kbw+1),min(kb2,kbe)
          ! Flops  (2*min/max, 2*8 Mult/Add/Sub); 18
          ! Memory (2 + 4 + 4 +2s): 12 / 14
          ! Intensity: 18/12 or 18/14
          mi = mi0 + k
          ee = me0 + k
          mu = mi - 1
          eu = ee - 1
          v1 = min(u(mu),zero)
          v3 = min(u(mi),zero)
          w1(k) = t1d01(mu)-t1d01(mi) + facx*( v1*(t1d01(eu)-t1d01(mu))        &
                                              -v3*(t1d01(ee)-t1d01(mi)))
          w2(k) = t1d02(mu)-t1d02(mi) + facx*( v1*(t1d02(eu)-t1d02(mu))        &
                                              -v3*(t1d02(ee)-t1d02(mi)))
        enddo
        do k=max(3,kbe+1),kb2
          ! Flops  (2 Mult/Add/Sub); 2
          ! Memory (2 + 2 +2s): 6 / 8
          ! Intensity: 2/6 or 2/8
          mi = mi0 + k
          mu = mi - 1
          w1(k) = t1d01(mu)-t1d01(mi)
          w2(k) = t1d02(mu)-t1d02(mi)
        enddo
      elseif (kbe <= min(kb2,kbw)) then
        do k=max(3,kbe+1),min(kb2,kbw)
          ! Flops  (2*min/max, 2*8 Mult/Add/Sub); 18
          ! Memory (2 + 4 + 4 +2s): 12 / 14
          ! Intensity: 18/12 or 18/14
          mi = mi0 + k
          ww = mw0 + k
          mu = mi - 1
          wu = ww - 1
          v2 = max(u(wu),zero)
          v4 = max(u(ww),zero)
          w1(k) = t1d01(mu)-t1d01(mi) + facx*( v2*(t1d01(mu)-t1d01(wu))        &
                                              -v4*(t1d01(mi)-t1d01(ww)))
          w2(k) = t1d02(mu)-t1d02(mi) + facx*( v2*(t1d02(mu)-t1d02(wu))        &
                                              -v4*(t1d02(mi)-t1d02(ww)))
        enddo
        do k=max(3,kbw+1),kb2
          ! Flops  (2 Mult/Add/Sub); 2
          ! Memory (2 + 2 +2s): 6 / 8
          ! Intensity: 2/6 or 2/8
          mi = mi0 + k
          mu = mi - 1
          w1(k) = t1d01(mu)-t1d01(mi)
          w2(k) = t1d02(mu)-t1d02(mi)
        enddo
      endif

      !## PART I:

      ! k=2 unrolled:
      k = 2
      mi = mi0 + k
      ee = me0 + k
      ww = mw0 + k
      mu = n
      eu = msrf(i,j+1)
      wu = msrf(i,j-1)
      md = 0
      ed = 0
      wd = 0
      if (kb  > 2) md = mi + 1
      if (kbe > 2) ed = ee + 1
      if (kbw > 2) wd = ww + 1
      if (kb2 > k) then
        wd1 = w1(k+1)
        wd2 = w2(k+1)
      else
        wd1 = zero
        wd2 = zero
      endif
      f1 = mhdtddx*min(u(mi),zero)
      f2 = mhdtddx*max(u(ww),zero)
      f3 = -onethird*dt/h_new(mi)
      f4 =-hdt/h_new(mi)
      f5 = (f2-f1)*f3
      v1 = min(w(ee),zero)
      v2 = max(w(ed),zero)
      v3 = min(w(mi),zero)
      v4 = max(w(md),zero)
      v5 = min(w(ww),zero)
      v6 = max(w(wd),zero)
      t201(mi) = t1d01(mi) + f4*(v3*w1srf + v4*wd1)                            &
               + f1*(t1d01(ee)-t1d01(mi) + f3*(  v1*(t1d01(eu)-t1d01(ee))      &
                                               + v2*(t1d01(ee)-t1d01(ed)) ))   &
               + f2*(t1d01(mi)-t1d01(ww) + f3*(- v5*(t1d01(wu)-t1d01(ww))      &
                                               - v6*(t1d01(ww)-t1d01(wd)) ))   &
               + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
      t202(mi) = t1d02(mi) + f4*(v3*w2srf + v4*wd2)                            &
               + f1*(t1d02(ee)-t1d02(mi) + f3*(  v1*(t1d02(eu)-t1d02(ee))      &
                                               + v2*(t1d02(ee)-t1d02(ed)) ))   &
               + f2*(t1d02(mi)-t1d02(ww) + f3*(- v5*(t1d02(wu)-t1d02(ww))      &
                                               - v6*(t1d02(ww)-t1d02(wd)) ))   &
               + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))

      ! main loop:
      do k=3,min(kbb,kbe,kbw)-1
        ! Flops  (8*min/max, 7+2*33+ *Div/Mult/Add/Sub); 81
        ! Memory (10 + 2*13 + 2s): 38 / 40
        ! Intensity: 81/38 or 81/40
        mi = mi0 + k
        ee = me0 + k
        ww = mw0 + k
        mu = mi - 1
        md = mi + 1
        eu = ee - 1
        ed = ee + 1
        wu = ww - 1
        wd = ww + 1
        f1 = mhdtddx*min(u(mi),zero)
        f2 = mhdtddx*max(u(ww),zero)
        f3 = -onethird*dt/h_new(mi)
        f4 =-hdt/h_new(mi)
        f5 = (f2-f1)*f3
        v1 = min(w(ee),zero)
        v2 = max(w(ed),zero)
        v3 = min(w(mi),zero)
        v4 = max(w(md),zero)
        v5 = min(w(ww),zero)
        v6 = max(w(wd),zero)
        t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                      &
                 + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee))     &
                                                 +v2*(t1d01(ee)-t1d01(ed))))   &
                 + f2*(t1d01(mi)-t1d01(ww) + f3*(-v5*(t1d01(wu)-t1d01(ww))     &
                                                 -v6*(t1d01(ww)-t1d01(wd))))   &
                 + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
        t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                      &
                 + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee))     &
                                                 +v2*(t1d02(ee)-t1d02(ed))))   &
                 + f2*(t1d02(mi)-t1d02(ww) + f3*(-v5*(t1d02(wu)-t1d02(ww))     &
                                                 -v6*(t1d02(ww)-t1d02(wd))))   &
                 + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
      enddo

      ! remainder stuff:
      if (kbb <= min(kbe,kbw)) then
        if (kbb >= 3) then
          k = kbb
          mi = mi0 + k
          ee = me0 + k
          ww = mw0 + k
          mu = mi - 1
          eu = ee - 1
          wu = ww - 1
          md = 0
          ed = 0
          wd = 0
          if (kb  > kbb) md = mi + 1
          if (kbe > kbb) ed = ee + 1
          if (kbw > kbb) wd = ww + 1
          if (kb2 > k) then
            wd1 = w1(k+1)
            wd2 = w2(k+1)
          else
            wd1 = zero
            wd2 = zero
          endif
          f1 = mhdtddx*min(u(mi),zero)
          f2 = mhdtddx*max(u(ww),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          f5 = (f2-f1)*f3
          v1 = min(w(ee),zero)
          v2 = max(w(ed),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ww),zero)
          v6 = max(w(wd),zero)
          t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                        &
                   + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee))   &
                                                   +v2*(t1d01(ee)-t1d01(ed)))) &
                   + f2*(t1d01(mi)-t1d01(ww) + f3*(-v5*(t1d01(wu)-t1d01(ww))   &
                                                   -v6*(t1d01(ww)-t1d01(wd)))) &
                   + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
          t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                        &
                   + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee))   &
                                                   +v2*(t1d02(ee)-t1d02(ed)))) &
                   + f2*(t1d02(mi)-t1d02(ww) + f3*(-v5*(t1d02(wu)-t1d02(ww))   &
                                                   -v6*(t1d02(ww)-t1d02(wd)))) &
                   + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
        endif
      elseif (kbe <= min(kbb,kbw)) then
        if (kbe >= 3) then
          k = kbe
          mi = mi0 + k
          ee = me0 + k
          ww = mw0 + k
          mu = mi - 1
          eu = ee - 1
          wu = ww - 1
          md = 0
          wd = 0
          if (kb  > kbe) md = mi + 1
          if (kbw > kbe) wd = ww + 1
          if (kb2 > k) then
            wd1 = w1(k+1)
            wd2 = w2(k+1)
          else
            wd1 = zero
            wd2 = zero
          endif
          f1 = mhdtddx*min(u(mi),zero)
          f2 = mhdtddx*max(u(ww),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          f5 = (f2-f1)*f3
          v1 = min(w(ee),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ww),zero)
          v6 = max(w(wd),zero)
          t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                        &
                   + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee)))) &
                   + f2*(t1d01(mi)-t1d01(ww) + f3*(-v5*(t1d01(wu)-t1d01(ww))   &
                                                   -v6*(t1d01(ww)-t1d01(wd)))) &
                   + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
          t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                        &
                   + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee)))) &
                   + f2*(t1d02(mi)-t1d02(ww) + f3*(-v5*(t1d02(wu)-t1d02(ww))   &
                                                   -v6*(t1d02(ww)-t1d02(wd)))) &
                   + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
        endif
        do k=max(3,kbe+1),min(kbb,kbw)-1
          mi = mi0 + k
          ww = mw0 + k
          mu = mi - 1
          md = mi + 1
          wu = ww - 1
          wd = ww + 1
          f2 = mhdtddx*max(u(ww),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ww),zero)
          v6 = max(w(wd),zero)
          t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                    &
                   + f2*(t1d01(mi)-t1d01(ww) + f3*( v3*(t1d01(mu)-t1d01(mi))   &
                                                   +v4*(t1d01(mi)-t1d01(md))   &
                                                   -v5*(t1d01(wu)-t1d01(ww))   &
                                                   -v6*(t1d01(ww)-t1d01(wd))))
          t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                    &
                   + f2*(t1d02(mi)-t1d02(ww) + f3*( v3*(t1d02(mu)-t1d02(mi))   &
                                                   +v4*(t1d02(mi)-t1d02(md))   &
                                                   -v5*(t1d02(wu)-t1d02(ww))   &
                                                   -v6*(t1d02(ww)-t1d02(wd))))
        enddo
        if (kbb <= kbw) then
          if (kbb >= 3 .and. kbe < kbb) then
            k = kbb
            mi = mi0 + k
            ww = mw0 + k
            mu = mi - 1
            wu = ww - 1
            md = 0
            wd = 0
            if (kb  > kbb) md = mi + 1
            if (kbw > kbb) wd = ww + 1
            if (kb2 > k) then
              wd1 = w1(k+1)
              wd2 = w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f2 = mhdtddx*max(u(ww),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            v5 = min(w(ww),zero)
            v6 = max(w(wd),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                      &
                     + f2*(t1d01(mi)-t1d01(ww) + f3*( v3*(t1d01(mu)-t1d01(mi)) &
                                                     +v4*(t1d01(mi)-t1d01(md)) &
                                                     -v5*(t1d01(wu)-t1d01(ww)) &
                                                     -v6*(t1d01(ww)-t1d01(wd))))
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                      &
                     + f2*(t1d02(mi)-t1d02(ww) + f3*( v3*(t1d02(mu)-t1d02(mi)) &
                                                     +v4*(t1d02(mi)-t1d02(md)) &
                                                     -v5*(t1d02(wu)-t1d02(ww)) &
                                                     -v6*(t1d02(ww)-t1d02(wd))))
          endif
        else ! kbw < kbb
          if (kbw >= 3 .and. kbe < kbw) then
            k = kbw
            mi = mi0 + k
            ww = mw0 + k
            mu = mi - 1
            wu = ww - 1
            md = 0
            wd = 0
            if (kb > kbw) md = mi + 1
            f2 = mhdtddx*max(u(ww),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            v5 = min(w(ww),zero)
            v6 = max(w(wd),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                  &
                     + f2*(t1d01(mi)-t1d01(ww) + f3*( v3*(t1d01(mu)-t1d01(mi)) &
                                                     +v4*(t1d01(mi)-t1d01(md)) &
                                                     -v5*(t1d01(wu)-t1d01(ww)) &
                                                     -v6*(t1d01(ww)-t1d01(wd))))
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                  &
                     + f2*(t1d02(mi)-t1d02(ww) + f3*( v3*(t1d02(mu)-t1d02(mi)) &
                                                     +v4*(t1d02(mi)-t1d02(md)) &
                                                     -v5*(t1d02(wu)-t1d02(ww)) &
                                                     -v6*(t1d02(ww)-t1d02(wd))))
          endif
          do k=max(3,kbw+1),kbb-1
            mi = mi0 + k
            md = mi + 1
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))
          enddo
          if (kbb >= 3) then
            k = kbb
            mi = mi0 + k
            if (kb2 > k) then
              md = mi + 1
              wd1 = max(w(md),zero)*w1(k+1)
              wd2 = max(w(md),zero)*w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + wd1)
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + wd2)
          endif
        endif

      elseif (kbw <= min(kbb,kbe)) then
        if (kbw >= 3) then
          k = kbw
          mi = mi0 + k
          ee = me0 + k
          ww = mw0 + k
          mu = mi - 1
          eu = ee - 1
          wu = ww - 1
          md = 0
          ed = 0
          if (kb  > kbw) md = mi + 1
          if (kbe > kbw) ed = ee + 1
          if (kb2 > k) then
            wd1 = w1(k+1)
            wd2 = w2(k+1)
          else
            wd1 = zero
            wd2 = zero
          endif
          f1 = mhdtddx*min(u(mi),zero)
          f2 = mhdtddx*max(u(ww),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          f5 = (f2-f1)*f3
          v1 = min(w(ee),zero)
          v2 = min(w(ed),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          v5 = min(w(ww),zero)
          t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                        &
                   + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee))   &
                                                   +v2*(t1d01(ee)-t1d01(ed)))) &
                   + f2*(t1d01(mi)-t1d01(ww) + f3*(-v5*(t1d01(wu)-t1d01(ww)))) &
                   + f5*(v3*(t1d01(mu)-t1d01(mi)) - v4*(t1d01(mi)-t1d01(md)))
          t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                        &
                   + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee))   &
                                                   +v2*(t1d02(ee)-t1d02(ed)))) &
                   + f2*(t1d02(mi)-t1d02(ww) + f3*(-v5*(t1d02(wu)-t1d02(ww)))) &
                   + f5*(v3*(t1d02(mu)-t1d02(mi)) - v4*(t1d02(mi)-t1d02(md)))
        endif
        do k=max(3,kbw+1),min(kbb,kbe)-1
          mi = mi0 + k
          ee = me0 + k
          mu = mi - 1
          md = mi + 1
          eu = ee - 1
          ed = ee + 1
          f1 = mhdtddx*min(u(mi),zero)
          f3 = -onethird*dt/h_new(mi)
          f4 =-hdt/h_new(mi)
          v1 = min(w(ee),zero)
          v2 = min(w(ed),zero)
          v3 = min(w(mi),zero)
          v4 = max(w(md),zero)
          t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                    &
                   + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee))   &
                                                   +v2*(t1d01(ee)-t1d01(ed))   &
                                                   -v3*(t1d01(mu)-t1d01(mi))   &
                                                   -v4*(t1d01(mi)-t1d01(md))))
          t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                    &
                   + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee))   &
                                                   +v2*(t1d02(ee)-t1d02(ed))   &
                                                   -v3*(t1d02(mu)-t1d02(mi))   &
                                                   -v4*(t1d02(mi)-t1d02(md))))
        enddo
        if (kbb <= kbe) then
          if (kbb >= 3 .and. kbw < kbb) then
            k = kbb
            mi = mi0 + k
            ee = me0 + k
            mu = mi - 1
            eu = ee - 1
            md = 0
            ed = 0
            if (kb  > kbb) md = mi + 1
            if (kbe > kbb) ed = ee + 1
            if (kb2 > k) then
              wd1 = w1(k+1)
              wd2 = w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f1 = mhdtddx*max(u(mi),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v1 = min(w(ee),zero)
            v2 = min(w(ed),zero)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*wd1)                      &
                     + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee)) &
                                                     +v2*(t1d01(ee)-t1d01(ed)) &
                                                     -v3*(t1d01(mu)-t1d01(mi)) &
                                                     -v4*(t1d01(mi)-t1d01(md))))
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*wd2)                      &
                     + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee)) &
                                                     +v2*(t1d02(ee)-t1d02(ed)) &
                                                     -v3*(t1d02(mu)-t1d02(mi)) &
                                                     -v4*(t1d02(mi)-t1d02(md))))
          endif
        else ! kbe < kbb
          if (kbe >= 3 .and. kbw < kbe) then
            k = kbe
            mi = mi0 + k
            ee = me0 + k
            mu = mi - 1
            eu = ee - 1
            md = 0
            ed = 0
            if (kb > kbe) md = mi + 1
            f1 = mhdtddx*max(u(mi),zero)
            f3 = -onethird*dt/h_new(mi)
            f4 =-hdt/h_new(mi)
            v1 = min(w(ee),zero)
            v2 = min(w(ed),zero)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))                  &
                     + f1*(t1d01(ee)-t1d01(mi) + f3*( v1*(t1d01(eu)-t1d01(ee)) &
                                                     +v2*(t1d01(ee)-t1d01(ed)) &
                                                     -v3*(t1d01(mu)-t1d01(mi)) &
                                                     -v4*(t1d01(mi)-t1d01(md))))
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))                  &
                     + f1*(t1d02(ee)-t1d02(mi) + f3*( v1*(t1d02(eu)-t1d02(ee)) &
                                                     +v2*(t1d02(ee)-t1d02(ed)) &
                                                     -v3*(t1d02(mu)-t1d02(mi)) &
                                                     -v4*(t1d02(mi)-t1d02(md))))
          endif
          do k=max(3,kbe+1),kbb-1
            mi = mi0 + k
            md = mi + 1
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            v4 = max(w(md),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + v4*w1(k+1))
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + v4*w2(k+1))
          enddo
          if (kbb >= 3) then
            k = kbb
            mi = mi0 + k
            if (kb2 > k) then
              md = mi + 1
              wd1 = max(w(md),zero)*w1(k+1)
              wd2 = max(w(md),zero)*w2(k+1)
            else
              wd1 = zero
              wd2 = zero
            endif
            f4 =-hdt/h_new(mi)
            v3 = min(w(mi),zero)
            t201(mi) = t1d01(mi) + f4*(v3*w1(k) + wd1)
            t202(mi) = t1d02(mi) + f4*(v3*w2(k) + wd2)
          endif
        endif

      endif

    enddo ! n loop
#ifdef _OPENACC
    enddo ! streams
!$ACC WAIT
!$ACC END DATA
#endif

  end subroutine c_tv

  subroutine advection(msrf,mcol,ind,n2d,kh,khu,khv,h_new,cosphi,idx,          &
                       t1d01,t1d02)
    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain
    
    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, msrf(0:,0:), mcol(0:), ind(:,:)
    integer(4), intent(in)    :: kh(0:), khu(0:), khv(0:)
    real(8),    intent(in)    :: h_new(0:), cosphi(:,0:)
    real(8),    intent(inout) :: t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_advection.inc'
    
    !- local vars --------------------------------------------------------------
    !  simple vars:
    integer(4) :: n, i, j, k, kb, mi, mn, mw, md, mi0, mn0, mw0
    integer(4) :: n2dl, n2du, khuw, khvn
    real(8)    :: fac0, fac, fx00, fy00, fy01, fy02

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC    present(ind,kh,h_new,msrf,mcol,kh,cosphi,khu,khv)                     &
!$ACC    pcreate(t1d01,t1d02,t401,t402,t501,t502,t601,t602,tt01,tt02) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (mi,n,kb,i,fac0,j,fy02,mw,mn,fy01,fac,md,fx00,fy00)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif
    ! --------------------------------------------------------------------------
    ! compute t
    do n=n2dl,n2du
      kb = kh(n)
      if (kb < 1) cycle
      i = ind(1,n)
      j = ind(2,n)
      fac0 = one/cosphi(1,i)
      fx00 = fac0/dx
      fy00 = fac0/dy
      fy01 = fy00*cosphi(2,i-1)
      fy02 = fy00*cosphi(2,i  )

      ! k=1 unrolled:
      mi  = n
      mw  = msrf(i,j-1)
      mn  = msrf(i-1,j)
      fac = dt/h_new(mi) 
      t1d01(mi) = tt01(mi)                                                     &
                + fac*(  fx00*(t401(mw) -      t401(mi))                       &
                       + fy02* t501(mi) - fy01*t501(mn))
      t1d02(mi) = tt02(mi)                                                     &
                + fac*(  fx00*(t402(mw) -      t402(mi))                       &
                       + fy02* t502(mi) - fy01*t502(mn))

      if (kb > 1) then
        md  = mcol(n)
        fac = dt/h_new(mi)
        t1d01(mi) = t1d01(mi) + fac*t601(md)
        t1d02(mi) = t1d02(mi) + fac*t602(md)
      endif
    enddo

#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (fy00,k,md,mi,mn,mw,n,kb,i,fac0,j,mi0,fy02,fx00,fac,           &
!$ACC            khuw,mw0,khvn,mn0,fy01)                                
#endif
    do n=n2dl,n2du
      kb = kh(n)
      if (kb <= 1) cycle

      i = ind(1,n)
      j = ind(2,n)
      fac0 = one/cosphi(1,i)
      fx00 = fac0/dx
      fy00 = fac0/dy
      fy01 = fy00*cosphi(2,i-1)
      fy02 = fy00*cosphi(2,i  )

      khuw = khu(msrf(i,j-1))
      khvn = khv(msrf(i-1,j))

      mi0 = mcol(n) - 2
      mn0 = mcol(msrf(i-1,j  )) - 2
      mw0 = mcol(msrf(i  ,j-1)) - 2
      do k=2,kb-1
        mi  = mi0 + k
        md  = mi + 1
        fac = dt/h_new(mi) 
        t1d01(mi) = tt01(mi)                                                   &
                   + fac*(                - fx00*t401(mi)                      &
                          + fy02*t501(mi)                                      &
                          +      t601(md) -      t601(mi))
        t1d02(mi) = tt02(mi)                                                   &
                   + fac*(                - fx00*t402(mi)                      &
                          + fy02*t502(mi)                                      &
                          +      t602(md) -      t602(mi))
      enddo
      do k=2,min(kb-1,khuw)
        mi  = mi0 + k
        mw  = mw0 + k
        fac = dt*fx00/h_new(mi)
        t1d01(mi) = t1d01(mi) + fac*t401(mw)
        t1d02(mi) = t1d02(mi) + fac*t402(mw)
      enddo
      do k=2,min(kb-1,khvn)
        mi  = mi0 + k
        mn  = mn0 + k
        fac = dt*fy01/h_new(mi)
        t1d01(mi) = t1d01(mi) - fac*t501(mn)
        t1d02(mi) = t1d02(mi) - fac*t502(mn)
      enddo

      ! k=kb unrolled:
      mi = mi0 + kb
      if (kb <= kh(msrf(i,j-1))) then
        mw = mw0 + kb
      else
        mw = 0
      endif
      if (kb <= kh(msrf(i-1,j))) then
        mn = mn0 + kb
      else
        mn = 0
      endif
      fac = dt/h_new(mi) 
      t1d01(mi) = tt01(mi)                                                     &
                + fac*(  fx00*(t401(mw) -      t401(mi))                       &
                       + fy02* t501(mi) - fy01*t501(mn)                        &
                                        -      t601(mi))
      t1d02(mi) = tt02(mi)                                                     &
                + fac*(  fx00*(t402(mw) -      t402(mi))                       &
                       + fy02* t502(mi) - fy01*t502(mn)                        &
                                       -       t602(mi))
    enddo
#ifdef _OPENACC
    enddo ! stream
!$ACC WAIT
!$ACC END DATA
#endif
  end subroutine advection 

  subroutine c_rin_rout (kmax,msrf,mcol,ind,n2d,kh,h_new,cosphi,idx,t1d01,t1d02)

    ! NOTE: it modifies the global module variables: t7, t8  
    ! and it uses tt, t1, t2, t3 

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: kmax, n2d
    integer(4), intent(in)    :: msrf(0:,0:), mcol(0:), ind(:,:), kh(0:)
    real(8),    intent(in)    :: cosphi(:,0:), h_new(0:), t1d01(0:), t1d02(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_rin.inc'

    !- local variables ---------------------------------------------------------
    integer(4) :: n, k, kb, i, j, mi, me, ms, mu, mw, mn, md, n2dl, n2du
    integer(4) :: mi0, mw0, mn0, me0, ms0, khuw, khvn, khue, khvs
    integer(4) :: klow1, klow2, klow3, klow4
    real(8)    :: pin1, qin1, pout1, qout1, pin2, qin2, pout2, qout2
    real(8)    :: fac0
    real(8)    :: facx, facy1, facy2, dh, tmin1, tmax1, tmin2, tmax2
    real(8)    :: tmint1(2:max(2,kmax)), tmaxt1(2:max(2,kmax))
    real(8)    :: tmint2(2:max(2,kmax)), tmaxt2(2:max(2,kmax))

    real(8), parameter :: small = 1.e-12_8

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,kh,msrf,mcol,h_new,cosphi)                                 &
!$ACC   pcreate(t101,t102,t201,t202,t301,t302,t701,t702,t801,t802,tt01,        &
!$ACC           tt02,t1d01,t1d02) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif
    ! k=1 unrolled:
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (mi,n,kb,i,fac0,j,dh,facy1,facx,mw,facy2,mn,md)         
#endif
    do n = n2dl,n2du
      kb = kh(n)
      i = ind(1,n)
      j = ind(2,n)

      fac0  = one/cosphi(1,i)
      facx  = fac0/dx
      facy1 = fac0*cosphi(2,i  )/dy
      facy2 = fac0*cosphi(2,i-1)/dy

      mi = n
      mw = msrf(i,  j-1)
      mn = msrf(i-1,j  )
      if (kb > 1) then
        md = mcol(n)
      else
        md = 0
      endif
      dh = dt/h_new(mi)
      t801(mi) = ((  max(t101(mi),zero)       - min(t101(mw),zero))*facx       &
                   - min(t201(mi),zero)*facy1 + max(t201(mn),zero)*facy2       &
                   + max(t301(mi),zero)       - min(t301(md),zero))*dh
      t802(mi) = ((  max(t102(mi),zero)       - min(t102(mw),zero))*facx       &
                   - min(t202(mi),zero)*facy1 + max(t202(mn),zero)*facy2       &
                   + max(t302(mi),zero)       - min(t302(md),zero))*dh
      t701(mi) = ((- min(t101(mi),zero)       + max(t101(mw),zero))*facx       &
                   + max(t201(mi),zero)*facy1 - min(t201(mn),zero)*facy2       &
                   - min(t301(mi),zero)       + max(t301(md),zero))*dh
      t702(mi) = ((- min(t102(mi),zero)       + max(t102(mw),zero))*facx       &
                   + max(t202(mi),zero)*facy1 - min(t202(mn),zero)*facy2       &
                   - min(t302(mi),zero)       + max(t302(md),zero))*dh
    enddo

#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,klow1,klow2,klow3,md,mi,mn,mw,n,kb,i,fac0,j,khvn,           &
!$ACC            khuw,mi0,mw0,mn0,dh,facy2,facx,facy1)                  
#endif
    do n = n2dl,n2du
      kb = kh(n)
      if (kb <= 1) cycle
      i = ind(1,n)
      j = ind(2,n)

      fac0  = one/cosphi(1,i)
      facx  = fac0/dx
      facy1 = fac0*cosphi(2,i  )/dy
      facy2 = fac0*cosphi(2,i-1)/dy

      khuw = kh(msrf(i,j-1))
      khvn = kh(msrf(i-1,j))

      !  column offsets:
      mi0 = mcol(n) - 2
      mw0 = mcol(msrf(i,  j-1)) - 2
      mn0 = mcol(msrf(i-1,j  )) - 2

      do k=2,min(kb-1,khuw,khvn)
        mi = mi0 + k
        md = mi + 1
        mw = mw0 + k
        mn = mn0 + k
        dh = dt/h_new(mi)
        t801(mi) = ( (max(t101(mi),zero)       - min(t101(mw),zero))*facx      &
                    + max(t201(mn),zero)*facy2 - min(t201(mi),zero) *facy1     &
                    + max(t301(mi),zero)       - min(t301(md),zero)       )*dh
        t802(mi) = ( (max(t102(mi),zero)       - min(t102(mw),zero))*facx      &
                    + max(t202(mn),zero)*facy2 - min(t202(mi),zero) *facy1     &
                    + max(t302(mi),zero)       - min(t302(md),zero)       )*dh
        t701(mi) = ( (max(t101(mw),zero)       - min(t101(mi),zero))*facx      &
                    + max(t201(mi),zero)*facy1 - min(t201(mn),zero) *facy2     &
                    + max(t301(md),zero)       - min(t301(mi),zero)       )*dh
        t702(mi) = ( (max(t102(mw),zero)       - min(t102(mi),zero))*facx      &
                    + max(t202(mi),zero)*facy1 - min(t202(mn),zero) *facy2     &
                    + max(t302(md),zero)       - min(t302(mi),zero)       )*dh
      enddo

      klow1 = max(2, min(khuw,khvn)+1)
      do k=klow1,kb-1
        mi = mi0 + k
        md = mi + 1
        dh = dt/h_new(mi)
        t801(mi) = (  max(t101(mi),zero)*facx                                  &
                    - min(t201(mi),zero)*facy1                                 &
                    + max(t301(mi),zero) - min(t301(md),zero))*dh
        t802(mi) = (  max(t102(mi),zero)*facx                                  &
                    - min(t202(mi),zero)*facy1                                 &
                    + max(t302(mi),zero) - min(t302(md),zero))*dh
        t701(mi) = (- min(t101(mi),zero)*facx                                  &
                    + max(t201(mi),zero)*facy1                                 &
                    - min(t301(mi),zero) + max(t301(md),zero))*dh
        t702(mi) = (- min(t102(mi),zero)*facx                                  &
                    + max(t202(mi),zero)*facy1                                 &
                    - min(t302(mi),zero) + max(t302(md),zero))*dh
      enddo
      klow2 = max(2, khvn+1)
      do k=klow2,min(kb-1,khuw)
        mi = mi0 + k
        mw = mw0 + k
        dh = facx*dt/h_new(mi)
        t801(mi) = t801(mi) - min(t101(mw),zero)*dh 
        t802(mi) = t802(mi) - min(t102(mw),zero)*dh 
        t701(mi) = t701(mi) + max(t101(mw),zero)*dh
        t702(mi) = t702(mi) + max(t102(mw),zero)*dh
      enddo
      klow3 = max(2, khuw+1)
      do k=klow3,min(kb-1,khvn)
        mi = mi0 + k
        mn = mn0 + k
        dh = facy2*dt/h_new(mi)
        t801(mi) = t801(mi) + max(t201(mn),zero)*dh
        t802(mi) = t802(mi) + max(t202(mn),zero)*dh
        t701(mi) = t701(mi) - min(t201(mn),zero)*dh
        t702(mi) = t702(mi) - min(t202(mn),zero)*dh
      enddo

      ! k=kb
      mi = mi0 + kb
      if (kb <= khuw) then
        mw = mw0 + kb
      else
        mw = 0
      endif
      if (kb <= khvn) then
        mn = mn0 + kb
      else
        mn = 0
      endif
      dh   = dt/h_new(mi)
      t801(mi) = ((  max(t101(mi),zero)       - min(t101(mw),zero))*facx       &
                   - min(t201(mi),zero)*facy1 + max(t201(mn),zero)*facy2       &
                   + max(t301(mi),zero)                            )*dh
      t802(mi) = ((  max(t102(mi),zero)       - min(t102(mw),zero))*facx       &
                   - min(t202(mi),zero)*facy1 + max(t202(mn),zero)*facy2       &
                   + max(t302(mi),zero)                            )*dh
      t701(mi) = ((- min(t101(mi),zero)       + max(t101(mw),zero))*facx       &
                   + max(t201(mi),zero)*facy1 - min(t201(mn),zero)*facy2       &
                   - min(t301(mi),zero)                            )*dh
      t702(mi) = ((- min(t102(mi),zero)       + max(t102(mw),zero))*facx       &
                   + max(t202(mi),zero)*facy1 - min(t202(mn),zero)*facy2       &
                   - min(t302(mi),zero)                            )*dh
    enddo
    
    ! limiter loop:
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,klow1,klow2,klow3,klow4,mi,mu,n,tmax1,tmax2,tmaxt1,         &
!$ACC            tmaxt2,tmin1,tmin2,tmint1,tmint2,kb,i,j,pin1,qin1,            &
!$ACC            pout1,qout1,pin2,qin2,pout2,qout2,mi0,khvs,khue,khvn,         &
!$ACC            khuw,me0,mw0,ms0,mn0,md,me,mn,ms,mw) 
#endif
    do n = n2dl,n2du
      kb = kh(n)
      i = ind(1,n)
      j = ind(2,n)

      ! k=1
      mi = n
      me = msrf(i  ,j+1)
      mw = msrf(i  ,j-1)
      ms = msrf(i+1,j  )
      mn = msrf(i-1,j  )
      tmax1 = max(t1d01(mi),tt01(mi))
      tmax2 = max(t1d02(mi),tt02(mi))
      tmin1 = min(t1d01(mi),tt01(mi))
      tmin2 = min(t1d02(mi),tt02(mi))
      if (mn > 0) then
        tmax1 = max(tmax1,t1d01(mn),tt01(mn))
        tmax2 = max(tmax2,t1d02(mn),tt02(mn))
        tmin1 = min(tmin1,t1d01(mn),tt01(mn))
        tmin2 = min(tmin2,t1d02(mn),tt02(mn))
      endif
      if (ms > 0) then
        tmax1 = max(tmax1,t1d01(ms),tt01(ms))
        tmax2 = max(tmax2,t1d02(ms),tt02(ms))
        tmin1 = min(tmin1,t1d01(ms),tt01(ms))
        tmin2 = min(tmin2,t1d02(ms),tt02(ms))
      endif
      if (mw > 0) then
        tmax1 = max(tmax1,t1d01(mw),tt01(mw))
        tmax2 = max(tmax2,t1d02(mw),tt02(mw))
        tmin1 = min(tmin1,t1d01(mw),tt01(mw))
        tmin2 = min(tmin2,t1d02(mw),tt02(mw))
      endif
      if (me > 0) then
        tmax1 = max(tmax1,t1d01(me),tt01(me))
        tmax2 = max(tmax2,t1d02(me),tt02(me))
        tmin1 = min(tmin1,t1d01(me),tt01(me))
        tmin2 = min(tmin2,t1d02(me),tt02(me))
      endif
      if (kb > 1) then
        md = mcol(n)
        tmax1 = max(tmax1,t1d01(md),tt01(md))
        tmax2 = max(tmax2,t1d02(md),tt02(md))
        tmin1 = min(tmin1,t1d01(md),tt01(md))
        tmin2 = min(tmin2,t1d02(md),tt02(md))
      endif
      pin1 = t701(mi)
      !     pin>=0
      if (pin1 > small) then
        qin1 = tmax1-tt01(mi)
      !     qin>=0
        t701(mi) = min(one,qin1/pin1)
      else
        t701(mi) = one
      endif
      pout1 = t801(mi)
      !     pout>=0
      if (pout1 > small) then
        qout1 = tt01(mi)-tmin1
      !     qout>=0
        t801(mi) = min(one,qout1/pout1)
      else
        t801(mi) = one
      endif

      pin2 = t702(mi)
      !     pin>=0
      if (pin2 > small) then
        qin2 = tmax2-tt02(mi)
      !     qin>=0
        t702(mi) = min(one,qin2/pin2)
      else
        t702(mi) = one
      endif
      pout2 = t802(mi)
      !     pout>=0
      if (pout2 > small) then
        qout2 = tt02(mi)-tmin2
      !     qout>=0
        t802(mi) = min(one,qout2/pout2)
      else
        t802(mi) = one
      endif

      if (kb <= 1) cycle

      mi0 = mcol(n) - 2

      ! unroll k=2:
      k = 2
      mi = mi0 + k
      mu = n
      if (kb > 2) then
        md = mi + 1
        tmaxt1(k) =max(t1d01(mi),tt01(mi),t1d01(md),tt01(md),t1d01(mu),tt01(mu))
        tmaxt2(k) =max(t1d02(mi),tt02(mi),t1d02(md),tt02(md),t1d02(mu),tt02(mu))
        tmint1(k) =min(t1d01(mi),tt01(mi),t1d01(md),tt01(md),t1d01(mu),tt01(mu))
        tmint2(k) =min(t1d02(mi),tt02(mi),t1d02(md),tt02(md),t1d02(mu),tt02(mu))
      else
        tmaxt1(k) = max(t1d01(mi),tt01(mi),t1d01(mu),tt01(mu))
        tmaxt2(k) = max(t1d02(mi),tt02(mi),t1d02(mu),tt02(mu))
        tmint1(k) = min(t1d01(mi),tt01(mi),t1d01(mu),tt01(mu))
        tmint2(k) = min(t1d02(mi),tt02(mi),t1d02(mu),tt02(mu))
      endif

      do k=3,kb-1
        mi = mi0 + k
        md = mi + 1
        mu = mi - 1
        tmaxt1(k) =max(t1d01(mi),tt01(mi),t1d01(md),tt01(md),t1d01(mu),tt01(mu))
        tmaxt2(k) =max(t1d02(mi),tt02(mi),t1d02(md),tt02(md),t1d02(mu),tt02(mu))
        tmint1(k) =min(t1d01(mi),tt01(mi),t1d01(md),tt01(md),t1d01(mu),tt01(mu))
        tmint2(k) =min(t1d02(mi),tt02(mi),t1d02(md),tt02(md),t1d02(mu),tt02(mu))
      enddo

      ! unroll k=kb:
      if (kb > 2) then
        k = kb
        mi = mi0 + k
        mu = mi - 1
        tmaxt1(k) = max(t1d01(mi),tt01(mi),t1d01(mu),tt01(mu))
        tmaxt2(k) = max(t1d02(mi),tt02(mi),t1d02(mu),tt02(mu))
        tmint1(k) = min(t1d01(mi),tt01(mi),t1d01(mu),tt01(mu))
        tmint2(k) = min(t1d02(mi),tt02(mi),t1d02(mu),tt02(mu))
      endif


      me0 = mcol(msrf(i,  j+1)) - 2
      mw0 = mcol(msrf(i,  j-1)) - 2
      ms0 = mcol(msrf(i+1,j  )) - 2
      mn0 = mcol(msrf(i-1,j  )) - 2

      khuw = kh(msrf(i,j-1))
      khue = kh(msrf(i,j+1))
      khvn = kh(msrf(i-1,j))
      khvs = kh(msrf(i+1,j))

      do k=2,min(kb,khue,khuw,khvn,khvs)
        me = me0 + k
        mw = mw0 + k
        ms = ms0 + k
        mn = mn0 + k

        tmaxt1(k) = max(tmaxt1(k),t1d01(me),tt01(me),t1d01(mw),tt01(mw),       &
                                  t1d01(ms),tt01(ms),t1d01(mn),tt01(mn))
        tmaxt2(k) = max(tmaxt2(k),t1d02(me),tt02(me),t1d02(me),tt02(me),       &
                                  t1d02(ms),tt02(ms),t1d02(mn),tt02(mn))
        tmint1(k) = min(tmint1(k),t1d01(me),tt01(me),t1d01(mw),tt01(mw),       &
                                  t1d01(ms),tt01(ms),t1d01(mn),tt01(mn))
        tmint2(k) = min(tmint2(k),t1d02(me),tt02(me),t1d02(mw),tt02(mw),       &
                                  t1d02(ms),tt02(ms),t1d02(mn),tt02(mn))
      enddo

      klow1 = max(2, min(khuw,khvn,khvs)+1)
      do k=klow1,min(kb,khue)
        me = me0 + k
        tmaxt1(k) = max(tmaxt1(k),t1d01(me),tt01(me))
        tmaxt2(k) = max(tmaxt2(k),t1d02(me),tt02(me))
        tmint1(k) = min(tmint1(k),t1d01(me),tt01(me))
        tmint2(k) = min(tmint2(k),t1d02(me),tt02(me))
      enddo
      klow2 = max(2, min(khue,khvn,khvs)+1)
      do k=klow2,min(kb,khuw)
        mw = mw0 + k
        tmaxt1(k) = max(tmaxt1(k),t1d01(mw),tt01(mw))
        tmaxt2(k) = max(tmaxt2(k),t1d02(mw),tt02(mw))
        tmint1(k) = min(tmint1(k),t1d01(mw),tt01(mw))
        tmint2(k) = min(tmint2(k),t1d02(mw),tt02(mw))
      enddo
      klow3 = max(2, min(khue,khuw,khvn)+1)
      do k=klow3,min(kb,khvs)
        ms = ms0 + k
        tmaxt1(k) = max(tmaxt1(k),t1d01(ms),tt01(ms))
        tmaxt2(k) = max(tmaxt2(k),t1d02(ms),tt02(ms))
        tmint1(k) = min(tmint1(k),t1d01(ms),tt01(ms))
        tmint2(k) = min(tmint2(k),t1d02(ms),tt02(ms))
      enddo
      klow4 = max(2, min(khue,khuw,khvs)+1)
      do k=klow4,min(kb,khvn)
        mn = mn0 + k
        tmaxt1(k) = max(tmaxt1(k),t1d01(mn),tt01(mn))
        tmaxt2(k) = max(tmaxt2(k),t1d02(mn),tt02(mn))
        tmint1(k) = min(tmint1(k),t1d01(mn),tt01(mn))
        tmint2(k) = min(tmint2(k),t1d02(mn),tt02(mn))
      enddo

      do k=2,kb
        mi = mi0 + k

        pin1 = t701(mi)
        !     pin>=0
        if (pin1 > small) then
          qin1 = tmaxt1(k) - tt01(mi)
          !     qin>=0
          t701(mi) = min(one,qin1/pin1)
        else
          t701(mi) = one
        endif
        pout1 = t801(mi)
        !     pout>=0
        if (pout1 > small) then
          qout1 = tt01(mi) - tmint1(k)
          !     qout>=0
          t801(mi) = min(one,qout1/pout1)
        else
          t801(mi) = one
        endif

        pin2 = t702(mi)
        !     pin>=0
        if (pin2 > small) then
          qin2 = tmaxt2(k) - tt02(mi)
          !     qin>=0
          t702(mi) = min(one,qin2/pin2)
        else
          t702(mi) = one
        endif
        pout2 = t802(mi)
        !     pout>=0
        if (pout2 > small) then
          qout2 = tt02(mi) - tmint2(k)
          !     qout>=0
          t802(mi) = min(one,qout2/pout2)
        else
          t802(mi) = one
        endif
      enddo

    enddo
#ifdef _OPENACC
    enddo ! stream
!$ACC WAIT
!$ACC END DATA
#endif

  end subroutine c_rin_rout

  subroutine c_cx_cy_cz (msrf,mcol,ind,n2d,kh,idx)

    ! NOTE: it modifies the global module variables: t4, t5, t6
    ! and it uses t7, t8, t1, t2, t3 

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d
    integer(4), intent(in)    :: msrf(0:,0:), mcol(0:), ind(:,:), kh(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow.c_cx_cy_cz.inc'

    !- local variables ---------------------------------------------------------
    integer(4) :: n, k, kb, i, j, mi, me, ms, md
    integer(4) :: n2dl, n2du, mi0, me0, ms0, kbe, kbs

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,kh,msrf,mcol)                                              &
!$ACC   pcreate(t101,t102,t201,t202,t401,t402,t501,                            &
!$ACC           t502,t601,t602,t701,t702,t801,t802,t301,t302) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif

    ! compute t4, t5, t6

    ! k=1 unrolled
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (mi,n,i,j,me,ms)                                        
#endif
    do n = n2dl,n2du
      i = ind(1,n)
      j = ind(2,n)

      mi = n
      me = msrf(i,  j+1)
      ms = msrf(i+1,j  )

      t601(mi) = zero
      t602(mi) = zero
        
      t401(mi) = min(t701(me),t801(mi))*max(t101(mi),zero)                     &
               + min(t701(mi),t801(me))*min(t101(mi),zero)
      t402(mi) = min(t702(me),t802(mi))*max(t102(mi),zero)                     &
               + min(t702(mi),t802(me))*min(t102(mi),zero)

      t501(mi) = min(t701(ms),t801(mi))*min(t201(mi),zero)                     &
               + min(t701(mi),t801(ms))*max(t201(mi),zero)
      t502(mi) = min(t702(ms),t802(mi))*min(t202(mi),zero)                     &
               + min(t702(mi),t802(ms))*max(t202(mi),zero)
    enddo

#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,md,me,mi,ms,n,kb,i,j,mi0,kbe,me0,kbs,ms0)            
#endif
    do n = n2dl,n2du
      kb = kh(n)
      if (kb <= 1) cycle
      i = ind(1,n)
      j = ind(2,n)

      me = msrf(i,j+1)
      ms = msrf(i+1,j)

      mi0 = mcol(n) - 2
      me0 = mcol(me) - 2
      ms0 = mcol(ms) - 2
      kbe = kh(me)
      kbs = kh(ms)

      mi = n
      md = mcol(n)
      t601(md) = min(t701(md),t801(mi))*min(t301(md),zero)                     &
               + min(t701(mi),t801(md))*max(t301(md),zero)
      t602(md) = min(t702(md),t802(mi))*min(t302(md),zero)                     &
               + min(t702(mi),t802(md))*max(t302(md),zero)

      do k=2,kb-1
        mi = mi0 + k
        md = mi + 1
        t601(md) = min(t701(md),t801(mi))*min(t301(md),zero)                   &
                 + min(t701(mi),t801(md))*max(t301(md),zero)
        t602(md) = min(t702(md),t802(mi))*min(t302(md),zero)                   &
                 + min(t702(mi),t802(md))*max(t302(md),zero)
      enddo

      do k=2,min(kb,kbe)
        mi = mi0 + k
        me = me0 + k
        t401(mi) = min(t701(me),t801(mi))*max(t101(mi),zero)                   &
                 + min(t701(mi),t801(me))*min(t101(mi),zero)
        t402(mi) = min(t702(me),t802(mi))*max(t102(mi),zero)                   &
                 + min(t702(mi),t802(me))*min(t102(mi),zero)
      enddo
      if (kb > kbe) then
        do k=max(2,kbe+1),kb
          mi = mi0 + k
          t401(mi) = min(zero,t801(mi))*max(t101(mi),zero)                     &
                   + min(t701(mi),zero)*min(t101(mi),zero)
          t402(mi) = min(zero,t802(mi))*max(t102(mi),zero)                     &
                   + min(t702(mi),zero)*min(t102(mi),zero)
        enddo
      endif

      do k=2,min(kb,kbs)
        mi = mi0 + k
        ms = ms0 + k
        t501(mi) = min(t701(ms),t801(mi))*min(t201(mi),zero)                   &
                 + min(t701(mi),t801(ms))*max(t201(mi),zero)
        t502(mi) = min(t702(ms),t802(mi))*min(t202(mi),zero)                   &
                 + min(t702(mi),t802(ms))*max(t202(mi),zero)
      enddo
      if (kb > kbs) then
        do k=max(2,kbs+1),kb
          mi = mi0 + k
          t501(mi) = min(zero,t801(mi))*min(t201(mi),zero)                     &
                   + min(t701(mi),zero)*max(t201(mi),zero)
          t502(mi) = min(zero,t802(mi))*min(t202(mi),zero)                     &
                   + min(t702(mi),zero)*max(t202(mi),zero)
        enddo
      endif

    enddo
#ifdef _OPENACC
     enddo
!$ACC WAIT
!$ACC END DATA
#endif

  end subroutine c_cx_cy_cz

  subroutine c_dtz (msrf,mcol,ind,n2d,advecstab,h_new,kh,u,v,w,cosphi,idx)

    ! NOTE: it modifies the global module variables: t3 and it uses t4, t5 ,t6 

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, msrf(0:,0:), mcol(0:), ind(:,:), kh(0:)
    real(8),    intent(in)    :: advecstab, h_new(0:), cosphi(:,0:)
    real(8),    intent(in)    :: u(0:), v(0:), w(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_dtz.inc'

    !- local variables ---------------------------------------------------------
    integer(4) :: n, k, i, j, mi, md, me, mw, ms, mn
    integer(4) :: mi0, me0, ms0, mw0, mn0, mde, mdw, mds, mdn
    integer(4) :: klow, kupper, kb, kbm1, kbnm1, kbsm1, kbwm1, kbem1
    real(8)    :: sx, sy, sz, dtddx, dtddy, hdt, dtddx0, fac
    real(8)    :: s1, s2, s3, s4, s5, s6, s7, s8
    real(8)    :: sx1, sx2, sx3, sx4, sy1, sy2, sy3, sy4 
    real(8)    :: wmn, wmx, wmd
    real(8)    :: szz, szp, szm
    integer(4) :: n2dl, n2du

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,kh,msrf,mcol,w,h_new,u,v,cosphi)                           &
!$ACC   pcreate(t301,t302,t401,t402,t501,t502,t601,t602) 
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif

    hdt    = half*dt
    dtddy  = hdt/dy
    dtddx0 = hdt/dx

    ! compute t1  --------------------------------------------------------------
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   private (k,kbm1,klow,kupper,mde,mdn,mds,mdw,mi,n,s3,s4,s5,s6,          &
!$ACC            s7,s8,sz,wmn,wmx,kb,i,j,dtddx,sx,szp,sy,szm,kbsm1,            &
!$ACC            kbnm1,kbem1,kbwm1,me0,mw0,ms0,mn0,wmd,sx1,sx3,sy1,sy3,        &
!$ACC            sx2,sx4,sy2,sy4,MI0,fac,md,me,mn,ms,mw,s1,s2,szz)
#endif
    do n=n2dl, n2du
      kb = kh(n)

      ! make sure to have nice values throughout
      t301(n) = zero
      t302(n) = zero

      if (kb <= 1) cycle
      i = ind(1,n)
      j = ind(2,n)
      kbm1 = kb-1

      dtddx = dtddx0/cosphi(1,i)

      ! unrol at surface:
      mi = n
      md = mcol(n)

      if (w(md) <= zero) then
        sz  = w(md)*hdt/(h_new(mi) - advecstab)
        szp = -(half+fourthird*sz)
        fac = -half-sz
        t301(md) = fac*t601(mi)
        t302(md) = fac*t602(mi)

        mw = msrf(i,j-1)
        if (u(mw) > zero) then
          sx = dtddx*u(mw)
          s1 = szp*sx
          s2 = (half - twothird*sx)*sx
          t301(md) = t301(md) - s1*(t601(mi)-t601(mw)) - s2*(t401(mi)-t401(mw))
          t302(md) = t302(md) - s1*(t602(mi)-t602(mw)) - s2*(t402(mi)-t402(mw))
        endif

        if (u(mi) < zero) then
          me = msrf(i,j+1)
          sx = dtddx*u(mi)
          s1 = szp*sx
          s2 =-(half + twothird*sx)*sx
          t301(md) = t301(md) - s1*(t601(me)-t601(mi)) - s2*(t401(me)-t401(mi))
          t302(md) = t302(md) - s1*(t602(me)-t602(mi)) - s2*(t402(me)-t402(mi))
        endif

        mn = msrf(i-1,j)
        if (v(mn) < zero) then
          sy = dtddy*v(mn)
          s1 = szp*sy
          s2 =-(half + twothird*sy)*sy
          t301(md) = t301(md) - s1*(t601(mn)-t601(mi)) - s2*(t501(mn)-t501(mi))
          t302(md) = t302(md) - s1*(t602(mn)-t602(mi)) - s2*(t502(mn)-t502(mi))
        endif

        if (v(mi) > zero) then
          ms = msrf(i+1,j)
          sy = dtddy*v(mi)
          s1 = szp*sy
          s2 = (half - twothird*sy)*sy
          t301(md) = t301(md) - s1*(t601(mi)-t601(ms)) - s2*(t501(mi)-t501(ms))
          t302(md) = t302(md) - s1*(t602(mi)-t602(ms)) - s2*(t502(mi)-t502(ms))
        endif

      else ! w > 0.0
        sz  = w(md)*hdt/(h_new(mi) - advecstab)
        szm = (half-fourthird*sz)
        fac = half-sz
        t301(md) = fac*t601(md)
        t302(md) = fac*t602(md)

        mw = mcol(msrf(i,j-1))
        if (u(mw) > zero) then
          sx = dtddx*u(mw)
          s1 = szm*sx
          s2 = (half - twothird*sx)*sx
          t301(md) = t301(md) - s1*(t601(md)-t601(mw)) - s2*(t401(md)-t401(mw))
          t302(md) = t302(md) - s1*(t602(md)-t602(mw)) - s2*(t402(md)-t402(mw))
        endif

        if (u(md) < zero) then
          me = mcol(msrf(i,j+1))
          sx = dtddx*u(md)
          s1 = szm*sx
          s2 =-(half + twothird*sx)*sx
          t301(md) = t301(md) - s1*(t601(me)-t601(md)) - s2*(t401(me)-t401(md))
          t302(md) = t302(md) - s1*(t602(me)-t602(md)) - s2*(t402(me)-t402(md))
        endif

        mn = mcol(msrf(i-1,j))
        if (v(mn) < zero) then
          sy = dtddy*v(mn)
          s1 = szm*sy
          s2 =-(half + twothird*sy)*sy
          t301(md) = t301(md) - s1*(t601(mn)-t601(md)) - s2*(t501(mn)-t501(md))
          t302(md) = t302(md) - s1*(t602(mn)-t602(md)) - s2*(t502(mn)-t502(md))
        endif

        if (v(md) > zero) then
          ms = mcol(msrf(i+1,j))
          sy = dtddy*v(md)
          s1 = szm*sy
          s2 = (half - twothird*sy)*sy
          t301(md) = t301(md) - s1*(t601(md)-t601(ms)) - s2*(t501(md)-t501(ms))
          t302(md) = t302(md) - s1*(t602(md)-t602(ms)) - s2*(t502(md)-t502(ms))
        endif

      endif
      t301(md) = t301(md)*w(md)
      t302(md) = t302(md)*w(md)

      ! k lims:
      kbwm1 = kh(msrf(i,j-1)) - 1
      kbem1 = kh(msrf(i,j+1)) - 1
      kbnm1 = kh(msrf(i-1,j)) - 1
      kbsm1 = kh(msrf(i+1,j)) - 1

      ! index off-sets:
      mi0 = mcol(n) - 2
      me0 = mcol(msrf(i,  j+1)) - 2
      mw0 = mcol(msrf(i  ,j-1)) - 2
      ms0 = mcol(msrf(i+1,j  )) - 2
      mn0 = mcol(msrf(i-1,j  )) - 2

      kupper = min(kbm1,kbwm1,kbem1,kbnm1,kbsm1)

      ! main loop #1, t60[12] terms:
      ! ld/st=10+2*(2*1+5+5)=34, flop=24+2*(2*18+1)=98, CI=98/34 ~ 2.9
      do k=2,kupper
        mi  = mi0 + k
        me  = me0 + k
        mw  = mw0 + k
        ms  = ms0 + k
        mn  = mn0 + k
        md  = mi + 1
        mde = me + 1
        mdw = mw + 1
        mds = ms + 1
        mdn = mn + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sx1 = dtddx*max(u(mw),zero)

        sx2 = dtddx*max(u(mdw),zero)

        sx3 = dtddx*min(u(mi),zero)

        sx4 = dtddx*min(u(md),zero)

        sy1 = dtddy*min(v(mn),zero)

        sy2 = dtddy*min(v(mdn),zero)

        sy3 = dtddy*max(v(mi),zero)

        sy4 = dtddy*max(v(md),zero)

        t301(md) = wmx*(  fac*t601(md)                                         &
                       -szz*sx2*(t601(md)-t601(mdw))                           &
                       -szz*sx4*(t601(mde)-t601(md))                           &
                       -szz*sy2*(t601(mdn)-t601(md))                           &
                       -szz*sy4*(t601(md)-t601(mds)))                          &
                 - wmn*(  fac*t601(mi)                                         &
                       -szz*sx1*(t601(mi)-t601(mw))                            &
                       -szz*sx3*(t601(me)-t601(mi))                            &
                       -szz*sy1*(t601(mn)-t601(mi))                            &
                       -szz*sy3*(t601(mi)-t601(ms)) )
        t302(md) = wmx*(  fac*t602(md)                                         &
                       -szz*sx2*(t602(md)-t602(mdw))                           &
                       -szz*sx4*(t602(mde)-t602(md))                           &
                       -szz*sy2*(t602(mdn)-t602(md))                           &
                       -szz*sy4*(t602(md)-t602(mds)))                          &
                 - wmn*(  fac*t602(mi)                                         &
                       -szz*sx1*(t602(mi)-t602(mw))                            &
                       -szz*sx3*(t602(me)-t602(mi))                            &
                       -szz*sy1*(t602(mn)-t602(mi))                            &
                       -szz*sy3*(t602(mi)-t602(ms)) )
      enddo

      ! main loop #2, t40[12] and t50[12] terms:
      ! ld/st=9+2*(2*1+5+5)=33, flop=5*8+2*(2*14)=96, CI=96/33 ~ 2.9
      do k=2,kupper
        mi  = mi0 + k
        me  = me0 + k
        mw  = mw0 + k
        ms  = ms0 + k
        mn  = mn0 + k
        md  = mi + 1
        mde = me + 1
        mdw = mw + 1
        mds = ms + 1
        mdn = mn + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sx1 = dtddx*max(u(mw),zero)
        s1  = (half - twothird*sx1)*sx1

        sx2 = dtddx*max(u(mdw),zero)
        s2  = (half - twothird*sx2)*sx2

        sx3 = dtddx*min(u(mi),zero)
        s3  = (half + twothird*sx3)*sx3

        sx4 = dtddx*min(u(md),zero)
        s4  = (half + twothird*sx4)*sx4

        sy1 = dtddy*min(v(mn),zero)
        s5  = (half + twothird*sy1)*sy1

        sy2 = dtddy*min(v(mdn),zero)
        s6  = (half + twothird*sy2)*sy2

        sy3 = dtddy*max(v(mi),zero)
        s7  = (half - twothird*sy3)*sy3

        sy4 = dtddy*max(v(md),zero)
        s8  = (half- twothird*sy4)*sy4

        t301(md) = t301(md)                                                    &
                 + wmx*(                                                       &
                       -s2*(t401(md)-t401(mdw))                                &
                       +s4*(t401(mde)-t401(md))                                &
                       +s6*(t501(mdn)-t501(md))                                &
                       -s8*(t501(md)-t501(mds)))                               &
                 - wmn*(                                                       &
                       +s1*(t401(mi)-t401(mw))                                 &
                       -s3*(t401(me)-t401(mi))                                 &
                       -s5*(t501(mn)-t501(mi))                                 &
                       +s7*(t501(mi)-t501(ms)))
        t302(md) = t302(md)                                                    &
                 + wmx*(                                                       &
                       -s2*(t402(md)-t402(mdw))                                &
                       +s4*(t402(mde)-t402(md))                                &
                       +s6*(t502(mdn)-t502(md))                                &
                       -s8*(t502(md)-t502(mds)))                               &
                 - wmn*(                                                       &
                       +s1*(t402(mi)-t402(mw))                                 &
                       -s3*(t402(me)-t402(mi))                                 &
                       -s5*(t502(mn)-t502(mi))                                 &
                       -s7*(t502(mi)-t502(ms)))
      enddo

      ! remainder loops:
      klow = max(2, kupper+1)
      do k=klow,kbm1
        mi  = mi0 + k
        md  = mi + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        t301(md) = wmx*(  fac*t601(md) )                                       &
                 - wmn*(  fac*t601(mi) )
        t302(md) = wmx*(  fac*t602(md) )                                       &
                 - wmn*(  fac*t602(mi) )
      enddo
      do k=klow,min(kbm1,kbwm1)
        mi  = mi0 + k
        mw  = mw0 + k
        md  = mi + 1
        mdw = mw + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sx1 = dtddx*max(u(mw),zero)
        s1  = half - twothird*sx1

        sx2 = dtddx*max(u(mdw),zero)
        s2  = half - twothird*sx2

        t301(md) = t301(md)                                                    &
                 + wmx*(                                                       &
                       -sx2*(szz*(t601(md)-t601(mdw))+s2*(t401(md)-t401(mdw))))&
                 - wmn*(                                                       &
                       -sx1*(szz*(t601(mi)-t601(mw))-s1*(t401(mi)-t401(mw))))
        t302(md) = t302(md)                                                    &
                 + wmx*(                                                       &
                       -sx2*(szz*(t602(md)-t602(mdw))+s2*(t402(md)-t402(mdw))))&
                 - wmn*(                                                       &
                       -sx1*(szz*(t602(mi)-t602(mw))-s1*(t402(mi)-t402(mw))))
      enddo
      ! unroll k=kbw if klow=<kbw<=kbm1
      if (klow <= kbwm1+1 .and. kbwm1+1 <= kbm1) then
        k = kbwm1+1
        mi  = mi0 + k
        mw  = mw0 + k
        md  = mi + 1

        wmd = w(md)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sx1 = dtddx*max(u(mw),zero)
        s1  = half - twothird*sx1

        t301(md) = t301(md)                                                    &
                 - wmn*(                                                       &
                       -sx1*(szz*(t601(mi)-t601(mw))-s1*(t401(mi)-t401(mw))))
        t302(md) = t302(md)                                                    &
                 - wmn*(                                                       &
                       -sx1*(szz*(t602(mi)-t602(mw))-s1*(t402(mi)-t402(mw))))
      endif
      do k=klow,min(kbm1,kbem1)
        mi  = mi0 + k
        me  = me0 + k
        md  = mi + 1
        mde = me + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sx3 = dtddx*min(u(mi),zero)
        s3  = half + twothird*sx3

        sx4 = dtddx*min(u(md),zero)
        s4  = half + twothird*sx4

        t301(md) = t301(md)                                                    &
                 + wmx*(                                                       &
                       -sx4*(szz*(t601(mde)-t601(md))-s4*(t401(mde)-t401(md))))&
                 - wmn*(                                                       &
                       -sx3*(szz*(t601(me)-t601(mi))+s3*(t401(me)-t401(mi))))
        t302(md) = t302(md)                                                    &
                 + wmx*(                                                       &
                       -sx4*(szz*(t602(mde)-t602(md))-s4*(t402(mde)-t402(md))))&
                 - wmn*(                                                       &
                       -sx3*(szz*(t602(me)-t602(mi))+s3*(t402(me)-t402(mi))))
      enddo
      ! unroll k=kbe if klow=<kbe<=kbm1
      if (klow <= kbem1+1 .and. kbem1+1 <= kbm1) then
        k = kbem1+1
        mi  = mi0 + k
        me  = me0 + k
        md  = mi + 1

        wmd = w(md)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sx3 = dtddx*min(u(mi),zero)
        s3  = half + twothird*sx3

        t301(md) = t301(md)                                                    &
                 - wmn*(                                                       &
                       -sx3*(szz*(t601(me)-t601(mi))+s3*(t401(me)-t401(mi))))
        t302(md) = t302(md)                                                    &
                 - wmn*(                                                       &
                       -sx3*(szz*(t602(me)-t602(mi))+s3*(t402(me)-t402(mi))))
      endif
      do k=klow,min(kbm1,kbnm1)
        mi  = mi0 + k
        mn  = mn0 + k
        md  = mi + 1
        mdn = mn + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sy1 = dtddy*min(v(mn),zero)
        s5  = half + twothird*sy1

        sy2 = dtddy*min(v(mdn),zero)
        s6  = half + twothird*sy2

        t301(md) = t301(md)                                                    &
                 + wmx*(                                                       &
                       -sy2*(szz*(t601(mdn)-t601(md))-s6*(t501(mdn)-t501(md))))&
                 - wmn*(                                                       &
                       -sy1*(szz*(t601(mn)-t601(mi))+s5*(t501(mn)-t501(mi))))
        t302(md) = t302(md)                                                    &
                 + wmx*(                                                       &
                       -sy2*(szz*(t602(mdn)-t602(md))-s6*(t502(mdn)-t502(md))))&
                 - wmn*(                                                       &
                       -sy1*(szz*(t602(mn)-t602(mi))+s5*(t502(mn)-t502(mi))))
      enddo
      ! unroll k=kbn if klow=<kbn<=kbm1
      if (klow <= kbnm1+1 .and. kbnm1+1 <= kbm1) then
        k = kbnm1+1
        mi  = mi0 + k
        mn  = mn0 + k
        md  = mi + 1

        wmd = w(md)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sy1 = dtddy*min(v(mn),zero)
        s5  = half + twothird*sy1

        t301(md) = t301(md)                                                    &
                 - wmn*(                                                       &
                       -sy1*(szz*(t601(mn)-t601(mi))+s5*(t501(mn)-t501(mi))))
        t302(md) = t302(md)                                                    &
                 - wmn*(                                                       &
                       -sy1*(szz*(t602(mn)-t602(mi))+s5*(t502(mn)-t502(mi))))
      endif
      do k=klow,min(kbm1,kbsm1)
        mi  = mi0 + k
        ms  = ms0 + k
        md  = mi + 1
        mds = ms + 1

        wmd = w(md)
        wmx = max(wmd,zero)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sy3 = dtddy*max(v(mi),zero)
        s7  = half - twothird*sy3

        sy4 = dtddy*max(v(md),zero)
        s8  = half- twothird*sy4

        t301(md) = t301(md)                                                    &
                 + wmx*(                                                       &
                       -sy4*(szz*(t601(md)-t601(mds))+s8*(t501(md)-t501(mds))))&
                 - wmn*(                                                       &
                       -sy3*(szz*(t601(mi)-t601(ms))-s7*(t501(mi)-t501(ms))))
        t302(md) = t302(md)                                                    &
                 + wmx*(                                                       &
                       -sy4*(szz*(t602(md)-t602(mds))+s8*(t502(md)-t502(mds))))&
                 - wmn*(                                                       &
                       -sy3*(szz*(t602(mi)-t602(ms))-s7*(t502(mi)-t502(ms))))
      enddo
      ! unroll k=kbs if klow=<kbs<=kbm1
      if (klow <= kbsm1+1 .and. kbsm1+1 <= kbm1) then
        k = kbsm1+1
        mi  = mi0 + k
        ms  = ms0 + k
        md  = mi + 1

        wmd = w(md)
        wmn = min(wmd,zero)

        sz  = abs(wmd)*hdt/h_new(mi)
        szz = half - fourthird*sz
        fac = half - sz

        sy3 = dtddy*max(v(mi),zero)
        s7  = half - twothird*sy3

        t301(md) = t301(md)                                                    &
                 - wmn*(                                                       &
                       -sy3*(szz*(t601(mi)-t601(ms))-s7*(t501(mi)-t501(ms))))
        t302(md) = t302(md)                                                    &
                 - wmn*(                                                       &
                       -sy3*(szz*(t602(mi)-t602(ms))-s7*(t502(mi)-t502(ms))))
      endif

    enddo
#ifdef _OPENACC
    enddo ! stream
!$ACC WAIT
!$ACC END DATA
#endif
    
  end subroutine c_dtz

  subroutine c_dty (msrf,mcol,ind,n2d,advecstab,h_new,kh,khu,khv,u,v,w,hy,     &
                    cosphi,idx)

    ! NOTE: it modifies the global module variables: t2 and it uses t4, t5 ,t6 

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, msrf(0:,0:), mcol(0:), ind(:,:)
    integer(4), intent(in)    :: kh(0:), khu(0:), khv(0:)
    real(8),    intent(in)    :: advecstab, h_new(0:), cosphi(:,0:)
    real(8),    intent(in)    :: u(0:), v(0:), w(0:), hy(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_dty.inc'

    !- local variables ---------------------------------------------------------
    integer(4) :: n, k, kb, i, j, kbu, kbv, mi, md, me, mw, ms
    integer(4) :: mse, mu, msw, mds, mus
    integer(4) :: mi0, me0, ms0, mw0, msw0, mse0, klow, kupper
    integer(4) :: kbs, kbuw, kbus, kbusw
    real(8)    :: sx, sy, sz, dtddx, dtddy, dtddz, hdt, dtddx0, fac
    real(8)    :: s1, s2, s3, s4, s5, s6, s7, s8 
    real(8)    :: hyv, vmn, vmx, sx1, sx2, sx3, sx4
    real(8)    :: sz1, sz2, sz3, sz4, syy, syp, sym
    integer(4) :: n2dl, n2du

#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(khv,kh,ind,msrf,mcol,cosphi,v,u,w,h_new,hy,khu)                &
!$ACC   pcreate(t201,t202,t401,t402,t501,t502,t601,t602)
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif

    hdt    = half*dt
    dtddy  = hdt/dy
    dtddx0 = hdt/dx

    ! compute t1  --------------------------------------------------------------
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC  async(is)                                                               &
!$ACC   private (k,klow,kupper,mi,mu,n,sy,kb,kbv,i,j,mi0,dtddx,sx,syp,         &
!$ACC            dtddz,sz,sym,kbs,hyv,kbusw,kbus,kbuw,kbu,me0,mw0,ms0,         &
!$ACC            msw0,mse0,sx2,sx1,sz1,sz2,syy,sx4,sx3,sz3,sz4,s4,s3,          &
!$ACC            s7,s8,s5,s6,fac,md,mds,me,ms,mse,msw,mus,mw,s1,s2,vmn,vmx)
#endif
    do n=n2dl, n2du
      kbv = khv(n)
      kb  = kh(n)
      i   = ind(1,n)
      j   = ind(2,n)

      ! make sure to have nice values throughout
      if (kb > 1 .and. kb > kbv) then
        mi0 = mcol(n) - 2
        t201(mi0+max(2,kbv+1):mi0+kb) = zero
        t202(mi0+max(2,kbv+1):mi0+kb) = zero
      endif
      if (kbv == 0) then
        t201(n) = zero
        t202(n) = zero
        cycle
      endif

      kbs   = kh(msrf(i+1,j))

      dtddx = dtddx0/cosphi(1,i)

      ! unroll k=1:
      if (v(n) < zero) then

        sy  = dtddy*v(n)
        syp =-(half+fourthird*sy)
        fac = -half-sy
        t201(n) = fac*t501(n)
        t202(n) = fac*t502(n)

        mw = msrf(i,j-1)
        if (u(mw) > zero) then
          sx = dtddx*u(mw)
          s1 = syp*sx
          s2 = (half- twothird*sx)*sx
          t201(n) = t201(n) - s1*(t501(n)-t501(mw)) - s2*(t401(n)-t401(mw))
          t202(n) = t202(n) - s1*(t502(n)-t502(mw)) - s2*(t402(n)-t402(mw))
        endif

        if (u(n) < zero) then
          sx = dtddx*u(n)
          s1 = syp*sx
          s2 =-(half+ twothird*sx)*sx
          me = msrf(i,j+1)
          t201(n) = t201(n) - s1*(t501(me)-t501(n)) - s2*(t401(me)-t401(n))
          t202(n) = t202(n) - s1*(t502(me)-t502(n)) - s2*(t402(me)-t402(n))
        endif

        if (kb > 1) then
          md = mcol(n)
          if (w(md) > zero) then
            dtddz = hdt/(h_new(n) - advecstab)
            sz = dtddz*w(md)
            s1 = syp*sz
            s2 = (half- twothird*sz)*sz
            t201(n) = t201(n) - s1*(t501(n)-t501(md)) - s2*(t601(n)-t601(md))
            t202(n) = t202(n) - s1*(t502(n)-t502(md)) - s2*(t602(n)-t602(md))
          endif
        endif

      else ! v > 0.0

        sy  = dtddy*v(n)
        sym = (half-fourthird*sy)
        fac = half-sy
        ms = msrf(i+1,j)
        t201(n) = fac*t501(ms)
        t202(n) = fac*t502(ms)

        msw = msrf(i+1,j-1)
        if (u(msw) > zero) then
          sx = dtddx*u(msw)
          s1 = sym*sx
          s2 = (half- twothird*sx)*sx
          t201(n) = t201(n) - s1*(t501(ms)-t501(msw)) - s2*(t401(ms)-t401(msw))
          t202(n) = t202(n) - s1*(t502(ms)-t502(msw)) - s2*(t402(ms)-t402(msw))
        endif

        if (u(ms) < zero) then
          sx = dtddx*u(ms)
          s1 = sym*sx
          s2 =-(half+ twothird*sx)*sx
          mse = msrf(i+1,j+1)
          t201(n) = t201(n) - s1*(t501(mse)-t501(ms)) - s2*(t401(mse)-t401(ms))
          t202(n) = t202(n) - s1*(t502(mse)-t502(ms)) - s2*(t402(mse)-t402(ms))
        endif

        if (kbs > 1) then
          mds = mcol(msrf(i+1,j))
          if (w(mds) > zero) then
            dtddz = hdt/(h_new(n) - advecstab)
            sz = dtddz*w(mds)
            s1 = sym*sz
            s2 = (half- twothird*sz)*sz
            t201(n) =t201(n) - s1*(t501(ms)-t501(mds)) - s2*(t601(ms)-t601(mds))
            t202(n) =t202(n) - s1*(t502(ms)-t502(mds)) - s2*(t602(ms)-t602(mds))
          endif
        endif

      endif

      hyv = hy(n)*v(n)
      t201(n) = t201(n)*hyv
      t202(n) = t202(n)*hyv
      ! end of surface layer treatement.
      if (kbv < 2) cycle

      ! k lims:
      kbu   = khu(n)
      kbuw  = khu(msrf(i,j-1))
      kbus  = khu(msrf(i+1,j)) 
      kbusw = khu(msrf(i+1,j-1)) 

      ! index off-sets:
      mi0  = mcol(n) - 2
      me0  = mcol(msrf(i,  j+1)) - 2
      mw0  = mcol(msrf(i  ,j-1)) - 2
      ms0  = mcol(msrf(i+1,j  )) - 2
      mse0 = mcol(msrf(i+1,j+1)) - 2
      msw0 = mcol(msrf(i+1,j-1)) - 2

      ! k=2 unrolled:
      mi = mi0 + 2
      hyv = hy(mi)*v(mi)
      if (hyv < zero) then
        sy  = dtddy*v(mi)
        syp = -(half+fourthird*sy)
        fac = -half-sy
        t201(mi) = fac*t501(mi)
        t202(mi) = fac*t502(mi)
        mw = mcol(msrf(i,j-1))
        if (u(mw) > zero) then
          sx = dtddx*u(mw)
          s1 = syp*sx
          s2 = (half- twothird*sx)*sx
          t201(mi) = t201(mi) - s1*(t501(mi)-t501(mw)) - s2*(t401(mi)-t401(mw))
          t202(mi) = t202(mi) - s1*(t502(mi)-t502(mw)) - s2*(t402(mi)-t402(mw))
        endif
        if (u(mi) < zero) then
          me = mcol(msrf(i,j+1))
          sx = dtddx*u(mi)
          s1 = syp*sx
          s2 =-(half+ twothird*sx)*sx
          t201(mi) = t201(mi) - s1*(t501(me)-t501(mi)) - s2*(t401(me)-t401(mi))
          t202(mi) = t202(mi) - s1*(t502(me)-t502(mi)) - s2*(t402(me)-t402(mi))
        endif
        if (w(mi) < zero) then
          mu = n
          dtddz = hdt/h_new(mi)
          sz = dtddz*w(mi)
          s1 = syp*sz
          s2 =-(half+ twothird*sz)*sz
          t201(mi) = t201(mi) - s1*(t501(mu)-t501(mi)) - s2*(t601(mu)-t601(mi))
          t202(mi) = t202(mi) - s1*(t502(mu)-t502(mi)) - s2*(t602(mu)-t602(mi))
        endif
        if (kb > 2) then
          md = mi + 1
          if (w(md) > zero) then
            dtddz = hdt/h_new(mi)
            sz = dtddz*w(md)
            s1 = syp*sz
            s2 = (half- twothird*sz)*sz
            t201(mi) = t201(mi) - s1*(t501(mi)-t501(md))- s2*(t601(mi)-t601(md))
            t202(mi) = t202(mi) - s1*(t502(mi)-t502(md))- s2*(t602(mi)-t602(md))
          endif
        endif
      else ! v > 0.0
        sy  = dtddy*v(mi)
        sym = (half-fourthird*sy)
        fac = half-sy
        ms = mcol(msrf(i+1,j))
        t201(mi) = fac*t501(ms)
        t202(mi) = fac*t502(ms)
        msw = mcol(msrf(i+1,j-1))
        if (u(msw) > zero) then
          sx = dtddx*u(msw)
          s1 = sym*sx
          s2 = (half- twothird*sx)*sx
          t201(mi) = t201(mi) - s1*(t501(ms)-t501(msw))- s2*(t401(ms)-t401(msw))
          t202(mi) = t202(mi) - s1*(t502(ms)-t502(msw))- s2*(t402(ms)-t402(msw))
        endif
        if (u(ms) < zero) then
          sx = dtddx*u(ms)
          s1 = sym*sx
          s2 =-(half+ twothird*sx)*sx
          mse = mcol(msrf(i+1,j+1))
          t201(mi) = t201(mi) - s1*(t501(mse)-t501(ms))- s2*(t401(mse)-t401(ms))
          t202(mi) = t202(mi) - s1*(t502(mse)-t502(ms))- s2*(t402(mse)-t402(ms))
        endif
        if (w(ms) < zero) then
          dtddz = hdt/h_new(mi)
          sz = dtddz*w(ms)
          s1 = sym*sz
          s2 =-(half+ twothird*sz)*sz
          mus = msrf(i+1,j)
          t201(mi) = t201(mi) - s1*(t501(mus)-t501(ms))- s2*(t601(mus)-t601(ms))
          t202(mi) = t202(mi) - s1*(t502(mus)-t502(ms))- s2*(t602(mus)-t602(ms))
        endif
        if (kbs > 2) then
          mds = ms + 1
          if (w(mds) > zero) then
            dtddz = hdt/h_new(mi)
            sz = dtddz*w(mds)
            s1 = sym*sz
            s2 = (half- twothird*sz)*sz
            t201(mi)= t201(mi) - s1*(t501(ms)-t501(mds))-s2*(t601(ms)-t601(mds))
            t202(mi)= t202(mi) - s1*(t502(ms)-t502(mds))-s2*(t602(ms)-t602(mds))
          endif
        endif
      endif
      t201(mi) = t201(mi)*hyv
      t202(mi) = t202(mi)*hyv

      kupper = min(kbu,kbv,kbuw,kbus,kbusw,kb-1,kbs-1)

      ! main loop #1, t50[12]
      do k=3,kupper
        mi  = mi0 + k
        me  = me0 + k
        mw  = mw0 + k
        ms  = ms0 + k
        msw = msw0 + k
        mse = mse0 + k
        mu  = mi - 1
        md  = mi + 1
        mus = ms - 1
        mds = ms + 1
        
        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        sx1 = dtddx*max(u(mw),zero)

        sx2 = dtddx*min(u(mi),zero)

        sx3 = dtddy*max(u(msw),zero)

        sx4 = dtddy*min(u(ms),zero)

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)

        sz2 = dtddz*max(w(md),zero)

        sz3 = dtddz*min(w(ms),zero)

        sz4 = dtddz*max(w(mds),zero)

        t201(mi) = vmx*( fac*t501(ms)                                          &
                       -syy*( sx3*(t501(ms)-t501(msw))                         &
                             +sx4*(t501(mse)-t501(ms))                         &
                             +sz3*(t501(mus)-t501(ms))                         &
                             +sz4*(t501(ms)-t501(mds))))                       &
                 - vmn*( fac*t501(mi)                                          &
                       -syy*( sx1*(t501(mi)-t501(mw))                          &
                             +sx2*(t501(me)-t501(mi))                          &
                             +sz1*(t501(mu)-t501(mi))                          &
                             +sz2*(t501(mi)-t501(md))))
        t202(mi) = vmx*( fac*t502(ms)                                          &
                       -syy*( sx3*(t502(ms)-t502(msw))                         &
                             +sx4*(t502(mse)-t502(ms))                         &
                             +sz3*(t502(mus)-t502(ms))                         &
                             +sz4*(t502(ms)-t502(mds))))                       &
                 - vmn*( fac*t502(mi)                                          &
                       -syy*( sx1*(t502(mi)-t502(mw))                          &
                             +sx2*(t502(me)-t502(mi))                          &
                             +sz1*(t502(mu)-t502(mi))                          &
                             +sz2*(t502(mi)-t502(md))))
      enddo

      ! main loop #2, t40[12] and t60[12]
      do k=3,kupper
        mi  = mi0 + k
        me  = me0 + k
        mw  = mw0 + k
        ms  = ms0 + k
        msw = msw0 + k
        mse = mse0 + k
        mu  = mi - 1
        md  = mi + 1
        mus = ms - 1
        mds = ms + 1
        
        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sx1 = dtddx*max(u(mw),zero)
        s1  = (half - twothird*sx1)*sx1

        sx2 = dtddx*min(u(mi),zero)
        s2  = (half + twothird*sx2)*sx2

        sx3 = dtddy*max(u(msw),zero)
        s3  = (half - twothird*sx3)*sx3

        sx4 = dtddy*min(u(ms),zero)
        s4  = (half + twothird*sx4)*sx4

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)
        s5  = (half + twothird*sz1)*sz1

        sz2 = dtddz*max(w(md),zero)
        s6  = (half - twothird*sz2)*sz2

        sz3 = dtddz*min(w(ms),zero)
        s7  = (half + twothird*sz3)*sz3

        sz4 = dtddz*max(w(mds),zero)
        s8  = (half - twothird*sz4)*sz4

        t201(mi) = t201(mi)                                                    &
                 + vmx*(                                                       &
                       -s3*(t401(ms)+t401(msw))                                &
                       +s4*(t401(mse)-t401(ms))                                &
                       +s7*(t601(mus)-t601(ms))                                &
                       -s8*(t601(ms)-t601(mds)))                               &
                 - vmn*(                                                       &
                       +s1*(t401(mi)+t401(mw))                                 &
                       -s2*(t401(me)-t401(mi))                                 &
                       -s5*(t601(mu)-t601(mi))                                 &
                       +s6*(t601(mi)-t601(md)) )
        t202(mi) = t202(mi)                                                    &
                 + vmx*(                                                       &
                       -s3*(t402(ms)+t402(msw))                                &
                       +s4*(t402(mse)-t402(ms))                                &
                       +s7*(t602(mus)-t602(ms))                                &
                       -s8*(t602(ms)-t602(mds)))                               &
                 - vmn*(                                                       &
                       +s1*(t402(mi)+t402(mw))                                 &
                       -s2*(t402(me)-t402(mi))                                 &
                       -s5*(t602(mu)-t602(mi))                                 &
                       +s6*(t602(mi)-t602(md)) )
      enddo

      ! remainder loops:
      klow = max(3, kupper+1)
      do k=klow,min(kbv,kb-1)
        mi  = mi0 + k
        ms  = ms0 + k
        mu  = mi - 1
        md  = mi + 1
        mus = ms - 1
        
        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)
        s5  = half + twothird*sz1

        sz2 = dtddz*max(w(md),zero)
        s6  = half - twothird*sz2

        sz3 = dtddz*min(w(ms),zero)
        s7  = half + twothird*sz3

        t201(mi) = vmx*( fac*t501(ms)                                          &
                       -sz3*(syy*(t501(mus)-t501(ms))-s7*(t601(mus)-t601(ms))))&
                 - vmn*( fac*t501(mi)                                          &
                       -sz1*(syy*(t501(mu)-t501(mi)) +s5*(t601(mu)-t601(mi)) ) &
                       -sz2*(syy*(t501(mi)-t501(md)) -s6*(t601(mi)-t601(md)) ))
        t202(mi) = vmx*( fac*t502(ms)                                          &
                       -sz3*(syy*(t502(mus)-t502(ms))-s7*(t602(mus)-t602(ms))))&
                 - vmn*( fac*t502(mi)                                          &
                       -sz1*(syy*(t502(mu)-t502(mi)) +s5*(t602(mu)-t602(mi)) ) &
                       -sz2*(syy*(t502(mi)-t502(md)) -s6*(t602(mi)-t602(md)) ))
      enddo
      do k=klow,min(kbv,kb-1,kbs-1)
        mi  = mi0 + k
        ms  = ms0 + k
        mds = ms + 1

        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        dtddz = hdt/h_new(mi)

        sz4 = dtddz*max(w(mds),zero)
        s8  = half - twothird*sz4

        t201(mi) = t201(mi)                                                    &
                 + vmx*(                                                       &
                       -sz4*(syy*(t501(ms)-t501(mds))+s8*(t601(ms)-t601(mds))))
        t202(mi) = t202(mi)                                                    &
                 + vmx*(                                                       &
                       -sz4*(syy*(t502(ms)-t502(mds))+s8*(t602(ms)-t602(mds))))
      enddo
      do k=klow,min(kbv,kb-1,kbu)
        mi  = mi0 + k
        me  = me0 + k

        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        sx2 = dtddx*min(u(mi),zero)
        s2  = half + twothird*sx2

        t201(mi) = t201(mi)                                                    &
                 - vmn*(                                                       &
                       -sx2*(syy*(t501(me)-t501(mi)) +s2*(t401(me)-t401(mi)) ))
        t202(mi) = t202(mi)                                                    &
                 - vmn*(                                                       &
                       -sx2*(syy*(t502(me)-t502(mi)) +s2*(t402(me)-t402(mi)) ))
      enddo
      do k=klow,min(kbv,kb-1,kbuw)
        mi  = mi0 + k
        mw  = mw0 + k

        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        sx1 = dtddx*max(u(mw),zero)
        s1  = half - twothird*sx1

        t201(mi) = t201(mi)                                                    &
                 - vmn*(                                                       &
                       -sx1*(syy*(t501(mi)-t501(mw)) -s1*(t401(mi)+t401(mw)) ))
        t202(mi) = t202(mi)                                                    &
                 - vmn*(                                                       &
                       -sx1*(syy*(t502(mi)-t502(mw)) -s1*(t402(mi)+t402(mw)) ))
      enddo
      do k=klow,min(kbv,kb-1,kbus)
        mi  = mi0 + k
        ms  = ms0 + k
        mse = mse0 + k
 
        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        sx4 = dtddy*min(u(ms),zero)
        s4  = half + twothird*sx4

        t201(mi) = t201(mi)                                                    &
                 + vmx*( fac*t501(ms)                                          &
                       -sx4*(syy*(t501(mse)-t501(ms))-s4*(t401(mse)-t401(ms))))
        t202(mi) = t202(mi)                                                    &
                 + vmx*( fac*t502(ms)                                          &
                       -sx4*(syy*(t502(mse)-t502(ms))-s4*(t402(mse)-t402(ms))))
      enddo
      do k=klow,min(kbv,kb-1,kbusw)
        mi  = mi0 + k
        ms  = ms0 + k
        msw = msw0 + k
 
        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        sx3 = dtddy*max(u(msw),zero)
        s3  = half - twothird*sx3

        t201(mi) = t201(mi)                                                    &
                 + vmx*(                                                       &
                       -sx3*(syy*(t501(ms)-t501(msw))+s3*(t401(ms)+t401(msw))))
        t202(mi) = t202(mi)                                                    &
                 + vmx*(                                                       &
                       -sx3*(syy*(t502(ms)-t502(msw))+s3*(t402(ms)+t402(msw))))
      enddo


      ! bottom unrolled:
      if (kb == kbv .and. kb > 2) then
        k = kbv
        mi  = mi0 + k
        ms  = ms0 + k
        if (kh(msrf(i,j-1)) >= k) then
          mw  = mw0 + k
        else
          mw  = 0
        endif
        if (kh(msrf(i+1,j-1)) >= k) then
          msw = msw0 + k
        else
          msw = 0
        endif
        if (kh(msrf(i+1,j+1)) >= k) then
          mse = mse0 + k
        else
          mse = 0
        endif
        if (kh(msrf(i,j+1)) >= k) then
          me  = me0 + k
        else
          me  = 0
        endif
        mu  = mi - 1
        mus = ms - 1
        if (kh(msrf(i+1,j)) > k) then
          mds = ms + 1
        else
          mds = 0
        endif

        hyv = hy(mi)*v(mi)
        vmx = max(hyv,zero)
        vmn = min(hyv,zero)

        sy  = dtddy*abs(v(mi))
        syy = half - fourthird*sy
        fac = half - sy

        sx1 = dtddx*max(u(mw),zero)
        s1  = half - twothird*sx1

        sx2 = dtddx*min(u(mi),zero)
        s2  = half + twothird*sx2

        sx3 = dtddy*max(u(msw),zero)
        s3  = half - twothird*sx3

        sx4 = dtddy*min(u(ms),zero)
        s4  = half + twothird*sx4

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)
        s5  = half + twothird*sz1

        sz3 = dtddz*min(w(ms),zero)
        s7  = half + twothird*sz3

        sz4 = dtddz*max(w(mds),zero)
        s8  = half - twothird*sz4

        t201(mi) = vmx*( fac*t501(ms)                                          &
                       -sx3*(syy*(t501(ms)-t501(msw))+s3*(t401(ms)+t401(msw))) &
                       -sx4*(syy*(t501(mse)-t501(ms))-s4*(t401(mse)-t401(ms))) &
                       -sz3*(syy*(t501(mus)-t501(ms))-s7*(t601(mus)-t601(ms))) &
                       -sz4*(syy*(t501(ms)-t501(mds))+s8*(t601(ms)-t601(mds))))&
                 - vmn*( fac*t501(mi)                                          &
                       -sx1*(syy*(t501(mi)-t501(mw)) -s1*(t401(mi)+t401(mw)) ) &
                       -sx2*(syy*(t501(me)-t501(mi)) +s2*(t401(me)-t401(mi)) ) &
                       -sz1*(syy*(t501(mu)-t501(mi)) +s5*(t601(mu)-t601(mi)) ))
        t202(mi) = vmx*( fac*t502(ms)                                          &
                       -sx3*(syy*(t502(ms)-t502(msw))+s3*(t402(ms)+t402(msw))) &
                       -sx4*(syy*(t502(mse)-t502(ms))-s4*(t402(mse)-t402(ms))) &
                       -sz3*(syy*(t502(mus)-t502(ms))-s7*(t602(mus)-t602(ms))) &
                       -sz4*(syy*(t502(ms)-t502(mds))+s8*(t602(ms)-t602(mds))))&
                 - vmn*( fac*t502(mi)                                          &
                       -sx1*(syy*(t502(mi)-t502(mw)) -s1*(t402(mi)+t402(mw)) ) &
                       -sx2*(syy*(t502(me)-t502(mi)) +s2*(t402(me)-t402(mi)) ) &
                       -sz1*(syy*(t502(mu)-t502(mi)) +s5*(t602(mu)-t602(mi)) ))
      endif
    enddo
#ifdef _OPENACC
    enddo
!$ACC WAIT
!$ACC END DATA
#endif
    
  end subroutine c_dty

  subroutine c_dtx (msrf,mcol,ind,n2d,advecstab,h_new,kh,khu,khv,u,v,w,hx,     &
                    cosphi,idx)

    ! NOTE: it modifies the global module variables: t1 and it uses t4, t5 ,t6 

    !- modules -----------------------------------------------------------------
    use dmi_omp, only : domp_get_domain

    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: n2d, msrf(0:,0:), mcol(0:), ind(:,:)
    integer(4), intent(in)    :: kh(0:), khu(0:), khv(0:)
    real(8),    intent(in)    :: advecstab, h_new(0:), cosphi(:,0:)
    real(8),    intent(in)    :: u(0:), v(0:), w(0:), hx(0:)
    integer(4), intent(inout) :: idx(1:,1:)

    include 'tflow_dtx.inc'

    !- local variables ---------------------------------------------------------
    integer(4) :: n, k, kb, i, j, kbu, kbv, mi, md, me, mn, ms
    integer(4) :: mse, mne, mde, mu, mue
    integer(4) :: mi0, me0, ms0, mn0, mne0, mse0, klow, kupper
    integer(4) :: kbe, kbvn, kbve, kbvne
    real(8)    :: sx, sy, sz, dtddx, dtddy, dtddz, hdt, dtddx0, fac
    real(8)    :: s1, s2, s3, s4, s5,s6, s7, s8 
    real(8)    :: hxu, umn, umx, sy1, sy2, sy3, sy4
    real(8)    :: sz1, sz2, sz3, sz4, sxx, sxp, sxm
    integer(4) :: n2dl, n2du
#ifdef _OPENACC
    integer(4) :: is
!$ACC data                                                                     &
!$ACC   present(ind,kh,msrf,mcol,w,h_new,u,v,cosphi)                           &
!$ACC   pcreate(t101,t102,t401,t402,t501,t502,t601,t602)
    do is = 1, 16
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, s_idx,is)
#else
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)
#endif

    hdt    = half*dt
    dtddy  = hdt/dy
    dtddx0 = hdt/dx

    ! compute t1  --------------------------------------------------------------
#ifdef _OPENACC
!$ACC parallel loop vector_length(32)                                          &
!$ACC   async(is)                                                              &
!$ACC   present(khu,kh,ind,msrf,mcol,cosphi,v,u,w,h_new,hx,khv)                &
!$ACC   private (fac,k,klow,kupper,mi,mu,n,sx,kb,kbu,i,j,mi0,dtddx,sy,         &
!$ACC            dtddz,sz,sxp,kbe,hxu,kbvne,kbve,kbvn,kbv,me0,mn0,ms0,         &
!$ACC            mne0,mse0,sy4,sy2,sz2,sz4,sxx,umn,sy3,sy1,sz1,sz3,umx,        &
!$ACC            s3,s5,s7,s4,s6,s8,md,mde,me,mn,mne,ms,mse,mue,s1,s2,sxm)
#endif
    do n=n2dl, n2du
      kbu = khu(n)
      kb  = kh(n)
      i   = ind(1,n)
      j   = ind(2,n)

      ! make sure to have nice values throughout
      if (kb > 1 .and. kb > kbu) then
        mi0 = mcol(n) - 2
        t101(mi0+max(2,kbu+1):mi0+kb) = zero
        t102(mi0+max(2,kbu+1):mi0+kb) = zero
      endif
      if (kbu == 0) then
        t101(n) = zero
        t102(n) = zero
        cycle
      endif

      kbe   = kh(msrf(i,j+1))

      dtddx = dtddx0/cosphi(1,i)

      ! unroll k=1:
      if (u(n) >= zero) then
    
        sx  = dtddx*u(n)
        sxm = (half-fourthird*sx)
        fac = half-sx
        t101(n) = fac*t401(n)
        t102(n) = fac*t402(n)
    
        mn = msrf(i-1,j)
        if (v(mn) < zero) then
          sy = dtddy*v(mn)
          s1 = sxm*sy
          s2 = (half+ twothird*sy)*sy
          t101(n) = t101(n) - s1*(t401(mn)-t401(n)) + s2*(t501(mn)-t501(n))
          t102(n) = t102(n) - s1*(t402(mn)-t402(n)) + s2*(t502(mn)-t502(n))
        endif
    
        if (v(n) > zero) then
          ms = msrf(i+1,j)
          sy = dtddy*v(n)
          s1 = sxm*sy
          s2 = (half- twothird*sy)*sy
          t101(n) = t101(n) - s1*(t401(n)-t401(ms)) - s2*(t501(n)-t501(ms))
          t101(n) = t102(n) - s1*(t402(n)-t402(ms)) - s2*(t502(n)-t502(ms))
        endif
    
        if (kb > 1) then
          md = mcol(n)
          if (w(md) > zero) then
            dtddz = hdt/(h_new(n) - advecstab)
            sz = dtddz*w(md)
            s1 = sxm*sz
            s2 = (half- twothird*sz)*sz
            t101(n) = t101(n) - s1*(t401(n)-t401(md)) - s2*(t601(n)-t601(md))
            t102(n) = t102(n) - s1*(t402(n)-t402(md)) - s2*(t602(n)-t602(md))
          endif
        endif
    
      else ! u < 0.0
    
        sx  = dtddx*u(n)
        sxp = (half+fourthird*sx)
        fac = -half-sx
        me  = msrf(i,j+1)
        t101(n) = fac*t401(me)
        t102(n) = fac*t402(me)
    
        mne = msrf(i-1,j+1)
        if (v(mne) < zero) then
          sy = dtddy*v(mne)
          s1 = sxp*sy
          s2 = (half+ twothird*sy)*sy
          t101(n) = t101(n) + s1*(t401(mne)-t401(me)) - s2*(t501(mne)-t501(me))
          t102(n) = t102(n) + s1*(t402(mne)-t402(me)) - s2*(t502(mne)-t502(me))
        endif
    
        if (v(me) > zero) then
          mse = msrf(i+1,j+1)
          sy = dtddy*v(me)
          s1 = sxp*sy
          s2 = (half- twothird*sy)*sy
          t101(n) = t101(n) + s1*(t401(me)-t401(mse)) + s2*(t501(me)-t501(mse))
          t102(n) = t102(n) + s1*(t402(me)-t402(mse)) + s2*(t502(me)-t502(mse))
        endif
    
        if (kbe > 1) then
          mde = mcol(msrf(i,j+1))
          if (w(mde) > zero) then
            dtddz = hdt/(h_new(n) - advecstab)
            sz = dtddz*w(mde)
            s1 = sxp*sz
            s2 = (half- twothird*sz)*sz
            t101(n) =t101(n) + s1*(t401(me)-t401(mde)) + s2*(t601(me)-t601(mde))
            t102(n) =t102(n) + s1*(t402(me)-t402(mde)) + s2*(t602(me)-t602(mde))
          endif
        endif

      endif

      hxu = hx(n)*u(n)
      t101(n) = t101(n)*hxu
      t102(n) = t102(n)*hxu
      ! end of surface layer treatement.
      if (kbu < 2) cycle

      ! k lims:
      kbv   = khv(n)
      kbvn  = khv(msrf(i-1,j))
      kbve  = khv(msrf(i,j+1)) 
      kbvne = khv(msrf(i-1,j+1)) 

      ! index off-sets:
      mi0  = mcol(n) - 2
      me0  = mcol(msrf(i,  j+1)) - 2
      mn0  = mcol(msrf(i-1,j  )) - 2
      ms0  = mcol(msrf(i+1,j  )) - 2
      mne0 = mcol(msrf(i-1,j+1)) - 2
      mse0 = mcol(msrf(i+1,j+1)) - 2

      ! k=2 unrolled:
      mi = mi0 + 2
      hxu = hx(mi)*u(mi)
      if (hxu > zero) then
        sx  = dtddx*u(mi)
        sxm = (half-fourthird*sx)
        fac = half-sx
        t101(mi) = fac*t401(mi)
        t102(mi) = fac*t402(mi)
        mn = mcol(msrf(i-1,j))
        if (v(mn) < zero) then
          sy = dtddy*v(mn)
          s1 = sxm*sy
          s2 = (half + twothird*sy)*sy
          t101(mi) = t101(mi) - s1*(t401(mn)-t401(mi)) + s2*(t501(mn)-t501(mi))
          t102(mi) = t102(mi) - s1*(t402(mn)-t402(mi)) + s2*(t502(mn)-t502(mi))
        endif
        if (v(mi) > zero) then
          ms = mcol(msrf(i+1,j))
          sy = dtddy*v(mi)
          s1 = sxm*sy
          s2 = (half- twothird*sy)*sy
          t101(mi) = t101(mi) - s1*(t401(mi)-t401(ms)) - s2*(t501(mi)-t501(ms))
          t102(mi) = t102(mi) - s1*(t402(mi)-t402(ms)) - s2*(t502(mi)-t502(ms))
        endif
        if (w(mi) < zero) then
          mu = n
          sz = w(mi)*hdt/h_new(mi)
          s1 = sxm*sz
          s2 = (half + twothird*sz)*sz
          t101(mi) = t101(mi) - s1*(t401(mu)-t401(mi)) + s2*(t601(mu)-t601(mi))
          t102(mi) = t102(mi) - s1*(t402(mu)-t402(mi)) + s2*(t602(mu)-t602(mi))
        endif
        if (kb > 2) then
          md = mi + 1
          if (w(md) > zero) then
            sz = w(md)*hdt/h_new(mi)
            s1 = sxm*sz
            s2 = (half - twothird*sz)*sz
            t101(mi) = t101(mi) - s1*(t401(mi)-t401(md))- s2*(t601(mi)-t601(md))
            t102(mi) = t102(mi) - s1*(t402(mi)-t402(md))- s2*(t602(mi)-t602(md))
          endif
        endif
      else ! u < 0.0
        sx  = dtddx*u(mi)
        sxp = (half+fourthird*sx)
        fac = -half-sx
        me  = mcol(msrf(i,j+1))
        t101(mi) = fac*t401(me)
        t102(mi) = fac*t402(me)
        mne = mcol(msrf(i-1,j+1))
        if (v(mne) < zero) then
          sy = dtddy*v(mne)
          s1 = sxp*sy
          s2 = (half + twothird*sy)*sy
          t101(mi) = t101(mi) + s1*(t401(mne)-t401(me))- s2*(t501(mne)-t501(me))
          t102(mi) = t102(mi) + s1*(t402(mne)-t402(me))- s2*(t502(mne)-t502(me))
        endif
        if (v(me) > zero) then
          sy = dtddy*v(me)
          s1 = sxp*sy
          s2 = (half - twothird*sy)*sy
          mse = mcol(msrf(i+1,j+1))
          t101(mi) = t101(mi) + s1*(t401(me)-t401(mse))+ s2*(t501(me)-t501(mse))
        endif
        if (w(me) < zero) then
          mue = msrf(i,j+1)
          sz = w(me)*hdt/h_new(mi)
          s1 = sxp*sz
          s2 = (half + twothird*sz)*sz
          t101(mi) = t101(mi) + s1*(t401(mue)-t401(me))- s2*(t601(mue)-t601(me))
          t102(mi) = t102(mi) + s1*(t402(mue)-t402(me))- s2*(t602(mue)-t602(me))
        endif
        if (kbe > 2) then
          mde = me + 1
          if (w(mde) > zero) then
            dtddz = hdt/h_new(mi)
            sz = dtddz*w(mde)
            s1 = sxp*sz
            s2 = (half- twothird*sz)*sz
            t101(mi)= t101(mi) + s1*(t401(me)-t401(mde))+s2*(t601(me)-t601(mde))
            t102(mi)= t102(mi) + s1*(t402(me)-t402(mde))+s2*(t602(me)-t602(mde))
          endif
        endif
      endif
      t101(mi) = t101(mi)*hxu
      t102(mi) = t102(mi)*hxu

      ! main loops:
      kupper = min(kbu,kbv,kbvn,kbve,kbvne,kb-1,kbe-1)

      ! main loop #1, t40[12]
      do k=3,kupper
        mi  = mi0 + k
        me  = me0 + k
        mn  = mn0 + k
        ms  = ms0 + k
        mne = mne0 + k
        mse = mse0 + k
        mu  = mi - 1
        md  = mi + 1
        mue = me - 1
        mde = me + 1
        
        hxu = hx(mi)*u(mi)
        umx = max(hxu,zero)
        umn = min(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx
        fac = half - sx

        sy1 = dtddy*min(v(mn),zero)

        sy2 = dtddy*min(v(mne),zero)

        sy3 = dtddy*max(v(mi),zero)

        sy4 = dtddy*max(v(me),zero)

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)

        sz2 = dtddz*min(w(me),zero)

        sz3 = dtddz*max(w(md),zero)

        sz4 = dtddz*max(w(mde),zero)

        t101(mi) = umx*( fac*t401(mi)                                          &
                       -sxx*( sy1*(t401(mn)-t401(mi))                          &
                             +sy3*(t401(mi)-t401(ms))                          &
                             +sz1*(t401(mu)-t401(mi))                          &
                             +sz3*(t401(mi)-t401(md)) ))                       &
                 - umn*( fac*t401(me)                                          &
                       -sxx*( sy2*(t401(mne)-t401(me))                         &
                             +sy4*(t401(me)-t401(mse))                         &
                             +sz2*(t401(mue)-t401(me))                         &
                             +sz4*(t401(me)-t401(mde)) ))
        t102(mi) = umx*( fac*t402(mi)                                          &
                       -sxx*( sy1*(t402(mn)-t402(mi))                          &
                             +sy3*(t402(mi)-t402(ms))                          &
                             +sz1*(t402(mu)-t402(mi))                          &
                             +sz3*(t402(mi)-t402(md)) ))                       &
                 - umn*( fac*t402(me)                                          &
                       -sxx*( sy2*(t402(mne)-t402(me))                         &
                             +sy4*(t402(me)-t402(mse))                         &
                             +sz2*(t402(mue)-t402(me))                         &
                             +sz4*(t402(me)-t402(mde)) ))
      enddo

      ! main loop #2, t50[12] and t60[12]
      do k=3,kupper
        mi  = mi0 + k
        me  = me0 + k
        mn  = mn0 + k
        ms  = ms0 + k
        mne = mne0 + k
        mse = mse0 + k
        mu  = mi - 1
        md  = mi + 1
        mue = me - 1
        mde = me + 1
        
        hxu = hx(mi)*u(mi)
        umx = max(hxu,zero)
        umn = min(hxu,zero)

        sy1 = dtddy*min(v(mn),zero)
        s1  = (half + twothird*sy1)*sy1

        sy2 = dtddy*min(v(mne),zero)
        s2  = (half + twothird*sy2)*sy2

        sy3 = dtddy*max(v(mi),zero)
        s3  = (half - twothird*sy3)*sy3

        sy4 = dtddy*max(v(me),zero)
        s4  = (half - twothird*sy4)*sy4

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)
        s5  = (half + twothird*sz1)*sz1

        sz2 = dtddz*min(w(me),zero)
        s6  = (half + twothird*sz2)*sz2

        sz3 = dtddz*max(w(md),zero)
        s7  = (half - twothird*sz3)*sz3

        sz4 = dtddz*max(w(mde),zero)
        s8  = (half - twothird*sz4)*sz4

        t101(mi) = t101(mi)                                                    &
                 + umx*(                                                       &
                       +s1*(t501(mn)-t501(mi))                                 &
                       -s3*(t501(mi)-t501(ms))                                 &
                       +s5*(t601(mu)-t601(mi))                                 &
                       -s7*(t601(mi)-t601(md)) )                               &
                 - umn*(                                                       &
                       +s2*(t501(mne)-t501(me))                                &
                       -s4*(t501(me)-t501(mse))                                &
                       +s6*(t601(mue)-t601(me))                                &
                       -s8*(t601(me)-t601(mde)))
        t102(mi) = t102(mi)                                                    &
                 + umx*(                                                       &
                       +s1*(t502(mn)-t502(mi))                                 &
                       -s3*(t502(mi)-t502(ms))                                 &
                       +s5*(t602(mu)-t602(mi))                                 &
                       -s7*(t602(mi)-t602(md)) )                               &
                 - umn*(                                                       &
                       +s2*(t502(mne)-t502(me))                                &
                       -s4*(t502(me)-t502(mse))                                &
                       +s6*(t602(mue)-t602(me))                                &
                       -s8*(t602(me)-t602(mde)))
      enddo

      ! remainder loops:
      klow = max(3, kupper+1)
      do k=klow,min(kbu,kb-1)
        mi  = mi0 + k
        me  = me0 + k
        mu  = mi - 1
        md  = mi + 1
        mue = me - 1
        
        hxu = hx(mi)*u(mi)
        umx = max(hxu,zero)
        umn = min(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx
        fac = half - sx

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)
        s5  = half + twothird*sz1

        sz2 = dtddz*min(w(me),zero)
        s6  = half + twothird*sz2

        sz3 = dtddz*max(w(md),zero)
        s7  = half - twothird*sz3

        t101(mi) = umx*( fac*t401(mi)                                          &
                       -sz1*(sxx*(t401(mu)-t401(mi)) -s5*(t601(mu)-t601(mi)) ) &
                       -sz3*(sxx*(t401(mi)-t401(md)) +s7*(t601(mi)-t601(md)) ))&
                 - umn*( fac*t401(me)                                          &
                       -sz2*(sxx*(t401(mue)-t401(me))-s6*(t601(mue)-t601(me))))
        t102(mi) = umx*( fac*t402(mi)                                          &
                       -sz1*(sxx*(t402(mu)-t402(mi)) -s5*(t602(mu)-t602(mi)) ) &
                       -sz3*(sxx*(t402(mi)-t402(md)) +s7*(t602(mi)-t602(md)) ))&
                 - umn*( fac*t402(me)                                          &
                       -sz2*(sxx*(t402(mue)-t402(me))-s6*(t602(mue)-t602(me))))
      enddo
      do k=klow,min(kbu,kb-1,kbe-1)
        mi  = mi0 + k
        me  = me0 + k
        mde = me + 1

        hxu = hx(mi)*u(mi)
        umn = min(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx

        dtddz = hdt/h_new(mi)

        sz4 = dtddz*max(w(mde),zero)
        s8  = half - twothird*sz4
        
        t101(mi) = t101(mi)                                                    &
                 - umn*(                                                       &
                       -sz4*(sxx*(t401(me)-t401(mde))+s8*(t601(me)-t601(mde))))
        t102(mi) = t102(mi)                                                    &
                 - umn*(                                                       &
                       -sz4*(sxx*(t402(me)-t402(mde))+s8*(t602(me)-t602(mde))))
      enddo
      do k=klow,min(kbu,kb-1,kbv)
        mi  = mi0 + k
        ms  = ms0 + k

        hxu = hx(mi)*u(mi)
        umx = max(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx

        sy3 = dtddy*max(v(mi),zero)
        s3  = half - twothird*sy3

        t101(mi) = t101(mi)                                                    &
                 + umx*(                                                       &
                       -sy3*(sxx*(t401(mi)-t401(ms)) +s3*(t501(mi)-t501(ms)) ))
        t102(mi) = t102(mi)                                                    &
                 + umx*(                                                       &
                       -sy3*(sxx*(t402(mi)-t402(ms)) +s3*(t502(mi)-t502(ms)) ))
      enddo
      do k=klow,min(kbu,kb-1,kbvn)
        mi  = mi0 + k
        mn  = mn0 + k

        hxu = hx(mi)*u(mi)
        umx = max(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx

        sy1 = dtddy*min(v(mn),zero)
        s1  = half + twothird*sy1

        t101(mi) = t101(mi)                                                    &
                 + umx*(                                                       &
                       -sy1*(sxx*(t401(mn)-t401(mi)) -s1*(t501(mn)-t501(mi)) ))
        t102(mi) = t102(mi)                                                    &
                 + umx*(                                                       &
                       -sy1*(sxx*(t402(mn)-t402(mi)) -s1*(t502(mn)-t502(mi)) ))
      enddo
      do k=klow,min(kbu,kb-1,kbve)
        mi  = mi0 + k
        me  = me0 + k
        mse = mse0 + k

        hxu = hx(mi)*u(mi)
        umn = min(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx

        sy4 = dtddy*max(v(me),zero)
        s4  = half - twothird*sy4

        t101(mi) = t101(mi)                                                    &
                 - umn*(                                                       &
                       -sy4*(sxx*(t401(me)-t401(mse))+s4*(t501(me)-t501(mse))))
        t102(mi) = t102(mi)                                                    &
                 - umn*(                                                       &
                       -sy4*(sxx*(t402(me)-t402(mse))+s4*(t502(me)-t502(mse))))
      enddo
      do k=klow,min(kbu,kb-1,kbvne)
        mi  = mi0 + k
        me  = me0 + k
        mne = mne0 + k

        hxu = hx(mi)*u(mi)
        umn = min(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx

        sy2 = dtddy*min(v(mne),zero)
        s2  = half + twothird*sy2

        t101(mi) = t101(mi)                                                    &
                 - umn*(                                                       &
                       -sy2*(sxx*(t401(mne)-t401(me))-s2*(t501(mne)-t501(me))))
        t102(mi) = t102(mi)                                                    &
                 - umn*(                                                       &
                       -sy2*(sxx*(t402(mne)-t402(me))-s2*(t502(mne)-t502(me))))
      enddo


      ! bottom unrolled:
      if (kb == kbu .and. kb > 2) then
        k = kbu
        mi  = mi0 + k
        me  = me0 + k
        if (kh(msrf(i-1,j)) >= k) then
          mn  = mn0 + k
        else
          mn  = 0
        endif
        if (kh(msrf(i+1,j)) >= k) then
          ms  = ms0 + k
        else
          ms  = 0
        endif
        if (kh(msrf(i-1,j+1)) >= k) then
          mne = mne0 + k
        else
          mne = 0
        endif
        if (kh(msrf(i+1,j+1)) >= k) then
          mse = mse0 + k
        else
          mse = 0
        endif
        mu  = mi - 1
        mue = me - 1
        if (kh(msrf(i,j+1)) > k) then
          mde = me + 1
        else
          mde = 0
        endif
        
        hxu = hx(mi)*u(mi)
        umx = max(hxu,zero)
        umn = min(hxu,zero)

        sx  = dtddx*abs(u(mi))
        sxx = half - fourthird*sx
        fac = half - sx

        sy1 = dtddy*min(v(mn),zero)
        s1  = half + twothird*sy1

        sy2 = dtddy*min(v(mne),zero)
        s2  = half + twothird*sy2

        sy3 = dtddy*max(v(mi),zero)
        s3  = half - twothird*sy3

        sy4 = dtddy*max(v(me),zero)
        s4  = half - twothird*sy4

        dtddz = hdt/h_new(mi)

        sz1 = dtddz*min(w(mi),zero)
        s5  = half + twothird*sz1

        sz2 = dtddz*min(w(me),zero)
        s6  = half + twothird*sz2

        sz4 = dtddz*max(w(mde),zero)
        s8  = half - twothird*sz4

        t101(mi) = umx*( fac*t401(mi)                                          &
                       -sy1*(sxx*(t401(mn)-t401(mi)) -s1*(t501(mn)-t501(mi)) ) &
                       -sy3*(sxx*(t401(mi)-t401(ms)) +s3*(t501(mi)-t501(ms)) ) &
                       -sz1*(sxx*(t401(mu)-t401(mi)) -s5*(t601(mu)-t601(mi)) ))&
                 - umn*( fac*t401(me)                                          &
                       -sy2*(sxx*(t401(mne)-t401(me))-s2*(t501(mne)-t501(me))) &
                       -sy4*(sxx*(t401(me)-t401(mse))+s4*(t501(me)-t501(mse))) &
                       -sz2*(sxx*(t401(mue)-t401(me))-s6*(t601(mue)-t601(me))) &
                       -sz4*(sxx*(t401(me)-t401(mde))+s8*(t601(me)-t601(mde))))
        t102(mi) = umx*( fac*t402(mi)                                          &
                       -sy1*(sxx*(t402(mn)-t402(mi)) -s1*(t501(mn)-t501(mi)) ) &
                       -sy3*(sxx*(t402(mi)-t402(ms)) +s3*(t501(mi)-t501(ms)) ) &
                       -sz1*(sxx*(t402(mu)-t402(mi)) -s5*(t601(mu)-t601(mi)) ))&
                 - umn*( fac*t402(me)                                          &
                       -sy2*(sxx*(t402(mne)-t402(me))-s2*(t502(mne)-t502(me))) &
                       -sy4*(sxx*(t402(me)-t402(mse))+s4*(t502(me)-t502(mse))) &
                       -sz2*(sxx*(t402(mue)-t402(me))-s6*(t602(mue)-t602(me))) &
                       -sz4*(sxx*(t402(me)-t402(mde))+s8*(t602(me)-t602(mde))))
      endif
    enddo
#ifdef _OPENACC
    enddo
!$ACC WAIT
!$ACC END DATA
#endif
    
  end subroutine c_dtx

  ! ----------------------------------------------------------------------------
    
  subroutine tflow_up_ext (ia,msrf,mcol,n2d,ind,kh,khu,khv,idx,u,v,w,hx,hy,    &
                           h_old,h_new,kmax,cosphi,to01,to02,ti01,ti02,dt,dx,dy)
  ! Simple 3D advection scheme using 
  !    control volume with upwind centred cell face tracers.
  !    uses t1, t2, tt
  ! External version, used for calls outside the tflow_int framework, i.e. with
  ! an input tracers array, ti, and output to a separate array, to. Also, we 
  ! need to define some constants that are module vars defined in tflow_int.

!FIXME: This is a first draft, still loads of stuff to do, e.g. to:
!       * secure stride-1 k-loops by splitting into main loops and 
!         remainer-loops as we do elsewhere,
!       * split the iR/iS/W calls for the halo swap; as is we do all in one call
!       * handle river inflow

    !- modules -----------------------------------------------------------------
    use dmi_omp,    only : domp_get_domain
    use dmi_mpi,    only : dmpi_halo, mpi_size
    use e_timer, only : timer

    !- directives --------------------------------------------------------------
    implicit none

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: ia, n2d, kmax
    integer(4), intent(in)    :: msrf(0:,0:), mcol(0:), ind(1:,1:)
    integer(4), intent(in)    :: kh(0:), khu(0:), khv(0:)
    integer(4), intent(inout) :: idx(1:,1:)
    real(8),    intent(in)    :: u(0:), v(0:), w(0:), hx(0:), hy(0:)
    real(8),    intent(in)    :: h_old(0:), h_new(0:), cosphi(1:,0:)
    real(8),    intent(inout) :: to01(0:), to02(0:)
    real(8),    intent(in)    :: ti01(0:), ti02(0:)
    real(8),    intent(in)    :: dt, dx, dy

    include 'tflow_up_ext.inc'

    !- local vars --------------------------------------------------------------
    integer(4) :: n, i, j, k, nl, nu, mmi, mme, mmw, mms, mmn, mmd, kb, mm2
    integer(4) :: mme2, mms2, mmw2, mmn2
    real(8)    :: fac, dtdxi, dtdxidy, dv, t1w, t2w, t1n, t2n

    !- some initialisations ----------------------------------------------------
    call domp_get_domain(kh, 1, n2d, nl, nu, idx)
!$OMP MASTER
    t101(0) = zero
    t201(0) = zero
    tt01(0) = zero
    t102(0) = zero
    t202(0) = zero
    tt02(0) = zero
!$OMP END MASTER
!$OMP BARRIER

    dtdx = dt*dx
    dtdy = dt*dy
    dxdy = dx*dy
    dtdxdy = dt*dx*dy

    !  X-transports ------------------------------------------------------------
    do n=nl,nu  !  x-transports, k=1
      i = ind(1,n)
      j = ind(2,n)
  
      mme = msrf(i,j+1)
      if (mme == 0) then
        t101(n) = zero
        t102(n) = zero
      else
        fac = dtdy*hx(n)*u(n)
        if (fac >= zero) then
          t101(n) = ti01(n)*fac
          t102(n) = ti02(n)*fac
        else
          t101(n) = ti01(mme)*fac
          t102(n) = ti02(mme)*fac
        endif
      endif
    enddo
    do n=nl,nu  !  x-transports, k>1
      if (kh(n) < 2) cycle

      i = ind(1,n)
      j = ind(2,n)

      mm2 = mcol(n) - 2

      do mmi=mm2+max(khu(n),1)+1,mm2+kh(n)
        t101(mmi) = zero
        t102(mmi) = zero
      enddo

      if (khu(n) < 2) cycle
      mme2 = mcol(msrf(i,j+1)) - 2

      do k=2,khu(n)
        fac = dtdy*hx(mm2+k)*u(mm2+k)
        t101(mm2+k) = ti01(mm2 +k)*max(fac,zero) + ti01(mme2+k)*min(fac,zero)
        t102(mm2+k) = ti02(mm2 +k)*max(fac,zero) + ti02(mme2+k)*min(fac,zero)
      enddo
    enddo

    !  Y-transports ------------------------------------------------------------
    do n=nl,nu  !  y-transports, k=1
      i = ind(1,n)
      j = ind(2,n)

      mms = msrf(i+1,j)
      if (mms == 0) then
        t201(n) = zero
        t202(n) = zero
      else
        fac = dtdx*cosphi(2,i)*hy(n)*v(n)
        if (fac <= zero) then
          t201(n) = ti01(n)*fac
          t202(n) = ti02(n)*fac
        else
          t201(n) = ti01(mms)*fac
          t202(n) = ti02(mms)*fac
        endif
      endif
    enddo
    do n=nl,nu  !  y-transports, k>1
      if (kh(n) < 2) cycle

      i = ind(1,n)
      j = ind(2,n)

      dtdxi = dtdx*cosphi(2,i)
      mm2   = mcol(n) - 2

      do mmi=mm2+max(khv(n),1)+1,mm2+kh(n)
        t201(mmi) = zero
        t202(mmi) = zero
      enddo

      if (khv(n) < 2) cycle
      mms2 = mcol(msrf(i+1,j)) - 2

      do k=2,khv(n)
        fac = dtdxi*hy(mm2+k)*v(mm2+k)
        t201(mm2+k) = ti01(mm2 +k)*min(fac,zero) + ti01(mms2+k)*max(fac,zero)
        t202(mm2+k) = ti02(mm2 +k)*min(fac,zero) + ti02(mms2+k)*max(fac,zero)
      enddo
    enddo

    !  Z-transports ------------------------------------------------------------
    tt01(nl:nu) = zero   !  z-transports, k=1
    tt02(nl:nu) = zero   !  z-transports, k=1
    do n=nl,nu           !  z-transports, k>1
      kb = kh(n)
      if (kb < 2) cycle

      dtdxidy = dtdxdy*cosphi(1,ind(1,n))
      mm2     = mcol(n) - 2

      ! unroll k=2:
      mmi = mm2 + 2
      fac = dtdxidy*w(mmi)
      if (fac >= zero) then
        tt01(mmi) = ti01(mmi)*fac
        tt02(mmi) = ti02(mmi)*fac
      else
        tt01(mmi) = ti01(n)*fac
        tt02(mmi) = ti02(n)*fac
      endif

      do mmi=mm2+3,mm2+kb
        fac = dtdxidy*w(mmi)
        tt01(mmi) = ti01(mmi)*max(fac,zero) + ti01(mmi-1)*min(fac,zero)
        tt02(mmi) = ti02(mmi)*max(fac,zero) + ti02(mmi-1)*min(fac,zero)
      enddo
    enddo

!$OMP BARRIER

    if (mpi_size > 1) then
!FIXME: consider pre-post recv earlier (here FLAG=6 takes all comm in one go)
      call timer(1,'MPI_haloTF')
      call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,6,t101,t201,t102,t202)
!$OMP BARRIER
      call timer(2,'MPI_haloTF')
    endif

    !  Update tracers  ---------------------------------------------------------
    do n=nl,nu  ! k=1
      i = ind(1,n)
      j = ind(2,n)

      dv  = cosphi(1,i)*dxdy

      mmw = msrf(i,  j-1)
      mmn = msrf(i-1,j  )
      if (kh(n) > 1) then
        mmd = mcol(n)
      else
        mmd = 0
      endif

      to01(n) = (  ti01(n)*h_old(n)                                            &
                    - (  (t101(n)  -t101(mmw))                                 &
                       + (t201(mmn)-t201(n)  )                                 &
                       + (tt01(n)  -tt01(mmd)) )/dv                            &
                )/h_new(n)
      to02(n) = (  ti02(n)*h_old(n)                                            &
                    - (  (t102(n)  -t102(mmw))                                 &
                       + (t202(mmn)-t202(n)  )                                 &
                       + (tt02(n)  -tt02(mmd)) )/dv                            &
                )/h_new(n)
    enddo
    do n=nl,nu  ! k>1
      kb = kh(n)
      if (kb < 2) cycle

      i = ind(1,n)
      j = ind(2,n)

      dv  = cosphi(1,i)*dxdy
      mm2 = mcol(n) - 2

      if (kh(msrf(i,j-1)) > 1) mmw2 = mcol(msrf(i,j-1)) - 2
      if (kh(msrf(i-1,j)) > 1) mmn2 = mcol(msrf(i-1,j)) - 2
 
      ! This is cheating: Instead of doing a main-loop plus a couple of
      ! remiander loops, we do three (almost) full loop.
      do k=2,kb-1
        to01(mm2+k) = (  ti01(mm2+k)*h_old(mm2+k)                              &
                          - (   t101(mm2 +k)                                   &
                                            -t201(mm2 +k)                      &
                             + (tt01(mm2 +k)-tt01(mm2 +k+1)) )/dv              &
                      )/h_new(mm2+k)
        to02(mm2+k) = (  ti02(mm2+k)*h_old(mm2+k)                              &
                          - (   t102(mm2 +k)                                   &
                                            -t202(mm2 +k)                      &
                             + (tt02(mm2 +k)-tt02(mm2 +k+1)) )/dv              &
                      )/h_new(mm2+k)
      enddo
      do k=2,min(kb-1,kh(msrf(i,j-1)))
        to01(mm2+k) = to01(mm2+k)                                              &
                      - ((                  -t101(mmw2+k) )/dv                 &
                        )/h_new(mm2+k)
        to02(mm2+k) = to02(mm2+k)                                              &
                      - ((                  -t102(mmw2+k) )/dv                 &
                        )/h_new(mm2+k)
      enddo
      do k=2,min(kb-1,kh(msrf(i-1,j)))
        to01(mm2+k) = to01(mm2+k)                                              &
                      - ((   t201(mmn2+k)                 )/dv                 &
                        )/h_new(mm2+k)
        to02(mm2+k) = to02(mm2+k)                                              &
                      - ((   t202(mmn2+k)                 )/dv                 &
                        )/h_new(mm2+k)
      enddo

      k = kb
      mmi = mm2 + k
      if (kh(msrf(i,j-1)) >= kb) then
        t1w = t101(mmw2+k) 
        t2w = t102(mmw2+k) 
      else
        t1w = zero
        t2w = zero
      endif
      if (kh(msrf(i-1,j)) >= kb) then
        t1n = t201(mmn2+k)
        t2n = t202(mmn2+k)
      else
        t1n = zero
        t2n = zero
      endif

      fac = one/(h_new(mmi)*dv)

      to01(mmi) = (  ti01(mmi)*h_old(mmi)                                      &
                      - (  (t101(mmi)-t1w      )                               &
                         + (t1n      -t201(mmi))                               &
                         +  tt01(mmi)               )/dv                       &
                  )/h_new(mmi)
      to02(mmi) = (  ti02(mmi)*h_old(mmi)                                      &
                      - (  (t102(mmi)-t2w      )                               &
                         + (t2n      -t202(mmi))                               &
                         +  tt02(mmi)               )/dv                       &
                  )/h_new(mmi)
    enddo

  end subroutine tflow_up_ext
    
  ! ----------------------------------------------------------------------------

  subroutine tflow_int_simd(kmax,msrf,mcol,ind,n2d,dx_j,dy_i,advecstab,dt_in,  &
                            kh,khu,khv,h,hx,hy,h_old,h_new,u,v,w,t,s,          &
                            avt,avs,eddyh,nbpz,krz,nbpq,krq,rwqk,              &
                            nbpu,kru,nbpv,krv,rwzkt,rwqkt,cosphi,              &
                            doadvec,dodiff,idx,halo2d,ia,dyn,uvdam,vdiff,q)

    !- modules -----------------------------------------------------------------
    use dmi_mpi,          only : dmpi_halo, mpi_size
    use e_timer,       only : timer
    use masscorr_srf_col, only : masscorr_init, masscorr
    use constants,        only : msscrr
    
    !- directives --------------------------------------------------------------
    implicit none
    
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)           :: ia,kmax,n2d,halo2d
    integer(4), intent(in)           :: ind(:,:)
    integer(4), intent(in)           :: msrf(0:,0:), mcol(0:)
    integer(4), intent(in)           :: kh(0:),khu(0:),khv(0:)
    integer(4), intent(in)           :: nbpz,nbpq,nbpu,nbpv
    integer(4), intent(in)           :: krq(:,0:),krz(:,0:),kru(:,0:),krv(:,0:)
    real(8),    intent(in)           :: dt_in, advecstab
    real(8),    intent(in)           :: dx_j,dy_i
    real(8),    intent(in)           :: h(0:),hx(0:),hy(0:)
    real(8),    intent(in)           :: h_old(0:),h_new(0:)
    real(8),    intent(in)           :: avt(0:), avs(0:), uvdam(:,0:)
    real(8),    intent(in)           :: u(0:),v(0:),w(0:),eddyh(0:)
    real(8),    intent(in)           :: rwqk(0:,:),rwzkt(:,:,0:)
    real(8),    intent(in)           :: rwqkt(:,:,0:),cosphi(:,0:)
    logical,    intent(in)           :: doadvec, dodiff, dyn
    integer(4), intent(inout)        :: idx(1:,1:)
    logical,    intent(in), optional :: vdiff
    real(8),    intent(in), optional :: q(0:)
    real(8),    intent(inout)        :: t(0:), s(0:)

    include 'tflow_int.inc' 
    !- local vars --------------------------------------------------------------
    integer(4) :: ndry_act

    !- obtain constant dx, dy, local dt and derived quantities: ----------------
    dx     = dx_j
    dy     = dy_i
    dt     = dt_in
    dxdy   = dx*dy
    ddxdy  = one/dxdy
    dtdx   = dt*dx
    dtdy   = dt*dy
    dtdxdy = dt*dx*dy

    ndry_act = 0

    if (doadvec) then ! do the original advection scheme

#ifdef _OPENACC
!$acc data pcopyin(kmax,msrf,mcol,ind,n2d,advecstab,h_new,khu,v,w,kh,idx,      &
!$acc&             khv,u,cosphi,h,nbpz,krz,nbpu,kru,nbpv,krv,                  &
!$acc&             hx,hy,h_old,nbpq,krq,rwqk,rwzkt,rwqkt)
#endif
    
      if (mpi_size > 1) then
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,1,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,5,2,1,t401,t501,t601,            &
                                                    t402,t502,t602)
        call timer(2,'MPI_haloTF')
      endif

      if ((.not.dyn) .and. msscrr) then
        ! prepared mass-corrective step ----------------------------------------
        call masscorr_init (n2d, h_new, h_old, hx, hy, advecstab, u, v, dt,    &
                            dx, dy, cosphi, ind, msrf, kh, idx, t, s, ndry_act)
      endif

      ! ------------------------------------------------------------------------
      ! begin advection
      ! ------------------------------------------------------------------------

      ! compute t1 -------------------------------------------------------------
      call c_tu(kmax,msrf,mcol,ind,n2d,advecstab,h_new,khu,v,w,kh,idx,t,s)

      ! compute t2 -------------------------------------------------------------
      call c_tv(kmax,msrf,mcol,ind,n2d,advecstab,h_new,khv,u,w,cosphi,kh,idx,  &
                t,s)

      if (mpi_size > 1) then
!$OMP BARRIER  ! for the halo swap
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,3,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,4,t101,t201,t102,t202)
        call timer(2,'MPI_haloTF')
      endif

      ! compute t3 -------------------------------------------------------------
      call c_tw(msrf,mcol,ind,n2d,kh,u,v,cosphi,idx,t,s)

      ! compute t4, t5, t6  ----------------------------------------------------
      call c_delta(msrf,mcol,ind,n2d,kh,khu,khv,h,nbpz,krz,nbpu,kru,nbpv,krv,  &
                   idx,t,s)
!$OMP BARRIER ! A must even with single node runs

      if (mpi_size > 1) then
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,5,2,3,t401,t501,t601,            &
                                                    t402,t502,t602)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,5,2,4,t401,t501,t601,            &
                                                    t402,t502,t602)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,2,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,5,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,3,1,1,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,11,1,3,1,tt01,tt02)
        call timer(2,'MPI_haloTF')
      endif

      ! compute tt -------------------------------------------------------------
      ! this computation relies on t1, t2 and t3, i.e. barrier here
      call c_tt (kmax,msrf,mcol,ind,n2d,kh,khu,khv,hx,hy,h_old,h_new,          &
                 u,v,w,nbpz,krz,nbpq,krq,rwqk,rwzkt,rwqkt,cosphi,idx,t,s)
!$OMP BARRIER ! A must even with single node runs

      if (mpi_size > 1) then
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia,11,1,3,3,tt01,tt02)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,11,1,3,4,tt01,tt02)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,5,2,2,t401,t501,t601,            &
                                                    t402,t502,t602)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,5,2,5,t401,t501,t601,            &
                                                    t402,t502,t602)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,1,2,1,t701,t801,t702,t802)
        call timer(2,'MPI_haloTF')
      endif

      ! compute t1, t2, t3  ----------------------------------------------------
      ! this computation relies on t4, t5, t6, i.e. barrier here
      call c_dtx(msrf,mcol,ind,n2d,advecstab,h_new,kh,khu,khv,u,v,w,hx,cosphi, &
                 idx)

      call c_dty(msrf,mcol,ind,n2d,advecstab,h_new,kh,khu,khv,u,v,w,hy,cosphi, &
                 idx)

      call c_dtz(msrf,mcol,ind,n2d,advecstab,h_new,kh,u,v,w,cosphi,idx)

!$OMP BARRIER ! A must even with single node runs

      if (mpi_size > 1) then
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia, 1,3,1,3,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia, 1,3,1,4,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia, 1,3,1,2,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia, 1,3,1,5,t101,t201,t102,t202)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,11,1,3,2,tt01,tt02)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,11,1,3,5,tt01,tt02)
        call dmpi_halo(kmax,msrf,mcol,kh,ia, 1,1,1,1,t401,t501,t402,t502)
        call timer(2,'MPI_haloTF')
      endif

      ! compute t8, t7 ---------------------------------------------------------
      ! this computation relies on tt, t1, t2, t3, i.e. barrier here
      call c_rin_rout(kmax,msrf,mcol,ind,n2d,kh,h_new,cosphi,idx,t,s)
!$OMP BARRIER ! A must even with single node runs

      if (mpi_size > 1) then
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,1,2,3,t701,t801,t702,t802)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,1,2,4,t701,t801,t702,t802)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,1,2,2,t701,t801,t702,t802)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,5,1,2,5,t701,t801,t702,t802)
        call timer(2,'MPI_haloTF')
      endif

      ! compute t4, t5, t6 ----------------------------------------------------
      ! this computation relies on t7, t8, t1, t2, t3, i.e. barrier here
      call c_cx_cy_cz (msrf,mcol,ind,n2d,kh,idx)
!$OMP BARRIER ! A must even with single node runs

      if (mpi_size > 1) then
        call timer(1,'MPI_haloTF')
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,3,t401,t501,t402,t502)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,4,t401,t501,t402,t502)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,2,t401,t501,t402,t502)
        call dmpi_halo(kmax,msrf,mcol,kh,ia,1,1,1,5,t401,t501,t402,t502)
        call timer(2,'MPI_haloTF')
      endif

      ! Finally, compute t -----------------------------------------------------
      ! this computation relies on tt, t1, t2, t3, t4, t5, t6, i.e. barrier here
      call advection(msrf,mcol,ind,n2d,kh,khu,khv,h_new,cosphi,idx,t,s)
!$OMP BARRIER ! A must even with single node runs
!$acc end data

      if ((.not.dyn) .and. ndry_act > 0) then
        ! do mass-correction on t and s ----------------------------------------
        call masscorr (t,s)
      endif
    endif

    if (dodiff .and. doadvec .and. ndry_act > 0) then
      ! barrier here due to fusion of two OMP parallel sections in case both
      ! advection and diffusion is activated in the same call
!$OMP BARRIER
    endif

    if (dodiff .and. present(vdiff)) then
      if (vdiff) then
        ! inlined version
        call diffusion_vi(kmax,mcol,n2d,kh,h,h_new,avt,avs,idx,t,s,q)
      else
        ! avoid data race for t&s in diffusion() by making a copy:
        call copy_t2tt(mcol,n2d,halo2d,kh,idx,t,s)
!$OMP BARRIER
        call diffusion_hx(msrf,mcol,ind,n2d,kh,h_new,hx,hy,eddyh,cosphi,idx,   &
                          uvdam,t,s)
      endif

    elseif (dodiff) then
! MPI: Halo swap already done outside tflow-int for: t&s
! MPI: Halo of tt also assigned in copy_t2tt(), i.e. no need for halo swap
      ! avoid data race for t1d in diffusion() by making a copy:
      call copy_t2tt(mcol,n2d,halo2d,kh,idx,t,s)
!$OMP BARRIER

! MPI: Halo swap already done outside tflow-int for: eddyh
      call diffusion_hx(msrf,mcol,ind,n2d,kh,h_new,hx,hy,eddyh,cosphi,idx,     &
                        uvdam,t,s)
      ! inlined version
      call diffusion_vi(kmax,mcol,n2d,kh,h,h_new,avt,avs,idx,t,s)
    endif

  end subroutine tflow_int_simd

end module tflow_simd_srf_col
