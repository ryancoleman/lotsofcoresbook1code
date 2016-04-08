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

module masscorr_srf_col
  !  Wrapper for the mass correction to be applied in with tracer advection
  !  when drying is present.
  !  The present implementation attemps to keep the increased memory usage 
  !  and computational load low by using temporary compressed work arrays.
  !  This means that OpenMP is not supported here, but MPI is if supported by
  !  the caller. I.e. it _MUST_ be called outside OpenMP parallel sections, and
  !  and masscorr_init() should be called on the same MPI-task where the 
  !  masscorr() is later applied but then there is no need for explicit 
  !  MPI-communication (gather/scatter/haloswap).

  implicit none
  private

  public :: masscorr_firstentry, masscorr_init, masscorr

  integer(4), save,              private :: ntr, ndry_max, ndry
  integer(4), save, allocatable, private :: idry(:), tndry(:,:)
  real(8),    save, allocatable, private :: buf(:,:), cor(:)

contains

  !-----------------------------------------------------------------------------

  subroutine masscorr_firstentry( nthr, ierr, kh, iw2, msscrr )
    implicit none

    integer, intent(in)    :: nthr, iw2, kh(0:)
    logical, intent(inout) :: msscrr
    integer, intent(inout) :: ierr

    integer(4) :: n
    
    if (iw2 > 0 .and. (.not.msscrr)) then
      do n=1,iw2
        if (kh(n) == 1) then
          msscrr = .true.
          return
        endif
      enddo

    elseif (iw2 <= 0 .and. msscrr) then
      ntr = nthr 
      allocate ( tndry(1:nthr,1:2), stat=ierr ) 
      if (ierr /= 0) return
      ndry_max = 0
    endif

  end subroutine masscorr_firstentry

  !-----------------------------------------------------------------------------

  subroutine masscorr_init(n2d, hnew, hold, hx, hy, hstab, u, v, dt, dx, dy,   &
                           cosphi, ind, msrf, kh, idx, t, s, ndry_act)

    use constants, only : hgr, heps, half
    use dmi_omp,   only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: n2d, ind(1:,1:), msrf(0:,0:), kh(0:)
    integer(4), intent(inout) :: idx(1:,1:), ndry_act
    real(8),    intent(in)    :: hnew(0:), hold(0:), u(0:), v(0:), t(0:), s(0:)
    real(8),    intent(in)    :: cosphi(1:,0:), hx(0:), hy(0:)
    real(8),    intent(in)    :: hstab, dt, dx, dy
    
    include 'masscorr_init.inc'

    real(8),    parameter :: hdry = hgr + heps
    integer(4), parameter :: NTS = 2
    integer(4)            :: ms, mw, mn, ii, i, j, tnum, itr, n2dl, n2du

    call domp_get_thread_no( tnum )
    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    !  find No. of dried points in this area/task:
    tndry(tnum,1) = 0
    do ms=n2dl,n2du
      if (hnew(ms) - hstab <= hdry) tndry(tnum,1) = tndry(tnum,1) + 1
    enddo

!$OMP BARRIER
   if (maxval(tndry(:,1)) > 0) then
!$OMP MASTER
      ndry       = tndry(1,1)
      tndry(1,2) = 0
      do itr=2,ntr
        ndry         = ndry + tndry(itr,1)
        tndry(itr,2) = tndry(itr-1,2) + tndry(itr-1,1)
      enddo
      !  poor man's re-allocator:
      if (ndry > ndry_max) then
        if (ndry_max > 0) deallocate( buf, cor, idry )
        allocate( buf(1:NTS,1:ndry), cor(1:ndry), idry(1:ndry) )
        ndry_max = ndry
      endif
!$OMP END MASTER
!$OMP BARRIER
      ndry_act = ndry
    else
      ndry_act = 0
    endif

    !  now, fill the work arrays:
    if (tndry(tnum,1) > 0) then
      ii = tndry(tnum,2)
      do ms=n2dl,n2du
        if (hnew(ms) - hstab > hdry) cycle
  
        i  = ind(1,ms)
        j  = ind(2,ms)
        mw = msrf(i,  j-1)
        mn = msrf(i-1,j  )

        ! save index:
        ii = ii + 1
        idry(ii) = ms

        ! save components:
        buf(1,ii) = t(ms)
        buf(2,ii) = s(ms)
       
        ! save correction:
        cor(ii) = (  (hnew(ms) - hold(ms))/dt                                  &
                   + ( ( hx(ms)*u(ms)                                          &
                        -hx(mw)*u(mw))/dx                                      &
                      +( hy(mn)*v(mn)*cosphi(2,i-1)                            &
                        -hy(ms)*v(ms)*cosphi(2,i)  )/dy)/cosphi(1,i) )         &
                  *half*dt/hnew(ms)
      enddo 
    endif

    if (ndry_act > 0) then
!$OMP BARRIER
    endif
 
  end subroutine masscorr_init

  !-----------------------------------------------------------------------------

  subroutine masscorr (t, s)

    use constants, only : one
    use dmi_omp,   only : domp_get_thread_no

    implicit none

    real(8),    intent(inout) :: t(0:), s(0:)

    include 'masscorr.inc'

    integer(4) :: id, mm, tnum
    real(8)    :: mcor

    call domp_get_thread_no( tnum )

    !  do the mass correction:
    if (tndry(tnum,1) > 0) then
      do id=tndry(tnum,2)+1,tndry(tnum,2)+tndry(tnum,1)
        mm   = idry(id)
        mcor = cor(id)
        t(mm) = (t(mm) + mcor*buf(1,id))/(one - mcor)
        s(mm) = (s(mm) + mcor*buf(2,id))/(one - mcor)
      enddo
    endif
  end subroutine masscorr

  !-----------------------------------------------------------------------------

end module masscorr_srf_col

