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

module io_miscsubs
  implicit none

  !- Public module vars and params
  !- 3d IO fields
  integer(4), parameter, public :: f_temp3 = 301  ! 3d temp field
  integer(4), parameter, public :: f_salt3 = 302  ! 3d salt field
  integer(4), parameter, public :: f_velu3 = 303  ! 3d u field
  integer(4), parameter, public :: f_velv3 = 304  ! 3d v field
  integer(4), parameter, public :: f_velw3 = 305  ! 3d w field
  integer(4), parameter, public :: f_tke3  = 306  ! 3d tke field
  integer(4), parameter, public :: f_diss3 = 307  ! 3d diss field
  integer(4), parameter, public :: f_avv3  = 308  ! 3d avv field

  !- 2d IO fields
  integer(4), parameter, public :: f_temp2 = 201  ! 2d temp field
  integer(4), parameter, public :: f_salt2 = 202  ! 2d salt field
  integer(4), parameter, public :: f_velu2 = 203  ! 2d u field
  integer(4), parameter, public :: f_velv2 = 204  ! 2d v field
  integer(4), parameter, public :: f_velw2 = 205  ! 2d w field
  integer(4), parameter, public :: f_tke2  = 206  ! 2d tke field
  integer(4), parameter, public :: f_diss2 = 207  ! 2d diss field
  integer(4), parameter, public :: f_avv2  = 208  ! 2d avv field
  integer(4), parameter, public :: f_sl    = 209  ! 2d sealevel
  integer(4), parameter, public :: f_wu    = 210  ! 2d zonal wind field
  integer(4), parameter, public :: f_wv    = 211  ! 2d meridional wind field
  integer(4), parameter, public :: f_ice1  = 212  ! 2d ice 
  integer(4), parameter, public :: f_ice2  = 213  ! 2d ice
  integer(4), parameter, public :: f_ice3  = 214  ! 2d ice
  integer(4), parameter, public :: f_tsnow = 215  ! 2d snow temp
  integer(4), parameter, public :: f_uice  = 216  ! 2d zonal ice vel
  integer(4), parameter, public :: f_vice  = 217  ! 2d meridional ice vel
  integer(4), parameter, public :: f_tsoil1= 218  ! 2d soil temp
  integer(4), parameter, public :: f_tsoil2= 219  ! 2d soil temp
  integer(4), parameter, public :: f_tsoil3= 220  ! 2d soil temp
  integer(4), parameter, public :: f_casus = 221  ! 2d ice cover

  !- collective IO fields
  integer(4), parameter, public :: f_tempdat  = 1001  ! all tempdat fields
  integer(4), parameter, public :: f_tempdat3 = 1002  ! all 3D tempdat fields
  integer(4), parameter, public :: f_tempdat2 = 1003  ! all 2D tempdat fields
  integer(4), parameter, public :: f_coupler  = 2001  ! all coupler fields

  !- io variables 
  integer(4), parameter, public :: max_io_fields_def = 50
  integer(4), public, save :: num_rfields3, num_rfields2
  integer(4), public, save :: num_rfields3_coupler = 0
  integer(4), public, save :: num_rfields2_coupler = 0

  interface sum_of_values
    module procedure priest_sum
    module procedure priest_sum_nc
  end interface


  !- private vars
  real(8), parameter :: sec2day = 86400.0_8
 
  public :: set_mo_n,calc_fileview, select_scatter_var, localTracerStats,      &
            depermute_output, set_predefined, set_restart, select_coupler_vars,&
            set_output_header, get_header_len, set_commtags 

contains

!===============================================================================

  subroutine get_header_len(tima, time, fields, freqout, dtmain, num_writes,   &
                            lheader)

    implicit none

    !- args
    real(8),    intent(in)  :: time, tima, dtmain
    integer(4), intent(in)  :: freqout, fields(:)
    integer(4), intent(out) :: lheader, num_writes

    !- locals
    real(8)    :: outtime, dt
    integer(4) :: num_iovars

    dt = dtmain / sec2day

    !- how many times do we write to the file?
    num_writes = 0
    outtime = tima + freqout*dt
    do while (outtime<=time)
      num_writes = num_writes + 1
      outtime    = outtime + freqout*dt
    enddo

    !- get number of output vars
    do num_iovars=1,max_io_fields_def
      if (fields(num_iovars) <= 0) exit
    enddo
    num_iovars = num_iovars - 1
    
    !- header length = (1 + output times) + (num_iovars + 1)
    lheader = num_writes + num_iovars + 2
  end subroutine get_header_len

!===============================================================================

  subroutine set_output_header(tima, time, dtmain, freqout, num_writes, iovars,&
                               headerbuff)
    implicit none

    !- args
!FIXME: time is never used for anything
    real(8),    intent(in)  :: tima, time, dtmain
    integer(4), intent(in)  :: iovars(:), freqout, num_writes
    real(8),    intent(out) :: headerbuff(:)

    !- locals
    real(8)              :: dt
    integer(4)           :: i, ibuf, num_iovars

    dt = dtmain / sec2day

    !- get number of output vars
    do num_iovars=1,max_io_fields_def
      if (iovars(num_iovars) <= 0) exit
    enddo
    num_iovars = num_iovars - 1
 
    !- write variable info
    headerbuff(1)              = real(num_iovars,8)
    headerbuff(2:num_iovars+1) = real(iovars(1:num_iovars),8)
    !- write number of timestamps
    ibuf = num_iovars+2
    headerbuff(ibuf) = num_writes
    ibuf = num_iovars+3

    !- save the actual output times
    do i = 1,num_writes
      headerbuff(ibuf)  = tima + i*freqout*dt
      ibuf = ibuf + 1
    enddo
  end subroutine set_output_header

!===============================================================================

  subroutine set_predefined(f, ofields, odims, coupler2, coupler3)
    use exits, only : exitme

    implicit none

    !- arguments
    integer(4), intent(in)           :: f
    integer(4), intent(inout)        :: ofields(:), odims(:)
    integer(4), intent(in), optional :: coupler2, coupler3

    !- local
    integer(4) :: i, up, low

    select case (f)
      case (f_tempdat)
        ofields(1:8) = (/f_wu, f_wv, f_velu3, f_velv3, f_velw3, f_sl,          &
                                f_temp3,f_salt3/)
        odims(1:8)   = (/ 2,2,3,3,3,2,3,3 /) 
      case (f_tempdat2)
        ofields(1:8) = (/f_wu, f_wv, f_velu2, f_velv2, f_velw2, f_sl,          &
                                f_temp2,f_salt2/)
        odims(1:8)   = (/ 2,2,2,2,2,2,2,2 /) 
      case (f_tempdat3)
        ofields(1:5) = (/f_velu3, f_velv3, f_velw3, f_temp3,f_salt3/)
        odims(1:5)   = (/ 3,3,3,3,3 /) 
      case (f_coupler)
        if (present(coupler2) .and. present(coupler3) ) then
          low = 1
          up  = coupler3
          do i=low,up
            ofields(i) = i
            odims(i) = 3
          enddo
  
          low = coupler3+1
          up  = coupler3+coupler2
          do i=low,up
            ofields(i) = i
            odims(i) = 2
          enddo
        else
          call exitme(1,'IO-Error: incorrect nr of args to set_predefined')
        endif
      end select
  end subroutine

!===============================================================================

  subroutine set_restart(ofields, odims, coupler2, coupler3)
    implicit none

    !- arguments
    integer(4), intent(inout)        :: ofields(:), odims(:)
    integer(4), intent(in), optional :: coupler2, coupler3

    !- local
    integer(4) :: i, up, low
    logical    :: coupler

    coupler = .false.
    if (present(coupler2) .and. present(coupler3)) coupler = .true.
    
    if (coupler) then
      num_rfields3_coupler = coupler3
      num_rfields2_coupler = coupler2

      low = 1
      up  = coupler3
      do i=low,up
        ofields(i) = i
        odims(i) = 3
      enddo

      low = coupler3+1
      up  = coupler3+coupler2
      do i=low,up
        ofields(i) = i
        odims(i) = 2
      enddo
    else
      num_rfields3 = 7
      num_rfields2 = 11
      up = num_rfields2 + num_rfields3
      ofields(1:up) = (/f_velu3,f_velv3,f_sl,f_salt3,f_temp3,f_casus,          &
                               f_ice1,f_ice2,f_ice3,f_tsnow,f_uice,f_vice,     &
                               f_tsoil1,f_tsoil2,f_tsoil3,f_tke3,f_diss3,      &
                               f_avv3/) 
      odims(1:up) = (/3,3,2,3,3,2,2,2,2,2,2,2,2,2,2,3,3,3/)
    endif

  end subroutine set_restart

!===============================================================================

  subroutine set_mo_n(ia, msrf, mcol, kh, mo, mop, iw2io, iw2g, hbm_lb, hbm_ub)
    !---------------------------------------------------------------------------
    ! Routine for calculating the permuted indicies. If we only have 1 
    ! IO-tasks, the mmo array is simply the permuted m1. However, do we have
    ! multiple IO-tasks, mmo will only contain the indicies for the HBM-tasks
    ! related to the specific IO-taks.
    !---------------------------------------------------------------------------
    use dmi_mpi,  only : mpi_tt
    use params_n, only : mmx,nmx,kmx
    
    implicit none

    !- args
    integer(4), intent(in)  :: ia, iw2io, iw2g, msrf(0:,0:), mcol(0:),         &
                               hbm_lb, hbm_ub
    integer(4), intent(in)  :: kh(0:)
    integer(4), intent(out) :: mo(1:,1:,1:), mop(1:,1:,1:)

    !- local vars
    integer(4) :: it, low_i, low_j, up_i, up_j, i, j, k, ms, ns, ms_all, ns_all

    ! init to nice values 
    mo(:,:,:)  = 0
    mop(:,:,:) = 0

    !- How many surface wet points are linked to this IO tasks
    !- NOTE:
    !    we also set mo(:,:,1) to 1 if this (i,j) is included in this io task
    !    so we dont have to look it up again in the second loop
    do j=1,nmx(ia)
      do i=1,mmx(ia) 
        ! wet point?
        if (msrf(i,j) == 0) cycle

        ! Is this grid point part of any HBM task linked to this IO task?
        do it=hbm_lb,hbm_ub
          low_i = mpi_tt(ia,it+1)%low_i
          low_j = mpi_tt(ia,it+1)%low_j
          up_i  = mpi_tt(ia,it+1)%up_i
          up_j  = mpi_tt(ia,it+1)%up_j
          
          !- cycle if it is not included
          if (i>=low_i .and. i<=up_i .and. j>=low_j .and. j<=up_j) then
            mop(i,j,1) = 1
            exit
          endif
        enddo !hbm
        
      enddo!i
    enddo!j

    !- Create permuted mm1 = mop and permuted IO-tasks mm1 = mo
    !- mop is overritten here  
    ns     = iw2io       ! subsurface index counter for the specific IO-tasks
    ms     = 0           ! surface index counter for the specific IO-task
    ms_all = 0           ! global counters for permuting mm1
    ns_all = iw2g        !- 
    do j=1,nmx(ia)
      do i=1,mmx(ia) 
        ! also cycle if water point
        if (msrf(i,j) == 0) cycle

        ms_all = ms_all + 1
        mo(i,j,1) = ms_all
        !- do subsurface layers
        
        do k=2,kh(msrf(i,j))
          ns_all = ns_all + 1
          mo(i,j,k) = ns_all
        enddo

        ! Is this grid point part of any HBM task linked to this IO task?
        if (mop(i,j,1) == 0) cycle 

        ms = ms + 1
        mop(i,j,1) = ms
        !- do subsurface layers
        do k=2,kh(msrf(i,j))
          ns = ns + 1
          mop(i,j,k) = ns
        enddo

      enddo!i
    enddo!j
  end subroutine set_mo_n

!===============================================================================

  subroutine depermute_output(ia,nf,mo,mop,vardims,n2dio,n3dio, hbm_lb, hbm_ub,&
                              rbuf, wbuf )
    !---------------------------------------------------------------------------
    ! Subroutine for depermuting the input buffer to the order it should be 
    ! writte to file
    !---------------------------------------------------------------------------
    !- modules
    use dmi_mpi,     only : mpi_tt
    use params_n, only : kmx
 
    implicit none

    !- arguments
    integer(4), intent(in)  :: ia, nf, mo(1:,1:,1:), mop(1:,1:,1:),            &
                               vardims(1:), n2dio, n3dio, hbm_lb, hbm_ub
    real(8), intent(in)     :: rbuf(1:)
    real(8), intent(out)    :: wbuf(1:)

    !- locals
    integer(4) :: ibuf, it, low_i, low_j, up_i, up_j, i, j, k, mip, ivar,      &
                  sumvaroffset
    !FIXME: it can be threaded by e.g. splitting the hbm tasks or variables
    ibuf = 0 
    do it=hbm_lb,hbm_ub
      if (mpi_tt(ia,it+1)%nwet3 <= 0) cycle
      
      low_i = mpi_tt(ia,it+1)%low_i
      low_j = mpi_tt(ia,it+1)%low_j
      up_i  = mpi_tt(ia,it+1)%up_i
      up_j  = mpi_tt(ia,it+1)%up_j

      sumvaroffset = 0
      do ivar = 1,nf
        do j = low_j,up_j
          do i = low_i,up_i
            if ( mo(i,j,1) == 0 ) cycle
            ibuf = ibuf + 1 
            mip  = mop(i,j,1) + sumvaroffset
   
            wbuf(mip) = rbuf(ibuf)
          enddo
        enddo
          
        if (vardims(ivar) == 2) then
          sumvaroffset = sumvaroffset + n2dio
          cycle
        endif
 
        do j = low_j,up_j
          do i = low_i,up_i
            do k=2,kmx(ia) !kh not available
              if ( mo(i,j,k) == 0 ) exit
              ibuf = ibuf + 1 
              mip = mop(i,j,k) + sumvaroffset

              wbuf(mip) = rbuf(ibuf)
  
            enddo !k
          enddo !i
        enddo !j
        sumvaroffset = sumvaroffset + n3dio
      enddo ! ivar
    enddo !it

  end subroutine depermute_output

!===============================================================================

  subroutine calc_fileview(ia,nf,mo,mop,blen,filedispl,vardims,n2dio, n3dio,   &
                            n2dg, n3dg, hbm_lb, hbm_ub, iibloc)
    !---------------------------------------------------------------------------
    ! Subroutine for calculating the permutation of the reciev buffer for 
    ! each IO tasks.
    !---------------------------------------------------------------------------
    !- modules
    use dmi_mpi,     only : mpi_tt
    use params_n, only : kmx
    implicit none

    !- arguments
    integer, intent(in)  :: ia, nf, mo(1:,1:,1:), mop(1:,1:,1:), vardims(1:),  &
                            n2dio, n3dio, n2dg, n3dg, hbm_lb, hbm_ub
    integer, intent(out) :: filedispl(1:), iibloc, blen(1:)

    !- locals
    integer(4) :: low_i, low_j, up_i, up_j, i, j, k, skip, mip, mio, hbmlen,   &
                  ivar, sumvaroffset, sumvaroffset_glob, prevbloclen, ibloc
    logical    :: first
    
    ibloc = 1
    filedispl(ibloc) = 0
    prevbloclen = 0
    skip   = 0 ! only matters if we have multiple io tasks

    !- need to take care of dry tasks so cannot simply take minval/maxval
    hbmlen = hbm_ub-hbm_lb+1
    first = .true.
    do i=1,hbmlen
      j = hbm_lb + i
      if (first .and. mpi_tt(ia,j)%nwet3 > 0) then
        low_i = mpi_tt(ia,j)%low_i
        low_j = mpi_tt(ia,j)%low_j
        up_i  = mpi_tt(ia,j)%up_i
        up_j  = mpi_tt(ia,j)%up_j
        first = .false.
        cycle
      endif
      
      if (mpi_tt(ia,j)%nwet3 > 0) then
        low_i = min(low_i,mpi_tt(ia,j)%low_i)
        low_j = min(low_j,mpi_tt(ia,j)%low_j)
        up_i  = max(up_i,mpi_tt(ia,j)%up_i)
        up_j  = max(up_j,mpi_tt(ia,j)%up_j)
      endif
    enddo

    !- how many surf points should we skip before we start?
    JLOOP: do j = low_j,up_j
      do i = low_i,up_i
        if ( mo(i,j,1) /= 0 ) then
          filedispl(1) = mo(i,j,1)-1
          skip = mo(i,j,1)-mop(i,j,1)
          exit JLOOP
        endif
      enddo
    enddo JLOOP

    sumvaroffset      = 0
    sumvaroffset_glob = 0
    do ivar = 1,nf
      do j = low_j,up_j
        do i = low_i,up_i
          if ( mo(i,j,1) == 0 ) cycle
          ! what is the unpermuted grid point at this (i,j)
          mip = mop(i,j,1) + sumvaroffset
          mio = mo(i,j,1)  + sumvaroffset_glob
 
          if (mip+skip<mio) then
            ibloc            = ibloc + 1
            filedispl(ibloc) = mio-1
            blen(ibloc-1)    = prevbloclen
            skip             = mio - mip
            prevbloclen      = 1
          else
            prevbloclen = prevbloclen+1 ! index in reciev buffer
          endif
        enddo
      enddo
        
      if (vardims(ivar) == 2) then
        sumvaroffset      = sumvaroffset + n2dio
        sumvaroffset_glob = sumvaroffset_glob + n2dg
        cycle
      endif

      do j = low_j,up_j
        do i = low_i,up_i
          do k=2,kmx(ia) !kh not available
            if ( mo(i,j,k) == 0 ) exit

            mip = mop(i,j,k) + sumvaroffset
            mio = mo(i,j,k) + sumvaroffset_glob

            if (mip+skip<mio) then 
              ibloc            = ibloc + 1
              filedispl(ibloc) = mio-1
              blen(ibloc-1)    = prevbloclen
              skip             = mio - mip
              prevbloclen      = 1
            else
              prevbloclen  = prevbloclen+1 ! index in reciev buffer
            endif 
          enddo !k
        enddo !i
      enddo !j
      sumvaroffset      = sumvaroffset + n3dio
      sumvaroffset_glob = sumvaroffset_glob + n3dg
    enddo ! ivar

    blen(ibloc)    = prevbloclen
    if (prevbloclen == 0) then
      iibloc = ibloc-1
    else
      iibloc = ibloc
    endif
  end subroutine calc_fileview

!===============================================================================
  
  subroutine select_coupler_vars(tmpbuf, selvar, buflow, bufup, varlow,        &
                                 varup, coupler_var_3, coupler_var_2)
    !- modules
    use dmi_omp,            only : domp_get_domain

    implicit none

    !- args
    integer(4), intent(in)    :: buflow, bufup, varlow, varup, selvar
    real(8),    intent(in)    :: coupler_var_3(:,0:), coupler_var_2(:,0:)
    real(8),    intent(inout) :: tmpbuf(1:)


    !- locals
    integer(4) :: i, blow, bup, vlow, vup, nt3 !, nt2
!FIXME: nt2 set but never referenced

    nt3 = ubound(coupler_var_3,1)
!   nt2 = ubound(coupler_var_2,1)

    call domp_get_domain(buflow, bufup, blow, bup)
    call domp_get_domain(varlow, varup, vlow, vup)

    if (selvar <= nt3 ) then
      i = selvar
      tmpbuf(blow:bup) = coupler_var_3(i,vlow:vup)
    else
      i = selvar - nt3
      tmpbuf(blow:bup) = coupler_var_2(i,vlow:vup)
    endif
  end subroutine select_coupler_vars

!===============================================================================

  subroutine select_scatter_var(ia, tmpbuf,selvar, buflow, bufup, varlow, varup)
    !- modules 
    use local_arrays_n,  only : s_l, t_l
    use dmi_omp,         only : domp_get_domain
    implicit none

    !- args
    integer(4), intent(in)    :: ia, buflow, bufup, varlow, varup, selvar
    real(8),    intent(inout) :: tmpbuf(1:)

    !- locals
    real(8), parameter :: one = 1.0_8, zero = 0.0_8
    integer(4)         :: i, slicelen, blow, bup, vlow, vup

    call domp_get_domain(buflow, bufup, blow, bup)
    call domp_get_domain(varlow, varup, vlow, vup)
    
    select case(selvar)
      case (f_temp3, f_temp2)
        tmpbuf(blow:bup) = t_l(ia)%p(vlow:vup)
      case (f_salt3, f_salt2)
        tmpbuf(blow:bup) = s_l(ia)%p(vlow:vup)
    end select
  end subroutine select_scatter_var


!===============================================================================
 subroutine fast2sum(a,b,s,t)
    implicit none
    real(8), intent(in)  :: a,b
    real(8), intent(out) :: s,t
    real(8) :: z
    if (b>a) then
      s=a+b
      z=s-b
      t=a-z
    else
      s=a+b
      z=s-a
      t=b-z
    endif
  end subroutine fast2sum
  subroutine priest_sum(nsize,x,sumout)
    ! double compensated sum
    use constants, only : zero
    implicit none
    integer(4), intent(in)  :: nsize
    real(8),    intent(in)  :: x(0:)
    real(8),    intent(out) :: sumout
    real(8)    :: c, u, y, t, s, v, z
    integer(4) :: i
    s = x(1)
    c = zero
    do i=2,nsize
      call fast2sum(c,x(i),y,u)
      call fast2sum(s,y,t,v)
      z = u + v 
      call fast2sum(t,z,s,c)
    enddo
    sumout = s
  end subroutine priest_sum

  subroutine priest_sum_nc(nsize,nc,x,sumout)
    ! double compensated sum
    use constants, only : zero
    implicit none
    integer(4), intent(in)  :: nsize, nc
    real(8),    intent(in)  :: x(:,0:)
    real(8),    intent(out) :: sumout
    real(8)    :: c, u, y, t, s, v, z
    integer(4) :: i
    s = x(nc,1)
    c = zero
    do i=2,nsize
      call fast2sum(c,x(nc,i),y,u)
      call fast2sum(s,y,t,v)
      z = u + v 
      call fast2sum(t,z,s,c)
    enddo
    sumout = s
  end subroutine priest_sum_nc

  subroutine localTracerStats(ia,var,varmin, varmax, do3d, w_varmean)
    use params_n,   only : iw3, iw2
    use dmi_mpi,       only : dd
    implicit none

    !- arguments
    real(8),    intent(in)            :: var(0:)
    integer(4), intent(in)            :: ia
    logical,    intent(in)            :: do3d
    real(8),    intent(out)           :: varmin, varmax
    real(8),    intent(out), optional :: w_varmean

    !- locals
    real(8)    :: partsum(0:3), frac, tmp1, tmp2
    integer(4) :: n3dg, n2dg, n3d, n2d, n3dlb, n3dub, n3dlen
    logical    :: do_mean

    ! we must take care of the halo points after 2d, so sum is splitted
    ! into surf and subsurf points. The are then added together and divided
    ! with iw3. Last, it is weighted with n3d_local/n3d_global
 
    n3dg   = iw3(ia)
    n2dg   = iw2(ia)
    n3d    = dd(ia)%nwet3
    n2d    = dd(ia)%up_ws - dd(ia)%low_ws + 1
    n3dlb  = n2d+1 + dd(ia)%halo2
    n3dlen = n3d - n2d
    n3dub  = n3dlb + n3dlen - 1
 
    if (present(w_varmean)) then
      do_mean = .true.
    else
      do_mean = .false.
    endif

    !- MIN/MAX
    if (do3d) then
      tmp1   = maxval( var(1:n2d) )!srf
      tmp2   = maxval( var(n3dlb:n3dub) ) !subsrf
      varmax = max(tmp1,tmp2)
      tmp1   = minval( var(1:n2d) )
      tmp2   = minval( var(n3dlb:n3dub) )
      varmin = min(tmp1,tmp2)
    else
      varmin = minval( var(1:n2d) )
      varmax = maxval( var(1:n2d) )
    endif
 
    if (do_mean) then
      if (do3d) then
        call sum_of_values(n2d,var(0:n2d),partsum(1))
        call sum_of_values(n3dlen,var(n3dlb-1:n3dub),partsum(2)) 
        call sum_of_values(2,partsum,w_varmean)
        frac = 1/real(n3dg,8)
        w_varmean = w_varmean*frac
      else
        call sum_of_values(n2d,var(0:n2d),w_varmean)
        frac = 1/real(n2dg,8)
        w_varmean = w_varmean*frac
      endif
    endif
  end subroutine localTracerStats

!===============================================================================

  subroutine set_commtags(commtags, narea, maxfh, nbufs, base)
    implicit none

    !- args
    integer(4), intent(in)  :: narea, maxfh, nbufs, base
    integer(4), intent(out) :: commtags(:,:,:)

    !- locals
    integer(4) :: ia, ifh, buf, itag

    itag = base
    do ia=1,narea
      do ifh=1,maxfh
        do buf=1,nbufs
          itag = itag + 1
          commtags(ia,ifh,buf) = itag
        enddo
      enddo
    enddo

  end subroutine set_commtags

!===============================================================================

end module io_miscsubs
