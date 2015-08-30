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

module params_n
  !- Modules -------------------------------------------------------------------
  use exits,     only : exitme
  use cmod_mem,  only : cmi1, cmi2, cmr1
  use constants, only : nproci,nprocj,decomp_version,only_islices,             &
                        decomp_coeff,ionml,iot

  !- Implicit Directives -------------------------------------------------------
  implicit none
  private
  integer(4),    save, public :: i2, izeite
  integer(4), save, public :: narea=1 
  integer(4), save, public :: ntotal=1 
  integer(4), save, public :: niter=500 
  logical,    save, public :: io=.false.
  real(8),    save, public :: dx=1836.0_8 
  real(8),    save, public :: dy=1836.0_8 
  real(8),    save, public :: ui=0.0_8 
  real(8),    save, public :: vi=0.0_8 
  real(8),    save, public :: wi=0.0_8 
  logical,    save, public :: ldiff=.false.
  logical,    save, public :: ladv=.true.
  logical,    save, public :: cube=.false.
  type(cmi2), allocatable, save, public :: msrf(:)
  type(cmi1), allocatable, save, public :: kh(:),  mcol(:)
  type(cmr1), allocatable, save, public :: hz(:), hz_f(:)
  integer(4), allocatable, save, public :: mmx(:), nmx(:), kmx(:), iw2(:),iw3(:)
  integer(4), allocatable, save, public :: mmxp(:), nmxp(:) ! mxp = mx + 1
  integer(4),              save, public :: kmax, kmaxp      ! kmaxp=kmax+1
  integer(4), parameter, public :: mainarea = 1  ! No. corresponding to the main
  !- OMP params
  integer(4), save, public :: nthreads    ! #threads
  !- Subroutines/functions -----------------------------------------------------
  public  :: getspecs
  private :: PrintParams
contains

  subroutine getspecs()
    !- Modules -----------------------------------------------------------------
    use io_subs,        only : io_new_unit
    use dmi_mpi_global, only : iu06

    !- local vars --------------------------------------------------------------
    integer(4)         :: lun, ios, i, mi,mj,mk, iw, j, iwc
    character(256)     :: fnam
    integer(4), parameter :: ia=1
    character(5), parameter :: w = '(a57)'

    namelist /optionlist/ mi,mj,mk,dx,dy,cube,ionml,niter,iot, io,             &
                          nproci, nprocj, only_islices, decomp_version,        &
                          ui, vi, wi, ldiff, ladv

   decomp_coeff   = (/0.0_8, 1.0_8, 0.0_8/) 

    !-----Defaults -------------------------------------------------------------
  
    lun = io_new_unit()
    open (unit=lun, file='options.nml', status='old', iostat=ios)
    if (ios /= 0) then
      write(iu06,w)  'Warning: Cannot open file options.nml ...                '
      write(iu06,w)  '         ... continue with default options.              '
      write(iu06,'(a25,i5)') 'I/O problem - iostat is: ',ios
    else
      read (unit=lun, nml=optionlist, iostat=ios)
      if (ios /= 0) then
        write(iu06,w)'Warning: Cannot read optionlist in options.nml file ...  '
        write(iu06,w)'         ... continue without                            '
        write(iu06,'(a25,i5)') 'I/O problem - iostat is: ',ios
      endif
      close (unit=lun)
    endif
    write(iu06,*) '  kernel iterations (niter) : ', niter
    !  total No. of domains:
    allocate( mmx(ntotal), nmx(ntotal),kmx(ntotal),iw2(ntotal),iw3(ntotal),    &
              mmxp(ntotal), nmxp(ntotal),                                      &
              hz(ntotal), hz_f(ntotal), mcol(ntotal), msrf(ntotal),kh(ntotal), &
              stat=ios)
    if (ios /= 0) then
      write(iu06,*) 'Error: Failed to allocate params, status was: ', ios
      call exitme(1)
    endif
    mmx(ia)=mi
    nmx(ia)=mj
    kmx(ia)=mk

    allocate( msrf(ia)%p(0:mmx(ia)+1,0:nmx(ia)+1), stat=ios )
    if (ios /= 0) then
      write(iu06,*) 'Error: Allocate failure was:: ', ios
      call exitme(1)
    endif
    if (cube) then
      msrf(ia)%p(0:,0:) = 0
      iw=1
      do j=1,nmx(ia)
        do i=1,mmx(ia)
          msrf(ia)%p(i,j) = iw
          iw = iw + 1
        enddo
      enddo
      iw=iw-1
      iw2(ia) = maxval(msrf(ia)%p(:,:))
      allocate( kh(ia)%p(0:iw2(ia)), mcol(ia)%p(0:iw2(ia)), stat=ios )
      if (ios /= 0) then
        write(iu06,*) 'Error: Allocating kh, mcol will exit'
        write(iu06,*) 'Error: The allocation failure was:: ', ios
        call exitme(1)
      endif
      mcol(ia)%p(0:iw2(ia))=0
      kh(ia)%p(0) = 0
      do iw=1,iw2(ia) 
        kh(ia)%p(iw) = mk
      enddo
      iwc=iw2(ia)+1
      do iw=1,iw2(ia) 
        mcol(ia)%p(iw) = iwc
        iwc=iwc+mk-1
      enddo
      if (mk .ge. 2) then
        iw3(ia) = maxval(mcol(ia)%p(:))+mk-2
      else
        iw3(ia) = max(maxval(mcol(ia)%p(:)),iw2(ia))
      endif
      call PrintParams (iw2(ia), iw3(ia))
      allocate( hz(ia)%p(0:iw3(ia)), hz_f(ia)%p(0:iw3(ia)) )
      hz(ia)%p(0) = 0
      hz(ia)%p(1:iw3(ia)) = 1.0_8
    else
      lun = io_new_unit()
      fnam = 'nidx.bin'  ! FIXME avoid hardcoded name
      open(lun, file=trim(fnam), form='unformatted', action='read',            &
         status='old', access='stream', iostat=ios)
      if (ios /= 0) then
        write(iu06,*) 'Error: Cannot open binary file '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      read(lun,iostat=ios) msrf(ia)%p(0:mmx(ia)+1,0:nmx(ia)+1)
      if (ios /= 0) then
        write(iu06,*) 'Error: Cannot read msrf from '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      iw2(ia) = maxval(msrf(ia)%p(:,:))
      allocate( kh(ia)%p(0:iw2(ia)), mcol(ia)%p(0:iw2(ia)), stat=ios )
      if (ios /= 0) then
        write(iu06,*) 'Error: Allocate failure was:: ', ios
        call exitme(1)
      endif
      read(lun,iostat=ios) mcol(ia)%p(0:iw2(ia))
      if (ios /= 0) then
        write(iu06,*) 'Error: Cannot read mcol from '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      read(lun,iostat=ios) kh(ia)%p(0:iw2(ia))
      if (ios /= 0) then
        write(iu06,*) 'Error: Cannot read kh from '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      close(lun)
      iw3(ia) = 0
      do i=1,iw2(ia)
        iw3(ia) = iw3(ia)+kh(ia)%p(i)
      enddo
      call PrintParams (iw2(ia), iw3(ia))
      lun = io_new_unit()
      fnam  = 'hz.bin'  ! FIXME avoid hardcoded name
      open(lun, file=trim(fnam), form='unformatted', action='read',            &
         status='old', access='stream', iostat=ios)
      if (ios /= 0) then
        write(iu06,*) 'Error: Cannot open binary file '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      allocate( hz(ia)%p(0:iw3(ia)), hz_f(ia)%p(0:iw3(ia)) )
      read(lun,iostat=ios) hz(ia)%p
      if (ios /= 0) then
          write(iu06,*) 'Error: read hz from binary file '//trim(fnam)//' exit'
          write(iu06,*) 'Error: The iostat failure was:: ', ios
          call exitme(1)
      endif
      close(lun)
      endif
  end subroutine

  subroutine PrintParams (n2d, n3d)
    use dmi_mpi_global, only : iu06

    !- Directives --------------------------------------------------------------
    implicit none

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in) :: n2d, n3d

    !- local vars --------------------------------------------------------------
    !  none, so far

    !- echo parameters ---------------------------------------------------------
    write(iu06,*) '  surface wet points:       ', n2d
    write(iu06,*) '  wet points:               ', n3d

    !- Do some sanity checking -------------------------------------------------
    if (n3d <= 0)  call exitme(1, 'Invalid number of wet points')
    if (n2d <= 0)  call exitme(1, 'Invalid number of surface wet points')
    if (n2d > n3d) call exitme(1, 'Invalid model setup (hmm?)')
  end subroutine PrintParams
end module params_n

module local_arrays_n
  !- Modules -------------------------------------------------------------------
  use cmod_mem, only : cmi1, cmi2, cmr1, cmr2

  !- Directives ----------------------------------------------------------------
  implicit none
  private

  !- Dimensions ----------------------------------------------------------------
  integer(4), allocatable, save, public :: iw2_l(:), iw3_l(:)
  real(8),    allocatable, save, public :: dx_l(:), dy_l(:)

  type(cmi2), allocatable, save, public :: msrf_l(:), msrf_lf(:)
  type(cmi2), allocatable, save, public :: ind_l(:), ind_lf(:)
  type(cmi1), allocatable, save, public :: kh_l(:),kh_lf(:),mcol_l(:),mcol_lf(:)
  type(cmi1), allocatable, save, public :: khu_l(:),khv_l(:)

  type(cmr1), allocatable, save, public :: h_l(:)
  type(cmr1), allocatable, save, public :: hu_l(:), hv_l(:), hz_l(:)
  type(cmr1), allocatable, save, public :: h_old_l(:), h_new_l(:)

  type(cmr1), allocatable, save, public :: ui_l(:), vi_l(:), wi_l(:)
  type(cmr1), allocatable, save, public :: eddyh_l(:)

  type(cmr1), allocatable, save, public :: dispt_l(:), disps_l(:)
  type(cmr1), allocatable, save, public :: t_l(:), s_l(:)
  type(cmr2), allocatable, save, public :: cosphi_l(:)
  type(cmr2), allocatable, save, public :: uvdam_l(:)

  !- Subroutines/functions -----------------------------------------------------
  public  :: Alloc_Local_Arrays
  private :: Alloc_Local

contains

  subroutine Alloc_Local_Arrays (narea, nthr)

    !- modules -----------------------------------------------------------------
    use dmi_mpi,       only : dd
    use exits,         only : exitme
    use tflow_simd_srf_col,     only : tflow_alloc_simd
    use masscorr_srf_col,only : masscorr_firstentry
    use constants,only : msscrr
 
    !- dircetives --------------------------------------------------------------
    implicit none

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in) :: narea, nthr

    !- local vars --------------------------------------------------------------
    integer(4) :: mnp, nnp, mn, nn, n2d, n3d, ia, n2s, n3s, n2h, n3h
    integer(4) :: maxn3s, maxn2s, ierr

    !- initialise --------------------------------------------------------------
    ierr = 0

    !- Allocate all the pointers -----------------------------------------------
    call Alloc_Local( narea, ierr )
    if (ierr /= 0) call exitme(-1,'Alloc_Local_Arrays failed, Alloc_Local')

    !- Allocate the local arrays -----------------------------------------------
    maxn3s = -1
    maxn2s = -1
    do ia=1,narea

      !- Full-size local arrays ------------------------------------------------
      !
      !   dim nos for the global quantities:
      mn  = dd(ia)%dim1
      nn  = dd(ia)%dim2
      mnp = mn + 1
      nnp = nn + 1
      n2d = dd(ia)%iwet2
      n3d = dd(ia)%iwet3
      !
      !   allocate full-size local arrays: - FIXME - should be local size
      allocate( msrf_l(ia)%p(0:mnp,0:nnp), msrf_lf(ia)%p(0:mnp,0:nnp),         &
                cosphi_l(ia)%p(1:2,0:mnp),stat=ierr)
      if (ierr /= 0) call exitme(-1,'Alloc_Local_Arrays failed, mm arrays')

      !- Task-size local arrays ------------------------------------------------
      !
      !   dim nos for the local 2D and 3D quantities:
      n2h = dd(ia)%halo2
      n3h = dd(ia)%halo3
      n3d = dd(ia)%nwet3
      if (n3d == 0) then
        if (dd(ia)%up_ws /= 0 .or. dd(ia)%low_ws /= 0 .or.                     &
            dd(ia)%halo2 /= 0 .or. dd(ia)%halo3  /= 0      )                   &
            call exitme(1,'Cannot figure out how to set up local arrays')
        n2d  = 0 
      else
        n2d  = dd(ia)%up_ws - dd(ia)%low_ws + 1
      endif
      n2s    = n2d + n2h
      n3s    = n3d + n3h
      maxn3s = max(n3s,maxn3s)
      maxn2s = max(n2s,maxn2s)
 
      allocate( hz_l(ia)%p(0:n3s), h_l(ia)%p(0:n3s),         &
                ui_l(ia)%p(0:n3s), vi_l(ia)%p(0:n3s), wi_l(ia)%p(0:n3s),       &
                hu_l(ia)%p(0:n3s), hv_l(ia)%p(0:n3s),                          &
                h_old_l(ia)%p(0:n3s),h_new_l(ia)%p(0:n3s),                     &
                kh_l(ia)%p(0:n2s),  kh_lf(ia)%p(0:n2s),khu_l(ia)%p(0:n2s),     &
                khv_l(ia)%p(0:n2s), t_l(ia)%p(0:n3s), s_l(ia)%p(0:n3s),        &
                eddyh_l(ia)%p(0:n3s),dispt_l(ia)%p(0:n3s),disps_l(ia)%p(0:n3s),&
                ind_l(ia)%p(1:2,1:n2s), uvdam_l(ia)%p(1:2,0:n2s),              &
                ind_lf(ia)%p(1:2,1:n2s),                                       &
                mcol_l(ia)%p(0:n2s), mcol_lf(ia)%p(0:n2s), stat=ierr )
      if (ierr /= 0) then
        call exitme(-1,'Alloc_Local_Arrays failed, task local arrays')
      endif

      !   transfer values:
      iw2_l(ia) = n2d  
      iw3_l(ia) = n3d  
    enddo
    !  tflow tmp arrays:
!    if (msscrr) call masscorr_firstentry( nthr, 0, kh_l(ia)%p, 0, msscrr )
    if (ierr /= 0) call exitme(-1,'Alloc_Local_Arrays failed, mass corr')
    call tflow_alloc_simd(maxn3s, ierr)
    if (ierr /= 0) call exitme(-1,'Alloc_Local_Arrays failed, tflow_alloc')
  end subroutine Alloc_Local_Arrays

  subroutine Alloc_Local (n, ierr)
    use params_n, only : dy, dx
    implicit none
    integer(4), intent(in)    :: n
    integer(4), intent(inout) :: ierr
    allocate( msrf_l(n), mcol_l(n), iw2_l(n), iw3_l(n), hz_l(n),h_l(n),        &
              hu_l(n), hv_l(n), kh_l(n), ind_l(n),                             &
              khu_l(n), khv_l(n), h_old_l(n), h_new_l(n),                      &
              dispt_l(n), disps_l(n), t_l(n), s_l(n),                          &
              ui_l(n), vi_l(n), wi_l(n), cosphi_l(n),                          &
              eddyh_l(n), dx_l(n), dy_l(n), uvdam_l(n), kh_lf(n),              &
              msrf_lf(n), ind_lf(n), mcol_lf(n), stat=ierr )
    dy_l(1) = dy
    dx_l(1) = dx
  end subroutine Alloc_Local

end module local_arrays_n

module arrays_n

  !- Modules -------------------------------------------------------------------
  use cmod_mem, only : cmi2, cmr2, cmr3

  !- Directives ----------------------------------------------------------------
  implicit none
  private
  type(cmi2), allocatable, save, public :: idx(:)
  type(cmi2), allocatable, save, public :: kru(:), krv(:)
  type(cmi2), allocatable, save, public :: krq(:), krz(:)
  type(cmr3), allocatable, save, public :: bndz(:), bndq(:)
  type(cmr2), allocatable, save, public :: rwqk(:)

  !- Subroutines/functions -----------------------------------------------------
  public  :: AllocArrays
  private :: AllocPointers 

contains

  subroutine AllocArrays (ia, nthreads) 
    !- modules -----------------------------------------------------------------
    use params_n, only : kmx,ntotal
    use exits,         only : exitme
 
    !- dircetives --------------------------------------------------------------
    implicit none
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in) :: ia, nthreads
    !- local vars --------------------------------------------------------------
    logical, save :: FirstEntry = .true.
    integer(4)    :: kn, ierr
    ierr = 0
    if (FirstEntry) then
      FirstEntry = .false.
      call AllocPointers( ntotal, ierr )
      if (ierr /= 0) call exitme(-1,'AllocArrays failed, AllocPointers')
      kn =  kmx(ia) 
      allocate( krz(ia)%p(1:3,0:0), krq(ia)%p(2,0:0), kru(ia)%p(4,0:0),        &
                krv(ia)%p(4,0:0), bndz(ia)%p(1:2,1:kn,0:0),                    &
                bndq(ia)%p(1:2,1:kn,0:0), rwqk(ia)%p(0:0,1:kn),                &
                idx(ia)%p(1:2,1:nthreads), stat=ierr )
    endif
  end subroutine AllocArrays 

  subroutine AllocPointers (n, ierr)
    implicit none
    integer(4), intent(in)    :: n
    integer(4), intent(inout) :: ierr
    allocate( krz(n), krq(n), kru(n), krv(n),bndz(n),bndq(n),rwqk(n),          &
              idx(n),stat=ierr )
  end subroutine AllocPointers


end module arrays_n

module init_n
 !- Directives ----------------------------------------------------------------
  implicit none
  private
  interface init_local_vars
    module procedure init_local_vars_ser
    module procedure init_local_vars_par
  end interface 
  !- Subroutines/functions -----------------------------------------------------
  public  :: init_local_vars, Validate
  private :: init_local_vars_ser, init_local_vars_par

  interface sum_of_values
    module procedure priest_sum
    module procedure priest_sum_nc
  end interface
  interface sum_of_squares
    module procedure sum_of_squares_1
  end interface
  

contains

  subroutine init_local_vars_ser(ia,kh_g,kh,hz_g,hz,mcol_g,mcol,msrf_g,        &
                                 msrf,ind,cosphi)
    use dmi_mpi,   only : dd
    implicit none
    !- args
    integer(4), intent(in)    :: ia
    integer(4), intent(in)    :: msrf_g(0:,0:),mcol_g(0:),kh_g(0:)
    real(8),    intent(in)    :: hz_g(0:)
    integer(4), intent(out)   :: ind(1:,1:),msrf(0:,0:),mcol(0:),kh(0:)
    real(8), intent(out)      :: hz(0:)
    real(8), intent(out)      :: cosphi(1:,0:)
    !- locals
    integer(4) :: kb
    integer(4) :: il, iu, jl, ju, ilh, iuh, jlh, juh, i, j, k
    integer(4) :: ns, nc, ms, mc, nwet2, m0

    il = dd(ia)%low_i
    iu = dd(ia)%up_i
    jl = dd(ia)%low_j
    ju = dd(ia)%up_j
    if (max(il, iu, jl, ju) == 0) return

! STEP1: setup temporary local msrf, mcol, ind, hz, kh
    !- Return if this task has nothing to do in this area ----------------------
    !  incl. halo zone:
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj

    ind(1:2,1:)    = 0
    mcol(0:)       = 0
    kh(0)          = 0
    msrf(0:,0:)    = 0
    cosphi(1:2,0:) = 0.5_8  ! FIXME
    cosphi(1:2,0)  = 0

    !- do this task's region ---------------------------------------------------
    !    ns: surface index according to local permuation
    !    ms: surface index according to global permuation
    !    nc: index in the water column below surface, local permutation
    !    mc: index in the water column below surface, global permutation
    nwet2 = dd(ia)%up_ws - dd(ia)%low_ws + 1
    ns    = 0
    nc    = nwet2 + dd(ia)%halo2
    do j=jl,ju
      do i=il,iu
        ms = msrf_g(i,j) 
        if (ms <= 0) cycle
        ! handle surface points:
        ns = ns + 1
        kb = kh_g(ms)
        ind(1,ns)  = i
        ind(2,ns)  = j
        kh(ns)       = kb
        msrf(i,j)    = ns
        hz(ns)       = hz_g(ms)
        ! handle sub-surface points:
        m0 =  mcol_g(ns)-2
        if (kb > 1) mcol(ns) = nc + 1
        do k=2,kb
          nc = nc + 1
          mc = m0+k
          hz(nc)     = hz_g(mc)
        enddo
      enddo
    enddo

    !- add the halo at the tail of this task's data ----------------------------
    ns = nwet2
    nc = dd(ia)%nwet3 + dd(ia)%halo2
    !    west + n/w + s/w:
    do j=jlh,jl-1
      do i=ilh,iuh
        ms = msrf_g(i,j) 
        if (ms <= 0) cycle
        ns = ns + 1
        kb = kh_g(ms)
        ind(1,ns)  = i
        ind(2,ns)  = j
        kh(ns)     = kb
        msrf(i,j)  = ns
        hz(ns)     = hz_g(ms)

        m0 =  mcol_g(ns)-2
        if (kb > 1) mcol(ns) = nc + 1
        do k=2,kb
          nc = nc + 1
          mc = m0+k
          hz(nc)     = hz_g(mc)
        enddo
      enddo
    enddo
    !    east + n/e + s/e:
    do j=ju+1,juh
      do i=ilh,iuh
        ms = msrf_g(i,j) 
        if (ms <= 0) cycle
        ns = ns + 1
        kb = kh_g(ms)
        ind(1,ns)  = i
        ind(2,ns)  = j
        kh(ns)     = kb
        msrf(i,j)  = ns
        hz(ns)     = hz_g(ms)
        
        m0 =  mcol_g(ns)-2
        if (kb > 1) mcol(ns) = nc + 1
        do k=2,kb
          nc = nc + 1
          mc = m0+k
          hz(nc)     = hz_g(mc)
        enddo
      enddo
    enddo
    !    north + south:
    do j=jl,ju
      !  north:
      do i=ilh,il-1
        ms = msrf_g(i,j) 
        if (ms <= 0) cycle
        ns = ns + 1
        kb = kh_g(ms)
        ind(1,ns)  = i
        ind(2,ns)  = j
        kh(ns)     = kb
        msrf(i,j)  = ns
        hz(ns)     = hz_g(ms)

        m0 =  mcol_g(ns)-2
        if (kb > 1) mcol(ns) = nc + 1
        do k=2,kb
          nc = nc + 1
          mc = m0+k
          hz(nc)     = hz_g(mc)
        enddo
      enddo
      !  south:
      do i=iu+1,iuh
        ms = msrf_g(i,j) 
        if (ms <= 0) cycle
        ns = ns + 1
        kb = kh_g(ms)
        ind(1,ns)  = i
        ind(2,ns)  = j
        kh(ns)     = kb
        msrf(i,j)  = ns
        hz(ns)     = hz_g(ms)

        m0 =  mcol_g(ns)-2
        if (kb > 1) mcol(ns) = nc + 1
        do k=2,kb
          nc = nc + 1
          mc = m0+k
          hz(nc)     = hz_g(mc)
        enddo
      enddo
    enddo
  end subroutine init_local_vars_ser

  subroutine init_local_vars_par(ia,n2d,kh,mcol,msrf,ind,hz,kh_n,mcol_n,msrf_n, &
                        ind_n,hz_n,khu,khv,h,h_old,h_new,t,s,u,v,w,hu,hv,       &
                        dispt,disps,eddyh,uvdam,idx)
    use dmi_mpi,       only : dd
    use dmi_omp,  only : domp_get_domain
    use params_n, only : ui,vi,wi
    implicit none
    !- args
    integer(4), intent(in)    :: n2d,ia
    integer(4), intent(in)    :: ind(1:,1:),msrf(0:,0:),mcol(0:),kh(0:)
    real(8), intent(in)       :: hz(0:)
    integer(4), intent(out)   :: ind_n(1:,1:),msrf_n(0:,0:),mcol_n(0:)
    integer(4), intent(out)   :: kh_n(0:),khu(0:),khv(0:)
    real(8), intent(out)      :: u(0:),v(0:),w(0:),hu(0:),hv(0:)
    real(8), intent(out)      :: dispt(0:),disps(0:),eddyh(0:)
    real(8), intent(out)      :: t(0:),s(0:)
    real(8), intent(out)      :: h(0:),h_old(0:),h_new(0:),hz_n(0:)
    real(8), intent(out)      :: uvdam(1:,0:)
    integer(4), intent(inout) :: idx(1:,1:)
    !- locals
    integer(4) :: n, n2dl, n2du, kb, ml, mu
    real(8), parameter :: FIXMEr = 0.0_8
    real(8), parameter :: zero = 0.0_8
    real(8), parameter :: land = 0.0_8 
    integer(4) :: i, j, mn, nn, mnp, nnp
    real(8), parameter :: s0 = 35.0_8, t0 = 10.0_8

! STEP2: MPI local kh, ind, msrf, mcol are now well-defined and we can
!        NUMA initialize all the MPI local vars properly

    call domp_get_domain(kh, 1, n2d, n2dl, n2du, idx)

    ! NUMA first-touch in surface (2D):
    do n=n2dl, n2du
      i=ind(1,n)
      j=ind(2,n)
      khu(n)=min(kh(msrf(i,j)),kh(msrf(i,j+1)))
      khv(n)=min(kh(msrf(i,j)),kh(msrf(i+1,j)))
      ind_n(1,n)=i
      ind_n(2,n)=j
      msrf_n(i,j)=n
    enddo
    kh_n(n2dl:n2du)   = kh(n2dl:n2du)
    mcol_n(n2dl:n2du) = mcol(n2dl:n2du)

    uvdam(1:2,n2dl:n2du)  = zero

    ! NUMA first-touch in surface (3D fields):
    s(n2dl:n2du)     = s0
    t(n2dl:n2du)     = t0
    u(n2dl:n2du)     = ui
    v(n2dl:n2du)     = vi
    w(n2dl:n2du)     = wi
    hz_n(n2dl:n2du)  = hz(n2dl:n2du)
    h(n2dl:n2du)     = hz_n(n2dl:n2du)
    hu(n2dl:n2du)    = h(n2dl:n2du)
    hv(n2dl:n2du)    = h(n2dl:n2du)
    h_old(n2dl:n2du) = h(n2dl:n2du)
    h_new(n2dl:n2du) = h(n2dl:n2du)
    ! ONLY for diffusion
    dispt(n2dl:n2du) = FIXMEr
    disps(n2dl:n2du) = FIXMEr
    eddyh(n2dl:n2du) = FIXMEr

    ! NUMA first-touch in sub-surface:
    do n=n2dl,n2du
      kb = kh(n)
      if (kb > 1) then
        ml = mcol(n)
        mu = ml + kb - 2
        s(ml:mu) = s0
        t(ml:mu) = t0
        u(ml:mu) = ui
        v(ml:mu) = vi
        w(ml:mu) = wi
        hz_n(ml:mu)  = hz(ml:mu)
        h(ml:mu)     = hz_n(ml:mu)
        hu(ml:mu)    =  h(ml:mu)
        hv(ml:mu)    =  h(ml:mu)
        h_old(ml:mu) = hz(ml:mu) 
        h_new(ml:mu) = hz(ml:mu)
        ! ONLY for diffusion
        dispt(ml:mu) =  FIXMEr
        disps(ml:mu) =  FIXMEr
        eddyh(ml:mu) =  FIXMEr
      endif
    enddo

!$OMP MASTER
    mn  = dd(ia)%dim1
    nn  = dd(ia)%dim2
    mnp = mn + 1
    nnp = nn + 1
    !- land values (2D fields)
    kh_n(0)    = 0
    khu(0)   = 0 
    khv(0)   = 0
    uvdam(1:2,0)  = land
    !- land values (3D fields)
    hz_n(0)    = land
    s(0)     = land
    t(0)     = land
    u(0)     = land
    v(0)     = land
    w(0)     = land
    hz_n(0)    = land
    hu(0)    = land
    hv(0)    = land
    dispt(0) = land 
    disps(0) = land
    eddyh(0) = land
    h(0)     = land
    h_old(0) = land
    h_new(0) = land
    msrf_n(0,0:)= 0
    msrf_n(0:,0)= 0
    mcol_n(0)= 0
    msrf_n(0:,nnp)= 0
    msrf_n(mnp,0:)= 0
!$OMP END MASTER
  end subroutine init_local_vars_par

  subroutine Validate ()

    use local_arrays_n,   only : t => t_l, s => s_l
    use params_n, only : iw3
    use dmi_mpi_global, only : iu06
    implicit none

    integer(4) :: ia
    real(8)    :: smean, smin, smax, tmean, tmin, tmax
    real(8)    :: srms, sstd, trms, tstd
    character(8),  parameter :: fmti1 = '(a33,i5)'
    character(12), parameter :: fmtr1 = '(a33,f16.13)'
    character(22), parameter :: fmtr2 = '(a33,f16.13,5x,f16.13)'
    real(8), parameter :: one = 1.0_8
    real(8), parameter :: zero = 0.0_8

    write(iu06,'(a18)') 'Validation prints:'
    do ia=1,1 
      write(iu06,'(a72)')                                                      &
      '------------------------------------------------------------------------'
      write(iu06,fmti1) 'Statistics for domain:           ', ia

      call sum_of_values(iw3(ia),s(ia)%p,smean)
      smean = smean/real(iw3(ia),8)
      call sum_of_squares(iw3(ia),s(ia)%p,zero,srms)
      srms = sqrt( srms/real(iw3(ia),8) )
      if (iw3(ia) > 1) then
        call sum_of_squares(iw3(ia),s(ia)%p,smean,sstd)
        sstd = sqrt( sstd/(real(iw3(ia)-1,8)) )
      else
        sstd = -one
      endif
      smin =  minval(s(ia)%p(1:))
      smax =  maxval(s(ia)%p(1:))

      call sum_of_values(iw3(ia),t(ia)%p,tmean)
      tmean = tmean/real(iw3(ia),8)
      call sum_of_squares(iw3(ia),t(ia)%p,zero,trms)
      trms = sqrt( trms/real(iw3(ia),8) )
      if (iw3(ia) > 1) then
        call sum_of_squares(iw3(ia),t(ia)%p,tmean,tstd)
        tstd = sqrt( tstd/(real(iw3(ia)-1,8)) )
      else
        tstd = -one
      endif
      tmin =  minval(t(ia)%p(1:))
      tmax =  maxval(t(ia)%p(1:))

      write(iu06,fmtr1) 'Average value for salinity:      ', smean
      write(iu06,fmtr2) 'RMS and STD for salinity:        ', srms, sstd
      write(iu06,fmtr1) 'Average value for temperature:   ', tmean
      write(iu06,fmtr2) 'RMS and STD for temperature:     ', trms, tstd

      write(iu06,fmtr2) 'Min and max for salinity:        ', smin, smax
      write(iu06,fmtr2) 'Min and max for temperature:     ', tmin, tmax
    enddo

    write(iu06,'(a72)')                                                        &
      '------------------------------------------------------------------------'

  end subroutine Validate

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

!==============================================================================+

  subroutine sum_of_squares_1(nsize,x,mean,sumout)
    ! double compensated sum-of-squares
    use constants, only : zero
    implicit none
    integer(4), intent(in)  :: nsize
    real(8),    intent(in)  :: x(0:), mean
    real(8),    intent(out) :: sumout
    real(8)    :: c, u, y, t, s, v, z
    integer(4) :: i
    s = (x(1) - mean)**2
    c = zero
    do i=2,nsize
      call fast2sum(c,(x(i)-mean)**2,y,u)
      call fast2sum(s,y,t,v)
      z = u + v 
      call fast2sum(t,z,s,c)
    enddo
    sumout = s
  end subroutine sum_of_squares_1


end module init_n

