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

module dmi_mpi

  !- modules: ------------------------------------------------------------------
  use dmi_mpi_global, only : iu06, debug_print, mpi_io_rank, mpi_decomp_file,  &
                             mpi_comm_model
  use exits,          only : exitme
  use cmod_mem,       only : cmr1, cmi1, cmi2
#if defined (MPI)
  use mpi
#endif

  !- implicit directives: ------------------------------------------------------
  implicit none
  private

  !- some simple vars: ---------------------------------------------------------
  integer(4), parameter, private :: halo_width=1 ! present implementation can 
                                                 ! only handle a halo width of 1
  integer(4), parameter, private :: max_nb = 8   ! max No. of neighbours
  integer(4), save,      public  :: mpi_rank
  logical,    save,      public  :: mpi_io_serial
  logical,    save,      public  :: ice_serial
  integer(4), save,      public  :: mpi_size
  integer(4), save,      private :: mpi_nc       ! max No. of components
  integer(4), save,      private :: mpi_ni       ! max No. of ice cmps 
  integer(4), save,      private :: mpi_nz       ! max No. of open bndry points
  integer(4), save,      private :: mpi_ma       ! store main area number
  integer(4), save,      private :: mpi_nt       ! No. of threads
  integer(4), save, allocatable, private :: nbnd(:)  ! # boundary points on task
  integer(4), parameter, public  :: ice_rank=0 ! FIXME hardcoding, hard to avoid
  integer(4), save,      private :: ierr
  integer(4), save,      private :: NextFreeItag, thr_b_size
  logical,    save,      private :: do_comm
  real(8),    parameter, private :: dummy_value = -999.99_8
  integer(4), save, allocatable, private :: lstatus(:,:), itagreq(:,:,:)
  type(cmi1), save, allocatable, private :: nestitags(:), rsmask(:,:)
  type(cmr1), save, allocatable, private :: rbufnb(:), sbufnb(:)
  integer(4), save, allocatable, private :: ireq(:), rbs(:), rbs1(:)
  integer(4), save, allocatable, private :: irqs(:)
  integer(4), save, allocatable, private :: lstat2(:,:), itagreq2(:,:,:,:)
  integer(4), save, allocatable, private :: itagreq7(:,:,:)
  type(cmi2), save, allocatable, private :: nestitags2(:), rsmask2(:,:)
  type(cmi1), save, allocatable, private :: nestitags7(:), rsmask7(:,:)
  type(cmr1), save, allocatable, private :: rbufnb2(:), sbufnb2(:)
  type(cmr1), save, allocatable, private :: rbufnb3(:), sbufnb3(:)
  type(cmr1), save, allocatable, private :: rbufnb7(:), sbufnb7(:)
  integer(4), save, allocatable, private :: ireq2(:), irqs2(:)
  integer(4), save, allocatable, private :: ireq3(:), irqs3(:), lstat3(:,:)
  integer(4), save, allocatable, private :: ireq7(:), irqs7(:), lstat7(:,:)
  integer(4), save, allocatable, private :: bsiz2(:,:,:,:), bsiz3(:,:,:,:)
  integer(4),       allocatable, private :: offset(:,:), offsethalo(:,:,:)
  logical,    save, allocatable, private :: lbuf(:)   ! logical buffer
  real(8),    save, allocatable, private :: gbuf(:)   ! large gather buffer
  real(8),    save, allocatable, private :: sbuf(:), rbuf(:)  ! halo S/R buffers
#if defined (MPI)
  integer(4),                    private :: istatus(MPI_STATUS_SIZE)
  type(cmr1), save, allocatable, private :: rbuf2(:)  !non-blocking halo buffers
  type(cmr1), save, allocatable, private :: tfrbuf(:), tfsbuf(:)  ! tracer buffs
  integer(4), save,              private :: ireqrh(3*max_nb), ireqsh(3*max_nb)
  integer(4), save,              private :: lstatush(MPI_STATUS_SIZE,3*max_nb)
#else
  integer(4), parameter, private :: MPI_STATUS_SIZE = 1
#endif

  !- type defs and array vars: -------------------------------------------------

  !- some task decompo info is stored in a distribution data set for each task:
  type distribution_data
    integer(4) :: dim1, dim2, dim3, iwet2, iwet3, ndims
    integer(4) :: nwet3                        ! no. of wet points
    integer(4) :: low_w3                       ! start idx of subsurface wetpnts
    integer(4) :: low_i, low_j, low_ws         ! active data areas
    integer(4) ::  up_i,  up_j,  up_ws         ! active data areas
    integer(4) :: low_hi, low_hj               ! halo region lower bound
    integer(4) ::  up_hi,  up_hj               ! halo region upper bound
    integer(4) :: nt_w, nt_n, nt_e, nt_s       ! task id for neighbours 
    integer(4) :: nt_nw, nt_ne, nt_se, nt_sw   ! task id for corners 
    integer(4) :: recv_bufsize(max_nb)         ! # active halopoints in 8 dir.
    integer(4) :: recv_buf_srf(max_nb)         ! as above, but only surface
    integer(4) :: send_bufsize(max_nb)         ! # active halopoints in 8 dir.
    integer(4) :: send_buf_srf(max_nb)         ! as above, but only surface
    integer(4) :: halo2, halo3                 ! total halo/buffer sizes
  end type distribution_data
  public :: distribution_data
  type (distribution_data), allocatable, save, public :: dd(:)

  !- save some decompo info for all tasks in a task table:
  type task_table
    integer(4) :: nwet3                        ! no. of wet points
    integer(4) :: low_w3                       ! start idx of subsurface wetpnts
    integer(4) :: low_i,  low_j, low_ws        ! active data areas, lower bounds
    integer(4) ::  up_i,   up_j,  up_ws        ! active data areas, upper bounds
    integer(4) :: low_hi, low_hj               ! halo region lower bound
    integer(4) ::  up_hi,  up_hj               ! halo region upper bound
    integer(4) :: hnb(max_nb)                  ! halo neighbours
  end type task_table
  public :: task_table
  type (task_table), allocatable, save, public :: mpi_tt(:,:)

  !- neighbour-to-neighbour table:
  integer(4), parameter, private :: n2n(max_nb) = (/ 3,4,1,2,7,8,5,6 /)

  !- scatter table:
  logical, allocatable, public, save :: tscat(:,:,:)
  logical, allocatable, public, save :: nctable(:,:,:), nftable(:,:,:)
  logical, allocatable, public, save :: uctable(:,:,:), vctable(:,:,:)
  logical, allocatable, public, save :: uftable(:,:,:), vftable(:,:,:)
  logical, allocatable, public, save :: cctable(:,:), cftable(:,:)
  logical, allocatable, public, save :: uhtable(:,:), vhtable(:,:)
  logical, allocatable, public, save :: mctable(:,:)

#if defined (MPI)
  !- masks for halo communication:
  type comm_mask
    integer(4), allocatable :: mask(:)
  end type comm_mask
  private comm_mask
  type (comm_mask), allocatable, private, save :: encode(:,:), decode(:,:)
#endif

  !- Interfaces ----------------------------------------------------------------
  interface dmpi_gather
    module procedure dmpi_gather_log
    module procedure dmpi_gather_srf
    module procedure dmpi_gather_srf_nc
    module procedure dmpi_gather_iota
    module procedure dmpi_gather_all
    module procedure dmpi_gather_all_nc
    module procedure dmpi_gather_global_srf
    module procedure dmpi_gather_global_all
    module procedure dmpi_gather_bnd_data_out
    module procedure dmpi_gather_bnd_data
    module procedure dmpi_gather_bnd_data_2
    module procedure dmpi_gather_bnd_3
    module procedure dmpi_gather_bnd_2
    module procedure dmpi_gather_uv_bdr_nb
    module procedure dmpi_gather_uv_bdr_cmp_nb
    module procedure dmpi_gather_copy_cmp
    module procedure dmpi_gather_copy
    module procedure dmpi_gather_mcf_nb
  end interface

  interface dmpi_decompose
    module procedure decompose_read ! starting point 
  end interface

  interface dmpi_halo
    module procedure dmpi_distribute_halo_nb
    module procedure dmpi_distribute_halo_nb_col
    module procedure dmpi_distribute_halo_TF_nb
    module procedure dmpi_distribute_halo_TF_nb2
    module procedure dmpi_distribute_halo_log
  end interface

  interface dmpi_broadcast
    module procedure dmpi_broadcast_met_info
    module procedure dmpi_broadcast_met_data
    module procedure dmpi_broadcast_bnd_data
    module procedure dmpi_broadcast_bnd_data_2
    module procedure dmpi_broadcast_bnd_1d
    module procedure dmpi_broadcast_logical
  end interface

  interface dmpi_scatter
    module procedure dmpi_scatter_ice
    module procedure dmpi_scatter_ice2
    module procedure dmpi_scatter_cmp_table
    module procedure dmpi_scatter_cmp
    module procedure dmpi_scatter_arr_table
    module procedure dmpi_scatter_arr
    module procedure dmpi_scatter_bnd_data
    module procedure dmpi_scatter_bnd_data_2
  end interface

  interface dmpi_reinit
    module procedure dmpi_reinit1
    module procedure dmpi_reinit2
  end interface

  interface dmpi_init
    module procedure dmpi_init1
    module procedure dmpi_init2
    module procedure dmpi_init3
    module procedure dmpi_init4
    module procedure dmpi_init5
  end interface

  interface dmpi_validate
    module procedure dmpi_validate_min_max
    module procedure dmpi_validate_pos
    module procedure dmpi_validate_r8
    module procedure dmpi_validate_log
  end interface

  interface dmpi_send
    module procedure dmpi_send_real8
    module procedure dmpi_send_integer4
    module procedure dmpi_send_logical
  end interface

  interface dmpi_recv
    module procedure dmpi_recv_real8
    module procedure dmpi_recv_integer4
    module procedure dmpi_recv_logical
  end interface

  interface dmpi_irecv
    module procedure dmpi_irecv_real8
  end interface

  interface dmpi_isend
    module procedure dmpi_isend_real8
    module procedure dmpi_isend_integer4
  end interface

  interface dmpi_bcast
    module procedure dmpi_bcast_real8
    module procedure dmpi_bcast_integer4
    module procedure dmpi_bcast_logical
    module procedure dmpi_bcast_char
  end interface

  interface dmpi_reduce
    module procedure dmpi_allreduce_real8
    module procedure dmpi_allreduce_integer4
    module procedure dmpi_allreduce_real8_scalar
    module procedure dmpi_allreduce_integer4_scalar
  end interface

#if defined (MPI)
  interface dmpi_decode_buf_nc
    module procedure dmpi_decode_buf_nc1
    module procedure dmpi_decode_buf_nc2
    module procedure dmpi_decode_buf_nc3
  end interface

  interface dmpi_encode_buf_nc
    module procedure dmpi_encode_buf_nc1
    module procedure dmpi_encode_buf_nc2
    module procedure dmpi_encode_buf_nc3
  end interface

  interface dmpi_decode_buf
    module procedure dmpi_decode_buf_col
    module procedure dmpi_decode_buf_mmk
  end interface

  interface dmpi_encode_buf
    module procedure dmpi_encode_buf_col
    module procedure dmpi_encode_buf_mmk
  end interface

#endif

  !- private methods -----------------------------------------------------------
#if defined (MPI)
  private :: dmpi_bufsize, dmpi_roff, dmpi_soff, dmpi_halosize
  private :: dmpi_encode_buf, dmpi_decode_buf
  private :: dmpi_encode_buf_mmk, dmpi_decode_buf_mmk
  private :: dmpi_encode_buf_col, dmpi_decode_buf_col
  private :: dmpi_encode_buf_nc, dmpi_decode_buf_nc, dmpi_comm_error
  private :: dmpi_decode_buf_nc1, dmpi_decode_buf_nc2, dmpi_decode_buf_nc3
  private :: dmpi_encode_buf_nc1, dmpi_encode_buf_nc2, dmpi_encode_buf_nc3
  private :: dmpi_encode_buf_log, dmpi_decode_buf_log
  private :: dmpi_init_comm_masks, dmpi_set_comm_masks
  private :: decompose_plane_set_bounds
#endif
  private :: dmpi_send_real8, dmpi_send_integer4, dmpi_send_logical
  private :: dmpi_recv_real8, dmpi_recv_integer4, dmpi_recv_logical
  private :: dmpi_irecv_real8, dmpi_isend_real8
  private :: dmpi_bcast_real8, dmpi_bcast_integer4, dmpi_bcast_logical
  private :: dmpi_bcast_char, dmpi_mcf, decompose_read
  private :: dmpi_broadcast_met_info, dmpi_broadcast_met_data
  private :: dmpi_broadcast_bnd_data, dmpi_broadcast_bnd_1d
  private :: dmpi_broadcast_logical, dmpi_broadcast_bnd_data_2
  private :: dmpi_gather_srf, dmpi_gather_srf_nc, dmpi_gather_log
  private :: dmpi_gather_all, dmpi_gather_all_nc, dmpi_gather_bnd_data
  private :: dmpi_gather_bnd_data_2, dmpi_gather_bnd_3, dmpi_gather_bnd_data_out
  private :: dmpi_gather_bnd_2
  private :: dmpi_gather_global_srf, dmpi_gather_global_all, dmpi_gather_iota
  private :: dmpi_distribute_halo_log
  private :: dmpi_scatter_bnd_data, dmpi_scatter_bnd_data_2
  private :: dmpi_gather_uv_bdr_nb, dmpi_gather_uv_bdr_cmp_nb
  private :: dmpi_scatter_ice, dmpi_scatter_cmp_table, dmpi_scatter_cmp
  private :: dmpi_scatter_arr_table, dmpi_scatter_arr, dmpi_scatter_ice2
  private :: dmpi_init1, dmpi_init2, dmpi_validate_min_max, dmpi_validate_pos
  private :: dmpi_init3, dmpi_init4, dmpi_init5
  private :: dmpi_validate_r8, dmpi_validate_log, dmpi_distribute_halo_TF_nb
  private :: dmpi_distribute_halo_TF_nb2, dmpi_distribute_halo_nb
  private :: dmpi_distribute_halo_nb_col
  private :: dmpi_gather_copy, dmpi_gather_copy_cmp, dmpi_gather_mcf_nb
  private :: dmpi_allreduce_real8, dmpi_allreduce_integer4
  private :: dmpi_allreduce_real8_scalar, dmpi_allreduce_integer4_scalar
  !- public methods ------------------------------------------------------------
  public  :: dmpi_init, dmpi_decompose, dmpi_debug, dmpi_halo
  public  :: dmpi_barrier, dmpi_validate, dmpi_reinit
  public  :: dmpi_broadcast, dmpi_gather, dmpi_scatter, dmpi_scatter_table
  public  :: dmpi_reset_size, dmpi_recv, dmpi_send, dmpi_irecv, dmpi_isend
  public  :: bathy2mpibin, bathyread_mm, bathyread_hz
  public  :: dmpi_bcast, dmpi_reduce, decompose_gen

contains

!===============================================================================

  subroutine dmpi_init1(ldebug)

    use io_subs, only : io_new_unit

    implicit none

    logical, intent(in), optional :: ldebug

    character(12) :: logfilename
    integer(4)    :: ios
#if defined (_OPENMP) && defined (MPI)
    integer(4)    :: provided
#endif
#if defined (MPI)
    integer(4)    :: ierr2
#endif
    logical, parameter :: stdout = .true.

    if (present(ldebug)) then
      debug_print = ldebug
    else
      debug_print = .false.
    endif


#if defined (_OPENMP) && defined (MPI)
    call mpi_init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,'(a)')    ' MPI_INIT failed.'
      write (*,'(a,i4)') ' MPI_INIT Error =  ', ierr
      call exitme(1,'Bailing out due to MPI errors',stdout)
    endif
    if (provided /= MPI_THREAD_MULTIPLE) then
      write (*,'(a)')    ' MPI_INIT failed.'
      write (*,'(a,i4)') ' Required =      ', MPI_THREAD_MULTIPLE
      write (*,'(a,i4)') ' Provided =      ', provided
      call exitme(1,'Bailing out due to MPI errors',stdout)
    endif
#elif defined (MPI)
    call mpi_init(ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,'(a)')    ' MPI_INIT failed.'
      write (*,'(a,i4)') ' MPI_INIT Error =  ', ierr
      call exitme(1,'Bailing out due to MPI errors',stdout)
    endif
#endif


#if defined (MPI)
    call mpi_comm_size(MPI_COMM_WORLD, mpi_size, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,'(a)')    ' MPI_COMM_SIZE failed.'
      write (*,'(a,i4)') ' MPI_COMM_SIZE Error =  ', ierr
      call mpi_abort(MPI_COMM_WORLD,ierr,ierr2)
    endif
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      write (*,'(a)')    ' mpi_comm_rank failed.'
      write (*,'(a,i4)') ' mpi_comm_rank Error =  ', ierr
      call mpi_abort(MPI_COMM_WORLD,ierr,ierr2)
    endif
#else
    mpi_size = 1 
    mpi_rank = 0
#endif

    iu06 = io_new_unit()
    write (logfilename,'(a,i4.4)') 'logfile.',mpi_rank
    open (iu06, file=logfilename,status='replace', iostat=ios)
    if (ios /= 0) then
      call exitme(1,'Unable to open '//trim(logfilename),stdout)
    endif

#if defined (MPI)
    if (mpi_size > 1) then
      write(iu06,'(a23,i5,a10)') 'Running with MPI using ',mpi_size,' MPI tasks'
    else
      write(iu06,'(a23)') 'Running with 1 MPI task'
    endif
#else
     write(iu06,'(a19)') 'Running without MPI'
#endif

    if (ice_rank /= mpi_io_rank) then
      call exitme(1,'Bailing out due to MPI errors: ice_rank /= mpi_io_rank')
    endif

    mpi_io_serial = ( mpi_rank == mpi_io_rank)
    ice_serial    = ( mpi_rank == ice_rank)

    NextFreeItag = 1
    do_comm      = .false.

  end subroutine dmpi_init1

!===============================================================================

  subroutine dmpi_reset_size(nsize)
    implicit none
    integer(4), intent(in) :: nsize
    mpi_size=nsize
  end subroutine dmpi_reset_size

!===============================================================================

  subroutine dmpi_reinit1(msize,mrank,mcm)
    implicit none
    integer(4), intent(in) :: msize,mrank,mcm
    mpi_rank=mrank
    mpi_size=msize
    mpi_comm_model=mcm
    mpi_io_serial = ( mpi_rank == mpi_io_rank)
    ice_serial    = ( mpi_rank == ice_rank)
  end subroutine dmpi_reinit1

  subroutine dmpi_reinit2(mcm)
    implicit none
    integer(4), intent(out) :: mcm
#if defined (MPI)
    call mpi_comm_dup(MPI_COMM_WORLD,mpi_comm_model,ierr) 
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    mcm=mpi_comm_model
#else
    mcm=0
#endif
  end subroutine dmpi_reinit2

!===============================================================================

  subroutine dmpi_init2(nc, nz, ma, ni, nt)

    implicit none

    integer(4), intent(in) :: nc, nz, ma, ni, nt

    mpi_nc = nc
    mpi_nz = nz
    mpi_ma = ma
    mpi_ni = ni
    mpi_nt = nt

    allocate( offsethalo(1:mpi_nt,1:2,1:2) )

  end subroutine dmpi_init2

!===============================================================================

  subroutine dmpi_init3( narea, nestingfrom, znesting, krz, kmx, iga, jga,     &
                         mm1k, kh, mm1k_l, kh_l )
                         
    ! prepare for non-blocking MPI-comm around nesting 1+4

    use cmod_mem, only : cmi1, cmi2, cmi3
    implicit none

    integer(4), intent(in) :: narea, kmx(:)
    type(cmi1), intent(in) :: nestingfrom(:), kh(:), kh_l(:)
    type(cmi2), intent(in) :: znesting(:), krz(:), iga(:), jga(:)
    type(cmi3), intent(in) :: mm1k(:), mm1k_l(:)

    integer(4) :: iam, MaxRreq, MaxSreq, Maxitag, ia, iia, iao, bs, inreq, it
    integer(4) :: ig1, ig2, jg1, jg2, nz1, nz2, ilf, iuf, jlf, juf, iff, jff
    integer(4) :: il, iu, jl, ju, iz, maxnestlev
    integer(4) :: nreq, irr, itt, irl, irreq, isreq, bs1
    integer(4), allocatable :: itmp(:,:), itmp1(:)

    ! some initial stuff:
    maxnestlev = 0
    allocate( nestitags(narea) )
    do ia=1,narea
      iao = nestingfrom(ia)%p(0)
      if (iao > 0) then  
        maxnestlev = max( maxnestlev, iao )
        allocate( nestitags(ia)%p(1:iao) )
        nestitags(ia)%p(:) = 0
      endif
    enddo

    !  initialized here, used in the followig
    NextFreeItag = 1

    ! in trivial cases, get out of here:
    if (narea == 1 .or. maxnestlev == 0) return

#if defined (MPI)
    do_comm = .true.  ! also for init4

    ! Initially, allow that all areas need to do the maximum number
    ! of MPI-communications around nesting:
    allocate( itmp(narea * maxnestlev * mpi_size, 2), rsmask(narea,mpi_size),  &
              itmp1(narea * maxnestlev * mpi_size) )
    itmp(:,:) = 0 ! temp size of r/s buffers
    itmp1(:)  = 0
    do ia=1,narea
      iao = nestingfrom(ia)%p(0)
      if (iao > 0) then
        do it=1,mpi_size
          allocate( rsmask(ia,it)%p(1:iao) )
          rsmask(ia,it)%p(:) = 0
        enddo
      endif
    enddo

    iam = mpi_rank + 1
    MaxRreq = 0 
    MaxSreq = 0 
    Maxitag = 0
    irreq   = 0
    isreq   = 0
    do ia=1,narea
      if (ia == mpi_ma) cycle         ! no nesting to main area
      do iao=1,nestingfrom(ia)%p(0)
        iia = nestingfrom(ia)%p(iao)  ! coarse grid

        ! define tags:
        Maxitag = Maxitag + 1
        nestitags(ia)%p(iao) = Maxitag

        ig1 = iga(ia)%p(1,iao)
        ig2 = iga(ia)%p(2,iao)
        jg1 = jga(ia)%p(1,iao)
        jg2 = jga(ia)%p(2,iao)
        nz1 = znesting(ia)%p(iao,1)
        nz2 = znesting(ia)%p(iao,2)

        ! define recv:
        if (nftable(iam,iia,ia)) then
          !  fine grid dims on task #iam:
          ilf = dd(ia)%low_i
          iuf = dd(ia)%up_i
          jlf = dd(ia)%low_j
          juf = dd(ia)%up_j

          do it=1,mpi_size
            irreq = irreq + 1

            !  coarse grid dims on task #it
            il = mpi_tt(iia,it)%low_i
            iu = mpi_tt(iia,it)%up_i
            jl = mpi_tt(iia,it)%low_j
            ju = mpi_tt(iia,it)%up_j

            ! receive but not from myself:
            if (it /= iam) then
              bs  = 0
              bs1 = 0
              do iz=nz1,nz2
                iff = krz(ia)%p(1,iz)
                jff = krz(ia)%p(2,iz)
                if (.not.(ilf <= iff .and. iff <= iuf .and.                    &
                          jlf <= jff .and. jff <= juf)      ) cycle
                if (mm1k(ia)%p(1,iff,jff) <= 0) cycle
                call dmpi_mcf(iz, nz2, krz(ia)%p, ig1, jg1, ig2, jg2,          &
                              il, iu, jl, ju,                                  &
                              kh(iia)%p, kmx(iia), (1), bs, mm1k(iia)%p )
                call dmpi_mcf(iz, nz2, krz(ia)%p, ig1, jg1, ig2, jg2,          &
                              il, iu, jl, ju,                                  &
                              kh(iia)%p, 1, (1), bs1, mm1k(iia)%p )
              enddo
              if (bs > 0) then
                MaxRreq = MaxRreq + 1
                itmp(irreq,1) = bs
                itmp1(irreq)  = bs1
                rsmask(ia,it)%p(iao) = rsmask(ia,it)%p(iao) + 1
              endif
            endif
          enddo
        endif

        ! define send:
        if (mctable(iia,ia)) then
          !  coarse grid dims on task #iam
          il = dd(iia)%low_i
          iu = dd(iia)%up_i
          jl = dd(iia)%low_j
          ju = dd(iia)%up_j
          do it=1,mpi_size
            isreq = isreq + 1

            !  fine grid dims on task #it:
            ilf = mpi_tt(ia,it)%low_i
            iuf = mpi_tt(ia,it)%up_i
            jlf = mpi_tt(ia,it)%low_j
            juf = mpi_tt(ia,it)%up_j

            ! send but not to myself:
            if (it /= iam) then
              bs = 0
              do iz=nz1,nz2
                iff = krz(ia)%p(1,iz)
                jff = krz(ia)%p(2,iz)
                if (.not.(ilf <= iff .and. iff <= iuf .and.                    &
                          jlf <= jff .and. jff <= juf)      ) cycle
                if (mm1k(ia)%p(1,iff,jff) <= 0) cycle
                call dmpi_mcf(iz, nz2, krz(ia)%p, ig1, jg1, ig2, jg2,          &
                              il, iu, jl, ju,                                  &
                              kh_l(iia)%p, kmx(iia), (1), bs, mm1k_l(iia)%p)
              enddo
              if (bs > 0) then
                MaxSreq = MaxSreq + 1
                itmp(isreq,2) = bs
                rsmask(ia,it)%p(iao) = rsmask(ia,it)%p(iao) + 10
              endif
            endif
          enddo
        endif
      enddo
    enddo

    NextFreeItag = Maxitag + 1   ! NextFreeItag value passed on to init4
    allocate( itagreq(Maxitag,1:2,1:2) )
    itagreq(:,:,:) = 0
    nreq = MaxRreq + MaxSreq
    if (nreq > 0) then
      allocate( lstatus(MPI_STATUS_SIZE,nreq) )
      if (MaxRreq > 0) then
        allocate( ireq(MaxRreq), rbs(MaxRreq), rbufnb(MaxRreq), rbs1(MaxRreq) )
      endif
      if (MaxSreq > 0) allocate( irqs(MaxSreq), sbufnb(MaxSreq) )

      ! now, alloc the recv buffers and define the request ranges for each itag
      ! as well as the sizes of the recv buffers:
      irr   = 0
      irreq = 0
      do ia=1,narea
        if (ia == mpi_ma) cycle         ! no nesting to main area
        do iao=1,nestingfrom(ia)%p(0)
          iia = nestingfrom(ia)%p(iao)  ! coarse grid
          itt = nestitags(ia)%p(iao)
          irl = 0
          if (nftable(iam,iia,ia)) then
            do it=1,mpi_size
              irreq = irreq + 1
              if (it /= iam) then
                bs = itmp(irreq,1)
                if (bs > 0) then
                  irr = irr + 1
                  rbs(irr)  = bs
                  rbs1(irr) = itmp1(irreq)
                  if (irl == 0) irl = irr
                  allocate( rbufnb(irr)%p(1:bs) )
                endif
              endif
            enddo
            if (irl > 0) then
              itagreq(itt,1,1) = irl
              itagreq(itt,1,2) = irr
            endif 
          endif
        enddo
      enddo

      ! then, alloc the send buffers and define the request ranges for each itag
      ! as well as the sizes of the send buffers:
      irr   = 0
      isreq = 0
      do ia=1,narea
        if (ia == mpi_ma) cycle         ! no nesting to main area
        do iao=1,nestingfrom(ia)%p(0)
          iia = nestingfrom(ia)%p(iao)  ! coarse grid
          itt = nestitags(ia)%p(iao)
          irl = 0
          if (mctable(iia,ia)) then
            do it=1,mpi_size
              isreq = isreq + 1
              if (it /= iam) then
                bs = itmp(isreq,2)
                if (bs > 0) then
                  irr = irr + 1
                  if (irl == 0) irl = irr
                  allocate( sbufnb(irr)%p(1:bs) )
                endif
              endif
            enddo
            if (irl > 0) then
              itagreq(itt,2,1) = irl
              itagreq(itt,2,2) = irr
            endif
          endif
        enddo
      enddo

    endif

    deallocate( itmp, itmp1 )
#else
    allocate( itagreq(1,1,1), rbs(1), rsmask(1,1), rbs1(1) )
    allocate( rsmask(1,1)%p(1) )
#endif

  end subroutine dmpi_init3

!===============================================================================

  subroutine dmpi_init4( narea, nestinglevels, nestingto, enclosing,           &
                         unesting, vnesting, kru, krv, iga, jga,               &
                         mm1k, kh, mm1k_l, kh_l )
                         
    ! prepare for non-blocking MPI-comm around nesting 2+8+b

    use cmod_mem, only : cmi1, cmi2, cmi3
    implicit none

    integer(4), intent(in) :: narea, nestinglevels(:), enclosing(:,:)
    type(cmi1), intent(in) :: nestingto(:), kh(:), kh_l(:)
    type(cmi2), intent(in) :: unesting(:), vnesting(:)
    type(cmi2), intent(in) :: kru(:), krv(:), iga(:), jga(:)
    type(cmi3), intent(in) :: mm1k(:), mm1k_l(:)

    integer(4) :: iam, MaxRreq, MaxSreq, Maxitag, ia, iia, iao, bs, inreq, it
    integer(4) :: ig1, ig2, jg1, jg2, ilf, iuf, jlf, juf, iff, jff, kb, bsc
    integer(4) :: il, iu, jl, ju, iz, maxnestlev, ii, nu1, nu2, nv1, nv2, ijpm
    integer(4) :: irr, itt, irl, irreq, isreq, ig, jg, midx, ill
    integer(4) :: iii, jjj, i, j
    integer(4), allocatable :: itmp(:,:,:), icmp(:,:,:)

    ! some initial stuff:
    maxnestlev = 0
    allocate( nestitags2(narea) )
    do ia=1,narea
      ii = nestinglevels(ia)
      if (ii > 0) then  
        maxnestlev = max( maxnestlev, ii )
        allocate( nestitags2(ia)%p(1:ii,1:4) )
        nestitags2(ia)%p(:,:) = 0
      endif
    enddo

    ! in trivial cases, get out of here:
    if (narea == 1 .or. maxnestlev == 0) return

#if defined (MPI)
    ! Initially, allow that all areas need to do the maximum number
    ! of MPI-communications around nesting:
    allocate( itmp(narea*maxnestlev*mpi_size, 2, 2), rsmask2(narea,mpi_size) )
    allocate( icmp(narea*maxnestlev*mpi_size, 2, 2) )
    itmp(:,:,:) = 0  ! temp size of r/s buffers
    icmp(:,:,:) = 0  ! temp size of r/s buffers
    do ia=1,narea
      ii = nestinglevels(ia)
      if (ii > 0) then
        do it=1,mpi_size
          allocate( rsmask2(ia,it)%p(1:ii,1:2) )
          rsmask2(ia,it)%p(:,:) = 0
        enddo
      endif
    enddo
    allocate( bsiz2(1:mpi_nt,1:mpi_size,1:2,1:2), offset(1:mpi_nt,1:2),        &
              bsiz3(1:mpi_nt,1:mpi_size,1:2,1:2) )

    iam = mpi_rank + 1
    MaxRreq = 0 
    MaxSreq = 0 
    Maxitag = NextFreeItag - 1 ! NextFreeItag first set in init3
    irreq   = 0
    isreq   = 0
    ijpm    = 1
    do ia=1,narea
      do ii=1,nestinglevels(ia)
        !  iia: fine grid
        !  ia:  coarse grid
        iia = nestingto(ia)%p(ii)
        iao = enclosing(iia,ia)
        nu1 = unesting(ia)%p(ii,1)
        nu2 = unesting(ia)%p(ii,2)
        nv1 = vnesting(ia)%p(ii,1)
        nv2 = vnesting(ia)%p(ii,2)

        ! define tags:
        nestitags2(ia)%p(ii,1) = Maxitag + 1   ! u-brd
        nestitags2(ia)%p(ii,2) = Maxitag + 2   ! v-brd
        nestitags2(ia)%p(ii,3) = Maxitag + 3   ! u-brd, cmp
        nestitags2(ia)%p(ii,4) = Maxitag + 4   ! v-brd, cmp
        Maxitag = Maxitag + 4

        ig1 = iga(iia)%p(1,iao)
        ig2 = iga(iia)%p(2,iao)
        jg1 = jga(iia)%p(1,iao)
        jg2 = jga(iia)%p(2,iao)

        ! define recv:
        if (uctable(iam,ia,iia)) then
          !  u-brd:
          do it=1,mpi_size
            irreq = irreq + 1
            if (it == iam) cycle
            bs  = 0
            bsc = 0
            do iz=nu1,nu2
              ig = kru(ia)%p(1,iz)
              jg = kru(ia)%p(2,iz)
              if (kru(ia)%p(3,iz) == 3) jg = jg+1
              if (.not. ((dd(ia)%low_i <= ig .and. ig <= dd(ia)%up_i) .and.    &
                         (dd(ia)%low_j <= jg .and. jg <= dd(ia)%up_j))  ) cycle
              if (mm1k(ia)%p(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=-ijpm,jg2-1+ijpm
                jff = j + jjj
                if (jff < mpi_tt(iia,it)%low_j .or.                            &
                    jff > mpi_tt(iia,it)%up_j       ) cycle
                do iii=-ijpm,ig2-1+ijpm
                  iff = i + iii
                  if (iff < mpi_tt(iia,it)%low_i .or.                          &
                      iff > mpi_tt(iia,it)%up_i       ) cycle
                  midx = mm1k(iia)%p(1,iff,jff)
                  if (midx > 0) bs = bs + kh(iia)%p(midx)
                enddo
              enddo
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < mpi_tt(iia,it)%low_j .or.                            &
                    jff > mpi_tt(iia,it)%up_j       ) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < mpi_tt(iia,it)%low_i .or.                          &
                      iff > mpi_tt(iia,it)%up_i       ) cycle
                  midx = mm1k(iia)%p(1,iff,jff)
                  if (midx > 0) bsc = bsc + kh(iia)%p(midx)
                enddo
              enddo
            enddo
            if (bs > 0) then
              MaxRreq = MaxRreq + 1
              itmp(irreq,1,1) = 4*bs
              icmp(irreq,1,1) = max(mpi_nc*bsc,1)
              rsmask2(ia,it)%p(ii,1) = rsmask2(ia,it)%p(ii,1) + 1
            endif
          enddo
        endif
        if (vctable(iam,ia,iia)) then
          !  v-brd:
          do it=1,mpi_size
            irreq = irreq + 1
            if (it == iam) cycle
            bs  = 0
            bsc = 0
            do iz=nv1,nv2
              ig = krv(ia)%p(1,iz)
              jg = krv(ia)%p(2,iz)
              if (krv(ia)%p(3,iz) == 4) ig = ig+1
              if (.not. ((dd(ia)%low_i <= ig .and. ig <= dd(ia)%up_i) .and.    &
                         (dd(ia)%low_j <= jg .and. jg <= dd(ia)%up_j))  ) cycle
              if (mm1k(ia)%p(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=-ijpm,jg2-1+ijpm
                jff = j + jjj
                if (jff < mpi_tt(iia,it)%low_j .or.                            &
                    jff > mpi_tt(iia,it)%up_j       ) cycle
                do iii=-ijpm,ig2-1+ijpm
                  iff = i + iii
                  if (iff < mpi_tt(iia,it)%low_i .or.                          &
                      iff > mpi_tt(iia,it)%up_i       ) cycle
                  midx = mm1k(iia)%p(1,iff,jff)
                  if (midx > 0) bs = bs + kh(iia)%p(midx)
                enddo
              enddo
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < mpi_tt(iia,it)%low_j .or.                            &
                    jff > mpi_tt(iia,it)%up_j       ) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < mpi_tt(iia,it)%low_i .or.                          &
                      iff > mpi_tt(iia,it)%up_i       ) cycle
                  midx = mm1k(iia)%p(1,iff,jff)
                  if (midx > 0) bsc = bsc + kh(iia)%p(midx)
                enddo
              enddo
            enddo
            if (bs > 0) then
              MaxRreq = MaxRreq + 1
              itmp(irreq,1,2) = 4*bs
              icmp(irreq,1,2) = max(mpi_nc*bsc,1)
              rsmask2(ia,it)%p(ii,2) = rsmask2(ia,it)%p(ii,2) + 1
            endif
          enddo
        endif

        ! define send:
        if (uhtable(ia,iia)) then
          !  u-brd:
          do it=1,mpi_size
            isreq = isreq + 1
            if (it == iam) cycle
            !  coarse grid dims on task #it
            il = mpi_tt(ia,it)%low_i
            iu = mpi_tt(ia,it)%up_i
            jl = mpi_tt(ia,it)%low_j
            ju = mpi_tt(ia,it)%up_j
            bs  = 0
            bsc = 0
            do iz=nu1,nu2
              ig = kru(ia)%p(1,iz)
              jg = kru(ia)%p(2,iz)
              if (kru(ia)%p(3,iz) == 3) jg = jg+1
              if (ig < il .or. iu < ig) cycle
              if (jg < jl .or. ju < jg) cycle
              if (mm1k(ia)%p(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=-ijpm,jg2-1+ijpm
                jff = j + jjj
                if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                do iii=-ijpm,ig2-1+ijpm
                  iff = i + iii
                  if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                  midx = mm1k_l(iia)%p(1,iff,jff)
                  if (midx > 0) bs = bs + kh_l(iia)%p(midx)
                enddo
              enddo
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                  midx = mm1k_l(iia)%p(1,iff,jff)
                  if (midx > 0) bsc = bsc + kh_l(iia)%p(midx)
                enddo
              enddo
            enddo
            if (bs > 0) then
              MaxSreq = MaxSreq + 1
              itmp(isreq,2,1) = 4*bs
              icmp(isreq,2,1) = max(mpi_nc*bsc,1)
              rsmask2(ia,it)%p(ii,1) = rsmask2(ia,it)%p(ii,1) + 10
            endif
          enddo
        endif
        if (vhtable(ia,iia)) then
          !  v-brd:
          do it=1,mpi_size
            isreq = isreq + 1
            if (it == iam) cycle
            !  coarse grid dims on task #it
            il = mpi_tt(ia,it)%low_i
            iu = mpi_tt(ia,it)%up_i
            jl = mpi_tt(ia,it)%low_j
            ju = mpi_tt(ia,it)%up_j
            bs  = 0
            bsc = 0
            do iz=nv1,nv2
              ig = krv(ia)%p(1,iz)
              jg = krv(ia)%p(2,iz)
              if (krv(ia)%p(3,iz) == 4) ig = ig+1
              if (ig < il .or. iu < ig) cycle
              if (jg < jl .or. ju < jg) cycle
              if (mm1k(ia)%p(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=-ijpm,jg2-1+ijpm
                jff = j + jjj
                if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                do iii=-ijpm,ig2-1+ijpm
                  iff = i + iii
                  if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                  midx = mm1k_l(iia)%p(1,iff,jff)
                  if (midx > 0) bs = bs + kh_l(iia)%p(midx)
                enddo
              enddo
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                  midx = mm1k_l(iia)%p(1,iff,jff)
                  if (midx > 0) bsc = bsc + kh_l(iia)%p(midx)
                enddo
              enddo
            enddo
            if (bs > 0) then
              MaxSreq = MaxSreq + 1
              itmp(isreq,2,2) = 4*bs
              icmp(isreq,2,2) = max(mpi_nc*bsc,1)
              rsmask2(ia,it)%p(ii,2) = rsmask2(ia,it)%p(ii,2) + 10
            endif
          enddo
        endif
      enddo
    enddo

    allocate( itagreq2(NextFreeItag:Maxitag,1:2,1:2,1:2) )
    Maxitag      = Maxitag - NextFreeItag + 1
    NextFreeItag = NextFreeItag + Maxitag
    itagreq2(:,:,:,:) = 0
    if (MaxRreq + MaxSreq > 0) then
      allocate(lstat2(MPI_STATUS_SIZE,MaxRreq+MaxSreq),                        &
               lstat3(MPI_STATUS_SIZE,MaxRreq+MaxSreq))
      if (MaxRreq > 0) allocate(ireq2(MaxRreq), rbufnb2(MaxRreq),              &
                                ireq3(MaxRreq), rbufnb3(MaxRreq))
      if (MaxSreq > 0) allocate(irqs2(MaxSreq), sbufnb2(MaxSreq),              &
                                irqs3(MaxSreq), sbufnb3(MaxSreq))

      ! now, alloc the recv buffers and define the request ranges for each itag
      ! as well as the sizes of the recv buffers:
      irr   = 0
      irreq = 0
      do ia=1,narea
        do ii=1,nestinglevels(ia)
          !  iia: fine grid
          !  ia:  coarse grid
          iia = nestingto(ia)%p(ii)
          iao = enclosing(iia,ia)

          ! u-brd
          irl = 0
          if (uctable(iam,ia,iia)) then
            do it=1,mpi_size
              irreq = irreq + 1
              if (it == iam) cycle
              bs = itmp(irreq,1,1)
              if (bs > 0) then
                irr = irr + 1
                if (irl == 0) irl = irr
                bsc = icmp(irreq,1,1)
                allocate( rbufnb2(irr)%p(1:bs), rbufnb3(irr)%p(1:bsc) )
              endif
            enddo
            if (irl > 0) then
              itt = nestitags2(ia)%p(ii,1)
              itagreq2(itt,1,1,1) = irl
              itagreq2(itt,1,2,1) = irr
              itt = nestitags2(ia)%p(ii,3)
              itagreq2(itt,1,1,1) = irl
              itagreq2(itt,1,2,1) = irr
            endif 
          endif

          ! v-brd
          ill = irr
          irl = ill
          if (vctable(iam,ia,iia)) then
            do it=1,mpi_size
              irreq = irreq + 1
              if (it == iam) cycle
              bs = itmp(irreq,1,2)
              if (bs > 0) then
                irr = irr + 1
                if (irl == ill) irl = irr
                bsc = icmp(irreq,1,2)
                allocate( rbufnb2(irr)%p(1:bs), rbufnb3(irr)%p(1:bsc) )
              endif
            enddo
            if (irl > ill) then
              itt = nestitags2(ia)%p(ii,2)
              itagreq2(itt,1,1,2) = irl
              itagreq2(itt,1,2,2) = irr
              itt = nestitags2(ia)%p(ii,4)
              itagreq2(itt,1,1,2) = irl
              itagreq2(itt,1,2,2) = irr
            endif 
          endif
        enddo
      enddo

      ! then, alloc the send buffers and define the request ranges for each itag
      ! as well as the sizes of the send buffers:
      irr   = 0
      isreq = 0
      do ia=1,narea
        do ii=1,nestinglevels(ia)
          !  iia: fine grid
          !  ia:  coarse grid
          iia = nestingto(ia)%p(ii)
          iao = enclosing(iia,ia)

          ! u-brd
          irl = 0
          if (uhtable(ia,iia)) then
            do it=1,mpi_size
              isreq = isreq + 1
              if (it == iam) cycle
              bs = itmp(isreq,2,1)
              if (bs > 0) then
                irr = irr + 1
                if (irl == 0) irl = irr
                bsc = icmp(isreq,2,1)
                allocate( sbufnb2(irr)%p(1:bs), sbufnb3(irr)%p(1:bsc) )
              endif
            enddo
            if (irl > 0) then
              itt = nestitags2(ia)%p(ii,1)
              itagreq2(itt,2,1,1) = irl
              itagreq2(itt,2,2,1) = irr
              itt = nestitags2(ia)%p(ii,3)
              itagreq2(itt,2,1,1) = irl
              itagreq2(itt,2,2,1) = irr
            endif 
          endif

          ! v-brd
          ill = irr
          irl = ill
          if (vhtable(ia,iia)) then
            do it=1,mpi_size
              isreq = isreq + 1
              if (it == iam) cycle
              bs = itmp(isreq,2,2)
              if (bs > 0) then
                irr = irr + 1
                if (irl == ill) irl = irr
                bsc = icmp(isreq,2,2)
                allocate( sbufnb2(irr)%p(1:bs), sbufnb3(irr)%p(1:bsc) )
              endif
            enddo
            if (irl > ill) then
              itt = nestitags2(ia)%p(ii,2)
              itagreq2(itt,2,1,2) = irl
              itagreq2(itt,2,2,2) = irr
              itt = nestitags2(ia)%p(ii,4)
              itagreq2(itt,2,1,2) = irl
              itagreq2(itt,2,2,2) = irr
            endif 
          endif
        enddo
      enddo

    endif

    deallocate( itmp, icmp )
#else
    allocate( itagreq2(1,1,1,1), rsmask2(1,1) )
    allocate( rsmask2(1,1)%p(1,1) )
    allocate( bsiz3(1:mpi_nt,1:1,1:2,1:2) )
#endif

  end subroutine dmpi_init4

!===============================================================================

  subroutine dmpi_init5( narea, nestinglevels, nestingto, nestingfrom,         &
                         enclosing, kmx, znesting, znest_serial, krz, iga, jga )
                         
    ! prepare for non-blocking MPI-comm around nesting 7
    ! that is, for dmpi_gather_bnd_data and dmpi_scatter_bnd_data

    use cmod_mem, only : cmi1, cmi2
    implicit none

    integer(4), intent(in) :: narea, nestinglevels(:), enclosing(:,:), kmx(:)
    type(cmi1), intent(in) :: nestingto(:), nestingfrom(:)
    type(cmi2), intent(in) :: znesting(:), znest_serial(:)
    type(cmi2), intent(in) :: krz(:), iga(:), jga(:)

    integer(4) :: iam, MaxRreq, MaxSreq, Maxitag, ia, iia, iao, bs, inreq, it
    integer(4) :: ii, iif, maxnestlev, ig1, ig2, jg1, jg2, krz3, ibnd, jbnd
    integer(4) :: i, j, irreq, isreq, iz, irr, irl, itt, nz1, nz2
    integer(4), allocatable :: icmp(:,:)

    ! some initial stuff:
    maxnestlev = 0
    allocate( nestitags7(narea) )
    do ia=1,narea
      ii = nestinglevels(ia)
      if (ii > 0) then  
        maxnestlev = max( maxnestlev, ii )
        allocate( nestitags7(ia)%p(1:ii) )
        nestitags7(ia)%p(1:) = 0
      endif
    enddo


#if defined (MPI)
    ! Initially, allow that all areas need to do the maximum number
    ! of MPI-communications around nesting:
    allocate(icmp(narea*maxnestlev*mpi_size, 2), rsmask7(narea,mpi_size))
    icmp(:,:) = 0  ! temp size of r/s buffers
    do ia=1,narea
      ii = nestinglevels(ia)
      if (ii > 0) then
        do it=1,mpi_size
          allocate( rsmask7(ia,it)%p(1:ii) )
          rsmask7(ia,it)%p(1:) = 0
        enddo
      endif
    enddo

    iam = mpi_rank + 1
    MaxRreq = 0 
    MaxSreq = 0 
    Maxitag = NextFreeItag - 1 ! NextFreeItag first set in init3+init4
    irreq   = 0
    isreq   = 0
    do ia=1,narea
      do ii=1,nestinglevels(ia)
        !  iia: fine grid
        !  ia:  coarse grid
        iia = nestingto(ia)%p(ii)
        iao = enclosing(iia,ia)
        iif = 1
        do while(ia /= nestingfrom(iia)%p(iif))
          iif = iif + 1
        enddo

        ! define tags:
        nestitags7(ia)%p(ii) = Maxitag + 1   ! dmpi_gather_bnd_data
        Maxitag = Maxitag + 1

        ig1 = iga(iia)%p(1,iao)
        ig2 = iga(iia)%p(2,iao)
        jg1 = jga(iia)%p(1,iao)
        jg2 = jga(iia)%p(2,iao)

        ! first, set up for the gather, both recv and send ---------------------
        !
        if (mpi_io_serial) then
          nz1 = znest_serial(iia)%p(iif,1)
          nz2 = znest_serial(iia)%p(iif,2)
          do it=1,mpi_size
            irreq = irreq + 1
            !  I should gather but not from myself
            if (it == iam) cycle
            !  find recv buffer size:
            bs = 0
            do iz=nz1,nz2
              ibnd = krz(iia)%p(1,iz)
              jbnd = krz(iia)%p(2,iz)
              i = ig1+int((ibnd-1)/ig2,4)
              j = jg1+int((jbnd-1)/jg2,4)
              krz3 = krz(iia)%p(3,iz)
              if     (krz3 == 1) then
                j = j-1
              elseif (krz3 == 2) then
                i = i-1
              elseif (krz3 == 3) then
                j = j+1
              elseif (krz3 == 4) then
                i = i+1
              endif
              if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.     &
                  j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
              bs = bs + kmx(iia)*mpi_nc
            enddo
            if (bs > 0) then
              MaxRreq              = MaxRreq + 1
              icmp(irreq,1)        = max(mpi_nc*bs,1)
              rsmask7(ia,it)%p(ii) = rsmask7(ia,it)%p(ii) + 1
            endif
          enddo

        else

          !  should I contribute to the gathering by sending data ? ------------
          if (nctable(iam,ia,iia)) then
            nz1 = znesting(iia)%p(iif,1)
            nz2 = znesting(iia)%p(iif,2)
          else
            nz1 =  1
            nz2 = -1
          endif
          isreq = isreq + 1
          !  find send buffer size:
          bs = 0
          do iz=nz1,nz2
            ibnd = krz(iia)%p(1,iz)
            jbnd = krz(iia)%p(2,iz)
            i = ig1+int((ibnd-1)/ig2,4)
            j = jg1+int((jbnd-1)/jg2,4)
            krz3 = krz(iia)%p(3,iz)
            if     (krz3 == 1) then
              j = j-1
            elseif (krz3 == 2) then
              i = i-1
            elseif (krz3 == 3) then
              j = j+1
            elseif (krz3 == 4) then
              i = i+1
            endif
            if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                     &
                j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
            bs = bs + kmx(iia)*mpi_nc
          enddo
          if (bs > 0) then
            MaxSreq               = MaxSreq + 1
            icmp(isreq,2)         = max(mpi_nc*bs,1)
            rsmask7(ia,iam)%p(ii) = rsmask7(ia,iam)%p(ii) + 10
          endif

        endif

      enddo ! ii
    enddo   ! ia


    allocate( itagreq7(NextFreeItag:Maxitag,1:2,1:2) )
    Maxitag      = Maxitag - NextFreeItag + 1
    NextFreeItag = NextFreeItag + Maxitag
    itagreq7(:,:,:) = 0
    if (MaxRreq + MaxSreq > 0) then
      allocate(lstat7(MPI_STATUS_SIZE,MaxRreq+MaxSreq))
      if (MaxRreq > 0) allocate(ireq7(MaxRreq), rbufnb7(MaxRreq))
      if (MaxSreq > 0) allocate(irqs7(MaxSreq), sbufnb7(MaxSreq))

      ! now, alloc the recv buffers and define the request ranges for each itag
      ! as well as the sizes of the recv buffers:
      irr   = 0
      irreq = 0
      do ia=1,narea
        do ii=1,nestinglevels(ia)
          irl = 0
          if (mpi_io_serial) then
            do it=1,mpi_size
              irreq = irreq + 1
              if (it == iam) cycle
              bs = icmp(irreq,1)
              if (bs > 0) then
                irr = irr + 1
                if (irl == 0) irl = irr
                allocate( rbufnb7(irr)%p(1:bs) )
              endif
            enddo
            if (irl > 0) then
              itt = nestitags7(ia)%p(ii)
              itagreq7(itt,1,1) = irl
              itagreq7(itt,1,2) = irr
            endif 
          endif
        enddo
      enddo
      ! then, alloc the send buffers and define the request ranges for each itag
      ! as well as the sizes of the send buffers:
      irr   = 0
      isreq = 0
      do ia=1,narea
        do ii=1,nestinglevels(ia)
          irl = 0
          if (.not.mpi_io_serial) then
            isreq = isreq + 1
            bs = icmp(isreq,2)
            if (bs > 0) then
              irr = irr + 1
              if (irl == 0) irl = irr
              allocate( sbufnb7(irr)%p(1:bs) )
            endif
            if (irl > 0) then
              itt = nestitags7(ia)%p(ii)
              itagreq7(itt,2,1) = irl
              itagreq7(itt,2,2) = irr
            endif
          endif
        enddo
      enddo

    endif

    deallocate( icmp )
#else
    allocate( itagreq7(1,1,1), rsmask7(1,1) )
    allocate( rsmask7(1,1)%p(1) )
#endif

  end subroutine dmpi_init5

!===============================================================================

  subroutine bathy2mpibin(fnam,mm1,smm1,hz,shz)
    implicit none
    character(256), intent(in)  :: fnam
    integer(4),     intent(in)  :: mm1(:,0:,0:)
    integer(4),     intent(in)  :: smm1,shz
    real(8),        intent(in)  :: hz(0:)
#if defined (MPI)
    integer(4) :: ierr, fh
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fnam, MPI_MODE_WRONLY +                 &
       MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr)
    call MPI_FILE_WRITE(fh,mm1(:,0:,0:),smm1,MPI_INTEGER4,                     &
       MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_WRITE(fh,hz(0:),shz,MPI_REAL8,MPI_STATUS_IGNORE, ierr)
    call MPI_FILE_CLOSE(fh,ierr)
#endif
  end subroutine bathy2mpibin

#if defined (MPI)
  subroutine bathyread_mm(fnam,mm1,smm1,fh)
    implicit none
    character(256), intent(in)  :: fnam
    integer(4),     intent(out) :: mm1(:,0:,0:)
    integer(4),     intent(in)  :: smm1
    integer(4),     intent(out) :: fh
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fnam, MPI_MODE_RDONLY, MPI_INFO_NULL,   &
                       fh, ierr)
    call MPI_FILE_READ_ALL(fh,mm1(:,0:,0:),smm1,MPI_INTEGER4,MPI_STATUS_IGNORE,&
                           ierr)
  end subroutine bathyread_mm

  subroutine bathyread_hz(fh,hz,iw3)
    integer(4), intent(inout)  :: fh
    integer(4), intent(in)  :: iw3
    real(8),    intent(out) :: hz(0:)
    call MPI_FILE_READ_ALL(fh,hz(0:),iw3+1,MPI_REAL8,MPI_STATUS_IGNORE, ierr)
    call MPI_FILE_CLOSE(fh,ierr)
  end subroutine bathyread_hz
#else
  subroutine bathyread_mm(fnam,mm1,smm1,fh)
    implicit none
    character(256), intent(in) :: fnam
    integer(4),     intent(in) :: mm1(:,0:,0:)
    integer(4),     intent(in) :: smm1
    integer(4),     intent(in) :: fh
    return
  end subroutine bathyread_mm

  subroutine bathyread_hz(fh,hz,iw3)
    integer(4), intent(in) :: fh
    integer(4), intent(in) :: iw3
    real(8),    intent(in) :: hz(0:)
    return
  end subroutine bathyread_hz
#endif

!===============================================================================

  subroutine decompose_read(narea, mmx, nmx, kmx, iw2, iw3, msrf, kh)
    use io_subs,           only : io_new_unit
    use cmod_mem,          only : cmi1, cmi2

    implicit none

    integer(4),     intent(in) :: narea
    integer(4),     intent(in) :: mmx(:), nmx(:), kmx(:), iw2(:), iw3(:)
    type (cmi1),     intent(in) :: kh(:)
    type (cmi2),     intent(in) :: msrf(:)

    integer(4)              :: ia, lun, ios, ilc, iia, iil, iiu, ijl, iju, iam
    integer(4)              :: in_w, in_n, in_e, in_s, ns, nw, i, j, k, it
    integer(4)              :: in_nw, in_ne, in_se, in_sw
    integer(4)              :: ibnd, jbnd, n
    logical, allocatable    :: found(:)
#if defined (MPI)
    integer(4)              :: s_buf, r_buf, g_buf, l_buf, i_buf, nact
    integer(4), allocatable :: jtask(:,:), jlows(:,:), jhghs(:,:), jnsrf(:,:)
    logical                 :: empty(1)
#endif

    !- Set up the distribution data and task table for one task ----------------
    if (mpi_size == 1) then
      it = mpi_size
      allocate( dd(narea), mpi_tt(narea,it), nbnd(it), tscat(it,narea,narea),  &
                nctable(it,narea,narea), nftable(it,narea,narea),              &
                uctable(it,narea,narea), vctable(it,narea,narea),              &
                uftable(it,narea,narea), vftable(it,narea,narea),              &
                cctable(narea,narea), cftable(narea,narea),                    &
                uhtable(narea,narea), vhtable(narea,narea),                    &
                mctable(narea,narea) )

      do ia=1,narea
        !  distribution data:
        dd(ia)%iwet2  = iw2(ia)
        dd(ia)%iwet3  = iw3(ia)
        dd(ia)%ndims  = 3
        dd(ia)%dim1   = mmx(ia)
        dd(ia)%dim2   = nmx(ia)
        dd(ia)%dim3   = kmx(ia)
        dd(ia)%low_ws = 1
        dd(ia)%up_ws  = iw2(ia)
        dd(ia)%nwet3  = iw3(ia)
        dd(ia)%low_w3 = iw2(ia) + 1
        dd(ia)%low_i  = 1
        dd(ia)%up_i   = mmx(ia)
        dd(ia)%low_j  = 1
        dd(ia)%up_j   = nmx(ia)
        dd(ia)%low_hi = 1
        dd(ia)%up_hi  = mmx(ia)
        dd(ia)%low_hj = 1
        dd(ia)%up_hj  = nmx(ia)
        dd(ia)%nt_w   = -1
        dd(ia)%nt_n   = -1
        dd(ia)%nt_e   = -1
        dd(ia)%nt_s   = -1
        dd(ia)%nt_nw  = -1
        dd(ia)%nt_ne  = -1
        dd(ia)%nt_se  = -1
        dd(ia)%nt_sw  = -1
        dd(ia)%recv_bufsize(:) = 0
        dd(ia)%send_bufsize(:) = 0
        dd(ia)%recv_buf_srf(:) = 0
        dd(ia)%send_buf_srf(:) = 0
        dd(ia)%halo2 = 0
        dd(ia)%halo3 = 0

        !  task table:
        mpi_tt(ia,it)%nwet3  = iw3(ia)
        mpi_tt(ia,it)%low_w3 = iw2(ia) + 1
        mpi_tt(ia,it)%low_ws = 1
        mpi_tt(ia,it)%up_ws  = iw2(ia)
        mpi_tt(ia,it)%low_i  = 1
        mpi_tt(ia,it)%up_i   = mmx(ia)
        mpi_tt(ia,it)%low_j  = 1
        mpi_tt(ia,it)%up_j   = nmx(ia)
        mpi_tt(ia,it)%low_hi = 1
        mpi_tt(ia,it)%up_hi  = mmx(ia)
        mpi_tt(ia,it)%low_hj = 1
        mpi_tt(ia,it)%up_hj  = nmx(ia)
        mpi_tt(ia,it)%hnb(:) = -1

      enddo
      !  boundary points:
      nbnd(it) = 0

    !- Set up the distribution data and task table for more than one task ------
    else
      ! who am I ?
      iam = mpi_rank + 1

      ! open file
      lun = io_new_unit()
      open (unit=lun, file=trim(mpi_decomp_file), status='old', iostat=ios)
      if (ios /= 0) call exitme(1,'Cannot open '//trim(mpi_decomp_file))
      read(lun, '(I5)', iostat=ios) it
      if (ios /= 0) then
        call exitme(1,'Cannot read number of tasks from '                      &
                       //trim(mpi_decomp_file))
      endif
      if (it /= mpi_size) call exitme(1,'Invalid number of tasks specified')
      allocate( dd(narea), mpi_tt(narea,it), found(narea), nbnd(it),           &
                tscat(it,narea,narea),                                         &
                nctable(it,narea,narea), nftable(it,narea,narea),              &
                uctable(it,narea,narea), vctable(it,narea,narea),              &
                uftable(it,narea,narea), vftable(it,narea,narea),              &
                cctable(narea,narea), cftable(narea,narea),                    &
                uhtable(narea,narea), vhtable(narea,narea),                    &
                mctable(narea,narea) )
     
      ! search the file for decomp info for each area:
      found(:) = .false.
      taskloop: do it=1,mpi_size
        do ia=1,narea
          read(lun, '(6I5)', iostat=ios) ilc, iia, iil, iiu, ijl, iju
          if (ios /= 0) then
            call exitme(1,'Cannot read grid extent data from '                 &
                           //trim(mpi_decomp_file))
          endif
          read(lun, '(10x,4I5)', iostat=ios) in_w, in_n, in_e, in_s
          if (ios /= 0) then
            call exitme(1,'Cannot read neighbour data from '                   &
                           //trim(mpi_decomp_file))
          endif
          read(lun, '(10x,4I5)', iostat=ios) in_nw, in_ne, in_se, in_sw
          if (ios /= 0) then
            call exitme(1,'Cannot read neighbour data from '                   &
                           //trim(mpi_decomp_file))
          endif

          ! test for sanity:
          if (ilc /= it .or. iia /= ia) then
            write(iu06,*) 'Failed reading mpi decomp for task and area'
            ios = 1
          endif
          if ( ( iil > iiu .or. iiu > mmx(ia) .or.                             &
                (iil < 1 .and. iil /= 0 .and. iiu /= 0)) .or.                  &
               ( ijl > iju .or. iju > nmx(ia) .or.                             &
                (ijl < 1 .and. ijl /= 0 .and. iju /= 0))      ) then
            write(iu06,*) 'Invalid mpi decomp for task and area'
            ios = 1
          endif
          if (in_w > mpi_size .or. in_w == it .or. in_w == 0 .or.              &
              in_n > mpi_size .or. in_n == it .or. in_n == 0 .or.              &
              in_e > mpi_size .or. in_e == it .or. in_e == 0 .or.              &
              in_s > mpi_size .or. in_s == it .or. in_s == 0 ) then
            write(iu06,*) 'Invalid neighbour task number for task and area'
            ios = 1
          endif
          if (in_nw > mpi_size .or. in_nw == it .or. in_nw == 0 .or.           &
              in_ne > mpi_size .or. in_ne == it .or. in_ne == 0 .or.           &
              in_se > mpi_size .or. in_se == it .or. in_se == 0 .or.           &
              in_sw > mpi_size .or. in_sw == it .or. in_sw == 0 ) then
            write(iu06,*) 'Invalid corner task number for task and area'
            ios = 1
          endif

          ! Oh no, got trapped ...
          if (ios /= 0) then
            write(iu06,'(i5,a5,i5)') iam, ' and ', ia
            write(iu06,*) 'Got this set of input data:'
            write(iu06,'(6I5)') ilc, iia, iil, iiu, ijl, iju
            write(iu06,'(10x,4I5)') in_w, in_n, in_e, in_s
            write(iu06,'(10x,4I5)') in_nw, in_ne, in_se, in_sw
            call exitme(1)
          endif

          ! Set up part of the task table:
          mpi_tt(ia,it)%low_i = iil
          mpi_tt(ia,it)%up_i  = iiu
          mpi_tt(ia,it)%low_j = ijl
          mpi_tt(ia,it)%up_j  = iju
          if (in_n > 0) then
            mpi_tt(ia,it)%low_hi = max(1, iil - halo_width)
          else
            mpi_tt(ia,it)%low_hi = iil
          endif
          if (in_s > 0) then
            mpi_tt(ia,it)%up_hi  = min(mmx(ia), iiu + halo_width)
          else
            mpi_tt(ia,it)%up_hi  = iiu
          endif
          if (in_w > 0) then
            mpi_tt(ia,it)%low_hj = max(1, ijl - halo_width)
          else
            mpi_tt(ia,it)%low_hj = ijl
          endif
          if (in_e > 0) then
            mpi_tt(ia,it)%up_hj  = min(nmx(ia), iju + halo_width)
          else
            mpi_tt(ia,it)%up_hj  = iju
          endif
          mpi_tt(ia,it)%hnb(1) = in_w
          mpi_tt(ia,it)%hnb(2) = in_n
          mpi_tt(ia,it)%hnb(3) = in_e
          mpi_tt(ia,it)%hnb(4) = in_s
          mpi_tt(ia,it)%hnb(5) = in_nw
          mpi_tt(ia,it)%hnb(6) = in_ne
          mpi_tt(ia,it)%hnb(7) = in_se
          mpi_tt(ia,it)%hnb(8) = in_sw

          ! did we find me?
          if (ilc /= iam) cycle
          found(ia) = .true. 

          ! Ok, lets assign part of the distribution data:
          dd(ia)%iwet2 = iw2(ia)
          dd(ia)%iwet3 = iw3(ia)
          dd(ia)%ndims = 3
          dd(ia)%dim1  = mmx(ia)
          dd(ia)%dim2  = nmx(ia)
          dd(ia)%dim3  = kmx(ia)
          dd(ia)%low_i = iil
          dd(ia)%up_i  = iiu
          dd(ia)%low_j = ijl
          dd(ia)%up_j  = iju
          if (in_n > 0) then
            dd(ia)%low_hi = max(1, iil - halo_width)
          else
            dd(ia)%low_hi = iil
          endif
          if (in_s > 0) then
            dd(ia)%up_hi  = min(mmx(ia), iiu + halo_width)
          else
            dd(ia)%up_hi  = iiu
          endif
          if (in_w > 0) then
            dd(ia)%low_hj = max(1, ijl - halo_width)
          else
            dd(ia)%low_hj = ijl
          endif
          if (in_e > 0) then
            dd(ia)%up_hj  = min(nmx(ia), iju + halo_width)
          else
            dd(ia)%up_hj  = iju
          endif
          dd(ia)%nt_w  = in_w
          dd(ia)%nt_n  = in_n
          dd(ia)%nt_e  = in_e
          dd(ia)%nt_s  = in_s
          dd(ia)%nt_nw = in_nw
          dd(ia)%nt_ne = in_ne
          dd(ia)%nt_se = in_se
          dd(ia)%nt_sw = in_sw
        enddo

        !- for main area, set No. of open bnd points on each task:
        ia = mpi_ma
        nbnd(it) = 0
        if (mpi_io_serial) nbnd(mpi_io_rank+1) = mpi_nz

      enddo taskloop

      !- check if we succeeded:
      do ia=1,narea
        if (.not.found(ia)) then
          write(iu06,*) 'ERROR: Did not find task decomposition for' 
          write(iu06,'(a5,i5,a10,i5)') 'task ', iam, ' and area ', ia
          call exitme(1)
        endif
      enddo

      !- clean up, please:
      close(lun)
      deallocate( found )

      !- Calculate # of wet points and range of surface wet points for 
      !  each task:
      !
      !  No. of wet points in 3D and (temporary) No. of surface wet points:
      do it=1,mpi_size
        do ia=1,narea
          !  excl. the halo:
          ns = 0
          nw = 0
          do j=mpi_tt(ia,it)%low_j,mpi_tt(ia,it)%up_j
            iloop: do i=mpi_tt(ia,it)%low_i,mpi_tt(ia,it)%up_i
              if (msrf(ia)%p(i,j) == 0) cycle iloop

              ns = ns + 1
              nw = nw + kh(ia)%p(msrf(ia)%p(i,j))
            enddo iloop
          enddo
          mpi_tt(ia,it)%up_ws = ns  ! temp. storage, modified below
          mpi_tt(ia,it)%nwet3 = nw
        enddo
      enddo
      !
      ! Range of surface wet points; possibly correct 3D wet points:
      do ia=1,narea
        nw = 0
        do it=1,mpi_size
          ns = mpi_tt(ia,it)%up_ws
          if (ns > 0) then
            mpi_tt(ia,it)%low_ws = nw + 1
            nw = nw + ns
            mpi_tt(ia,it)%up_ws  = nw
          else
            mpi_tt(ia,it)%low_ws = 0
            mpi_tt(ia,it)%up_ws  = 0
            ! make sure to have zero wet points in task: 
            mpi_tt(ia,it)%nwet3  = 0
          endif
        enddo
      enddo
      !
      !- Obtain start index for subsurface wetpoints:
      do ia=1,narea
        nw = iw2(ia)
        do it=1,mpi_size
          if (mpi_tt(ia,it)%nwet3 > 0) then
            mpi_tt(ia,it)%low_w3 = nw + 1
            ns = mpi_tt(ia,it)%up_ws - mpi_tt(ia,it)%low_ws + 1
            nw = nw + mpi_tt(ia,it)%nwet3 - ns
          endif
        enddo
      enddo
      !
      !- Transfer data from task table to distribution data for this task:
      it = iam
      do ia=1,narea
        !  srf points:
        dd(ia)%low_ws = mpi_tt(ia,it)%low_ws
        dd(ia)%up_ws  = mpi_tt(ia,it)%up_ws
        !  wet points:
        dd(ia)%nwet3  = mpi_tt(ia,it)%nwet3
        dd(ia)%low_w3 = mpi_tt(ia,it)%low_w3
      enddo

#if defined (MPI)
      !- set up the halo-communication masks -----------------------------------
      call dmpi_init_comm_masks (narea)
      do ia=1,narea
        call dmpi_set_comm_masks_msrf (ia, msrf(ia)%p)
      enddo
#endif

      !- Obtain the halo sizes and halo-buffersizes ----------------------------
      do ia=1,narea
        dd(ia)%halo2           = 0
        dd(ia)%halo3           = 0
        dd(ia)%recv_bufsize(:) = 0
        dd(ia)%recv_buf_srf(:) = 0
        dd(ia)%send_bufsize(:) = 0
        dd(ia)%send_buf_srf(:) = 0
      enddo

      ia=1
#if defined (MPI)
      do ia=1,narea
        call dmpi_halosize_srf( kmx(ia), msrf(ia)%p, kh(ia)%p,                 &
                            dd(ia)%low_i, dd(ia)%up_i,                         &
                            dd(ia)%low_j, dd(ia)%up_j,                         &
                            dd(ia)%low_hi, dd(ia)%up_hi,                       &
                            dd(ia)%low_hj, dd(ia)%up_hj,                       &
                            dd(ia)%halo3, dd(ia)%halo2 )
        call dmpi_bufsize_srf( kmx(ia), msrf(ia)%p, kh(ia)%p, ia,              &
                           dd(ia)%recv_bufsize, dd(ia)%send_bufsize,           &
                           dd(ia)%recv_buf_srf, dd(ia)%send_buf_srf )
      enddo
#endif
      

#if defined (MPI)
      !- allocate the communication buffers ------------------------------------
      s_buf = 1
      r_buf = 1
      g_buf = 1
      l_buf = 1
      do ia=1,narea
        s_buf = max( s_buf, maxval(dd(ia)%send_bufsize(:)) )
        r_buf = max( r_buf, maxval(dd(ia)%recv_bufsize(:)) )
        l_buf = max( l_buf, dd(ia)%iwet2 + 1 )
        g_buf = max( g_buf, dd(ia)%iwet3 + 1 )
      enddo
      i_buf = l_buf
      l_buf = max( l_buf, s_buf, r_buf )
      s_buf = s_buf*mpi_nc*3
      r_buf = r_buf*mpi_nc*3
      g_buf = max(g_buf*mpi_nc, mpi_nz*kmx(mpi_ma)*mpi_nc, i_buf*(mpi_ni+1))
      allocate( sbuf(s_buf), rbuf(r_buf), gbuf(g_buf),                         &
                lbuf(l_buf) )

      ! trim the size of recv buffers used with non-blocking halo swaps:
      allocate( rbuf2(max_nb), tfrbuf(max_nb), tfsbuf(max_nb) )
      do i=1,max_nb
        r_buf = 1
        s_buf = 1
        do ia=1,narea
          r_buf = max( r_buf, 6*mpi_nc*dd(ia)%recv_bufsize(i),                 &
                          (2*mpi_ni+1)*dd(ia)%recv_buf_srf(i) )
          s_buf = max( s_buf, 6*mpi_nc*dd(ia)%send_bufsize(i),                 &
                          (2*mpi_ni+1)*dd(ia)%send_buf_srf(i) )
        enddo
        allocate( rbuf2(i)%p(4*r_buf), tfrbuf(i)%p(r_buf), tfsbuf(i)%p(s_buf) )
      enddo
#endif

    endif

    !- better safe than sorry --------------------------------------------------
#if defined (MPI)
    call dmpi_barrier(MPI_COMM_WORLD) 
#else
    call dmpi_barrier( ) 
#endif

  end subroutine decompose_read


  subroutine decompose_gen(narea, mmx, nmx, kmx, iw2, iw3, msrf, kh, nproci, nprocj ) 

    use io_subs,           only : io_new_unit
    use cmod_mem,          only : cmi1,cmi2
    use auto_decompo, only : autodecompose

    implicit none

    integer(4),     intent(in) :: narea, nproci, nprocj
    integer(4),     intent(in) :: mmx(:), nmx(:), kmx(:), iw2(:), iw3(:)
    type (cmi2),    intent(in) :: msrf(:)
    type (cmi1),    intent(in) :: kh(:)

    integer(4)              :: ia

    !- simple autodecompose ----------------------------------------------------
    if (nproci >= 1 .and. nprocj > 1) then
      if (nproci == 1 .and. nprocj > 1) then
        ! get forbidden i-lines if any:
        ! fixme: hard-coded to take max 99 forbidden border lines (haha) and 
        !        boundary lines of each type in each area. Consider making this
        !        dynamically adjustable if required for future setups.
        ia = narea
        ! in serial mode, make decompositions and terminate
        call autodecompose(nprocj, 1, mmx, nmx, iw3, msrf, kh, iu06)
        call exitme(0,'Auto-decompositions made. Exit.')
      endif
    endif

  end subroutine decompose_gen


!===============================================================================

  subroutine dmpi_scatter_table(narea, mmx, nmx, iga, jga,                     &
                                nestinglevels, enclosing, nestingto)

    use cmod_mem, only : cmi2, cmi1

    implicit none 

    integer(4),  intent(in) :: narea, mmx(:), nmx(:)
    integer(4),  intent(in) :: enclosing(:,:), nestinglevels(:)
    type (cmi2), intent(in) :: iga(:), jga(:)
    type (cmi1), intent(in) :: nestingto(:)

    integer(4) :: ia, iia, it, iao, ii
    integer(4) :: ifl, ifu, jfl, jfu, il, iu, jl, ju, ig2, jg2

    ! OBS: tscat(:,:,:) was allocated in decompose_read()

    if (mpi_size == 1) then
      !  scatter table:
      tscat(:,:,:) = .false.

    else
      do it=1,mpi_size
        !- set up scatter table:
        tscat(it,:,:) = .false.
        do ia=1,narea
          do ii=1,nestinglevels(ia)
            !  iia: fine grid
            !  ia:  coarse grid
            iia = nestingto(ia)%p(ii)
            iao = enclosing(iia,ia)

            !  Identity if coarse and fine grid overlap on this task:
            ig2 = iga(iia)%p(2,iao)
            ifl = iga(iia)%p(1,iao) + (ig2/2+iga(iia)%p(3,iao))/ig2
            ifu = iga(iia)%p(1,iao) + (ig2/2+iga(iia)%p(3,iao)+mmx(iia)-1)/ig2
            jg2 = jga(iia)%p(2,iao)
            jfl = jga(iia)%p(1,iao) + (jg2/2+jga(iia)%p(3,iao))/jg2
            jfu = jga(iia)%p(1,iao) + (jg2/2+jga(iia)%p(3,iao)+nmx(iia)-1)/jg2
            il  = mpi_tt(ia,it)%low_hi
            iu  = mpi_tt(ia,it)%up_hi
            jl  = mpi_tt(ia,it)%low_hj
            ju  = mpi_tt(ia,it)%up_hj
            tscat(it,ia,iia) = ((ifl <= il .and. il <= ifu) .or.               &
                                (ifl <= iu .and. iu <= ifu)      ) .and.       &
                               ((jfl <= jl .and. jl <= jfu) .or.               &
                                (jfl <= ju .and. ju <= jfu)      )
          enddo
        enddo
      enddo
    endif

  end subroutine dmpi_scatter_table

!===============================================================================

#if defined (MPI)

  subroutine decompose_plane_set_bounds(na, nproci, nprocj, mmx, nmx, kmx,     &
                                        iw2, iw3, mm1, empty)
    use io_subs, only : io_new_unit

    implicit none

    integer(4), intent(in)  :: na, nproci, nprocj
    integer(4), intent(in)  :: mmx(:), nmx(:), kmx(:), iw2(:), iw3(:)
    integer(4), intent(in)  :: mm1(0:,0:,1:)
    logical,    intent(out) :: empty

    integer(4) :: ios, lun, ia, it, il, iu, jl, ju, n, tmin, tave, tmax, tnul
    integer(4) :: n_w, n_n, n_e, n_s, n_nw, n_ne, n_se, n_sw, newmpi
    integer(4) :: i, j, k, i_pe_grid, j_pe_grid
    integer(4) :: i_slices(nproci), j_slices(nprocj)
    integer(4) :: pe_grid(nproci,nprocj), tmp(nproci*nprocj,1:2)
    integer(4) :: ilows(nproci), jlows(nprocj)
    integer(4) :: i_slice, j_slice, rem_i, rem_j

    ia = na

    i_slice = mmx(ia) / nproci
    rem_i   = mmx(ia) - nproci * i_slice
    i_slices(rem_i+1:nproci) = i_slice
    if (rem_i > 0) i_slices(1:rem_i) = i_slice + 1

    j_slice = nmx(ia) / nprocj
    rem_j   = nmx(ia) - nprocj * j_slice
    j_slices(rem_j+1:nprocj) = j_slice
    if (rem_j > 0) j_slices(1:rem_j) = j_slice + 1

    do i = 1, nproci
      if (i-1 < rem_i) then
        ilows(i) = 1 + (i-1)*(i_slices(i))
      else
        ilows(i) = 1 + (i-1)*(i_slices(i)) + rem_i
      endif
    enddo
    do j = 1, nprocj
      if (j-1 < rem_j) then
        jlows(j) = 1 + (j-1)*(j_slices(j))
      else
        jlows(j) = 1 + (j-1)*(j_slices(j)) + rem_j
      endif
    enddo

    lun = io_new_unit()
    open (unit=lun, file=trim(mpi_decomp_file), status='replace', iostat=ios)
    if (ios /= 0) then
      call exitme(1,'Cannot create '//trim(mpi_decomp_file))
    endif
    write(lun, '(I5)') mpi_size
    pe_grid(:,:) = -1
    do it=1,mpi_size
      ! define pe_grid:
      i_pe_grid = mod(it-1,nproci)
      j_pe_grid = int((it-1)/nproci,4)
      pe_grid(i_pe_grid+1,j_pe_grid+1) = it
      
      ! define index range:
      il = ilows(i_pe_grid+1)
      iu = ilows(i_pe_grid+1) + i_slices(i_pe_grid+1) - 1
      jl = jlows(j_pe_grid+1)
      ju = jlows(j_pe_grid+1) + j_slices(j_pe_grid+1) - 1
      write(lun, '(6I5)') it, ia, il, iu, jl, ju

      ! define neighbours:
      if (i_pe_grid < nproci-1) then
        n_s = i_pe_grid + 1 + j_pe_grid*nproci + 1
      else
        n_s = -1
      endif
      if (i_pe_grid > 0) then
        n_n = i_pe_grid - 1 + j_pe_grid*nproci + 1
      else
        n_n = -1
      endif
      if (j_pe_grid > 0) then
        n_w = i_pe_grid + (j_pe_grid-1)*nproci + 1
      else
        n_w = -1
      endif
      if (j_pe_grid < nprocj-1) then
        n_e = i_pe_grid + (j_pe_grid+1)*nproci + 1
      else
        n_e = -1
      endif
      write(lun, '(10x,4I5)') n_w, n_n, n_e, n_s

      ! define corner-neighbours:
      if (j_pe_grid > 0 .and. i_pe_grid > 0) then
        n_nw = (i_pe_grid-1) + (j_pe_grid-1)*nproci + 1
      else
        n_nw = -1
      endif
      if (j_pe_grid < nprocj-1 .and. i_pe_grid > 0) then
        n_ne = (i_pe_grid-1) + (j_pe_grid+1)*nproci + 1
      else
        n_ne = -1
      endif
      if (j_pe_grid < nprocj-1 .and. i_pe_grid < nproci-1) then
        n_se = (i_pe_grid+1) + (j_pe_grid+1)*nproci + 1
      else
        n_se = -1
      endif
      if (j_pe_grid > 0 .and. i_pe_grid < nproci-1) then
        n_sw = (i_pe_grid+1) + (j_pe_grid-1)*nproci + 1
      else
        n_sw = -1
      endif
      write(lun, '(10x,4I5)') n_nw, n_ne, n_se, n_sw
    enddo
    close(lun)

    ! print stats for this decompo:
    write(iu06,*) 'Statistics for ',nproci,' by ',nprocj,' auto-decomposition'
    write(iu06,*) '  Surface, No. wetpoints:      ',iw2(ia)
    write(iu06,*) '           even split:         ',iw2(ia)/mpi_size
    tmin = iw3(ia) + 1
    tave = 0
    tmax = -tmin
    tnul = 0
    tmp(1:mpi_size,1:2) = -1
    il = 0
    iu = 0
    do it=1,mpi_size
      i_pe_grid = mod(it-1,nproci)
      j_pe_grid = int((it-1)/nproci,4)
      n = 0
      do j=jlows(j_pe_grid+1),jlows(j_pe_grid+1)+j_slices(j_pe_grid+1)-1      
        do i=ilows(i_pe_grid+1),ilows(i_pe_grid+1)+i_slices(i_pe_grid+1)-1
          if (mm1(i,j,1) > 0) n = n + 1
        enddo
      enddo
      if (n == 0) then 
        tnul = tnul + 1
        pe_grid(i_pe_grid+1,j_pe_grid+1) = -1
      endif
      tmin = min(tmin, n)
      tmax = max(tmax, n)
      tave = tave + n
      if (n < 0.1*(iw2(ia)/mpi_size)) then
        il = il+1
        tmp(il,1) = it
      elseif (n > 3*iw2(ia)/mpi_size) then
        iu = iu+1
        tmp(iu,2) = it
      endif
    enddo
    empty = (tnul > 0)
    write(iu06,*) '           min, ave, max:      ',tmin,tave/mpi_size,tmax
    if (empty) write(iu06,*) '           No. of empty tasks: ',tnul
    if (il > 0)                                                                &
      write(iu06,*) '           < 10% of even split:',(tmp(jl,1),jl=1,il)
    if (iu > 0)                                                                &
      write(iu06,*) '           > 3x even split:    ',(tmp(ju,2),ju=1,iu)
    write(iu06,*) '  3D, No. wetpoints:           ',iw3(ia)
    write(iu06,*) '           even split:         ',iw3(ia)/mpi_size
    tmin = iw3(ia) + 1
    tave = 0
    tmax = -tmin
    tmp(1:max(il,iu),1:2) = -1
    il = 0
    iu = 0
    do it=1,mpi_size
      i_pe_grid = mod(it-1,nproci)
      j_pe_grid = int((it-1)/nproci,4)
      n = 0
      do k=1,kmx(ia)
        do j=jlows(j_pe_grid+1),jlows(j_pe_grid+1)+j_slices(j_pe_grid+1)-1      
          do i=ilows(i_pe_grid+1),ilows(i_pe_grid+1)+i_slices(i_pe_grid+1)-1
            if (mm1(i,j,k) > 0) n = n + 1
          enddo
        enddo
      enddo
      tmin = min(tmin, n)
      tmax = max(tmax, n)
      tave = tave + n
      if (n < 0.1*(iw3(ia)/mpi_size)) then
        il = il+1
        tmp(il,1) = it
      elseif (n > 3*iw3(ia)/mpi_size) then
        iu = iu+1
        tmp(iu,2) = it
      endif
    enddo
    write(iu06,*) '           min, ave, max:      ',tmin,tave/mpi_size,tmax
    if (il > 0)                                                                &
      write(iu06,*) '           < 10% of even split:',(tmp(jl,1),jl=1,il)
    if (iu > 0)                                                                &
      write(iu06,*) '           > 3x even split:    ',(tmp(ju,2),ju=1,iu)
    write(iu06,*) '------------------------------------------------------------'
 
    ! make the pruned_mpi_decomp_file:
    if (empty) then
      ! shuffle pe_grid(i,j) by skipping empty tasks:
      n = 0
      do it=1,mpi_size
        i_pe_grid = mod(it-1,nproci)
        j_pe_grid = int((it-1)/nproci,4)
        if (pe_grid(i_pe_grid+1,j_pe_grid+1) > 0) then
          n = n + 1
          pe_grid(i_pe_grid+1,j_pe_grid+1) = n
        endif
      enddo
      newmpi = n

      ! write pruned_mpi_decomp_file
      lun = io_new_unit()
      open (unit=lun, file='pruned_'//trim(mpi_decomp_file), status='replace', &
            iostat=ios)
      if (ios /= 0) then
        call exitme(1,'Cannot create pruned_'//trim(mpi_decomp_file))
      endif
      write(lun, '(I5)') newmpi
      do it=1,mpi_size
        ! find non-empty task:
        i_pe_grid = mod(it-1,nproci)
        j_pe_grid = int((it-1)/nproci,4)
        if (pe_grid(i_pe_grid+1,j_pe_grid+1) < 0) cycle

        ! define index range:
        il = ilows(i_pe_grid+1)
        iu = ilows(i_pe_grid+1) + i_slices(i_pe_grid+1) - 1
        jl = jlows(j_pe_grid+1)
        ju = jlows(j_pe_grid+1) + j_slices(j_pe_grid+1) - 1
        write(lun, '(6I5)') pe_grid(i_pe_grid+1,j_pe_grid+1), ia, il, iu, jl, ju

        ! define neighbours:
        if (i_pe_grid < nproci-1) then
          n_s = pe_grid(i_pe_grid+2,j_pe_grid+1)
        else
          n_s = -1
        endif
        if (i_pe_grid > 0) then
          n_n = pe_grid(i_pe_grid,j_pe_grid+1)
        else
          n_n = -1
        endif
        if (j_pe_grid > 0) then
          n_w = pe_grid(i_pe_grid+1,j_pe_grid)
        else
          n_w = -1
        endif
        if (j_pe_grid < nprocj-1) then
          n_e = pe_grid(i_pe_grid+1,j_pe_grid+2)
        else
          n_e = -1
        endif
        write(lun, '(10x,4I5)') n_w, n_n, n_e, n_s

        ! define corner-neighbours:
        if (j_pe_grid > 0 .and. i_pe_grid > 0) then
          n_nw = pe_grid(i_pe_grid,j_pe_grid)
        else
          n_nw = -1
        endif
        if (j_pe_grid < nprocj-1 .and. i_pe_grid > 0) then
          n_ne = pe_grid(i_pe_grid,j_pe_grid+2)
        else
          n_ne = -1
        endif
        if (j_pe_grid < nprocj-1 .and. i_pe_grid < nproci-1) then
          n_se = pe_grid(i_pe_grid+2,j_pe_grid+2)
        else
          n_se = -1
        endif
        if (j_pe_grid > 0 .and. i_pe_grid < nproci-1) then
          n_sw = pe_grid(i_pe_grid+2,j_pe_grid)
        else
          n_sw = -1
        endif
        write(lun, '(10x,4I5)') n_nw, n_ne, n_se, n_sw
      enddo
      close(lun)

      ! print stats for pruned decompo.
      write(iu06,*) 'Statistics for pruned ',nproci,' by ',nprocj,             &
                    ' auto-decomposition'
      write(iu06,*) '  Surface, No. wetpoints:      ',iw2(ia)
      write(iu06,*) '           even split:         ',iw2(ia)/newmpi
      tmin = iw3(ia) + 1
      tave = 0
      tmax = -tmin
      tmp(1:mpi_size,1:2) = -1
      il = 0
      iu = 0
      do it=1,mpi_size
        i_pe_grid = mod(it-1,nproci)
        j_pe_grid = int((it-1)/nproci,4)
        if (pe_grid(i_pe_grid+1,j_pe_grid+1) < 0) cycle
        n = 0
        do j=jlows(j_pe_grid+1),jlows(j_pe_grid+1)+j_slices(j_pe_grid+1)-1      
          do i=ilows(i_pe_grid+1),ilows(i_pe_grid+1)+i_slices(i_pe_grid+1)-1
            if (mm1(i,j,1) > 0) n = n + 1
          enddo
        enddo
        tmin = min(tmin, n)
        tmax = max(tmax, n)
        tave = tave + n
        if (n < 0.1*(iw2(ia)/newmpi)) then
          il = il+1
          tmp(il,1) = pe_grid(i_pe_grid+1,j_pe_grid+1)
        elseif (n > 3*iw2(ia)/newmpi) then
          iu = iu+1
          tmp(iu,2) = pe_grid(i_pe_grid+1,j_pe_grid+1)
        endif
      enddo
      write(iu06,*) '           min, ave, max:      ',tmin,tave/newmpi,tmax
      if (il > 0)                                                              &
        write(iu06,*) '           < 10% of even split:',(tmp(jl,1),jl=1,il)
      if (iu > 0)                                                              &
        write(iu06,*) '           > 3x even split:    ',(tmp(ju,2),ju=1,iu)
      write(iu06,*) '  3D, No. wetpoints:           ',iw3(ia)
      write(iu06,*) '           even split:         ',iw3(ia)/newmpi
      tmin = iw3(ia) + 1
      tave = 0
      tmax = -tmin
      tmp(1:max(il,iu),1:2) = -1
      il = 0
      iu = 0
      do it=1,mpi_size
        i_pe_grid = mod(it-1,nproci)
        j_pe_grid = int((it-1)/nproci,4)
        if (pe_grid(i_pe_grid+1,j_pe_grid+1) < 0) cycle
        n = 0
        do k=1,kmx(ia)
          do j=jlows(j_pe_grid+1),jlows(j_pe_grid+1)+j_slices(j_pe_grid+1)-1
            do i=ilows(i_pe_grid+1),ilows(i_pe_grid+1)+i_slices(i_pe_grid+1)-1
              if (mm1(i,j,k) > 0) n = n + 1
            enddo
          enddo
        enddo
        tmin = min(tmin, n)
        tmax = max(tmax, n)
        tave = tave + n
        if (n < 0.1*(iw3(ia)/newmpi)) then
          il = il+1
          tmp(il,1) = pe_grid(i_pe_grid+1,j_pe_grid+1)
        elseif (n > 3*iw3(ia)/newmpi) then
          iu = iu+1
          tmp(iu,2) = pe_grid(i_pe_grid+1,j_pe_grid+1)
        endif
      enddo
      write(iu06,*) '           min, ave, max:      ',tmin,tave/newmpi,tmax
      if (il > 0)                                                              &
        write(iu06,*) '           < 10% of even split:',(tmp(jl,1),jl=1,il)
      if (iu > 0)                                                              &
        write(iu06,*) '           > 3x even split:    ',(tmp(ju,2),ju=1,iu)
      write(iu06,*)'-----------------------------------------------------------'
    endif
  end subroutine decompose_plane_set_bounds

!===============================================================================

  subroutine dmpi_encode_buf_mmk(iam, kmx, mmk, kh, ia, buf, receiver, ijoff,  &
                                 n, argfactor, a, b, c, d, e, f)
 
    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)           :: kmx, ia, iam, receiver, ijoff
    integer(4), intent(in)           :: mmk(:,0:,0:), n, kh(0:), argfactor
    real(8),    intent(inout)        :: buf(:) !only part of buf is updated here
    real(8),    intent(in)           :: a(0:)
    real(8),    intent(in), optional :: b(0:)
    real(8),    intent(in), optional :: c(0:)
    real(8),    intent(in), optional :: d(0:)
    real(8),    intent(in), optional :: e(0:)
    real(8),    intent(in), optional :: f(0:)

    integer(4) :: wsend, i, j, ij, midx, kidx, kb, mm2, tnum, itr
    integer(4) :: il, iu, jl, ju, isl, isu, jsl, jsu, jompl, jompu, ij0

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    isl = mpi_tt(ia,receiver)%low_hi
    isu = mpi_tt(ia,receiver)%up_hi
    jsl = mpi_tt(ia,receiver)%low_hj
    jsu = mpi_tt(ia,receiver)%up_hj

    call domp_get_domain( max(jsl,jl), min(jsu,ju), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ijoff + max(0,jompl-max(jsl,jl))*max(0,min(isu,iu)-max(isl,il)+1)

    ij    = ij0
    wsend = 0
    do j=jompl,jompu
      do i=max(isl,il),min(isu,iu)
        ij = ij + 1
        if (encode(ia,n)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wsend = wsend + argfactor*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,2) = wsend
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,2) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,2) = offsethalo(itr-1,1,2) + offsethalo(itr-1,2,2)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,2) > 0) then
      ij    = ij0
      wsend = offsethalo(tnum,2,2)
      do j=jompl,jompu
        do i=max(isl,il),min(isu,iu)
          ij = ij + 1
          if (encode(ia,n)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          wsend = wsend + 1
          buf(wsend) = a(midx)
          if (present(b)) then
            wsend = wsend + 1
            buf(wsend) = b(midx)
          endif
          if (present(c)) then
            wsend = wsend + 1
            buf(wsend) = c(midx)
          endif
          if (present(d)) then
            wsend = wsend + 1
            buf(wsend) = d(midx)
          endif
          if (present(e)) then
            wsend = wsend + 1
            buf(wsend) = e(midx)
          endif
          if (present(f)) then
            wsend = wsend + 1
            buf(wsend) = f(midx)
          endif
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            buf(wsend+1:wsend+kb-1) = a(mm2:mm2+kb-2)
            wsend = wsend + kb - 1
            if (present(b)) then
              buf(wsend+1:wsend+kb-1) = b(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(c)) then
              buf(wsend+1:wsend+kb-1) = c(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(d)) then
              buf(wsend+1:wsend+kb-1) = d(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(e)) then
              buf(wsend+1:wsend+kb-1) = e(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(f)) then
              buf(wsend+1:wsend+kb-1) = f(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_encode_buf_mmk

!===============================================================================

  subroutine dmpi_encode_buf_col (iam, kmx, msrf, mcol, kh, ia, buf, receiver, &
                                  ijoff, n, argfactor, a, b, c, d, e, f)
 
    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)           :: kmx, ia, iam, receiver, ijoff, n
    integer(4), intent(in)           :: msrf(0:,0:), mcol(0:), kh(0:), argfactor
    real(8),    intent(inout)        :: buf(:) !only part of buf is updated here
    real(8),    intent(in)           :: a(0:)
    real(8),    intent(in), optional :: b(0:)
    real(8),    intent(in), optional :: c(0:)
    real(8),    intent(in), optional :: d(0:)
    real(8),    intent(in), optional :: e(0:)
    real(8),    intent(in), optional :: f(0:)

    integer(4) :: wsend, i, j, ij, midx, kidx, kb, mm2, tnum, itr
    integer(4) :: il, iu, jl, ju, isl, isu, jsl, jsu, jompl, jompu, ij0

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    isl = mpi_tt(ia,receiver)%low_hi
    isu = mpi_tt(ia,receiver)%up_hi
    jsl = mpi_tt(ia,receiver)%low_hj
    jsu = mpi_tt(ia,receiver)%up_hj

    call domp_get_domain( max(jsl,jl), min(jsu,ju), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ijoff + max(0,jompl-max(jsl,jl))*max(0,min(isu,iu)-max(isl,il)+1)

    ij    = ij0
    wsend = 0
    do j=jompl,jompu
      do i=max(isl,il),min(isu,iu)
        ij = ij + 1
        if (encode(ia,n)%mask(ij) == 0) cycle
        midx = msrf(i,j)
        if (midx == 0) cycle
        wsend = wsend + argfactor*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,2) = wsend
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,2) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,2) = offsethalo(itr-1,1,2) + offsethalo(itr-1,2,2)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,2) > 0) then
      ij    = ij0
      wsend = offsethalo(tnum,2,2)
      do j=jompl,jompu
        do i=max(isl,il),min(isu,iu)
          ij = ij + 1
          if (encode(ia,n)%mask(ij) == 0) cycle
          midx = msrf(i,j)
          if (midx == 0) cycle
          wsend = wsend + 1
          buf(wsend) = a(midx)
          if (present(b)) then
            wsend = wsend + 1
            buf(wsend) = b(midx)
          endif
          if (present(c)) then
            wsend = wsend + 1
            buf(wsend) = c(midx)
          endif
          if (present(d)) then
            wsend = wsend + 1
            buf(wsend) = d(midx)
          endif
          if (present(e)) then
            wsend = wsend + 1
            buf(wsend) = e(midx)
          endif
          if (present(f)) then
            wsend = wsend + 1
            buf(wsend) = f(midx)
          endif
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mcol(midx)
            buf(wsend+1:wsend+kb-1) = a(mm2:mm2+kb-2)
            wsend = wsend + kb - 1
            if (present(b)) then
              buf(wsend+1:wsend+kb-1) = b(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(c)) then
              buf(wsend+1:wsend+kb-1) = c(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(d)) then
              buf(wsend+1:wsend+kb-1) = d(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(e)) then
              buf(wsend+1:wsend+kb-1) = e(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
            if (present(f)) then
              buf(wsend+1:wsend+kb-1) = f(mm2:mm2+kb-2)
              wsend = wsend + kb - 1
            endif
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_encode_buf_col

!===============================================================================

  subroutine dmpi_encode_buf_nc1(iam, kmx, mmk, kh,ia,a,buf,receiver,ijoff,n,nc)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)  :: kmx, mmk(:,0:,0:), kh(0:), ia, iam, receiver
    integer(4), intent(in)  :: ijoff, n, nc
    real(8),    intent(in)  :: a(:,0:)
    real(8),    intent(inout) :: buf(:) ! only part of buf is updated here

    integer(4) :: wtmp, i, j, kidx, ij, midx, kb, mm2, tnum, itr 
    integer(4) :: il, iu, jl, ju, isl, isu, jsl, jsu, jompl, jompu, ij0

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    isl = mpi_tt(ia,receiver)%low_hi
    isu = mpi_tt(ia,receiver)%up_hi
    jsl = mpi_tt(ia,receiver)%low_hj
    jsu = mpi_tt(ia,receiver)%up_hj

    call domp_get_domain( max(jsl,jl), min(jsu,ju), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ijoff + max(0,jompl-max(jsl,jl))*max(0,min(isu,iu)-max(isl,il)+1)

    ij   = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(isl,il),min(isu,iu)
        ij = ij + 1
        if (encode(ia,n)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wtmp = wtmp + nc*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,2) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,2) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,2) = offsethalo(itr-1,1,2) + offsethalo(itr-1,2,2)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,2) > 0) then
      ij   = ij0
      wtmp = offsethalo(tnum,2,2)
      do j=jompl,jompu
        do i=max(isl,il),min(isu,iu)
          ij = ij + 1
          if (encode(ia,n)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          buf(wtmp+1:wtmp+nc) = a(1:nc,midx)
          wtmp = wtmp + nc
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            do kidx=mm2,mm2+kb-2
              buf(wtmp+1:wtmp+nc) = a(1:nc,kidx)
              wtmp = wtmp + nc
            enddo
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_encode_buf_nc1

!===============================================================================

  subroutine dmpi_encode_buf_nc2(iam,kmx,mmk,kh,ia,a,b,buf,receiver,ijoff,n,nc)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)  :: kmx, mmk(:,0:,0:), kh(0:), ia, iam, receiver
    integer(4), intent(in)  :: ijoff, n, nc
    real(8),    intent(in)  :: a(:,0:), b(:,0:)
    real(8),    intent(inout) :: buf(:) ! only part of buf is updated here

    integer(4) :: wtmp, i, j, kidx, ij, midx, kb, mm2, tnum, itr, bfac
    integer(4) :: il, iu, jl, ju, isl, isu, jsl, jsu, jompl, jompu, ij0

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    isl = mpi_tt(ia,receiver)%low_hi
    isu = mpi_tt(ia,receiver)%up_hi
    jsl = mpi_tt(ia,receiver)%low_hj
    jsu = mpi_tt(ia,receiver)%up_hj

    bfac = 2*nc

    call domp_get_domain( max(jsl,jl), min(jsu,ju), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ijoff + max(0,jompl-max(jsl,jl))*max(0,min(isu,iu)-max(isl,il)+1)

    ij   = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(isl,il),min(isu,iu)
        ij = ij + 1
        if (encode(ia,n)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wtmp = wtmp + bfac*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,2) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,2) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,2) = offsethalo(itr-1,1,2) + offsethalo(itr-1,2,2)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,2) > 0) then
      ij   = ij0
      wtmp = offsethalo(tnum,2,2)
      do j=jompl,jompu
        do i=max(isl,il),min(isu,iu)
          ij = ij + 1
          if (encode(ia,n)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          buf(wtmp   +1:wtmp  +nc) = a(1:nc,midx)
          buf(wtmp+nc+1:wtmp+2*nc) = b(1:nc,midx)
          wtmp = wtmp + bfac
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            do kidx=mm2,mm2+kb-2
              buf(wtmp   +1:wtmp  +nc) = a(1:nc,kidx)
              buf(wtmp+nc+1:wtmp+2*nc) = b(1:nc,kidx)
              wtmp = wtmp + bfac
            enddo
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_encode_buf_nc2

!===============================================================================

  subroutine dmpi_encode_buf_nc3(iam,kmx,mmk,kh,ia,a,b,c,buf,receiver,ijoff,n, &
                                 nc)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)  :: kmx, mmk(:,0:,0:), kh(0:), ia, iam, receiver
    integer(4), intent(in)  :: ijoff, n, nc
    real(8),    intent(in)  :: a(:,0:), b(:,0:), c(:,0:)
    real(8),    intent(inout) :: buf(:) ! only part of buf is updated here

    integer(4) :: wtmp, i, j, kidx, ij, midx, kb, mm2, tnum, itr, bfac
    integer(4) :: il, iu, jl, ju, isl, isu, jsl, jsu, jompl, jompu, ij0

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    isl = mpi_tt(ia,receiver)%low_hi
    isu = mpi_tt(ia,receiver)%up_hi
    jsl = mpi_tt(ia,receiver)%low_hj
    jsu = mpi_tt(ia,receiver)%up_hj

    bfac = 3*nc

    call domp_get_domain( max(jsl,jl), min(jsu,ju), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ijoff + max(0,jompl-max(jsl,jl))*max(0,min(isu,iu)-max(isl,il)+1)

    ij   = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(isl,il),min(isu,iu)
        ij = ij + 1
        if (encode(ia,n)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wtmp = wtmp + bfac*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,2) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,2) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,2) = offsethalo(itr-1,1,2) + offsethalo(itr-1,2,2)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,2) > 0) then
      ij   = ij0
      wtmp = offsethalo(tnum,2,2)
      do j=jompl,jompu
        do i=max(isl,il),min(isu,iu)
          ij = ij + 1
          if (encode(ia,n)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          buf(wtmp     +1:wtmp  +nc) = a(1:nc,midx)
          buf(wtmp  +nc+1:wtmp+2*nc) = b(1:nc,midx)
          buf(wtmp+2*nc+1:wtmp+3*nc) = c(1:nc,midx)
          wtmp = wtmp + bfac
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            do kidx=mm2,mm2+kb-2
              buf(wtmp     +1:wtmp  +nc) = a(1:nc,kidx)
              buf(wtmp  +nc+1:wtmp+2*nc) = b(1:nc,kidx)
              buf(wtmp+2*nc+1:wtmp+3*nc) = c(1:nc,kidx)
              wtmp = wtmp + bfac
            enddo
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_encode_buf_nc3

!===============================================================================

  subroutine dmpi_encode_buf_log(iam, mmk, ia, a, buf, receiver, ijoff, n)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: mmk(:,0:,0:), ia, iam, receiver, ijoff, n
    logical,    intent(in)    :: a(0:)
    logical,    intent(inout) :: buf(:) ! only part of buf is updated here

    integer(4) :: wtmp, i, j, ij, tnum, jompl, jompu, ij0, itr
    integer(4) :: il, iu, jl, ju, isl, isu, jsl, jsu

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    isl = mpi_tt(ia,receiver)%low_hi
    isu = mpi_tt(ia,receiver)%up_hi
    jsl = mpi_tt(ia,receiver)%low_hj
    jsu = mpi_tt(ia,receiver)%up_hj

    call domp_get_domain( max(jsl,jl), min(jsu,ju), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ijoff + max(0,jompl-max(jsl,jl))*max(0,min(isu,iu)-max(isl,il)+1)

    ij   = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(isl,il),min(isu,iu)
        ij = ij + 1
        if (encode(ia,n)%mask(ij) == 0) cycle
        wtmp = wtmp + 1
      enddo
    enddo
    offsethalo(tnum,1,2) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,2) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,2) = offsethalo(itr-1,1,2) + offsethalo(itr-1,2,2)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,2) > 0) then
      ij   = ij0
      wtmp = offsethalo(tnum,2,2)
      do j=jompl,jompu
        do i=max(isl,il),min(isu,iu)
          ij = ij + 1
          if (encode(ia,n)%mask(ij) == 0) cycle
          wtmp = wtmp + 1
          buf(wtmp) = a(mmk(1,i,j))
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_encode_buf_log

!===============================================================================

  subroutine dmpi_decode_buf_mmk(iam, id, kmx, mmk, kh, ia, buf, argfactor,    &
                                 a, b, c, d, e, f)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)              :: kmx, ia, id, iam, argfactor
    integer(4), intent(in)              :: mmk(:,0:,0:), kh(0:)
    real(8),    intent(in)              :: buf(:)
    ! only part of halo is updated here
    real(8),    intent(inout)           :: a(0:)
    real(8),    intent(inout), optional :: b(0:)
    real(8),    intent(inout), optional :: c(0:)
    real(8),    intent(inout), optional :: d(0:)
    real(8),    intent(inout), optional :: e(0:)
    real(8),    intent(inout), optional :: f(0:)

    integer(4) :: wrecv, i, j, ij, nb, midx, kidx, mm2, kb, tnum, itr, ij0
    integer(4) :: irl, iru, jrl, jru, jompl, jompu
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    nb  = mpi_tt(ia,iam)%hnb(id)
    il  = mpi_tt(ia,nb)%low_i
    iu  = mpi_tt(ia,nb)%up_i
    jl  = mpi_tt(ia,nb)%low_j
    ju  = mpi_tt(ia,nb)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    ! determine zone borders:
    call dmpi_soff(id,il, iu, jl, ju, ihl, ihu, jhl, jhu, irl, iru, jrl,jru,ij0)

    call domp_get_domain( max(jrl,jhl), min(jru,jhu), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ij0 + max(0,jompl-max(jrl,jhl))*max(0,min(iru,ihu)-max(irl,ihl)+1)

    ij    = ij0
    wrecv = 0
    do j=jompl,jompu
      do i=max(irl,ihl),min(iru,ihu)
        ij = ij + 1
        if (decode(ia,id)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wrecv = wrecv + argfactor*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,1) = wrecv
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,1) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,1) = offsethalo(itr-1,1,1) + offsethalo(itr-1,2,1)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,1) > 0) then
      ij    = ij0 
      wrecv = offsethalo(tnum,2,1)
      do j=jompl,jompu
        do i=max(irl,ihl),min(iru,ihu)
          ij = ij + 1
          if (decode(ia,id)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          wrecv = wrecv + 1
          a(midx) = buf(wrecv)
          if (present(b)) then
            wrecv = wrecv + 1
            b(midx) = buf(wrecv)
          endif
          if (present(c)) then
            wrecv = wrecv + 1
            c(midx) = buf(wrecv)
          endif
          if (present(d)) then
            wrecv = wrecv + 1
            d(midx) = buf(wrecv)
          endif
          if (present(e)) then
            wrecv = wrecv + 1
            e(midx) = buf(wrecv)
          endif
          if (present(f)) then
            wrecv = wrecv + 1
            f(midx) = buf(wrecv)
          endif
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            a(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
            wrecv = wrecv + kb - 1
            if (present(b)) then
              b(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(c)) then
              c(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(d)) then
              d(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(e)) then
              e(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(f)) then
              f(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_decode_buf_mmk

!===============================================================================

  subroutine dmpi_decode_buf_col(iam, id, kmx, msrf, mcol, kh, ia, buf,        &
                                 argfactor, a, b, c, d, e, f)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)              :: kmx, ia, id, iam, argfactor
    integer(4), intent(in)              :: msrf(0:,0:), mcol(0:), kh(0:)
    real(8),    intent(in)              :: buf(:)
    ! only part of halo is updated here
    real(8),    intent(inout)           :: a(0:)
    real(8),    intent(inout), optional :: b(0:)
    real(8),    intent(inout), optional :: c(0:)
    real(8),    intent(inout), optional :: d(0:)
    real(8),    intent(inout), optional :: e(0:)
    real(8),    intent(inout), optional :: f(0:)

    integer(4) :: wrecv, i, j, ij, nb, midx, kidx, mm2, kb, tnum, itr, ij0
    integer(4) :: irl, iru, jrl, jru, jompl, jompu
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    nb  = mpi_tt(ia,iam)%hnb(id)
    il  = mpi_tt(ia,nb)%low_i
    iu  = mpi_tt(ia,nb)%up_i
    jl  = mpi_tt(ia,nb)%low_j
    ju  = mpi_tt(ia,nb)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    ! determine zone borders:
    call dmpi_soff(id,il, iu, jl, ju, ihl, ihu, jhl, jhu, irl, iru, jrl,jru,ij0)

    call domp_get_domain( max(jrl,jhl), min(jru,jhu), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ij0 + max(0,jompl-max(jrl,jhl))*max(0,min(iru,ihu)-max(irl,ihl)+1)

    ij    = ij0
    wrecv = 0
    do j=jompl,jompu
      do i=max(irl,ihl),min(iru,ihu)
        ij = ij + 1
        if (decode(ia,id)%mask(ij) == 0) cycle
        midx = msrf(i,j)
        if (midx == 0) cycle
        wrecv = wrecv + argfactor*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,1) = wrecv
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,1) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,1) = offsethalo(itr-1,1,1) + offsethalo(itr-1,2,1)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,1) > 0) then
      ij    = ij0 
      wrecv = offsethalo(tnum,2,1)
      do j=jompl,jompu
        do i=max(irl,ihl),min(iru,ihu)
          ij = ij + 1
          if (decode(ia,id)%mask(ij) == 0) cycle
          midx = msrf(i,j)
          if (midx == 0) cycle
          wrecv = wrecv + 1
          a(midx) = buf(wrecv)
          if (present(b)) then
            wrecv = wrecv + 1
            b(midx) = buf(wrecv)
          endif
          if (present(c)) then
            wrecv = wrecv + 1
            c(midx) = buf(wrecv)
          endif
          if (present(d)) then
            wrecv = wrecv + 1
            d(midx) = buf(wrecv)
          endif
          if (present(e)) then
            wrecv = wrecv + 1
            e(midx) = buf(wrecv)
          endif
          if (present(f)) then
            wrecv = wrecv + 1
            f(midx) = buf(wrecv)
          endif
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mcol(midx)
            a(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
            wrecv = wrecv + kb - 1
            if (present(b)) then
              b(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(c)) then
              c(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(d)) then
              d(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(e)) then
              e(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
            if (present(f)) then
              f(mm2:mm2+kb-2) = buf(wrecv+1:wrecv+kb-1)
              wrecv = wrecv + kb - 1
            endif
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_decode_buf_col

!===============================================================================

  subroutine dmpi_decode_buf_nc1(iam, id, kmx, mmk, kh, ia, a, buf, nc)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: kmx, mmk(:,0:,0:), ia, id, iam, nc, kh(0:)
    real(8),    intent(in)    :: buf(:)
    real(8),    intent(inout) :: a(:,0:)   ! only part of halo is updated here

    integer(4) :: wtmp, i, j, kidx, ij, nb, midx, kb, mm2, tnum, itr
    integer(4) :: irl, iru, jrl, jru
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu, jompl, jompu, joff, ij0

    nb  = mpi_tt(ia,iam)%hnb(id)
    il  = mpi_tt(ia,nb)%low_i
    iu  = mpi_tt(ia,nb)%up_i
    jl  = mpi_tt(ia,nb)%low_j
    ju  = mpi_tt(ia,nb)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    ! determine zone borders:
    call dmpi_soff(id,il, iu, jl, ju, ihl, ihu, jhl, jhu, irl, iru, jrl,jru,ij0)

    call domp_get_domain( max(jrl,jhl), min(jru,jhu), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ij0 + max(0,jompl-max(jrl,jhl))*max(0,min(iru,ihu)-max(irl,ihl)+1)
 
    ij = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(irl,ihl),min(iru,ihu)
        ij = ij + 1
        if (decode(ia,id)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wtmp = wtmp + nc*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,1) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,1) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,1) = offsethalo(itr-1,1,1) + offsethalo(itr-1,2,1)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,1) > 0) then
      ij = ij0 
      wtmp = offsethalo(tnum,2,1)
      do j=jompl,jompu
        do i=max(irl,ihl),min(iru,ihu)
          ij = ij + 1
          if (decode(ia,id)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          a(1:nc,midx) = buf(wtmp+1:wtmp+nc)
          wtmp = wtmp + nc
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            do kidx=mm2,mm2+kb-2
              a(1:nc,kidx) = buf(wtmp+1:wtmp+nc)
              wtmp = wtmp + nc
            enddo
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_decode_buf_nc1

!===============================================================================

  subroutine dmpi_decode_buf_nc2(iam, id, kmx, mmk, kh, ia, a, b, buf, nc)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: kmx, mmk(:,0:,0:), ia, id, iam, nc, kh(0:)
    real(8),    intent(in)    :: buf(:)
    real(8),    intent(inout) :: a(:,0:), b(:,0:)  ! only part of halo updated

    integer(4) :: wtmp, i, j, kidx, ij, nb, midx, kb, mm2, tnum, itr, bfac
    integer(4) :: irl, iru, jrl, jru
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu, jompl, jompu, joff, ij0

    nb  = mpi_tt(ia,iam)%hnb(id)
    il  = mpi_tt(ia,nb)%low_i
    iu  = mpi_tt(ia,nb)%up_i
    jl  = mpi_tt(ia,nb)%low_j
    ju  = mpi_tt(ia,nb)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    bfac = 2*nc

    ! determine zone borders:
    call dmpi_soff(id,il, iu, jl, ju, ihl, ihu, jhl, jhu, irl, iru, jrl,jru,ij0)

    call domp_get_domain( max(jrl,jhl), min(jru,jhu), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ij0 + max(0,jompl-max(jrl,jhl))*max(0,min(iru,ihu)-max(irl,ihl)+1)
 
    ij = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(irl,ihl),min(iru,ihu)
        ij = ij + 1
        if (decode(ia,id)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wtmp = wtmp + bfac*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,1) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,1) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,1) = offsethalo(itr-1,1,1) + offsethalo(itr-1,2,1)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,1) > 0) then
      ij = ij0 
      wtmp = offsethalo(tnum,2,1)
      do j=jompl,jompu
        do i=max(irl,ihl),min(iru,ihu)
          ij = ij + 1
          if (decode(ia,id)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          a(1:nc,midx) = buf(wtmp   +1:wtmp  +nc)
          b(1:nc,midx) = buf(wtmp+nc+1:wtmp+2*nc)
          wtmp = wtmp + bfac
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            do kidx=mm2,mm2+kb-2
              a(1:nc,kidx) = buf(wtmp   +1:wtmp  +nc)
              b(1:nc,kidx) = buf(wtmp+nc+1:wtmp+2*nc)
              wtmp = wtmp + bfac
            enddo
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_decode_buf_nc2

!===============================================================================

  subroutine dmpi_decode_buf_nc3(iam, id, kmx, mmk, kh, ia, a, b, c, buf, nc)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: kmx, mmk(:,0:,0:), ia, id, iam, nc, kh(0:)
    real(8),    intent(in)    :: buf(:)
    real(8),    intent(inout) :: a(:,0:), b(:,0:), c(:,0:) 

    integer(4) :: wtmp, i, j, kidx, ij, nb, midx, kb, mm2, tnum, itr, bfac
    integer(4) :: irl, iru, jrl, jru
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu, jompl, jompu, joff, ij0

    nb  = mpi_tt(ia,iam)%hnb(id)
    il  = mpi_tt(ia,nb)%low_i
    iu  = mpi_tt(ia,nb)%up_i
    jl  = mpi_tt(ia,nb)%low_j
    ju  = mpi_tt(ia,nb)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    bfac = 3*nc

    ! determine zone borders:
    call dmpi_soff(id,il, iu, jl, ju, ihl, ihu, jhl, jhu, irl, iru, jrl,jru,ij0)

    call domp_get_domain( max(jrl,jhl), min(jru,jhu), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ij0 + max(0,jompl-max(jrl,jhl))*max(0,min(iru,ihu)-max(irl,ihl)+1)
 
    ij = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(irl,ihl),min(iru,ihu)
        ij = ij + 1
        if (decode(ia,id)%mask(ij) == 0) cycle
        midx = mmk(1,i,j)
        if (midx == 0) cycle
        wtmp = wtmp + bfac*min(kmx,kh(midx))
      enddo
    enddo
    offsethalo(tnum,1,1) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,1) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,1) = offsethalo(itr-1,1,1) + offsethalo(itr-1,2,1)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,1) > 0) then
      ij = ij0 
      wtmp = offsethalo(tnum,2,1)
      do j=jompl,jompu
        do i=max(irl,ihl),min(iru,ihu)
          ij = ij + 1
          if (decode(ia,id)%mask(ij) == 0) cycle
          midx = mmk(1,i,j)
          if (midx == 0) cycle
          a(1:nc,midx) = buf(wtmp     +1:wtmp  +nc)
          b(1:nc,midx) = buf(wtmp  +nc+1:wtmp+2*nc)
          c(1:nc,midx) = buf(wtmp+2*nc+1:wtmp+3*nc)
          wtmp = wtmp + bfac
          kb = min(kmx,kh(midx))
          if (kb >= 2) then
            mm2 = mmk(2,i,j)
            do kidx=mm2,mm2+kb-2
              a(1:nc,kidx) = buf(wtmp     +1:wtmp  +nc)
              b(1:nc,kidx) = buf(wtmp  +nc+1:wtmp+2*nc)
              c(1:nc,kidx) = buf(wtmp+2*nc+1:wtmp+3*nc)
              wtmp = wtmp + bfac
            enddo
          endif
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_decode_buf_nc3

!===============================================================================

  subroutine dmpi_decode_buf_log(iam, id, mmk, ia, a, buf)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: mmk(:,0:,0:), ia, id, iam
    logical,    intent(in)    :: buf(:)
    logical,    intent(inout) :: a(0:)   ! only part of halo is updated here

    integer(4) :: wtmp, i, j, ij, nb, jompl, jompu, tnum, itr, ij0
    integer(4) :: irl, iru, jrl, jru
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    nb  = mpi_tt(ia,iam)%hnb(id)
    il  = mpi_tt(ia,nb)%low_i
    iu  = mpi_tt(ia,nb)%up_i
    jl  = mpi_tt(ia,nb)%low_j
    ju  = mpi_tt(ia,nb)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    ! determine zone borders:
    call dmpi_soff(id,il, iu, jl, ju, ihl, ihu, jhl, jhu, irl, iru, jrl,jru,ij0)

    call domp_get_domain( max(jrl,jhl), min(jru,jhu), jompl, jompu )
    call domp_get_thread_no( tnum )

    ij0 = ij0 + max(0,jompl-max(jrl,jhl))*max(0,min(iru,ihu)-max(irl,ihl)+1)
 
    ij = ij0
    wtmp = 0
    do j=jompl,jompu
      do i=max(irl,ihl),min(iru,ihu)
        ij = ij + 1
        if (decode(ia,id)%mask(ij) == 0) cycle
        wtmp = wtmp + 1
      enddo
    enddo
    offsethalo(tnum,1,1) = wtmp
!$OMP BARRIER
!$OMP MASTER
    offsethalo(1,2,1) = 0
    do itr=2,mpi_nt
      offsethalo(itr,2,1) = offsethalo(itr-1,1,1) + offsethalo(itr-1,2,1)
    enddo
!$OMP END MASTER
!$OMP BARRIER
    if (offsethalo(tnum,1,1) > 0) then
      ij = ij0 
      wtmp = offsethalo(tnum,2,1)
      do j=jompl,jompu
        do i=max(irl,ihl),min(iru,ihu)
          ij = ij + 1
          if (decode(ia,id)%mask(ij) == 0) cycle
          wtmp = wtmp + 1
          a(mmk(1,i,j)) = buf(wtmp)
        enddo
      enddo
    endif

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_decode_buf_log

!===============================================================================

  subroutine dmpi_init_comm_masks (na)
    implicit none

    integer(4), intent(in) :: na

    integer(4)               :: iam, ia, il, iu, jl, ju, nb, id
    integer(4), dimension(4) :: iln, iun, jln, jun

    iam = mpi_rank + 1

    allocate( encode(1:na,1:max_nb), decode(1:na,1:max_nb) )

    do ia=1,na

      !- Allocate encode masks -------------------------------------------------
      !  Task's coord.extremes:
      il = mpi_tt(ia,iam)%low_i
      iu = mpi_tt(ia,iam)%up_i
      jl = mpi_tt(ia,iam)%low_j
      ju = mpi_tt(ia,iam)%up_j
      !  Neighbour's halo-coord.extremes:
      do id=1,4
        nb = mpi_tt(ia,iam)%hnb(id)
        if (nb /= -1) then
          iln(id) = mpi_tt(ia,nb)%low_hi
          iun(id) = mpi_tt(ia,nb)%up_hi
          jln(id) = mpi_tt(ia,nb)%low_hj
          jun(id) = mpi_tt(ia,nb)%up_hj
        else
          iln(id) = il
          iun(id) = iu
          jln(id) = jl
          jun(id) = ju
        endif
      enddo

      allocate( encode(ia,1)%mask(max(il,iln(1)):min(iu,iun(1))),              &
                encode(ia,2)%mask(max(jl,jln(2)):min(ju,jun(2))),              &
                encode(ia,3)%mask(max(il,iln(3)):min(iu,iun(3))),              &
                encode(ia,4)%mask(max(jl,jln(4)):min(ju,jun(4))),              &
                encode(ia,5)%mask(1), encode(ia,6)%mask(1),                    &
                encode(ia,7)%mask(1), encode(ia,8)%mask(1)       )


      !- Allocate encode masks -------------------------------------------------
      !  Task's halo-coord.extremes:
      il = mpi_tt(ia,iam)%low_hi
      iu = mpi_tt(ia,iam)%up_hi
      jl = mpi_tt(ia,iam)%low_hj
      ju = mpi_tt(ia,iam)%up_hj
      !  Neighbour's coord.extremes:
      do id=1,4
        nb = mpi_tt(ia,iam)%hnb(id)
        if (nb /= -1) then
          iln(id) = mpi_tt(ia,nb)%low_i
          iun(id) = mpi_tt(ia,nb)%up_i
          jln(id) = mpi_tt(ia,nb)%low_j
          jun(id) = mpi_tt(ia,nb)%up_j
        else
          iln(id) = il
          iun(id) = iu
          jln(id) = jl
          jun(id) = ju
        endif
      enddo

      allocate( decode(ia,1)%mask(max(il,iln(1)):min(iu,iun(1))),              &
                decode(ia,2)%mask(max(jl,jln(2)):min(ju,jun(2))),              &
                decode(ia,3)%mask(max(il,iln(3)):min(iu,iun(3))),              &
                decode(ia,4)%mask(max(jl,jln(4)):min(ju,jun(4))),              &
                decode(ia,5)%mask(1), decode(ia,6)%mask(1),                    &
                decode(ia,7)%mask(1), decode(ia,8)%mask(1)       )
    enddo

    !- Initialise encode/code masks:
    do id=1,max_nb
      do ia=1,na
        encode(ia,id)%mask(:) = 0
        decode(ia,id)%mask(:) = 0
      enddo
    enddo
  end subroutine dmpi_init_comm_masks

!===============================================================================
  subroutine dmpi_set_comm_masks_msrf (ia, mm)
    implicit none

    integer(4), intent(in) :: ia, mm(0:,0:)

    integer(4) :: iam, it_recv, itp1, id, nb, idi, nbi
    integer(4) :: i, j, ij, ii, jj, ip, jp, im, jm
    integer(4) :: irl, iru, jrl, jru, isl, isu, jsl, jsu
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    iam = mpi_rank + 1

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    do it_recv=0,mpi_size-1
      itp1 = it_recv + 1

      ! I should not send to myself, but receive from up to 8 of my neighbours:
      if (itp1 == iam) then
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (nb > 0) then
            ! recv-buffer coordinates:
            select case (id)
              case (1) ! west
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhl
                jru = jhl
                ij  = max(irl,ihl) - 1
                ii  =  0
                jj  = +1
                ip  = +1
                jp  =  0
                im  = -1
                jm  =  0
              case (2) ! north
                irl = ihl
                iru = ihl
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
                ii  = +1
                jj  =  0
                ip  =  0
                jp  = +1
                im  =  0
                jm  = -1
              case (3) ! east
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhu
                jru = jhu
                ij  = max(irl,ihl) - 1
                ii  =  0
                jj  = -1
                ip  = +1
                jp  =  0
                im  = -1
                jm  =  0
              case (4) ! south
                irl = ihu
                iru = ihu
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
                ii  = -1
                jj  =  0
                ip  =  0
                jp  = +1
                im  =  0
                jm  = -1
              case (5) ! north-west
                irl = ihl
                iru = ihl
                jrl = jhl
                jru = jhl
                ij  =  0
                ii  = +1 
                jj  = +1
                ip  =  0
                jp  = +1
                im  = +1
                jm  =  0
              case (6) ! north-east
                irl = ihl
                iru = ihl
                jrl = jhu
                jru = jhu
                ij  =  0
                ii  = +1 
                jj  = -1
                ip  =  0
                jp  = -1
                im  = +1
                jm  =  0
              case (7) ! south-east
                irl = ihu
                iru = ihu
                jrl = jhu
                jru = jhu
                ij  =  0
                ii  = -1 
                jj  = -1
                ip  =  0
                jp  = -1
                im  = -1
                jm  =  0
              case (8) ! south-west
                irl = ihu
                iru = ihu
                jrl = jhl
                jru = jhl
                ij  =  0
                ii  = -1 
                jj  = +1
                ip  =  0
                jp  = +1
                im  = -1
                jm  =  0
            end select

            ! assign the decode mask values:
            do j=max(jrl,jhl),min(jru,jhu)
              do i=max(irl,ihl),min(iru,ihu)
                ij = ij + 1
                if (mm(i,j) == 0) cycle
                ! id is neighbour seen from receiver iam.
                if (id <= 4) then
                  if (mm(i+ii,j+jj) == 0) then
                    if ( (mm(i+im,   j+jm   ) == 0  .or.                     &
                          mm(i+im+ii,j+jm+jj) == 0) .and.                    &
                         (mm(i+ip,   j+jp   ) == 0  .or.                     &
                          mm(i+ip+ii,j+jp+jj) == 0)       )  cycle
                  endif
                else
                  if (mm(i+ii,j+jj) == 0) cycle
                  if (mm(i+im,j+jm) == 0 .and. mm(i+ip,j+jp) == 0) cycle
                endif
                decode(ia,id)%mask(ij) = ij
              enddo
            enddo
          endif
        enddo  ! id

      ! I should not receive from myself, but send to up to 8 of my neighbours:
      else

        do id=1,max_nb
          nb = mpi_tt(ia,itp1)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #itp1 and I must send to him

            ! Find my neighbour which corresponds to task #itp1
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == itp1) then
              ! send-buffer coordinates:
              select case (idi)
                case (1) ! west
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = jl
                  jsu = jl
                  ij  = max(isl,il) - 1
                  ii  =  0
                  jj  = -1
                  ip  = +1
                  jp  =  0
                  im  = -1
                  jm  =  0
                case (2) ! north
                  isl = il
                  isu = il
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                  ii  = -1
                  jj  =  0
                  ip  =  0
                  jp  = +1
                  im  =  0
                  jm  = -1
                case (3) ! east
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = ju
                  jsu = ju
                  ij  = max(isl,il) - 1
                  ii  =  0
                  jj  = +1
                  ip  = +1
                  jp  =  0
                  im  = -1
                  jm  =  0
                case (4) ! south
                  isl = iu
                  isu = iu
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                  ii  = +1
                  jj  =  0
                  ip  =  0
                  jp  = +1
                  im  =  0
                  jm  = -1
                case (5) ! north-west
                  isl = il
                  isu = il
                  jsl = jl
                  jsu = jl
                  ij  =  0
                  ii  = -1 
                  jj  = -1
                  ip  = -1
                  jp  =  0
                  im  =  0
                  jm  = -1
                case (6) ! north-east
                  isl = il
                  isu = il
                  jsl = ju
                  jsu = ju
                  ij  =  0
                  ii  = -1 
                  jj  = +1
                  ip  =  0
                  jp  = +1
                  im  = -1
                  jm  =  0
                case (7) ! south-east
                  isl = iu
                  isu = iu
                  jsl = ju
                  jsu = ju
                  ij  =  0
                  ii  = +1 
                  jj  = +1
                  ip  =  0
                  jp  = +1
                  im  = +1
                  jm  =  0
                case (8) ! south-west
                  isl = iu
                  isu = iu
                  jsl = jl
                  jsu = jl
                  ij  =  0
                  ii  = +1 
                  jj  = -1
                  ip  = +1
                  jp  =  0
                  im  =  0
                  jm  = -1
              end select

              ! assign the encode mask values:
              do j=max(jsl,jl),min(jsu,ju)
                do i=max(isl,il),min(isu,iu)
                  ij = ij + 1
                  if (mm(i,j) == 0) cycle

                  ! id  is neighbour seen from receiver itp1.
                  ! idi is neighbour seen from sender iam.
                  if (idi <= 4) then
                    if (mm(i+ii,j+jj) == 0) then
                      if ( (mm(i+im,   j+jm   ) == 0  .or.                   &
                            mm(i+im+ii,j+jm+jj) == 0) .and.                  &
                           (mm(i+ip,   j+jp   ) == 0  .or.                   &
                            mm(i+ip+ii,j+jp+jj) == 0)       )  cycle
                    endif
                  else
                    if (mm(i+ii,j+jj) == 0) cycle
                    if (mm(i+im,j+jm) == 0 .and. mm(i+ip,j+jp) == 0) cycle
                  endif
                  encode(ia,idi)%mask(ij) = ij
                enddo
              enddo
            endif  ! nbi == itp1
          endif    ! nb == iam
        enddo      ! id

      endif
    enddo   ! it_recv, itp1


  end subroutine dmpi_set_comm_masks_msrf

  subroutine dmpi_set_comm_masks (ia, mm)
    implicit none

    integer(4), intent(in) :: ia, mm(0:,0:,:)

    integer(4) :: iam, it_recv, itp1, id, nb, idi, nbi
    integer(4) :: i, j, ij, ii, jj, ip, jp, im, jm
    integer(4) :: irl, iru, jrl, jru, isl, isu, jsl, jsu
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    iam = mpi_rank + 1

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    do it_recv=0,mpi_size-1
      itp1 = it_recv + 1

      ! I should not send to myself, but receive from up to 8 of my neighbours:
      if (itp1 == iam) then
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (nb > 0) then
            ! recv-buffer coordinates:
            select case (id)
              case (1) ! west
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhl
                jru = jhl
                ij  = max(irl,ihl) - 1
                ii  =  0
                jj  = +1
                ip  = +1
                jp  =  0
                im  = -1
                jm  =  0
              case (2) ! north
                irl = ihl
                iru = ihl
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
                ii  = +1
                jj  =  0
                ip  =  0
                jp  = +1
                im  =  0
                jm  = -1
              case (3) ! east
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhu
                jru = jhu
                ij  = max(irl,ihl) - 1
                ii  =  0
                jj  = -1
                ip  = +1
                jp  =  0
                im  = -1
                jm  =  0
              case (4) ! south
                irl = ihu
                iru = ihu
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
                ii  = -1
                jj  =  0
                ip  =  0
                jp  = +1
                im  =  0
                jm  = -1
              case (5) ! north-west
                irl = ihl
                iru = ihl
                jrl = jhl
                jru = jhl
                ij  =  0
                ii  = +1 
                jj  = +1
                ip  =  0
                jp  = +1
                im  = +1
                jm  =  0
              case (6) ! north-east
                irl = ihl
                iru = ihl
                jrl = jhu
                jru = jhu
                ij  =  0
                ii  = +1 
                jj  = -1
                ip  =  0
                jp  = -1
                im  = +1
                jm  =  0
              case (7) ! south-east
                irl = ihu
                iru = ihu
                jrl = jhu
                jru = jhu
                ij  =  0
                ii  = -1 
                jj  = -1
                ip  =  0
                jp  = -1
                im  = -1
                jm  =  0
              case (8) ! south-west
                irl = ihu
                iru = ihu
                jrl = jhl
                jru = jhl
                ij  =  0
                ii  = -1 
                jj  = +1
                ip  =  0
                jp  = +1
                im  = -1
                jm  =  0
            end select

            ! assign the decode mask values:
            do j=max(jrl,jhl),min(jru,jhu)
              do i=max(irl,ihl),min(iru,ihu)
                ij = ij + 1
                if (mm(i,j,1) == 0) cycle
                ! id is neighbour seen from receiver iam.
                if (id <= 4) then
                  if (mm(i+ii,j+jj,1) == 0) then
                    if ( (mm(i+im,   j+jm   ,1) == 0  .or.                     &
                          mm(i+im+ii,j+jm+jj,1) == 0) .and.                    &
                         (mm(i+ip,   j+jp   ,1) == 0  .or.                     &
                          mm(i+ip+ii,j+jp+jj,1) == 0)       )  cycle
                  endif
                else
                  if (mm(i+ii,j+jj,1) == 0) cycle
                  if (mm(i+im,j+jm,1) == 0 .and. mm(i+ip,j+jp,1) == 0) cycle
                endif
                decode(ia,id)%mask(ij) = ij
              enddo
            enddo
          endif
        enddo  ! id

      ! I should not receive from myself, but send to up to 8 of my neighbours:
      else

        do id=1,max_nb
          nb = mpi_tt(ia,itp1)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #itp1 and I must send to him

            ! Find my neighbour which corresponds to task #itp1
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == itp1) then
              ! send-buffer coordinates:
              select case (idi)
                case (1) ! west
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = jl
                  jsu = jl
                  ij  = max(isl,il) - 1
                  ii  =  0
                  jj  = -1
                  ip  = +1
                  jp  =  0
                  im  = -1
                  jm  =  0
                case (2) ! north
                  isl = il
                  isu = il
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                  ii  = -1
                  jj  =  0
                  ip  =  0
                  jp  = +1
                  im  =  0
                  jm  = -1
                case (3) ! east
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = ju
                  jsu = ju
                  ij  = max(isl,il) - 1
                  ii  =  0
                  jj  = +1
                  ip  = +1
                  jp  =  0
                  im  = -1
                  jm  =  0
                case (4) ! south
                  isl = iu
                  isu = iu
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                  ii  = +1
                  jj  =  0
                  ip  =  0
                  jp  = +1
                  im  =  0
                  jm  = -1
                case (5) ! north-west
                  isl = il
                  isu = il
                  jsl = jl
                  jsu = jl
                  ij  =  0
                  ii  = -1 
                  jj  = -1
                  ip  = -1
                  jp  =  0
                  im  =  0
                  jm  = -1
                case (6) ! north-east
                  isl = il
                  isu = il
                  jsl = ju
                  jsu = ju
                  ij  =  0
                  ii  = -1 
                  jj  = +1
                  ip  =  0
                  jp  = +1
                  im  = -1
                  jm  =  0
                case (7) ! south-east
                  isl = iu
                  isu = iu
                  jsl = ju
                  jsu = ju
                  ij  =  0
                  ii  = +1 
                  jj  = +1
                  ip  =  0
                  jp  = +1
                  im  = +1
                  jm  =  0
                case (8) ! south-west
                  isl = iu
                  isu = iu
                  jsl = jl
                  jsu = jl
                  ij  =  0
                  ii  = +1 
                  jj  = -1
                  ip  = +1
                  jp  =  0
                  im  =  0
                  jm  = -1
              end select

              ! assign the encode mask values:
              do j=max(jsl,jl),min(jsu,ju)
                do i=max(isl,il),min(isu,iu)
                  ij = ij + 1
                  if (mm(i,j,1) == 0) cycle

                  ! id  is neighbour seen from receiver itp1.
                  ! idi is neighbour seen from sender iam.
                  if (idi <= 4) then
                    if (mm(i+ii,j+jj,1) == 0) then
                      if ( (mm(i+im,   j+jm   ,1) == 0  .or.                   &
                            mm(i+im+ii,j+jm+jj,1) == 0) .and.                  &
                           (mm(i+ip,   j+jp   ,1) == 0  .or.                   &
                            mm(i+ip+ii,j+jp+jj,1) == 0)       )  cycle
                    endif
                  else
                    if (mm(i+ii,j+jj,1) == 0) cycle
                    if (mm(i+im,j+jm,1) == 0 .and. mm(i+ip,j+jp,1) == 0) cycle
                  endif
                  encode(ia,idi)%mask(ij) = ij
                enddo
              enddo
            endif  ! nbi == itp1
          endif    ! nb == iam
        enddo      ! id

      endif
    enddo   ! it_recv, itp1


  end subroutine dmpi_set_comm_masks

!===============================================================================

  subroutine dmpi_bufsize(kmx, mm, ia, recv, send, rsrf, ssrf )

    implicit none

    integer(4), intent(in)  :: kmx, mm(0:,0:,:), ia
    integer(4), intent(out) :: recv(:), send(:), rsrf(:), ssrf(:)

    integer(4) :: wrecv, wsend, wrsrf, wssrf
    integer(4) :: iam, it_recv, itp1, id, nb, idi, nbi
    integer(4) :: i, j, k, ij
    integer(4) :: irl, iru, jrl, jru, isl, isu, jsl, jsu
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    iam = mpi_rank + 1

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    do it_recv=0,mpi_size-1
      itp1 = it_recv + 1

      ! I should not send to myself, but receive from up to 8 of my neighbours:
      if (itp1 == iam) then
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (nb > 0) then
            ! recv-buffer coordinates:
            select case (id)
              case (1) ! west
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhl
                jru = jhl
                ij  = max(irl,ihl) - 1
              case (2) ! north
                irl = ihl
                iru = ihl
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
              case (3) ! east
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhu
                jru = jhu
                ij  = max(irl,ihl) - 1
              case (4) ! south
                irl = ihu
                iru = ihu
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
              case (5) ! north-west
                irl = ihl
                iru = ihl
                jrl = jhl
                jru = jhl
                ij  =  0
              case (6) ! north-east
                irl = ihl
                iru = ihl
                jrl = jhu
                jru = jhu
                ij  =  0
              case (7) ! south-east
                irl = ihu
                iru = ihu
                jrl = jhu
                jru = jhu
                ij  =  0
              case (8) ! south-west
                irl = ihu
                iru = ihu
                jrl = jhl
                jru = jhl
                ij  =  0
            end select

            ! calculate No. of wet points to recv/decode:
            wrsrf = 0
            wrecv = 0
            do j=max(jrl,jhl),min(jru,jhu)
              do i=max(irl,ihl),min(iru,ihu)
                ij = ij + 1
                if (decode(ia,id)%mask(ij) == 0) cycle
                wrsrf = wrsrf + 1
                do k=1,kmx
                  if (mm(i,j,k) == 0) exit
                  wrecv = wrecv + 1
                enddo
              enddo
            enddo
            if (wrecv > 0) then
              rsrf(id) = wrsrf
              recv(id) = wrecv
            endif

          endif
        enddo  ! id

      ! I should not receive from myself, but send to up to 8 of my neighbours:
      else

        do id=1,max_nb
          nb = mpi_tt(ia,itp1)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #itp1 and I must send to him

            ! Find my neighbour which corresponds to task #itp1
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == itp1) then
              ! send-buffer coordinates:
              select case (idi)
                case (1) ! west
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = jl
                  jsu = jl
                  ij  = max(isl,il) - 1
                case (2) ! north
                  isl = il
                  isu = il
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                case (3) ! east
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = ju
                  jsu = ju
                  ij  = max(isl,il) - 1
                case (4) ! south
                  isl = iu
                  isu = iu
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                case (5) ! north-west
                  isl = il
                  isu = il
                  jsl = jl
                  jsu = jl
                  ij  =  0
                case (6) ! north-east
                  isl = il
                  isu = il
                  jsl = ju
                  jsu = ju
                  ij  =  0
                case (7) ! south-east
                  isl = iu
                  isu = iu
                  jsl = ju
                  jsu = ju
                  ij  =  0
                case (8) ! south-west
                  isl = iu
                  isu = iu
                  jsl = jl
                  jsu = jl
                  ij  =  0
              end select

              ! calculate No. of wet points to send/encode:
              wssrf = 0
              wsend = 0
              do j=max(jsl,jl),min(jsu,ju)
                do i=max(isl,il),min(isu,iu)
                  ij = ij + 1
                  if (encode(ia,idi)%mask(ij) == 0) cycle
                  wssrf = wssrf + 1
                  do k=1,kmx
                    if (mm(i,j,k) == 0) exit
                    wsend = wsend + 1
                  enddo
                enddo
              enddo
              if (wsend > 0) then
                ssrf(idi) = wssrf
                send(idi) = wsend
              endif

            endif  ! nbi == itp1
          endif    ! nb == iam
        enddo      ! id

      endif
    enddo   ! it_recv, itp1

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_bufsize

!===============================================================================

  subroutine dmpi_roff (id, ia, iam, trecv, ij)
    implicit none

    integer(4), intent(in)  :: id, ia, iam, trecv
    integer(4), intent(out) :: ij

    integer(4) :: isl, jsl, il, jl

    isl = mpi_tt(ia,trecv)%low_hi
    jsl = mpi_tt(ia,trecv)%low_hj
    il  = mpi_tt(ia,iam)%low_i
    jl  = mpi_tt(ia,iam)%low_j

    if     (id == 1) then
      ij = max(isl,il) - 1
    elseif (id == 2) then
      ij = max(jsl,jl) - 1
    elseif (id == 3) then
      ij = max(isl,il) - 1
    elseif (id == 4) then
      ij = max(jsl,jl) - 1
    elseif (id >= 5) then
      ij =  0
    endif

  end subroutine dmpi_roff

!===============================================================================

  subroutine dmpi_soff (id, il, iu, jl, ju, ihl, ihu, jhl, jhu,                &
                        irl, iru, jrl, jru, ij)
    implicit none

    integer(4), intent(in)  :: id, il, iu, jl, ju, ihl, ihu, jhl, jhu
    integer(4), intent(out) :: irl, iru, jrl, jru, ij

    select case (id)
      case (1) ! west
        irl = il
        iru = iu
        jrl = jhl
        jru = jhl
        ij  = max(irl,ihl) - 1

      case (2) ! north
        irl = ihl
        iru = ihl
        jrl = jl
        jru = ju
        ij  = max(jrl,jhl) - 1

      case (3) ! east
        irl = il
        iru = iu
        jrl = jhu
        jru = jhu
        ij  = max(irl,ihl) - 1

      case (4) ! south
        irl = ihu
        iru = ihu
        jrl = jl
        jru = ju
        ij  = max(jrl,jhl) - 1

      case (5) ! north-west
        irl = ihl
        iru = ihl
        jrl = jhl
        jru = jhl
        ij  =  0

      case (6) ! north-east
        irl = ihl
        iru = ihl
        jrl = jhu
        jru = jhu
        ij  =  0

      case (7) ! south-east
        irl = ihu
        iru = ihu
        jrl = jhu
        jru = jhu
        ij  =  0

      case (8) ! south-west
        irl = ihu
        iru = ihu
        jrl = jhl
        jru = jhl
        ij  =  0
    end select

  end subroutine dmpi_soff

!===============================================================================

  subroutine dmpi_halosize(kmx, mm, il, iu, jl, ju, ihl, ihu, jhl, jhu, nh3,nh2)

    implicit none

    integer(4), intent(in)  :: kmx, mm(0:,0:,:)
    integer(4), intent(in)  :: il, iu, jl, ju, ihl, ihu, jhl, jhu
    integer(4), intent(out) :: nh3, nh2

    integer(4) :: i, j, k

    nh3 = 0
    nh2 = 0
    !  w, n/w, s/w
    do j=jhl,jl-1
      do i=ihl,ihu
        if (mm(i,j,1) == 0) cycle
        nh2 = nh2 + 1
        do k=1,kmx
          if (mm(i,j,k) == 0) exit
          nh3 = nh3 + 1
        enddo
      enddo
    enddo
    !  e, n/e, s/e
    do j=ju+1,jhu
      do i=ihl,ihu
        if (mm(i,j,1) == 0) cycle
        nh2 = nh2 + 1
        do k=1,kmx
          if (mm(i,j,k) == 0) exit
          nh3 = nh3 + 1
        enddo
      enddo
    enddo
    !  n
    do j=jl,ju
      do i=ihl,il-1
        if (mm(i,j,1) == 0) cycle
        nh2 = nh2 + 1
        do k=1,kmx
          if (mm(i,j,k) == 0) exit
          nh3 = nh3 + 1
        enddo
      enddo
    enddo
    !  s
    i = ihu
    do j=jl,ju
      do i=iu+1,ihu
        if (mm(i,j,1) == 0) cycle
        nh2 = nh2 + 1
        do k=1,kmx
          if (mm(i,j,k) == 0) exit
          nh3 = nh3 + 1
        enddo
      enddo
    enddo

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_halosize

!===============================================================================

#endif

!===============================================================================

  subroutine debug_info(ia)

    use io_subs, only : flush_unit

    implicit none

    integer(4), intent(in) :: ia

    integer(4) :: i, iam

    write(iu06,'(5X,A20,I8)')          'rank: ', mpi_rank
    write(iu06,'(5X,A20,5I10)')        'array config: ',                       &
          dd(ia)%dim1, dd(ia)%dim2, dd(ia)%dim3, dd(ia)%iwet2, dd(ia)%iwet3
    write(iu06,'(5X,A20,I8)')          'no of wetpoints: ',                    &
          dd(ia)%nwet3
    write(iu06,'(5X,A20,I8,I8)')       'surface wet range: ',                  &
          dd(ia)%low_ws, dd(ia)%up_ws
    write(iu06,'(5X,A20,I8,I8)')       'plane data range i: ',                 &
          dd(ia)%low_i, dd(ia)%up_i
    write(iu06,'(5X,A20,I8,I8)')       'plane data range j: ',                 &
          dd(ia)%low_j, dd(ia)%up_j
    write(iu06,'(5X,A20,I8,I8,I8,I8)') 'plane halo range i: ',                 &
          dd(ia)%low_hi, dd(ia)%up_hi
    write(iu06,'(5X,A20,I8,I8,I8,I8)') 'plane halo range j: ',                 &
          dd(ia)%low_hj, dd(ia)%up_hj

    iam = mpi_rank + 1
    do i=1,8
      if (dd(ia)%recv_bufsize(i) > 0) then 
        write(iu06,'(5X,A20,2I8)') 'must receive from: ',                      &
              mpi_tt(ia,iam)%hnb(i), dd(ia)%recv_bufsize(i)
      endif
    enddo
    do i=1,8
      if (dd(ia)%send_bufsize(i) > 0) then 
        write(iu06,'(5X,A20,2I8)') 'must send to:      ',                      &
              mpi_tt(ia,iam)%hnb(i), dd(ia)%send_bufsize(i)
      endif
    enddo

    call flush_unit(iu06)

  end subroutine debug_info

!===============================================================================

  subroutine dmpi_debug(narea)
    implicit none

    integer(4), intent(in) :: narea

    integer(4) :: ia

    do ia=1,narea
      call debug_info(ia)
    enddo
    
  end subroutine dmpi_debug

!===============================================================================

  subroutine dmpi_validate_pos(a,mmx,nmx,kmx,mmk,txt)
    use io_subs, only : flush_unit

    implicit none

    real(8),                intent(in) :: a(0:)
    integer(4),             intent(in) :: mmx,nmx,kmx,mmk(:,0:,0:)
    character(*), optional, intent(in) :: txt

    integer(4) :: i,j,k,imn,jmn,kmn,imx,jmx,kxx, mm
    real(8)    :: amn, amx

    if (size(a) > 1) then  
      amn = maxval(a(1:))
      amx = minval(a(1:))
      imn = -1
      jmn = -1
      kmn = -1
      imx = -1
      jmx = -1
      kxx = -1
      do j=1,nmx
        do i=1,mmx
          do k=1,kmx
            mm = mmk(k,i,j)
            if (mm <= 0) exit
            if (a(mm) < amn) then
              amn = a(mm)
              imn = i
              jmn = j
              kmn = k
            endif
            if (a(mm) > amx) then
              amx = a(mm)
              imx = i
              jmx = j
              kxx = k
            endif
          enddo
        enddo
      enddo
      if (present(txt)) write(iu06,*) 'MPI VALIDATE ', trim(txt)
      write(iu06,*) 'MPI VALIDATE: MINPOS:',imn,jmn,kmn,amn
      write(iu06,*) 'MPI VALIDATE: MAXPOS:',imx,jmx,kxx,amx
      call flush_unit(iu06)
    endif

  end subroutine dmpi_validate_pos

!===============================================================================

  subroutine dmpi_validate_min_max(a,iw2,halo,txt)
    use io_subs, only : flush_unit

    implicit none

    real(8),                intent(in) :: a(0:)
    integer(4),             intent(in) :: iw2, halo
    character(*), optional, intent(in) :: txt

    integer(4) :: nl, nu

    if (present(txt)) write(iu06,*) 'MPI VALIDATE ', trim(txt)
    write(iu06,*) 'MPI VALIDATE: Size of array: ', size(a)
    write(iu06,*) 'MPI VALIDATE: Point 0: ', a(0)
    if (size(a) > 1) then  
      write(iu06,*) 'MPI VALIDATE: MIN, MAX for inner domain: ',               &
                    minval(a(1:iw2)), maxval(a(1:iw2))
      write(iu06,*) 'MPI VALIDATE: MAXLOC for inner domain: ',maxloc(a(1:iw2))
      if (halo > 0) then
        nl = iw2 + 1
        nu = iw2 + halo
        write(iu06,*) 'MPI VALIDATE: MIN, MAX for halo: ',                     &
                      minval(a(nl:nu)), maxval(a(nl:nu))
        write(iu06,*) 'MPI VALIDATE: MAXLOC for halo: ',maxloc(a(nl:nu))
      endif
      call flush_unit(iu06)
    endif

  end subroutine dmpi_validate_min_max

!===============================================================================

  subroutine dmpi_validate_r8(i, j, k, mmk, a, txt)

    use io_subs, only : flush_unit

    implicit none

    real(8),                intent(in) :: a(0:)
    integer(4),             intent(in) :: i, j, k, mmk(:,0:,0:)
    character(*), optional, intent(in) :: txt

    write(iu06,'(a24,3i5)') 'MPI VALIDATE: position: ',i,j,k 
    if (present(txt)) then
      write(iu06,*) 'MPI VALIDATE '//trim(txt)//':',a(mmk(k,i,j))
    else
      write(iu06,*) 'MPI VALIDATE: value:     ',a(mmk(k,i,j)) 
    endif
    call flush_unit(iu06)

  end subroutine dmpi_validate_r8

!===============================================================================

  subroutine dmpi_validate_log(i, j, k, mmk, a, txt)

    use io_subs, only : flush_unit

    implicit none

    logical,                intent(in) :: a(0:)
    integer(4),             intent(in) :: i, j, k, mmk(:,0:,0:)
    character(*), optional, intent(in) :: txt

    write(iu06,'(a24,3i5)') ' MPI VALIDATE: position: ',i,j,k 
    if (present(txt)) then
      write(iu06,*) 'MPI VALIDATE '//trim(txt)//':',a(mmk(k,i,j))
    else
      write(iu06,*) 'MPI VALIDATE: value:     ',a(mmk(k,i,j)) 
    endif
    call flush_unit(iu06)

  end subroutine dmpi_validate_log

!===============================================================================

  subroutine dmpi_gather_srf_nc(ia,nc,low,ind_l,mmk_l,a_l,mmk,a,trecv)
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, low, nc, trecv
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), ind_l(:,:)
    real(8),    intent(inout) :: a(:,low:)
    real(8),    intent(in)    :: a_l(:,low:)

    integer(4) :: iw, itag, ii, it, bufsize, itp1, ic, nl, nu, mij, n2d

    itag = 1

    if ((trecv >= mpi_size) .or. (trecv <0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    if (mpi_rank /= trecv) then

      if (dd(ia)%up_ws /= 0) then
        bufsize = dd(ia)%up_ws-dd(ia)%low_ws+1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ii, ic, iw)
        call domp_get_domain(1, bufsize, nl, nu)
        ii = (nl-1)*nc
        do iw=nl,nu
          do ic=1,nc
            ii       = ii + 1
            gbuf(ii) = a_l(ic,iw) 
          enddo
        enddo
!$OMP END PARALLEL
        call dmpi_send( gbuf, nc*bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        itp1 = it + 1
        if (it==trecv .and. dd(ia)%nwet3 > 0) then 
          n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
          mij = dd(ia)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
          call domp_get_domain(1, n2d, nl, nu)
          a(1:nc,mij+nl:mij+nu) = a_l(1:nc,nl:nu)
!$OMP END PARALLEL
        else
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            bufsize = mpi_tt(ia,itp1)%up_ws-mpi_tt(ia,itp1)%low_ws+1
            call dmpi_recv( gbuf, nc*bufsize, it, itag )
            mij = mpi_tt(ia,itp1)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ii, ic, iw)
            call domp_get_domain(1, bufsize, nl, nu)
            ii = (nl-1)*nc
            do iw=nl,nu
              do ic=1,nc
                ii = ii + 1
                a(ic,mij+iw) = gbuf(ii)
              enddo
            enddo
!$OMP END PARALLEL
          endif
        endif
      enddo

    endif
  end subroutine dmpi_gather_srf_nc

!===============================================================================

  subroutine dmpi_gather_srf(ia,ind_l,mmk_l,a_l,mmk,a,trecv)
    use dmi_omp, only : domp_get_domain
    implicit none

    integer(4), intent(in)    :: ia, trecv
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), ind_l(:,:) 
    real(8),    intent(inout) :: a(0:)
    real(8),    intent(in)    :: a_l(0:)

    integer(4) :: iw, itag, ii, it, bufsize, itp1, nl, nu, n2d, mij

    itag = 1

    if ((trecv >= mpi_size) .or. (trecv <0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    if (mpi_rank /= trecv) then

      if (dd(ia)%up_ws /= 0) then
        bufsize = dd(ia)%up_ws-dd(ia)%low_ws+1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
        call domp_get_domain(1, bufsize, nl, nu)
        gbuf(nl:nu) = a_l(nl:nu)
!$OMP END PARALLEL
        call dmpi_send( gbuf, bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        itp1 = it + 1
        if (it==trecv .and. dd(ia)%nwet3 > 0) then 
          n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
          mij = dd(ia)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
          call domp_get_domain(1,n2d, nl, nu)
          a(mij+nl:mij+nu) = a_l(nl:nu)
!$OMP END PARALLEL
        else
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            bufsize = mpi_tt(ia,itp1)%up_ws-mpi_tt(ia,itp1)%low_ws+1
            call dmpi_recv( gbuf, bufsize, it, itag )
            mij = mpi_tt(ia,itp1)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
            call domp_get_domain(1, bufsize, nl, nu)
            a(mij+nl:mij+nu) = gbuf(nl:nu)
!$OMP END PARALLEL
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_srf

!===============================================================================

  subroutine dmpi_gather_log(ia,a_l,a,trecv)
    use dmi_omp, only : domp_get_domain
    implicit none

    integer(4), intent(in)    :: ia, trecv
    logical,    intent(inout) :: a(0:)
    logical,    intent(in)    :: a_l(0:)

    integer(4) :: itag, it, bufsize, itp1, n2d, nl, nu, ii

    itag = 1

    if ((trecv >= mpi_size) .or. (trecv <0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    if (mpi_rank /= trecv) then

      if (dd(ia)%nwet3 > 0) then
        bufsize = dd(ia)%up_ws - dd(ia)%low_ws + 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
        call domp_get_domain(1, bufsize, nl, nu)
        lbuf(nl:nu) = a_l(nl:nu)
!$OMP END PARALLEL
        call dmpi_send( lbuf, bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        itp1 = it + 1
        if (it==trecv .and. dd(ia)%nwet3 > 0) then
          n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
          ii  = dd(ia)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
          call domp_get_domain(1, n2d, nl, nu)
          a(ii+nl:ii+nu) = a_l(nl:nu)
!$OMP END PARALLEL
        else
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            bufsize = mpi_tt(ia,itp1)%up_ws - mpi_tt(ia,itp1)%low_ws + 1
            call dmpi_recv( lbuf, bufsize, it, itag )
            ii = mpi_tt(ia,itp1)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
            call domp_get_domain(1, bufsize, nl, nu)
            a(ii+nl:ii+nu) = lbuf(nl:nu)
!$OMP END PARALLEL
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_log

!===============================================================================

  subroutine dmpi_gather_iota(ia,i_l,a_l,a,trecv)
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, trecv
    real(8),    intent(inout) :: a(:,0:)
    real(8),    intent(in)    :: a_l(:,0:), i_l(0:)

    integer(4), parameter :: nc = 3
    integer(4) :: itag, it, itp1, bufsize, n2d, ib, nl, nu, nn, i

    if ((trecv >= mpi_size) .or. (trecv < 0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    itag = 1

    if (mpi_rank /= trecv) then

      if (dd(ia)%nwet3 > 0) then
        n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (ib, nl, nu)
        call domp_get_domain(1, n2d, nl, nu)
        do ib=nl,nu
          gbuf(1 + nc*(ib-1)) = i_l(ib)
          gbuf(2 + nc*(ib-1)) = a_l(1,ib)
          gbuf(3 + nc*(ib-1)) = a_l(2,ib)
        enddo
!$OMP END PARALLEL
        bufsize = n2d
        call dmpi_send( gbuf, nc*bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        itp1 = it + 1
        if (it==trecv .and. dd(ia)%nwet3 > 0) then
          n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
          nn  = dd(ia)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
          call domp_get_domain(1, n2d, nl, nu)
          a(1   ,nn+nl:nn+nu) = i_l(nl:nu)
          a(2:nc,nn+nl:nn+nu) = a_l(1:2,nl:nu)
!$OMP END PARALLEL
        else
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            n2d = mpi_tt(ia,itp1)%up_ws - mpi_tt(ia,itp1)%low_ws + 1
            bufsize = n2d
            call dmpi_recv( gbuf, nc*bufsize, it, itag )

            i  = mpi_tt(ia,itp1)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (ib, nl, nu)
            call domp_get_domain(1, n2d, nl, nu)
            do ib=nl,nu
              a(1,i+ib) = gbuf(1 + nc*(ib-1))
              a(2,i+ib) = gbuf(2 + nc*(ib-1))
              a(3,i+ib) = gbuf(3 + nc*(ib-1))
            enddo
!$OMP END PARALLEL
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_iota

!===============================================================================

  subroutine dmpi_gather_all_nc(ia,kmx,a_l,a,nc,trecv)
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, kmx, trecv, nc
    real(8),    intent(inout) :: a(:,0:)
    real(8),    intent(in)    :: a_l(:,0:)

    integer(4) :: itag, ii, it, itp1, bufsize, i, j, n2d, ic, ib, bs, nl, nu, nn

    if ((trecv >= mpi_size) .or. (trecv < 0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    itag = 1

    if (mpi_rank /= trecv) then

      if (dd(ia)%nwet3 > 0) then
        n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
        bs  = dd(ia)%nwet3
        ii  = (n2d+1) + dd(ia)%halo2
        j   = bs - (n2d+1)
!$OMP PARALLEL DEFAULT (shared) PRIVATE (ib, ic, nl, nu)
        call domp_get_domain(1, n2d, nl, nu)
        do ib=nl,nu
          do ic=1,nc
            gbuf(ic+nc*(ib-1)) = a_l(ic,ib)
          enddo
        enddo
        if (kmx > 1) then
          call domp_get_domain(ii, ii+j, nl, nu)
          do ib=nl,nu
            do ic=1,nc
              gbuf(ic+nc*(ib-dd(ia)%halo2-1)) = a_l(ic,ib)
            enddo
          enddo
        endif
!$OMP END PARALLEL
        if (kmx > 1) then
          bufsize = bs
        else
          bufsize = n2d
        endif
        call dmpi_send( gbuf, nc*bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        itp1 = it + 1
        if (it==trecv .and. dd(ia)%nwet3 > 0) then
          n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
          nn  = dd(ia)%low_ws - 1
          j   = dd(ia)%nwet3 - (n2d+1)
          ii  = (n2d+1) + dd(ia)%halo2
          i   = dd(ia)%low_w3 - ii
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
          call domp_get_domain(1, n2d, nl, nu)
          a(1:nc,nn+nl:nn+nu) = a_l(1:nc,nl:nu)
          if (kmx > 1) then
            call domp_get_domain(ii, ii+j, nl, nu)
            a(1:nc,i+nl:i+nu) = a_l(1:nc,nl:nu)
          endif
!$OMP END PARALLEL
        else
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            n2d = mpi_tt(ia,itp1)%up_ws - mpi_tt(ia,itp1)%low_ws + 1
            if (kmx > 1) then
              bufsize = mpi_tt(ia,itp1)%nwet3
            else
              bufsize = n2d
            endif
            call dmpi_recv( gbuf, nc*bufsize, it, itag )

            i  = mpi_tt(ia,itp1)%low_ws - 1
            j  = mpi_tt(ia,itp1)%nwet3 - (n2d+1)
            ii = (n2d+1)
            nn = mpi_tt(ia,itp1)%low_w3 - ii
!$OMP PARALLEL DEFAULT (shared) PRIVATE (ib, ic, nl, nu)
            call domp_get_domain(1, n2d, nl, nu)
            do ib=nl,nu
              do ic=1,nc
                a(ic,i+ib) = gbuf(ic+nc*(ib-1))
              enddo
            enddo
            if (kmx > 1) then
              call domp_get_domain(ii, ii+j, nl, nu)
              do ib=nl,nu
                do ic=1,nc
                  a(ic,nn+ib) = gbuf(ic+nc*(ib-1))
                enddo
              enddo
            endif
!$OMP END PARALLEL
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_all_nc

!===============================================================================

  subroutine dmpi_gather_all(ia,kmx,a_l,a,trecv,b_l,b)
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, kmx, trecv
    real(8),    intent(inout) :: a(0:)
    real(8),    intent(in)    :: a_l(0:)
    real(8),    intent(inout), optional :: b(0:)
    real(8),    intent(in)   , optional :: b_l(0:)

    integer(4) :: itag, ii, it, itp1, bufsize, i, j, n2d, bs, nn, nl, nu
    integer(4) :: ml, mu, boff
    logical    :: b2

    if ((trecv >= mpi_size) .or. (trecv < 0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    itag = 1

    b2 = (present(b_l) .and. present(b))

    if (mpi_rank /= trecv) then

      if (dd(ia)%nwet3 > 0) then
        n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
        bs  = dd(ia)%nwet3
        ii  = (n2d+1) + dd(ia)%halo2
        j   = bs - (n2d+1)
        i   = (n2d+1) - ii
        if (kmx > 1) then
          bufsize = bs
        else
          bufsize = n2d
        endif
        if (b2) then
          boff    = bufsize + 1
          bufsize = 2*bufsize
        endif
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ml, mu)
        call domp_get_domain(1, n2d, nl, nu)
        gbuf(nl:nu) = a_l(nl:nu)
        if (kmx > 1) then
          call domp_get_domain(ii, ii+j, ml, mu)
          gbuf(i+ml:i+mu) = a_l(ml:mu)
        endif
        if (b2) then
          gbuf(boff+nl:boff+nu) = b_l(nl:nu)
          if (kmx > 1) gbuf(boff+i+ml:boff+i+mu) = b_l(ml:mu)
        endif
!$OMP END PARALLEL
        call dmpi_send( gbuf, bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        itp1 = it + 1
        if (it==trecv .and. dd(ia)%nwet3 > 0) then
          n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
          nn  = dd(ia)%low_ws - 1
          j   = dd(ia)%nwet3 - (n2d+1)
          ii  = (n2d+1) + dd(ia)%halo2
          i   = dd(ia)%low_w3 - ii
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ml, mu)
          call domp_get_domain(1, n2d, nl, nu)
          a(nn+nl:nn+nu) = a_l(nl:nu)
          if (kmx > 1) then
            call domp_get_domain(ii, ii+j, ml, mu)
            a(i+ml:i+mu) = a_l(ml:mu)
          endif
          if (b2) then
            b(nn+nl:nn+nu) = b_l(nl:nu)
            if (kmx > 1) b(i+ml:i+mu) = b_l(ml:mu)
          endif
!$OMP END PARALLEL
        else
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            n2d = mpi_tt(ia,itp1)%up_ws - mpi_tt(ia,itp1)%low_ws + 1
            if (kmx > 1) then
              bufsize = mpi_tt(ia,itp1)%nwet3
            else
              bufsize = n2d
            endif
            if (b2) then
              boff    = bufsize + 1
              bufsize = 2*bufsize
            endif
            call dmpi_recv( gbuf, bufsize, it, itag )

            i  = mpi_tt(ia,itp1)%low_ws - 1
            j  = mpi_tt(ia,itp1)%nwet3 - (n2d+1)
            ii = (n2d+1)
            nn = mpi_tt(ia,itp1)%low_w3 - ii
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ml, mu)
            call domp_get_domain(1, n2d, nl, nu)
            a(i+nl:i+nu) = gbuf(nl:nu)
            if (kmx > 1) then
              call domp_get_domain(ii, ii+j, ml, mu)
              a(nn+ml:nn+mu) = gbuf(ml:mu)
            endif
            if (b2) then
              b(i+nl:i+nu) = gbuf(boff+nl:boff+nu)
              if (kmx > 1) b(nn+ml:nn+mu) = gbuf(boff+ml:boff+mu)
            endif
!$OMP END PARALLEL
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_all

!===============================================================================

  subroutine dmpi_gather_uv_bdr_cmp_nb(ia, iia, n1, n2, kr, recv, send, nc,    &
                                       kmx, mmk_l, mmk, mmk_c, kh_l, kh,       &
                                       iga, jga, uvdir, ii, FLAG, cmp_l, cmp)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    logical,      intent(in)    :: recv, send
    integer(4),   intent(in)    :: ia, iia, nc, kmx, n1, n2, kr(:,0:), uvdir, ii
    integer(4),   intent(in)    :: FLAG
    integer(4),   intent(in)    :: mmk_l(:,0:,0:), mmk(:,0:,0:), iga(:), jga(:)
    integer(4),   intent(in)    :: mmk_c(:,0:,0:)
    integer(4),   intent(in)    :: kh_l(0:), kh(0:)
    real(8),      intent(in)    :: cmp_l(:,0:)
    real(8),      intent(inout) :: cmp(:,0:)

    integer(4) :: it, iam, itag, iz, k, ic, nuv, ig, jg, il, iu, jl, ju
    integer(4) :: ig1, jg1, ig2, jg2, iii, jjj, iff, jff, midx, i, j
    integer(4) :: kb, ml2, mm2, btmp, tnum, itr
    integer(4) :: nzl, nzu, ireqi, irequ, ireql, imask, nnreq
    ! these must only be used from MASTER:
    integer(4) :: b_size

    !-  nothing to do here? ----------------------------------------------------
    if (n2 > 0) then
      nuv = n2 - n1 + 1
    else
      nuv = 0
    endif
    if (nuv < 1) return

    !-  who am I? --------------------------------------------------------------
    call domp_get_thread_no( tnum )

    !-  misc vars --------------------------------------------------------------
    if (do_comm) itag = nestitags2(ia)%p(ii,uvdir+2)
    ig1  = iga(1)
    ig2  = iga(2)
    jg1  = jga(1)
    jg2  = jga(2)
    call domp_get_domain(n1, n2, nzl, nzu, .true.)

    !-  Prepost non blocking receive call to iam -------------------------------
    if (recv  .and. (FLAG == 1 .or. FLAG == 4) .and. do_comm) then
      irequ = itagreq2(itag,1,2,uvdir)
      if (irequ > 0) then
        ireqi = itagreq2(itag,1,1,uvdir)

        do it=1,mpi_size
          bsiz3(tnum,it,1:2,uvdir) = 0

          ! receive but not from myself:
          imask = rsmask2(ia,it)%p(ii,uvdir)
          if (imask /= 0 .and. imask /= 10) then

            !  find buffer size:
            btmp = 0
            do iz=nzl,nzu
              ig = kr(1,iz)
              jg = kr(2,iz)
              if (kr(3,iz) == 3) then
                jg = jg+1
              elseif (kr(3,iz) == 4) then
                ig = ig+1
              endif
              if (.not. ((dd(ia)%low_i <= ig .and. ig <= dd(ia)%up_i) .and.    &
                         (dd(ia)%low_j <= jg .and. jg <= dd(ia)%up_j))  ) cycle
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < mpi_tt(iia,it)%low_j .or.                            &
                    jff > mpi_tt(iia,it)%up_j       ) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < mpi_tt(iia,it)%low_i .or.                          &
                      iff > mpi_tt(iia,it)%up_i       ) cycle
                  midx = mmk(1,iff,jff)
                  if (midx > 0) btmp = btmp + nc*min(kmx,kh(midx))
                enddo
              enddo
            enddo
            bsiz3(tnum,it,1,uvdir) = btmp
!$OMP BARRIER
!$OMP MASTER
            b_size = bsiz3(1,it,1,uvdir)
            bsiz3(1,it,2,uvdir) = 0
            do itr=2,mpi_nt
              b_size = b_size + bsiz3(itr,it,1,uvdir)
              bsiz3(itr,it,2,uvdir) = bsiz3(itr-1,it,1,uvdir)                  &
                                    + bsiz3(itr-1,it,2,uvdir)
            enddo

            !  recv a dummy in case of degenerated buffer, FIXME
            if (b_size == 0) b_size = 1

            !  prepost the receive:
            call dmpi_irecv(rbufnb3(ireqi)%p, b_size, it-1, itag, ireq3(ireqi))
            ireqi = ireqi + 1
!$OMP END MASTER
          endif
        enddo
      endif
    endif


    !-  if I have data, send it to the relevant task ---------------------------
    if (send .and. (FLAG == 2 .or. FLAG == 4) .and. do_comm) then
      irequ = itagreq2(itag,2,2,uvdir)
      if (irequ > 0) then
        ireqi = itagreq2(itag,2,1,uvdir)

        do it=1,mpi_size
          ! send but not to myself:
          imask = rsmask2(ia,it)%p(ii,uvdir)
          if (imask /= 0 .and. imask /= 1) then

            !  coarse grid dims on task #it
            il = mpi_tt(ia,it)%low_i
            iu = mpi_tt(ia,it)%up_i
            jl = mpi_tt(ia,it)%low_j
            ju = mpi_tt(ia,it)%up_j

            !  find thread offset:
            offset(tnum,1:2) = 0
            btmp = 0
            do iz=nzl,nzu
              ig = kr(1,iz)
              jg = kr(2,iz)
              if (kr(3,iz) == 3) then
                jg = jg+1
              elseif (kr(3,iz) == 4) then
                ig = ig+1
              endif
              if (ig < il .or. iu < ig) cycle
              if (jg < jl .or. ju < jg) cycle
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                  midx = mmk_l(1,iff,jff)
                  if (midx > 0) btmp = btmp + nc*min(kmx,kh_l(midx))
                enddo
              enddo
            enddo
            offset(tnum,1) = btmp
!$OMP BARRIER
!$OMP MASTER
            b_size = offset(1,1)
            offset(1,2) = 0
            do itr=2,mpi_nt
              b_size = b_size + offset(itr,1)
              offset(itr,2) = offset(itr-1,1) + offset(itr-1,2)
            enddo
!$OMP END MASTER
!$OMP BARRIER

            if (offset(tnum,1) > 0) then
              ! encode and send:
              btmp = offset(tnum,2)
              do iz=nzl,nzu
                ig = kr(1,iz)
                jg = kr(2,iz)
                if (kr(3,iz) == 3) then
                  jg = jg+1
                elseif (kr(3,iz) == 4) then
                  ig = ig+1
                endif
                if (ig < il .or. iu < ig) cycle
                if (jg < jl .or. ju < jg) cycle
                if (mmk_c(1,ig,jg) <= 0) cycle
                i = ig2*(ig-ig1) + 1
                j = jg2*(jg-jg1) + 1
                do jjj=0,jg2-1
                  jff = j + jjj
                  if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                  do iii=0,ig2-1
                    iff = i + iii
                    if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                    midx = mmk_l(1,iff,jff)
                    if (midx <= 0) cycle
                    sbufnb3(ireqi)%p(btmp+1:btmp+nc) = cmp_l(1:nc,midx)
                    btmp = btmp + nc
                    kb = min(kmx,kh_l(midx))
                    if (kb < 2) cycle
                    ml2 = mmk_l(2,iff,jff) - 2
                    do k=2,kb
                      sbufnb3(ireqi)%p(btmp+1:btmp+nc) = cmp_l(1:nc,ml2+k)
                      btmp = btmp + nc
                    enddo
                  enddo
                enddo
              enddo
            endif
!$OMP BARRIER
!$OMP MASTER
            !  send a dummy in case of degenerated buffer, FIXME
            if (b_size == 0) then
              b_size = 1
              sbufnb3(ireqi)%p(1) = dummy_value
            endif

            call dmpi_isend(sbufnb3(ireqi)%p, b_size, it-1, itag, irqs3(ireqi))
!$OMP END MASTER
            ireqi = ireqi + 1
          endif
        enddo
      endif
    endif

    !-  recv on this task from all tasks having any data -----------------------
    if (recv .and. (FLAG == 3 .or. FLAG == 4)) then
      ! copy from local to global array
      if (send .and. recv) then
        do iz=nzl,nzu
          ig = kr(1,iz)
          jg = kr(2,iz)
          if (kr(3,iz) == 3) then
            jg = jg+1
          elseif (kr(3,iz) == 4) then
            ig = ig+1
          endif
          if (ig < dd(ia)%low_i .or. dd(ia)%up_i < ig) cycle
          if (jg < dd(ia)%low_j .or. dd(ia)%up_j < jg) cycle
          if (mmk_c(1,ig,jg) <= 0) cycle
          i = ig2*(ig-ig1) + 1
          j = jg2*(jg-jg1) + 1
          do jjj=0,jg2-1
            jff = j + jjj
            if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
            do iii=0,ig2-1
              iff = i + iii
              if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
              midx = mmk_l(1,iff,jff)
              if (midx <= 0) cycle
              if (mmk(1,iff,jff) <= 0) cycle
              cmp(1:nc,mmk(1,iff,jff)) = cmp_l(1:nc,midx)
              kb = min(kmx,kh_l(midx))
              if (kb < 2) cycle
              ml2 = mmk_l(2,iff,jff) - 2
              mm2 = mmk(2,iff,jff) - 2
              do k=2,kb
                cmp(1:nc,mm2+k) = cmp_l(1:nc,ml2+k)
              enddo
            enddo
          enddo
        enddo
      endif

      !  wait for the buffers to get filled and then decode
      if (do_comm) then
        irequ = itagreq2(itag,1,2,uvdir)
        if (irequ > 0) then
          ireql = itagreq2(itag,1,1,uvdir)
          nnreq = irequ - ireql + 1
!$OMP MASTER
          call MPI_Waitall(nnreq,ireq3(ireql:irequ),lstat3(:,ireql:irequ),ierr)
!$OMP END MASTER
!$OMP BARRIER

          ireqi = ireql

          do it=1,mpi_size
            ! recieve but not from myself:
            imask = rsmask2(ia,it)%p(ii,uvdir)
            if (imask /= 0 .and. imask /= 10) then
              if (bsiz3(tnum,it,1,uvdir) > 0) then
                ! decode buffer:
                btmp = bsiz3(tnum,it,2,uvdir)
                do iz=nzl,nzu
                  ig = kr(1,iz)
                  jg = kr(2,iz)
                  if (kr(3,iz) == 3) then
                    jg = jg+1
                  elseif (kr(3,iz) == 4) then
                    ig = ig+1
                  endif
                  if (.not.((dd(ia)%low_i <= ig .and. ig <= dd(ia)%up_i) .and. &
                            (dd(ia)%low_j <= jg .and. jg <= dd(ia)%up_j))) cycle
                  if (mmk_c(1,ig,jg) <= 0) cycle
                  i = ig2*(ig-ig1) + 1
                  j = jg2*(jg-jg1) + 1
                  do jjj=0,jg2-1
                    jff = j + jjj
                    if (jff < mpi_tt(iia,it)%low_j .or.                        &
                        jff > mpi_tt(iia,it)%up_j       ) cycle
                    do iii=0,ig2-1
                      iff = i + iii
                      if (iff < mpi_tt(iia,it)%low_i .or.                      &
                          iff > mpi_tt(iia,it)%up_i       ) cycle
                      midx = mmk(1,iff,jff)
                      if (midx <= 0) cycle
                      cmp(1:nc,midx) = rbufnb3(ireqi)%p(btmp+1:btmp+nc)
                      btmp = btmp + nc
                      kb = min(kmx,kh(midx))
                      if (kb < 2) cycle
                      mm2 = mmk(2,iff,jff) - 2
                      do k=2,kb
                        cmp(1:nc,mm2+k) = rbufnb3(ireqi)%p(btmp+1:btmp+nc)
                        btmp = btmp + nc
                      enddo
                    enddo
                  enddo
                enddo
              endif
              ireqi = ireqi + 1
            endif
          enddo
        endif
      endif
    endif


    !-  complete the send from this task ---------------------------------------
    if (send .and. (FLAG == 3 .or. FLAG == 4) .and. do_comm) then
      !  wait for the send to finish
      irequ = itagreq2(itag,2,2,uvdir)
      if (irequ > 0) then
!$OMP MASTER
        ireql = itagreq2(itag,2,1,uvdir)
        nnreq = irequ - ireql + 1
        call MPI_Waitall(nnreq, irqs3(ireql:irequ), lstat3(:,ireql:irequ), ierr)
!$OMP END MASTER
      endif
    endif

  end subroutine dmpi_gather_uv_bdr_cmp_nb

!===============================================================================

  subroutine dmpi_gather_mcf_nb(ia, iia, n1, n2, kr, recv, send, kh_l, kh,     &
                                kmx, mmk_l, mmk, mmk_f, iga, jga,              &
                                u_l, u, v_l, v, iao, FLAG)

    ! FLAG present and == 1: Pre-post the receives (if recv, of course)
    !                  == 2: Fill send buffers and pre-post thes sends (if send)
    !                  == 3: Waitall and decode the receive buffers (if recv)
    !                                and possibly copy l-->g (if send and recv)
    !                  == 4: do all the above.
!
! Usage:
!   ia = ...
!     ...
!   if (ia == mainarea) then
!     ...
!   else
!     !- fix velocity by nesting:
!     do iao=1,nestingfrom(ia)%p(0)
!       !  ia:  fine grid
!       !  iia: coarse grid
!       iia = nestingfrom(ia)%p(iao)
!       nz1 = znesting(ia)%p(iao,1)
!       nz2 = znesting(ia)%p(iao,2)
!
!       !  gather coarse grid velocity to global arrays:
!       call dmpi_gather(ia, iia, nz1, nz2, krz(ia)%p,                         &
!                        nftable(mpi_rank+1,iia,ia), mctable(iia,ia),          &
!                        kh_l(iia)%p, kh(iia)%p, kmx(iia),                     &
!                        mm1k_l(iia)%p, mm1k(iia)%p, mm1k(ia)%p,               &
!                        iga(ia)%p(:,iao), jga(ia)%p(:,iao),                   &
!                        un_l(iia)%p, wrk(iia)%p, vn_l(iia)%p, wrk2(iia)%p,    &
!                        iao, FLAG)
!     enddo
!   endif
!

    use dmi_omp, only : domp_get_domain

    implicit none

    logical,      intent(in)    :: recv, send
    integer(4),   intent(in)    :: iao, FLAG

    integer(4),   intent(in)    :: ia, iia, n1, n2, kr(:,0:), iga(:), jga(:)
    integer(4),   intent(in)    :: kh_l(0:), kh(0:), kmx
    integer(4),   intent(in)    :: mmk_l(:,0:,0:), mmk(:,0:,0:), mmk_f(:,0:,0:)
    real(8),      intent(in)    :: u_l(0:), v_l(0:)
    real(8),      intent(inout) :: u(0:), v(0:)

    integer(4) :: it, iam, b_size, iz, nnreq, ireql, irequ, ireqi, imask, itag
    integer(4) :: ig1, jg1, ig2, jg2, il, iu, jl, ju, itl, itu, nzl, nzu
    integer(4) :: ilf, iuf, jlf, juf, iff, jff

    !-  some initializations ---------------------------------------------------
    iam  = mpi_rank + 1
    if (do_comm) itag = nestitags(ia)%p(iao)

    !-  Prepost non blocking receive call to iam -------------------------------
    if (recv .and. (FLAG == 1 .or. FLAG == 4) .and. do_comm) then
      irequ = itagreq(itag,1,2)
      if (irequ > 0) then
        call domp_get_domain(1, mpi_size, itl, itu)
        ireqi = itagreq(itag,1,1)
        do it=1,mpi_size
          ! receive but not from myself:
          imask = rsmask(ia,it)%p(iao)
          if (imask == 1 .or. imask == 11) then
            if (itl <= it .and. it <= itu) then
              if (kmx > 1) then
                b_size = rbs(ireqi)
              else
                b_size = rbs1(ireqi)
              endif
              !  prepost the receive:
              call dmpi_irecv(rbufnb(ireqi)%p, b_size, it-1, itag, ireq(ireqi))
            endif
            ireqi = ireqi + 1
          endif
        enddo
      endif
    endif


    !-  if I have data, send it to the relevant task ---------------------------
    if (send .and. (FLAG == 2 .or. FLAG == 4) .and. do_comm) then
      ig1 = iga(1)
      ig2 = iga(2)
      jg1 = jga(1)
      jg2 = jga(2)
      !  coarse grid dims on task #iam
      il = dd(iia)%low_i
      iu = dd(iia)%up_i
      jl = dd(iia)%low_j
      ju = dd(iia)%up_j

      irequ = itagreq(itag,2,2)
      if (irequ > 0) then
        call domp_get_domain(1, mpi_size, itl, itu, .true.)
        ireqi = itagreq(itag,2,1)
        do it=1,mpi_size
          ! send but not to myself:
          imask = rsmask(ia,it)%p(iao)
          if (imask == 10 .or. imask == 11) then
            if (itl <= it .and. it <= itu) then
              !  fine grid dims on task #it:
              ilf = mpi_tt(ia,it)%low_i
              iuf = mpi_tt(ia,it)%up_i
              jlf = mpi_tt(ia,it)%low_j
              juf = mpi_tt(ia,it)%up_j

              !  encode and send:
              b_size = 0
              do iz=n1,n2
                iff = kr(1,iz)
                jff = kr(2,iz)
                if (.not.(ilf <= iff .and. iff <= iuf .and.                    &
                          jlf <= jff .and. jff <= juf)      ) cycle
                if (mmk_f(1,iff,jff) <= 0) cycle
                call dmpi_mcf(iz, n2, kr, ig1, jg1, ig2, jg2, il, iu, jl, ju,  &
                              kh_l, kmx, (3), b_size, mmk_l, u_l, v_l,         &
                              sbufnb(ireqi)%p)
              enddo
              call dmpi_isend(sbufnb(ireqi)%p, b_size, it-1, itag, irqs(ireqi))
            endif
            ireqi = ireqi + 1
          endif
        enddo
      endif
    endif


    !-  recv on this task from all tasks having any data -----------------------
    if (recv .and. (FLAG == 3 .or. FLAG == 4)) then
      ig1 = iga(1)
      ig2 = iga(2)
      jg1 = jga(1)
      jg2 = jga(2)
      !  fine grid dims on task #iam:
      ilf = dd(ia)%low_i
      iuf = dd(ia)%up_i
      jlf = dd(ia)%low_j
      juf = dd(ia)%up_j

      ! copy from local to global array
      if (send .and. recv) then
        it = iam
        !  coarse grid dims on task #it
        il = mpi_tt(iia,it)%low_i
        iu = mpi_tt(iia,it)%up_i
        jl = mpi_tt(iia,it)%low_j
        ju = mpi_tt(iia,it)%up_j

        call domp_get_domain(n1, n2, nzl, nzu)
        b_size = 0
        do iz=nzl,nzu
          iff = kr(1,iz)
          jff = kr(2,iz)
          if (.not.(ilf <= iff .and. iff <= iuf .and.                          &
                    jlf <= jff .and. jff <= juf)      ) cycle
          if (mmk_f(1,iff,jff) <= 0) cycle
          call dmpi_mcf(iz, n2, kr, ig1, jg1, ig2, jg2, il, iu, jl, ju,        &
                     kh_l, kmx, (4), b_size, mmk_l, u_l, v_l, gbuf, mmk, u, v)
        enddo
      endif

      !  wait for the recv buffers to get filled and then decode
      if (do_comm) then
        irequ = itagreq(itag,1,2)
        if (irequ > 0) then
          ireql = itagreq(itag,1,1)
          nnreq = irequ - ireql + 1
!$OMP MASTER
          call MPI_Waitall(nnreq,ireq(ireql:irequ),lstatus(:,ireql:irequ),ierr)
!$OMP END MASTER
!$OMP BARRIER

          call domp_get_domain(1, mpi_size, itl, itu)
          ireqi = ireql
          do it=1,mpi_size
            imask = rsmask(ia,it)%p(iao)
            if (imask == 1 .or. imask == 11) then
              if (itl <= it .and. it <= itu) then
                !  coarse grid dims on task #it
                il = mpi_tt(iia,it)%low_i
                iu = mpi_tt(iia,it)%up_i
                jl = mpi_tt(iia,it)%low_j
                ju = mpi_tt(iia,it)%up_j

                ! decode buffer:
                b_size = 0
                do iz=n1,n2
                  iff = kr(1,iz)
                  jff = kr(2,iz)
                  if (.not.(ilf <= iff .and. iff <= iuf .and.                  &
                            jlf <= jff .and. jff <= juf)      ) cycle
                  if (mmk_f(1,iff,jff) <= 0) cycle
                  call dmpi_mcf(iz, n2, kr, ig1, jg1, ig2, jg2, il, iu, jl, ju,&
                                kh, kmx, (2), b_size, mmk, u_l, v_l,           &
                                rbufnb(ireqi)%p, mmk, u, v)
                enddo
                if ( (b_size /= rbs(ireqi)  .and. kmx > 1)  .or.               &
                     (b_size /= rbs1(ireqi) .and. kmx == 1)      )             &
                  call exitme(2,'MPI-mcf error, decode')
              endif
              ireqi = ireqi + 1
            endif
          enddo
        endif
      endif
    endif


    !-  complete the send from this task ---------------------------------------
    if (send .and. (FLAG == 3 .or. FLAG == 4) .and. do_comm) then
      !  wait for the send to finish
      irequ = itagreq(itag,2,2)
      if (irequ > 0) then
!$OMP MASTER
        ireql = itagreq(itag,2,1)
        nnreq = irequ - ireql + 1
        call MPI_Waitall(nnreq, irqs(ireql:irequ), lstatus(:,ireql:irequ), ierr)
!$OMP END MASTER
!$OMP BARRIER
      endif
    endif

  end subroutine dmpi_gather_mcf_nb

!===============================================================================

  subroutine dmpi_gather_copy_cmp(iia, ia, recv, send, noff, nc, mmx, nmx, kmx,&
                                  mmk_l, mmk, mmk_c, iga, jga, cmp_l, cmp)
    ! ia:  fine grid
    ! iia: coarse grid

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none
    logical,      intent(in)    :: recv, send
    integer(4),   intent(in)    :: ia, iia, noff, nc, mmx, nmx, kmx
    integer(4),   intent(in)    :: mmk_l(:,0:,0:), mmk(:,0:,0:), iga(:), jga(:)
    integer(4),   intent(in)    :: mmk_c(:,0:,0:)
    real(8),      intent(in)    :: cmp_l(:,0:)
    real(8),      intent(inout) :: cmp(:,0:)

    integer(4) :: it, iam, itag, b_size, ik, ic, ig, jg, jglo, jghi, tnum, boff
    integer(4) :: ig1, jg1, ig2, jg2, ig3, jg3, iii, jjj, iff, jff, midx, i, j
    integer(4) :: is, ie, js, je, il, iu, jl, ju, ittt, itr
!Fixme: we must find the max size needed and allocate bbbb only once
!       the present code is only proof-of-concept
    real(8), allocatable :: bbbb(:)

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1
    call domp_get_thread_no( tnum )

    !-  misc vars --------------------------------------------------------------
    itag = NextFreeItag + tnum - 1
    ig1  = iga(1)
    ig2  = iga(2)
    ig3  = iga(3)
    jg1  = jga(1)
    jg2  = jga(2)
    jg3  = jga(3)

    !  coarse grid bounds of the fine grid:
    is = ig1 + (ig2/2+ig3)/ig2
    ie = ig1 + (ig2/2+ig3+mmx-1)/ig2
    js = jg1 + (jg2/2+jg3)/jg2
    je = jg1 + (jg2/2+jg3+nmx-1)/jg2

    do ittt=1,mpi_size
      !-  recv on this task from all tasks having any data ---------------------
      if (ittt == iam .and. recv) then
        ! recieve but not from myself:
        !  coarse grid bounds:
        il = dd(iia)%low_i
        iu = dd(iia)%up_i
        jl = dd(iia)%low_j
        ju = dd(iia)%up_j
        call domp_get_domain(max(js,jl), min(je,ju), jglo, jghi)

        do it=1,mpi_size
          if (it == iam) cycle
          !  find buffer size:
          bsiz3(tnum,it,1:2,1) = 0
          do jg=jglo,jghi
            do ig=max(is,il),min(ie,iu)
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < mpi_tt(ia,it)%low_j .or.                             &
                    jff > mpi_tt(ia,it)%up_j       ) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < mpi_tt(ia,it)%low_i .or.                           &
                      iff > mpi_tt(ia,it)%up_i       ) cycle
                  do ik=1,kmx
                    if (mmk(ik,iff,jff) > 0) then
                      bsiz3(tnum,it,1,1) = bsiz3(tnum,it,1,1) + (nc-noff+1)
                    else
                      exit
                    endif
                  enddo  ! ki
                enddo    ! iii
              enddo      ! jjj
            enddo        ! ig
          enddo          ! jg

          !  decode:
          if (bsiz3(tnum,it,1,1) > 0) then
            allocate(bbbb(bsiz3(tnum,it,1,1)))
            call dmpi_recv(bbbb, bsiz3(tnum,it,1,1), it-1, itag)
            boff = 0
            do jg=jglo,jghi
              do ig=max(is,il),min(ie,iu)
                if (mmk_c(1,ig,jg) <= 0) cycle
                i = ig2*(ig-ig1) + 1
                j = jg2*(jg-jg1) + 1
                do jjj=0,jg2-1
                  jff = j + jjj
                  if (jff < mpi_tt(ia,it)%low_j .or.                           &
                      jff > mpi_tt(ia,it)%up_j       ) cycle
                  do iii=0,ig2-1
                    iff = i + iii
                    if (iff < mpi_tt(ia,it)%low_i .or.                         &
                        iff > mpi_tt(ia,it)%up_i       ) cycle
                    do ik=1,kmx
                      midx = mmk(ik,iff,jff)
                      if (midx > 0) then
                        do ic=noff,nc
                          boff         = boff + 1
                          cmp(ic,midx) = bbbb(boff)
                        enddo
                      else
                        exit
                      endif
                    enddo  ! ki
                  enddo    ! iii
                enddo      ! jjj
              enddo        ! ig
            enddo          ! jg
            deallocate(bbbb)
          endif
        enddo  ! it

        ! copy from local to global array
        if (send) then
          do jg=jglo,jghi
            do ig=max(is,il),min(ie,iu)
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(ia)%low_j .or. jff > dd(ia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(ia)%low_i .or. iff > dd(ia)%up_i) cycle
                  do ik=1,kmx
                    midx = mmk_l(ik,iff,jff)
                    if (midx > 0) then
                      cmp(noff:nc,mmk(ik,iff,jff)) = cmp_l(noff:nc,midx)
                    else
                      exit
                    endif
                  enddo  ! ki
                enddo    ! iii
              enddo      ! jjj
            enddo        ! ig
          enddo          ! jg
        endif ! send+recv


      !-  if I have data, send it to the relevant task -------------------------
      elseif (ittt /= iam .and. send) then
        ! send but not to myself:
        do it=ittt,ittt
          if (it == iam) cycle
          !  coarse grid bounds:
          il = mpi_tt(iia,it)%low_i
          iu = mpi_tt(iia,it)%up_i
          jl = mpi_tt(iia,it)%low_j
          ju = mpi_tt(iia,it)%up_j
          call domp_get_domain(max(js,jl), min(je,ju), jglo, jghi)

          !  encode and send:
          bsiz3(tnum,it,1:2,1) = 0
          do jg=jglo,jghi
            do ig=max(is,il),min(ie,iu)
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(ia)%low_j .or. jff > dd(ia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(ia)%low_i .or. iff > dd(ia)%up_i) cycle
                  do ik=1,kmx
                    midx = mmk_l(ik,iff,jff)
                    if (midx > 0) then
                      do ic=noff,nc
                        bsiz3(tnum,it,1,1) = bsiz3(tnum,it,1,1) + 1
                      enddo
                    else
                      exit
                    endif
                  enddo  ! ki
                enddo    ! iii
              enddo      ! jjj
            enddo        ! ig
          enddo          ! jg
          if (bsiz3(tnum,it,1,1) > 0) then
            allocate(bbbb(bsiz3(tnum,it,1,1)))
            boff = 0
            do jg=jglo,jghi
              do ig=max(is,il),min(ie,iu)
                if (mmk_c(1,ig,jg) <= 0) cycle
                i = ig2*(ig-ig1) + 1
                j = jg2*(jg-jg1) + 1
                do jjj=0,jg2-1
                  jff = j + jjj
                  if (jff < dd(ia)%low_j .or. jff > dd(ia)%up_j) cycle
                  do iii=0,ig2-1
                    iff = i + iii
                    if (iff < dd(ia)%low_i .or. iff > dd(ia)%up_i) cycle
                    do ik=1,kmx
                      midx = mmk_l(ik,iff,jff)
                      if (midx > 0) then
                        do ic=noff,nc
                          boff       = boff + 1
                          bbbb(boff) = cmp_l(ic,midx)
                        enddo
                      else
                        exit
                      endif
                    enddo  ! ki
                  enddo    ! iii
                enddo      ! jjj
              enddo        ! ig
            enddo          ! jg
            call dmpi_send(bbbb, boff, it-1, itag)
            deallocate(bbbb)
          endif
        enddo 

      endif 
    enddo   ! ittt

  end subroutine dmpi_gather_copy_cmp

!===============================================================================

  subroutine dmpi_gather_copy(iia, ia, recv, send, mmx, nmx, kmx, khg, khl,    &
                              mmk_l, mmk, mmk_c, iga, jga,                     &
                              a_l, a, b_l, b, c_l, c)
    ! ia:  fine grid
    ! iia: coarse grid

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none
    logical,      intent(in)    :: recv, send
    integer(4),   intent(in)    :: ia, iia, mmx, nmx, kmx, khg(0:), khl(0:)
    integer(4),   intent(in)    :: mmk_l(:,0:,0:), mmk(:,0:,0:), iga(:), jga(:)
    integer(4),   intent(in)    :: mmk_c(:,0:,0:)
    real(8),      intent(in)    :: a_l(0:)
    real(8),      intent(inout) :: a(0:)
    real(8),      intent(in),    optional :: b_l(0:)
    real(8),      intent(inout), optional :: b(0:)
    real(8),      intent(in),    optional :: c_l(0:)
    real(8),      intent(inout), optional :: c(0:)

    integer(4) :: it, iam, itag, b_size, ik, ig, jg, jglo, jghi, tnum, boff
    integer(4) :: ig1, jg1, ig2, jg2, ig3, jg3, iii, jjj, iff, jff, i, j
    integer(4) :: is, ie, js, je, il, iu, jl, ju, ittt, itr, mi0g, mi0l
    logical    :: b2, b3
!Fixme: we must find the max size needed and allocate bbbb only once
!       the present code is only proof-of-concept
    real(8), allocatable :: bbbb(:)

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1
    call domp_get_thread_no( tnum )

    !-  misc vars --------------------------------------------------------------
    itag = NextFreeItag + tnum - 1
    ig1  = iga(1)
    ig2  = iga(2)
    ig3  = iga(3)
    jg1  = jga(1)
    jg2  = jga(2)
    jg3  = jga(3)
    b2   = (present(b) .and. present(b_l))
    b3   = (present(c) .and. present(c_l)) .and. b2

    !  coarse grid bounds of the fine grid:
    is = ig1 + (ig2/2+ig3)/ig2
    ie = ig1 + (ig2/2+ig3+mmx-1)/ig2
    js = jg1 + (jg2/2+jg3)/jg2
    je = jg1 + (jg2/2+jg3+nmx-1)/jg2

    do ittt=1,mpi_size
      !-  recv on this task from all tasks having any data ---------------------
      if (ittt == iam .and. recv) then
        ! recieve but not from myself:
        !  coarse grid bounds:
        il = dd(iia)%low_i
        iu = dd(iia)%up_i
        jl = dd(iia)%low_j
        ju = dd(iia)%up_j
        call domp_get_domain(max(js,jl), min(je,ju), jglo, jghi)

        do it=1,mpi_size
          if (it == iam) cycle
          !  find buffer size:
          bsiz3(tnum,it,1:2,1) = 0
          do jg=jglo,jghi
            do ig=max(is,il),min(ie,iu)
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < mpi_tt(ia,it)%low_j .or.                             &
                    jff > mpi_tt(ia,it)%up_j       ) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < mpi_tt(ia,it)%low_i .or.                           &
                      iff > mpi_tt(ia,it)%up_i       ) cycle
                  if (khg(mmk(1,iff,jff)) >= 1)                                &
                    bsiz3(tnum,it,1,1) = bsiz3(tnum,it,1,1)                    &
                                       + min(kmx,khg(mmk(1,iff,jff)))
                enddo    ! iii
              enddo      ! jjj
            enddo        ! ig
          enddo          ! jg
          if (b2) then
            if (b3) then
              bsiz3(tnum,it,1,1) = 3*bsiz3(tnum,it,1,1)
            else
              bsiz3(tnum,it,1,1) = 2*bsiz3(tnum,it,1,1)
            endif
          endif

          !  receive and decode:
          if (bsiz3(tnum,it,1,1) > 0) then
            allocate(bbbb(bsiz3(tnum,it,1,1)))
            call dmpi_recv(bbbb, bsiz3(tnum,it,1,1), it-1, itag)
            boff = 0
            do jg=jglo,jghi
              do ig=max(is,il),min(ie,iu)
                if (mmk_c(1,ig,jg) <= 0) cycle
                i = ig2*(ig-ig1) + 1
                j = jg2*(jg-jg1) + 1
                do jjj=0,jg2-1
                  jff = j + jjj
                  if (jff < mpi_tt(ia,it)%low_j .or.                           &
                      jff > mpi_tt(ia,it)%up_j       ) cycle
                  do iii=0,ig2-1
                    iff = i + iii
                    if (iff < mpi_tt(ia,it)%low_i .or.                         &
                        iff > mpi_tt(ia,it)%up_i       ) cycle
                    if (mmk(1,iff,jff) < 1) cycle
                    boff = boff + 1
                    a(mmk(1,iff,jff)) = bbbb(boff)
                    if (kmx >= 2 .and. khg(mmk(1,iff,jff)) >= 2) then
                      mi0g = mmk(2,iff,jff) - 2
                      do ik=2,min(kmx,khg(mmk(1,iff,jff)))
                        boff = boff + 1
                        a(mi0g+ik) = bbbb(boff)
                      enddo  ! ki
                    endif
                    if (b2) then
                      boff = boff + 1
                      b(mmk(1,iff,jff)) = bbbb(boff)
                      if (kmx >= 2 .and. khg(mmk(1,iff,jff)) >= 2) then
                        mi0g = mmk(2,iff,jff) - 2
                        do ik=2,min(kmx,khg(mmk(1,iff,jff)))
                          boff = boff + 1
                          b(mi0g+ik) = bbbb(boff)
                        enddo  ! ki
                      endif
                      if (b3) then
                        boff = boff + 1
                        c(mmk(1,iff,jff)) = bbbb(boff)
                        if (kmx >= 2 .and. khg(mmk(1,iff,jff)) >= 2) then
                          mi0g = mmk(2,iff,jff) - 2
                          do ik=2,min(kmx,khg(mmk(1,iff,jff)))
                            boff = boff + 1
                            c(mi0g+ik) = bbbb(boff)
                          enddo  ! ki
                        endif
                      endif
                    endif
                  enddo    ! iii
                enddo      ! jjj
              enddo        ! ig
            enddo          ! jg
            deallocate(bbbb)
          endif
        enddo  ! it

        ! copy from local to global array
        if (send) then
          do jg=jglo,jghi
            do ig=max(is,il),min(ie,iu)
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(ia)%low_j .or. jff > dd(ia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(ia)%low_i .or. iff > dd(ia)%up_i) cycle
                  if (mmk(1,iff,jff) < 1) cycle
                  a(mmk(1,iff,jff)) = a_l(mmk_l(1,iff,jff))
                  if (kmx >= 2 .and. khg(mmk(1,iff,jff)) >= 2) then
                    mi0g = mmk(2,iff,jff)   - 2
                    mi0l = mmk_l(2,iff,jff) - 2
                    do ik=2,min(kmx,khg(mmk(1,iff,jff)))
                      a(mi0g+ik) = a_l(mi0l+ik)
                    enddo  ! ki
                  endif
                  if (b2) then
                    b(mmk(1,iff,jff)) = b_l(mmk_l(1,iff,jff))
                    if (kmx >= 2 .and. khg(mmk(1,iff,jff)) >= 2) then
                      mi0g = mmk(2,iff,jff)   - 2
                      mi0l = mmk_l(2,iff,jff) - 2
                      do ik=2,min(kmx,khg(mmk(1,iff,jff)))
                        b(mi0g+ik) = b_l(mi0l+ik)
                      enddo  ! ki
                    endif
                    if (b3) then
                      c(mmk(1,iff,jff)) = c_l(mmk_l(1,iff,jff))
                      if (kmx >= 2 .and. khg(mmk(1,iff,jff)) >= 2) then
                        mi0g = mmk(2,iff,jff)   - 2
                        mi0l = mmk_l(2,iff,jff) - 2
                        do ik=2,min(kmx,khg(mmk(1,iff,jff)))
                          c(mi0g+ik) = c_l(mi0l+ik)
                        enddo  ! ki
                      endif
                    endif
                  endif
                enddo    ! iii
              enddo      ! jjj
            enddo        ! ig
          enddo          ! jg
        endif ! send


      !-  if I have data, send it to the relevant task -------------------------
      elseif (ittt /= iam .and. send) then
        ! send but not to myself:
        do it=ittt,ittt
          if (it == iam) cycle
          !  coarse grid bounds:
          il = mpi_tt(iia,it)%low_i
          iu = mpi_tt(iia,it)%up_i
          jl = mpi_tt(iia,it)%low_j
          ju = mpi_tt(iia,it)%up_j
          call domp_get_domain(max(js,jl), min(je,ju), jglo, jghi)

          !  encode and send:
          bsiz3(tnum,it,1:2,1) = 0
          do jg=jglo,jghi
            do ig=max(is,il),min(ie,iu)
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=0,jg2-1
                jff = j + jjj
                if (jff < dd(ia)%low_j .or. jff > dd(ia)%up_j) cycle
                do iii=0,ig2-1
                  iff = i + iii
                  if (iff < dd(ia)%low_i .or. iff > dd(ia)%up_i) cycle
                  if (khl(mmk_l(1,iff,jff)) >= 1)                              &
                    bsiz3(tnum,it,1,1) = bsiz3(tnum,it,1,1)                    &
                                       + min(kmx,khl(mmk_l(1,iff,jff)))
                enddo    ! iii
              enddo      ! jjj
            enddo        ! ig
          enddo          ! jg
          if (b2) then
            if (b3) then
              bsiz3(tnum,it,1,1) = 3*bsiz3(tnum,it,1,1)
            else
              bsiz3(tnum,it,1,1) = 2*bsiz3(tnum,it,1,1)
            endif
          endif

          if (bsiz3(tnum,it,1,1) > 0) then
            allocate(bbbb(bsiz3(tnum,it,1,1)))
            boff = 0
            do jg=jglo,jghi
              do ig=max(is,il),min(ie,iu)
                if (mmk_c(1,ig,jg) <= 0) cycle
                i = ig2*(ig-ig1) + 1
                j = jg2*(jg-jg1) + 1
                do jjj=0,jg2-1
                  jff = j + jjj
                  if (jff < dd(ia)%low_j .or. jff > dd(ia)%up_j) cycle
                  do iii=0,ig2-1
                    iff = i + iii
                    if (iff < dd(ia)%low_i .or. iff > dd(ia)%up_i) cycle
                    if (mmk_l(1,iff,jff) < 1) cycle
                    boff       = boff + 1
                    bbbb(boff) = a_l(mmk_l(1,iff,jff))
                    if (kmx >= 2 .and. khl(mmk_l(1,iff,jff)) >= 2) then 
                      mi0l = mmk_l(2,iff,jff) - 2
                      do ik=2,min(kmx,khl(mmk_l(1,iff,jff)))
                        boff       = boff + 1
                        bbbb(boff) = a_l(mi0l+ik)
                      enddo  ! ki
                    endif
                    if (b2) then
                      boff       = boff + 1
                      bbbb(boff) = b_l(mmk_l(1,iff,jff))
                      if (kmx >= 2 .and. khl(mmk_l(1,iff,jff)) >= 2) then 
                        mi0l = mmk_l(2,iff,jff) - 2
                        do ik=2,min(kmx,khl(mmk_l(1,iff,jff)))
                          boff       = boff + 1
                          bbbb(boff) = b_l(mi0l+ik)
                        enddo  ! ki
                      endif
                      if (b3) then
                        boff       = boff + 1
                        bbbb(boff) = c_l(mmk_l(1,iff,jff))
                        if (kmx >= 2 .and. khl(mmk_l(1,iff,jff)) >= 2) then 
                          mi0l = mmk_l(2,iff,jff) - 2
                          do ik=2,min(kmx,khl(mmk_l(1,iff,jff)))
                            boff       = boff + 1
                            bbbb(boff) = c_l(mi0l+ik)
                          enddo  ! ki
                        endif
                      endif
                    endif
                  enddo    ! iii
                enddo      ! jjj
              enddo        ! ig
            enddo          ! jg
            call dmpi_send(bbbb, boff, it-1, itag)
            deallocate(bbbb)
          endif
        enddo 

      endif 
    enddo   ! ittt

  end subroutine dmpi_gather_copy

!===============================================================================

  subroutine dmpi_scatter_ice(ia,iw,mmk_l,a_l,b_l,c_l,mmk,a,b,c,nc,send,l_l,l)
    ! Do a scatter from a global array on task #send to local arrays on all 
    ! tasks using broadcasting of the global array.
    ! FIXME: consider re-coding this to send/recv only the portions needed for
    !        the local arrays instead of full global arrays.
    
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, mmk(:,0:,0:), mmk_l(:,0:,0:), nc, send, iw
    real(8),    intent(inout) :: a_l(0:), b_l(0:), c_l(:,0:)
    real(8),    intent(inout) :: a(0:), b(0:), c(:,0:)
    logical,    intent(inout), optional :: l_l(0:), l(0:)

    integer(4) :: ic, i, j, ilh, iuh, jlh, juh, ms, ns, ii, n2d, nl, nu
    logical    :: dol

    !- quit if there's nothing to do here --------------------------------------
    if (iw < 0) return
    if ((send >= mpi_size) .or. (send < 0)) then
      call exitme(1,'Called dmpi_scatter with invalid send parameter')
    endif
    dol = (present(l) .and. present(l_l)) 

    !- broadcast the global arrays ---------------------------------------------
    if (mpi_size > 1) then
      if (mpi_rank == send) then
        !  encode:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, ii, ic)
        call domp_get_domain(0, iw, nl, nu)
        if (dol) lbuf(nl+1:nu+1) = l(nl:nu)
        gbuf(nl   +1:nu   +1) = a(nl:nu)
        gbuf(nl+iw+2:nu+iw+2) = b(nl:nu)
        if (nc > 0) then
          ii = nl*nc + 2*iw + 2
          do i=nl,nu
            do ic=1,nc
              ii = ii + 1
              gbuf(ii) = c(ic,i)
            enddo
          enddo
        endif
!$OMP END PARALLEL
      endif
      call dmpi_bcast( gbuf, (nc+2)*(iw+1), send )
      if (dol) call dmpi_bcast( lbuf, (iw+1), send )
      if (mpi_rank /= send) then
        !  decode:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, ii, ic)
        call domp_get_domain(0, iw, nl, nu)
        if (dol) l(nl:nu) = lbuf(nl+1:nu+1)
        a(nl:nu) = gbuf(nl   +1:nu   +1)
        b(nl:nu) = gbuf(nl+iw+2:nu+iw+2)
        if (nc > 0) then
          ii = nl*nc + 2*iw + 2
          do i=nl,nu
            do ic=1,nc
              ii = ii + 1
              c(ic,i) = gbuf(ii)
            enddo
          enddo
        endif
!$OMP END PARALLEL
      endif
    endif

    !- transfer to local arrays ------------------------------------------------
    !
    !- quit, if we have no more jobs to do here --------------------------------
    if (dd(ia)%nwet3 == 0) return

    !- transfer to local arrays incl. halo zone --------------------------------
    !  Constraint: Must comply with the index order which was set up in
    !              task_local_arrays.
    !  We do that in this way: First the inner task region, then the halo zone
    !              in this order:    n/w, w, s/w, n/e, e, s/e, n, s
    !              For simplicity, the halo width is here assumed to be 1.
    !
    !  Most things below here is stride-1, except the halo-surface which is only
    !  stride-1  wrt a_l(), not for a(), but this is difficult to take advantage
    !  of so performance for those four loops is expected to suck.
    !
    !  We could benefit from doing a couple of loop-fusions.
    !
    !
    ! ... first, the task without halo:
    n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
    ii  = dd(ia)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
    call domp_get_domain(1, n2d, nl, nu)
    if (dol) l_l(nl:nu) = l(ii+nl:ii+nu)
    a_l(nl:nu) = a(ii+nl:ii+nu)
    b_l(nl:nu) = b(ii+nl:ii+nu)
    c_l(1:nc,nl:nu) = c(1:nc,ii+nl:ii+nu)
!$OMP END PARALLEL
    !
    ! ... then, the halo:
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj
    !
    !  ns: index according to local permutation
    !  ms: index according to global permutation
    !
    j = jlh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      if (dol) l_l(ns) = l(ms)
      a_l(ns) = a(ms)
      b_l(ns) = b(ms)
      c_l(1:nc,ns) = c(1:nc,ms)
    enddo
    j = juh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      if (dol) l_l(ns) = l(ms)
      a_l(ns) = a(ms)
      b_l(ns) = b(ms)
      c_l(1:nc,ns) = c(1:nc,ms)
    enddo
    i = ilh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      if (dol) l_l(ns) = l(ms)
      a_l(ns) = a(ms)
      b_l(ns) = b(ms)
      c_l(1:nc,ns) = c(1:nc,ms)
    enddo
    i = iuh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      if (dol) l_l(ns) = l(ms)
      a_l(ns) = a(ms)
      b_l(ns) = b(ms)
      c_l(1:nc,ns) = c(1:nc,ms)
    enddo

  end subroutine dmpi_scatter_ice

!===============================================================================

  subroutine dmpi_scatter_ice2(ia, iw, mmk_l, l_l, t_l, ts_l, c_l,             &
                               mmk, l, t, ts, c, nc, send)
    ! Scatter logical array and component array (ment for ice vars)
    ! Do a scatter from a global array on task #send to local arrays on all 
    ! tasks using broadcasting of the global array.
    ! FIXME: consider re-coding this to send/recv only the portions needed for
    !        the local arrays instead of full global arrays.
    
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, mmk(:,0:,0:), mmk_l(:,0:,0:), nc, send, iw
    logical,    intent(inout) :: l_l(0:)
    real(8),    intent(inout) :: c_l(:,0:), t_l(:,0:), ts_l(0:)
    logical,    intent(inout) :: l(0:)
    real(8),    intent(inout) :: c(:,0:), t(:,0:), ts(0:)

    integer(4) :: ic, i, j, ilh, iuh, jlh, juh, ms, ns, ii, n2d, nl, nu

    !- quit if there's nothing to do here --------------------------------------
    if (iw < 0) return
    if ((send >= mpi_size) .or. (send < 0)) then
      call exitme(1,'Called dmpi_scatter with invalid send parameter')
    endif

    !- broadcast the global arrays ---------------------------------------------
    if (mpi_size > 1) then
      if (mpi_rank == send) then
        !  encode:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, ii, ic)
        call domp_get_domain(0, iw, nl, nu)
        lbuf(nl   +1:nu   +1) = l(  nl:nu)
        gbuf(nl   +1:nu   +1) = ts( nl:nu)
        gbuf(nl+iw+2:nu+iw+2) = t(1,nl:nu)
        if (nc > 0) then
          ii = nl*nc + 2*iw + 2
          do i=nl,nu
            do ic=1,nc
              ii = ii + 1
              gbuf(ii) = c(ic,i)
            enddo
          enddo
        endif
!$OMP END PARALLEL
      endif
      call dmpi_bcast( gbuf, (nc+2)*(iw+1), send )
      call dmpi_bcast( lbuf, (iw+1), send )
      if (mpi_rank /= send) then
        !  decode:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, ii, ic)
        call domp_get_domain(0, iw, nl, nu)
        l(  nl:nu) = lbuf(nl   +1:nu   +1)
        ts( nl:nu) = gbuf(nl   +1:nu   +1)
        t(1,nl:nu) = gbuf(nl+iw+2:nu+iw+2)
        if (nc > 0) then
          ii = nl*nc + 2*iw + 2
          do i=nl,nu
            do ic=1,nc
              ii = ii + 1
              c(ic,i) = gbuf(ii)
            enddo
          enddo
        endif
!$OMP END PARALLEL
      endif
    endif

    !- transfer to local arrays ------------------------------------------------
    !
    !- quit, if we have no more jobs to do here --------------------------------
    if (dd(ia)%nwet3 == 0) return

    !- transfer to local arrays incl. halo zone --------------------------------
    !  Constraint: Must comply with the index order which was set up in
    !              task_local_arrays.
    !  We do that in this way: First the inner task region, then the halo zone
    !              in this order:    n/w, w, s/w, n/e, e, s/e, n, s
    !              For simplicity, the halo width is here assumed to be 1.
    !
    !  Most things below here is stride-1, except the halo-surface which is only
    !  stride-1  wrt a_l(), not for a(), but this is difficult to take advantage
    !  of so performance for those four loops is expected to suck.
    !
    !  We could benefit from doing a couple of loop-fusions.
    !
    !
    ! ... first, the task without halo:
    n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
    ii  = dd(ia)%low_ws - 1
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu)
    call domp_get_domain(1, n2d, nl, nu)
    l_l(nl:nu) = l(ii+nl:ii+nu)
    if (nc > 0) c_l(1:nc,nl:nu) = c(1:nc,ii+nl:ii+nu)
    t_l(1,nl:nu) = t(1,ii+nl:ii+nu)
    ts_l(nl:nu)  = ts(ii+nl:ii+nu)
!$OMP END PARALLEL
    !
    ! ... then, the halo:
    !             (we need not the halo of ts)
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj
    !
    !  ns: index according to local permutation
    !  ms: index according to global permutation
    !
    j = jlh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      l_l(ns) = l(ms)
      if (nc > 0) c_l(1:nc,ns) = c(1:nc,ms)
      t_l(1,ns) = t(1,ms)
    enddo
    j = juh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      l_l(ns) = l(ms)
      if (nc > 0) c_l(1:nc,ns) = c(1:nc,ms)
      t_l(1,ns) = t(1,ms)
    enddo
    i = ilh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      l_l(ns) = l(ms)
      if (nc > 0) c_l(1:nc,ns) = c(1:nc,ms)
      t_l(1,ns) = t(1,ms)
    enddo
    i = iuh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      ms = mmk(1,i,j)
      l_l(ns) = l(ms)
      if (nc > 0) c_l(1:nc,ns) = c(1:nc,ms)
      t_l(1,ns) = t(1,ms)
    enddo

  end subroutine dmpi_scatter_ice2

!===============================================================================

  subroutine dmpi_scatter_cmp_table(ia,iw,mmk_l,c_l,kh,mmk,c,nc,send,scat)
    ! Do a scatter from a global array on task #send to local arrays on all 
    ! tasks using send/recv of the global array, using the scatter table.
    ! FIXME: consider re-coding this to send/recv only the portions needed for
    !        the local arrays instead of full global arrays.
    
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, nc, send, iw
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), kh(0:)
    real(8),    intent(inout) :: c_l(:,0:)
    real(8),    intent(inout) :: c(:,0:)
    logical,    intent(in)    :: scat(:)

    integer(4) :: ic, i, j, ilh, iuh, jlh, juh, ms, ns, k, ii, iam, itag,it_recv
    integer(4) :: nl, nu

    !- quit if there's nothing to do here --------------------------------------
    if (iw < 0) return
    if ((send >= mpi_size) .or. (send < 0)) then
      call exitme(1,'Called dmpi_scatter with invalid send parameter')
    endif
    iam = mpi_rank + 1
    if ((.not.scat(iam)) .and. (mpi_rank /= send)) return

    !- scatter the global array ------------------------------------------------
    if (mpi_size > 1) then

      if (send /= mpi_rank) then
        !  I'm not the sender, I should possibly receive from the sender -------
        if (scat(iam)) then
          !  receive buffer:
          itag = 1
          call dmpi_recv( gbuf, nc*(iw+1), send, itag )

          !  decode buffer:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ii, ic, i)
          call domp_get_domain(0, iw, nl, nu)
          ii = nl*nc
          do i=nl,nu
            do ic=1,nc
              ii = ii + 1
              c(ic,i) = gbuf(ii)
            enddo
          enddo
!$OMP END PARALLEL
        endif

      else
        !  I am the sender and I must send -------------------------------------

        !  fill buffer:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, ii, ic, i)
        call domp_get_domain(0, iw, nl, nu)
        ii = nl*nc
        do i=nl,nu
          do ic=1,nc
            ii = ii + 1
            gbuf(ii) = c(ic,i)
          enddo
        enddo
!$OMP END PARALLEL

        do it_recv=0,mpi_size-1
          !  don't send to myself ...
          if (it_recv == send) cycle

          !  ... but to those tasks in the table:
          if (scat(it_recv+1)) then
            !  send buffer:
            itag = 1
            call dmpi_send( gbuf, nc*(iw+1), it_recv, itag )
          endif
        enddo
      endif

    endif

    !- transfer to local arrays ------------------------------------------------
    !
    !  do this task's region incl. halo zone:
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj
    !
    !  ns: index according to local permutation
    !  ms: index according to global permutation
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, j, k, ns, ms, ic)
    call domp_get_domain(jlh, juh, nl, nu)
    do j=nl,nu
      do i=ilh,iuh
        if (mmk_l(1,i,j) <= 0) cycle
        do k=1,kh(mmk(1,i,j))
          ns = mmk_l(k,i,j)
          ms = mmk(k,i,j)
          do ic=1,nc
            c_l(ic,ns) = c(ic,ms)
          enddo
        enddo
      enddo
    enddo
!$OMP END PARALLEL

  end subroutine dmpi_scatter_cmp_table

!===============================================================================

  subroutine dmpi_scatter_cmp(ia,iw,mmk_l,c_l,kh_l,mmk,c,nc,noff,send)
    ! Do a scatter from a global array on task #send to local arrays on all 
    ! tasks using send/recv of the global array, without using a scattertable.
    !
    ! FIXME: consider re-coding this to send/recv only the portions needed for
    !        the local arrays instead of full global arrays.
    
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, nc, noff, send, iw
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), kh_l(0:)
    real(8),    intent(inout) :: c_l(:,0:)
    real(8),    intent(inout) :: c(:,0:)

    integer(4) :: ic, i, j, ilh, iuh, jlh, juh, ns, ii, itag, it_recv, n2d, kb
    integer(4) :: nl, nu

    !- quit if there's nothing to do here --------------------------------------
    if (iw < 0) return
    if ((send >= mpi_size) .or. (send < 0)) then
      call exitme(1,'Called dmpi_scatter with invalid send parameter')
    endif

    !- scatter the global array ------------------------------------------------
    if (mpi_size > 1) then

      if (send /= mpi_rank) then
        !  I'm not the sender, I should receive from the sender ----------------

        !  receive buffer:
        itag = 1
        call dmpi_recv( gbuf, nc*(iw+1), send, itag )

        !  decode buffer:
!$OMP PARALLEL DEFAULT (none)                                                  &
!$OMP SHARED  (iw, nc, noff, c, gbuf)                                          &
!$OMP PRIVATE (ii, i, ic, nl, nu)
        call domp_get_domain(0, iw, nl, nu)
        ii = nl*nc
        do i=nl,nu
          do ic=1+noff,nc+noff
            ii = ii + 1
            c(ic,i) = gbuf(ii)
          enddo
        enddo
!$OMP END PARALLEL

      else
        !  I am the sender and I must send -------------------------------------

        !  fill buffer:
!$OMP PARALLEL DEFAULT (none)                                                  &
!$OMP SHARED  (iw, nc, noff, c, gbuf)                                          &
!$OMP PRIVATE (ii, i, ic, nl, nu)
        call domp_get_domain(0, iw, nl, nu)
        ii = nl*nc
        do i=nl,nu
          do ic=1+noff,nc+noff
            ii = ii + 1
            gbuf(ii) = c(ic,i)
          enddo
        enddo
!$OMP END PARALLEL

        do it_recv=0,mpi_size-1
          !  don't send to myself ...
          if (it_recv == send) cycle

          !  send buffer:
          itag = 1
          call dmpi_send( gbuf, nc*(iw+1), it_recv, itag )
        enddo
      endif

    endif

    !- transfer to local arrays ------------------------------------------------
    !
    !- quit, if we have no more jobs to do here --------------------------------
    if (dd(ia)%nwet3 == 0) return

    !- transfer to local arrays incl. halo zone --------------------------------
    !  Constraint: Must comply with the index order which was set up in
    !              task_local_arrays.
    !  We do that in this way: First the inner task region, then the halo zone
    !              in this order:    n/w, w, s/w, n/e, e, s/e, n, s
    !              For simplicity, the halo width is here assumed to be 1.
    !
    !  Most things below here is stride-1, except the halo-surface which is only
    !  stride-1  wrt a_l(), not for a(), but this is difficult to take advantage
    !  of so performance for those four loops is expected to suck.
    !
    !  We could benefit from doing a couple of loop-fusions.
    !
    !
    ! ... first, the task without halo:
    n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
    j   = dd(ia)%nwet3 - (n2d+1)
    ii  = (n2d+1) + dd(ia)%halo2
!$OMP PARALLEL DEFAULT (shared) PRIVATE (i, nl, nu)
    call domp_get_domain(1, n2d, nl, nu)
    i = dd(ia)%low_ws + nl - 1
    c_l(noff+1:noff+nc,nl:nu) = c(noff+1:noff+nc,i:i+nu-nl)
    call domp_get_domain(ii, ii+j, nl, nu)
    i   = dd(ia)%low_w3 + nl - ii
    c_l(noff+1:noff+nc,nl:nu) = c(noff+1:noff+nc,i:i+nu-nl)
!$OMP END PARALLEL
    !
    ! ... then, the halo:
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj
    !
    !  ns: index according to local permutation
    !
    !  do the surface array:
    j = jlh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      c_l(noff+1:noff+nc,ns) = c(noff+1:noff+nc,mmk(1,i,j))
    enddo
    j = juh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      c_l(noff+1:noff+nc,ns) = c(noff+1:noff+nc,mmk(1,i,j))
    enddo
    i = ilh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      c_l(noff+1:noff+nc,ns) = c(noff+1:noff+nc,mmk(1,i,j))
    enddo
    i = iuh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      c_l(noff+1:noff+nc,ns) = c(noff+1:noff+nc,mmk(1,i,j))
    enddo
    !
    !  do the sub-surface array:
    j = jlh
    do i=ilh,iuh
      kb = kh_l(mmk_l(1,i,j))
      if (kb < 2) cycle
      c_l(noff+1:noff+nc,mmk_l(2,i,j):mmk_l(kb,i,j)) =                         &
                                       c(noff+1:noff+nc,mmk(2,i,j):mmk(kb,i,j))
    enddo
    j = juh
    do i=ilh,iuh
      kb = kh_l(mmk_l(1,i,j))
      if (kb < 2) cycle
      c_l(noff+1:noff+nc,mmk_l(2,i,j):mmk_l(kb,i,j)) =                         &
                                       c(noff+1:noff+nc,mmk(2,i,j):mmk(kb,i,j))
    enddo
    i = ilh
    do j=jlh+1,juh-1
      kb = kh_l(mmk_l(1,i,j))
      if (kb < 2) cycle
      c_l(noff+1:noff+nc,mmk_l(2,i,j):mmk_l(kb,i,j)) =                         &
                                       c(noff+1:noff+nc,mmk(2,i,j):mmk(kb,i,j))
    enddo
    i = iuh
    do j=jlh+1,juh-1
      kb = kh_l(mmk_l(1,i,j))
      if (kb < 2) cycle
      c_l(noff+1:noff+nc,mmk_l(2,i,j):mmk_l(kb,i,j)) =                         &
                                       c(noff+1:noff+nc,mmk(2,i,j):mmk(kb,i,j))
    enddo

  end subroutine dmpi_scatter_cmp

!===============================================================================

  subroutine dmpi_scatter_arr_table(ia,kmx,iw,mmk_l,a_l,mmk,a,send,scat)
    ! Do a scatter from a global array on task #send to local arrays on all 
    ! tasks using send/recv of the global array, using the scatter table.
    ! FIXME: consider re-coding this to send/recv only the portions needed for
    !        the local arrays instead of full global arrays.
    
    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, kmx, send, iw
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:)
    real(8),    intent(inout) :: a_l(0:)
    real(8),    intent(inout) :: a(0:)
    logical,    intent(in)    :: scat(:)

    integer(4) :: i, j, ilh, iuh, jlh, juh, ms, ns, k, itag, iam, it_recv
    integer(4) :: nl, nu

    !- quit if there's nothing to do here --------------------------------------
    if (iw < 0) return
    if ((send >= mpi_size) .or. (send < 0)) then
      call exitme(1,'Called dmpi_scatter with invalid send parameter')
    endif
    iam = mpi_rank + 1
    if ((.not.scat(iam)) .and. (mpi_rank /= send)) return

    !- scatter the global array ------------------------------------------------
    if (mpi_size > 1) then

      if (send /= mpi_rank) then
        !  I'm not the sender, I should possibly receive from the sender -------
        if (scat(iam)) then
          !  receive buffer:
          itag = 1
          call dmpi_recv( a, (iw+1), send, itag )
        endif

      else
        !  I am the sender and I must send -------------------------------------
        do it_recv=0,mpi_size-1
          !  don't send to myself ...
          if (it_recv == send) cycle
          !  ... but to those tasks in the table:
          if (scat(it_recv+1)) then
            !  send buffer:
            itag = 1
            call dmpi_send( a, (iw+1), it_recv, itag )
          endif
        enddo
      endif

    endif


    !- transfer to local arrays ------------------------------------------------
    !
    !  do this task's region incl. halo zone:
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj
    !
    !  ns: index according to local permutation
    !  ms: index according to global permutation
    if (kmx <= 1) then
      !  do a surface array:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, j, ns)
      call domp_get_domain(jlh, juh, nl, nu)
      do j=nl,nu
        do i=ilh,iuh
          ns = mmk_l(1,i,j)
          if (ns <= 0) cycle
          a_l(ns) = a(mmk(1,i,j))
        enddo
      enddo
!$OMP END PARALLEL

    else
      !  do a full 3D array:
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i, j, k, ns, ms)
      call domp_get_domain(jlh, juh, nl, nu)
      do j=nl,nu
        do i=ilh,iuh
          if (mmk_l(1,i,j) <= 0) cycle
          do k=1,kmx
            ns = mmk_l(k,i,j)
            if (ns <= 0) exit
            ms = mmk(k,i,j)
            a_l(ns) = a(ms)
          enddo
        enddo
      enddo
!$OMP END PARALLEL
    endif

  end subroutine dmpi_scatter_arr_table

!===============================================================================

  subroutine dmpi_scatter_arr(ia,kmx,iw,mmk_l,a_l,mmk,a,send,kh_l)
    ! Do a scatter from a global array on task #send to local arrays on all
    ! tasks using send/recv of the global array
    ! FIXME: consider re-coding this to send/recv only the portions needed for
    !        the local arrays instead of full global arrays.

    use dmi_omp, only : domp_get_domain
    implicit none
    integer(4), intent(in)    :: ia, kmx, send, iw
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), kh_l(0:)
    real(8),    intent(inout) :: a_l(0:)
    real(8),    intent(inout) :: a(0:)

    integer(4) :: i, j, ilh, iuh, jlh, juh, ns, itag, it_recv, kb, ii, n2d
    integer(4) :: nl, nu

    !- quit if there's nothing to do here --------------------------------------
    if (iw < 0) return
    if ((send >= mpi_size) .or. (send < 0)) then
      call exitme(1,'Called dmpi_scatter with invalid send parameter')
    endif

    !- scatter the global array ------------------------------------------------
    if (mpi_size > 1) then

      if (send /= mpi_rank) then
        !  I'm not the sender, I should receive from the sender

        !  receive buffer:
        itag = 1
        call dmpi_recv( a, (iw+1), send, itag )

      else
        !  I am the sender and I must send
        do it_recv=0,mpi_size-1
          !  don't send to myself ...
          if (it_recv == send) cycle

          !  send buffer:
          itag = 1
          call dmpi_send( a, (iw+1), it_recv, itag )
        enddo
      endif

    endif

    !- quit, if we have no more jobs to do here --------------------------------
    if (dd(ia)%nwet3 == 0) return

    !- transfer to local arrays incl. halo zone --------------------------------
    !  Constraint: Must comply with the index order which was set up in
    !              task_local_arrays.
    !  We do that in this way: First the inner task region, then the halo zone
    !              in this order:    n/w, w, s/w, n/e, e, s/e, n, s
    !              For simplicity, the halo width is here assumed to be 1.
    !
    !  Most things below here is stride-1, except the halo-surface which is only
    !  stride-1  wrt a_l(), not for a(), but this is difficult to take advantage
    !  of so performance for those four loops is expected to suck.
    !
    !  We could benefit from doing a couple of loop-fusions.
    !
    !
    ! ... first, the task without halo:
    n2d = dd(ia)%up_ws - dd(ia)%low_ws + 1
    j   = dd(ia)%nwet3 - (n2d+1)
    ii  = (n2d+1) + dd(ia)%halo2
!$OMP PARALLEL DEFAULT (shared) PRIVATE (nl, nu, i)
    call domp_get_domain(1, n2d, nl, nu)
    i = dd(ia)%low_ws + nl - 1
    a_l(nl:nu) = a(i:i+nu-nl)

    if (kmx > 1) then
      call domp_get_domain(ii, ii+j, nl, nu)
      i = dd(ia)%low_w3 + nl - ii
      a_l(nl:nu) = a(i:i+nu-nl)
    endif
!$OMP END PARALLEL
    !
    ! ... then, the halo:
    ilh = dd(ia)%low_hi
    iuh = dd(ia)%up_hi
    jlh = dd(ia)%low_hj
    juh = dd(ia)%up_hj
    !
    !  ns: index according to local permutation
    !
    !  do the surface array:
    j = jlh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      a_l(ns) = a(mmk(1,i,j))  !<== actually, this IS stride-1
    enddo
    j = juh
    do i=ilh,iuh
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      a_l(ns) = a(mmk(1,i,j))  !<== actually, this IS stride-1
    enddo
    i = ilh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      a_l(ns) = a(mmk(1,i,j))  !<== bad stride
    enddo
    i = iuh
    do j=jlh+1,juh-1
      ns = mmk_l(1,i,j)
      if (ns <= 0) cycle
      a_l(ns) = a(mmk(1,i,j))  !<== bad stride
    enddo
    !
    !  do the sub-surface array:
    if (kmx > 1) then
      j = jlh
      do i=ilh,iuh
        kb = kh_l(mmk_l(1,i,j))
        if (kb < 2) cycle
        a_l(mmk_l(2,i,j):mmk_l(kb,i,j)) = a(mmk(2,i,j):mmk(kb,i,j))
      enddo
      j = juh
      do i=ilh,iuh
        kb = kh_l(mmk_l(1,i,j))
        if (kb < 2) cycle
        a_l(mmk_l(2,i,j):mmk_l(kb,i,j)) = a(mmk(2,i,j):mmk(kb,i,j))
      enddo
      i = ilh
      do j=jlh+1,juh-1
        kb = kh_l(mmk_l(1,i,j))
        if (kb < 2) cycle
        a_l(mmk_l(2,i,j):mmk_l(kb,i,j)) = a(mmk(2,i,j):mmk(kb,i,j))
      enddo
      i = iuh
      do j=jlh+1,juh-1
        kb = kh_l(mmk_l(1,i,j))
        if (kb < 2) cycle
        a_l(mmk_l(2,i,j):mmk_l(kb,i,j)) = a(mmk(2,i,j):mmk(kb,i,j))
      enddo
    endif

  end subroutine dmpi_scatter_arr

!===============================================================================

  subroutine dmpi_mcf( n, n2, kr, ig1, jg1, ig2, jg2, il, iu, jl, ju, kh, kmx, &
                       func, b_size, mmkc_l, uc_l, vc_l, buf, mmkc, uc, vc)

    implicit none

    integer(4), intent(in)    :: n, n2, kr(:,0:), ig1, jg1, ig2, jg2
    integer(4), intent(in)    :: il, iu, jl, ju, kh(0:), kmx
    integer(4), intent(in)    :: func   ! ==  1: calc buffer size
                                        !     2: decode buffer (buf --> uc/vc)
                                        !     3: encode buffer (uc/vc --> buf)
                                        !     4: copy local to global
    integer(4), intent(inout) :: b_size

    integer(4), intent(in)              :: mmkc_l(:,0:,0:)
    integer(4), intent(in),    optional :: mmkc(:,0:,0:)
    real(8),    intent(in),    optional :: uc_l(0:), vc_l(0:)
    real(8),    intent(inout), optional :: buf(:), uc(0:), vc(0:)
    

    integer(4) :: bordertype, i, j, ncorner, ig, jg, kc, ign, igs, jgw, jge

    ! further documentation, cf mom_c_f_default()

    ! border point index (i,j) in fine grid
    bordertype = kr(3,n)
    if (bordertype == 0) return
    i = kr(1,n)
    j = kr(2,n)

    ! only do corner points once
    if (n > 1) then
      if (kr(3,n-1) > 0 .and. i == kr(1,n-1) .and. j == kr(2,n-1)) return
    endif
    ncorner = 0
    if (n < n2) then
      if (kr(3,n+1) > 0 .and. i == kr(1,n+1) .and. j == kr(2,n+1)) then
        ncorner = bordertype
      endif
    endif
  
    ! border point index (ig,jg) in coarse grid 
    ig = ig1 + int((i-1)/ig2,4)
    jg = jg1 + int((j-1)/jg2,4)

    if (bordertype == 1) then
      ! West border 
      !   extract U(west) from coarse grid.
      jg  = jg - 1
      ign = ig1 + int((i-2)/ig2,4)
      igs = ig1 + int((i  )/ig2,4)
      if (jl <= jg .and. jg <= ju) then
        if (il <= ig  .and. ig  <= iu) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig, jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,ig, jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ig, jg))
            if (func == 4) uc(mmkc(kc,ig, jg)) = uc_l(mmkc_l(kc,ig, jg))
          enddo
        endif
        if (il <= ign .and. ign <= iu) then
          do kc=1,min(kmx,kh(mmkc_l(1,ign,jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,ign,jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ign,jg))
            if (func == 4) uc(mmkc(kc,ign,jg)) = uc_l(mmkc_l(kc,ign,jg))
          enddo
        endif
        if (il <= igs .and. igs <= iu) then
          do kc=1,min(kmx,kh(mmkc_l(1,igs,jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,igs,jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,igs,jg))
            if (func == 4) uc(mmkc(kc,igs,jg)) = uc_l(mmkc_l(kc,igs,jg))
          enddo
        endif
      endif
      if (ncorner == bordertype) then
        ! N/W corner point
        !   extract V(north) from coarse grid.
        if (il <= ig-1 .and. ig-1 <= iu) then
          if (jl <= jg   .and. jg   <= ju) then
            do kc=1,min(kmx,kh(mmkc_l(1,ig-1,jg  )))
              b_size = b_size + 1
              if (func == 2) vc(mmkc(kc,ig-1,jg  )) = buf(b_size)
              if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig-1,jg  ))
              if (func == 4) vc(mmkc(kc,ig-1,jg  )) = vc_l(mmkc_l(kc,ig-1,jg  ))
            enddo
          endif
          if (jl <= jg+1 .and. jg+1 <= ju) then
            do kc=1,min(kmx,kh(mmkc_l(1,ig-1,jg+1)))
              b_size = b_size + 1
              if (func == 2) vc(mmkc(kc,ig-1,jg+1)) = buf(b_size)
              if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig-1,jg+1))
              if (func == 4) vc(mmkc(kc,ig-1,jg+1)) = vc_l(mmkc_l(kc,ig-1,jg+1))
            enddo
          endif
        endif
      endif

    elseif (bordertype == 2) then
      ! North border
      !   extract V(north) from coarse grid.
      ig  = ig - 1
      jgw = jg1 + (j-2)/jg2
      jge = jg1 + (j  )/jg2
      if (il <= ig .and. ig <= iu) then
        if (jl <= jg  .and. jg  <= ju) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jg )))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jg )) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig, jg))
            if (func == 4) vc(mmkc(kc,ig,jg )) = vc_l(mmkc_l(kc,ig,jg ))
          enddo
        endif
        if (jl <= jgw .and. jgw <= ju) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jgw)))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jgw)) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig,jgw))
            if (func == 4) vc(mmkc(kc,ig,jgw)) = vc_l(mmkc_l(kc,ig,jgw))
          enddo
        endif
        if (jl <= jge .and. jge <= ju) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jge)))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jge)) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig,jge))
            if (func == 4) vc(mmkc(kc,ig,jge)) = vc_l(mmkc_l(kc,ig,jge))
          enddo
        endif
      endif
      if (ncorner == bordertype) then
        ! N/E corner point
        !   extract U(east) from coarse grid.
        if ((jl <= jg .and. jg <= ju) .and. (il <= ig+1 .and. ig+1 <= iu)) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig+1,jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,ig+1,jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ig+1,jg))
            if (func == 4) uc(mmkc(kc,ig+1,jg)) = uc_l(mmkc_l(kc,ig+1,jg))
          enddo
        endif
      endif
  
    elseif (bordertype == 3) then
      ! East border 
      !   extract U(east) from coarse grid.
      ign = ig1 + (i-2)/ig2
      igs = ig1 + (i  )/ig2
      if (jl <= jg .and. jg <= ju) then
        if (il <= ig  .and. ig  <= iu) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig, jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,ig, jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ig, jg))
            if (func == 4) uc(mmkc(kc,ig, jg)) = uc_l(mmkc_l(kc,ig, jg))
          enddo
        endif
        if (il <= ign .and. ign <= iu) then
          do kc=1,min(kmx,kh(mmkc_l(1,ign,jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,ign,jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ign,jg))
            if (func == 4) uc(mmkc(kc,ign,jg)) = uc_l(mmkc_l(kc,ign,jg))
          enddo
        endif
        if (il <= igs .and. igs <= iu) then
          do kc=1,min(kmx,kh(mmkc_l(1,igs,jg)))
            b_size = b_size + 1
            if (func == 2) uc(mmkc(kc,igs,jg)) = buf(b_size)
            if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,igs,jg))
            if (func == 4) uc(mmkc(kc,igs,jg)) = uc_l(mmkc_l(kc,igs,jg))
          enddo
        endif
      endif
      if (ncorner == bordertype) then
        ! S/E corner point
        !   extract V(south) from coarse grid.
        if ((il <= ig .and. ig <= iu) .and. (jl <= jg .and. jg <= ju)) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jg)))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jg)) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig,jg))
            if (func == 4) vc(mmkc(kc,ig,jg)) = vc_l(mmkc_l(kc,ig,jg))
          enddo
        endif
      endif
  
    elseif (bordertype == 4) then
      ! South border
      !   extract V(south) from coarse grid.
      jgw = jg1 + (j-2)/jg2
      jge = jg1 + (j  )/jg2
      if (il <= ig .and. ig <= iu) then
        if (jl <= jg  .and. jg  <= ju) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jg )))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jg )) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig,jg ))
            if (func == 4) vc(mmkc(kc,ig,jg )) = vc_l(mmkc_l(kc,ig,jg ))
          enddo
        endif
        if (jl <= jgw .and. jgw <= ju) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jgw)))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jgw)) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig,jgw))
            if (func == 4) vc(mmkc(kc,ig,jgw)) = vc_l(mmkc_l(kc,ig,jgw))
          enddo
        endif
        if (jl <= jge .and. jge <= ju) then
          do kc=1,min(kmx,kh(mmkc_l(1,ig,jge)))
            b_size = b_size + 1
            if (func == 2) vc(mmkc(kc,ig,jge)) = buf(b_size)
            if (func == 3) buf(b_size) = vc_l(mmkc_l(kc,ig,jge))
            if (func == 4) vc(mmkc(kc,ig,jge)) = vc_l(mmkc_l(kc,ig,jge))
          enddo
        endif
      endif
      if (ncorner == bordertype) then
        ! S/W corner 
        !   extract U(west) from coarse grid.
        if (jl <= jg-1 .and. jg-1 <= ju) then
          if (il <= ig   .and. ig   <= iu) then
            do kc=1,min(kmx,kh(mmkc_l(1,ig,  jg-1)))
              b_size = b_size + 1
              if (func == 2) uc(mmkc(kc,ig,  jg-1)) = buf(b_size)
              if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ig,  jg-1))
              if (func == 4) uc(mmkc(kc,ig,  jg  )) = uc_l(mmkc_l(kc,ig,  jg  ))
            enddo
          endif
          if (il <= ig+1 .and. ig+1 <= iu) then
            do kc=1,min(kmx,kh(mmkc_l(1,ig+1,jg-1)))
              b_size = b_size + 1
              if (func == 2) uc(mmkc(kc,ig+1,jg-1)) = buf(b_size)
              if (func == 3) buf(b_size) = uc_l(mmkc_l(kc,ig+1,jg-1))
              if (func == 4) uc(mmkc(kc,ig+1,jg-1)) = uc_l(mmkc_l(kc,ig+1,jg-1))
            enddo
          endif
        endif
      endif

    endif

  end subroutine dmpi_mcf

!===============================================================================

  subroutine dmpi_gather_bnd_3(ia, nz1, nz2, krz, kh, mmk, m_l, a, b, c, bnd)
    implicit none
    integer(4), intent(in)    :: ia, nz1, nz2, krz(:,0:), kh(0:)
    integer(4), intent(in)    :: mmk(1:,0:,0:), m_l(1:,0:,0:)
    real(8),    intent(in)    :: a(0:), b(0:), c(0:)
    real(8),    intent(inout) :: bnd(1:,0:,1:)

    integer(4) :: it, iam, itag, b_size, b_pos, iz, idx, idl, i, j, k
    integer(4), parameter :: na = 3

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  misc vars --------------------------------------------------------------
    itag = 1

    !-  recv on I/O task from all tasks having any data ------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam) then
          !  find buffer size:
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            b_size = b_size + kh(mmk(1,i,j))
          enddo

          !  receive and decode:
          if (b_size > 0) then
            call dmpi_recv(gbuf, na*b_size, it-1, itag)
            b_pos = 0
            do iz=nz1,nz2
              ! skip if this point is not on the relevant task:
              i = krz(1,iz)
              j = krz(2,iz)
              if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.     &
                  j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
              do k=1,kh(mmk(1,i,j))
                bnd(k,iz,1) = gbuf(b_pos+1)
                bnd(k,iz,2) = gbuf(b_pos+2)
                bnd(k,iz,3) = gbuf(b_pos+3)
                b_pos = b_pos + na
              enddo
            enddo
          endif

        else ! it == iam  and we must copy a,b,c,d to bnd-array
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            do k=1,kh(mmk(1,i,j))
              idl = m_l(k,i,j)
              bnd(k,iz,1) = a(idl)
              bnd(k,iz,2) = b(idl)
              bnd(k,iz,3) = c(idl)
            enddo
          enddo

        endif
      enddo

    !-  if I have data, send it to I/O task ------------------------------------
    else
      !  encode and send:
      b_size = 0
      do iz=nz1,nz2
        ! skip if this point is not on this task:
        i = krz(1,iz)
        j = krz(2,iz)
        if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                         &
            j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
        do k=1,kh(mmk(1,i,j))
          idl = m_l(k,i,j)
          gbuf(b_size+1) = a(idl)
          gbuf(b_size+2) = b(idl)
          gbuf(b_size+3) = c(idl)
          b_size = b_size + na
        enddo
      enddo
      if (b_size > 0) call dmpi_send(gbuf, b_size, mpi_io_rank, itag)

    endif

  end subroutine dmpi_gather_bnd_3
  
!===============================================================================

  subroutine dmpi_gather_bnd_2(ia, nz1, nz2, krz, kh, mmk, m_l, a, b, bnd)
    implicit none
    integer(4), intent(in)    :: ia, nz1, nz2, krz(:,0:), kh(0:)
    integer(4), intent(in)    :: mmk(1:,0:,0:), m_l(1:,0:,0:)
    real(8),    intent(in)    :: a(0:), b(0:)
    real(8),    intent(inout) :: bnd(1:,0:,1:)

    integer(4) :: it, iam, itag, b_size, b_pos, iz, idx, idl, i, j, k
    integer(4), parameter :: na = 2

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  misc vars --------------------------------------------------------------
    itag = 1

    !-  recv on I/O task from all tasks having any data ------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam) then
          !  find buffer size:
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            b_size = b_size + kh(mmk(1,i,j))
          enddo

          !  receive and decode:
          if (b_size > 0) then
            call dmpi_recv(gbuf, na*b_size, it-1, itag)
            b_pos = 0
            do iz=nz1,nz2
              ! skip if this point is not on the relevant task:
              i = krz(1,iz)
              j = krz(2,iz)
              if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.     &
                  j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
              do k=1,kh(mmk(1,i,j))
                bnd(k,iz,1) = gbuf(b_pos+1)
                bnd(k,iz,2) = gbuf(b_pos+2)
                b_pos = b_pos + na
              enddo
            enddo
          endif

        else ! it == iam  and we must copy a,b to bnd-array
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            do k=1,kh(mmk(1,i,j))
              idl = m_l(k,i,j)
              bnd(k,iz,1) = a(idl)
              bnd(k,iz,2) = b(idl)
            enddo
          enddo

        endif
      enddo

    !-  if I have data, send it to I/O task ------------------------------------
    else
      !  encode and send:
      b_size = 0
      do iz=nz1,nz2
        ! skip if this point is not on this task:
        i = krz(1,iz)
        j = krz(2,iz)
        if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                         &
            j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
        do k=1,kh(mmk(1,i,j))
          idl = m_l(k,i,j)
          gbuf(b_size+1) = a(idl)
          gbuf(b_size+2) = b(idl)
          b_size = b_size + na
        enddo
      enddo
      if (b_size > 0) call dmpi_send(gbuf, b_size, mpi_io_rank, itag)

    endif

  end subroutine dmpi_gather_bnd_2

!===============================================================================

#if defined (MPI)
  subroutine dmpi_distribute_halo_nb(kmx, mmk, kh, ia,                         &
                                     irbuf, isbuf, itagoff, FLAG,              &
                                     a, b, c, d, e, f)
    implicit none

    integer(4), intent(in)              :: ia, kmx, mmk(:,0:,0:), kh(0:)
    integer(4), intent(in)              :: irbuf, isbuf, itagoff, FLAG
    real(8),    intent(inout)           :: a(0:)
    real(8),    intent(inout), optional :: b(0:), c(0:), d(0:), e(0:), f(0:)


    integer(4) :: itag, it, nb, b_size, id, iam, nbi, idi, ij, afac
    integer(4) :: roff, soff, nreqr, rreq, nreqs, sreq

    iam  = mpi_rank + 1
    itag = NextFreeItag + (itagoff-1)

    afac = 1
    if (present(b)) then
      afac = 2
      if (present(c)) then
        afac = 3
        if (present(d)) then
          afac = 4
          if (present(e)) then
            afac = 5
            if (present(f)) afac = 6
          endif
        endif
      endif
    endif

    ! RECV:
    if (FLAG == 1 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! nb is the sender, iam is the receiver
        nb = mpi_tt(ia,iam)%hnb(id)
        if (kmx > 1) then
          b_size = dd(ia)%recv_bufsize(id)
        else
          b_size = dd(ia)%recv_buf_srf(id)
        endif
        if (b_size > 0) then
          ! recv buffer
!$OMP MASTER
          roff  = (irbuf-1)*b_size + 1
          nreqr = nreqr + 1
          rreq  = (itagoff-1)*max_nb + nreqr
          call dmpi_irecv( tfrbuf(id)%p(roff:), afac*b_size, nb-1, itag,       &
                           ireqrh(rreq))
!$OMP END MASTER
!$OMP BARRIER
        endif
      enddo
    endif

    ! ENCODE:
    if (FLAG == 3 .or. FLAG == 6) then
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                soff = (isbuf-1)*b_size + 1
                ! fill buffer
                if (afac == 1) then
                  call dmpi_encode_buf(iam,kmx,mmk,kh,ia,tfsbuf(idi)%p(soff:), &
                                       it,ij,idi,afac,a)
                elseif (afac == 2) then
                  call dmpi_encode_buf(iam,kmx,mmk,kh,ia,tfsbuf(idi)%p(soff:), &
                                       it,ij,idi,afac,a,b)
                elseif (afac == 3) then
                  call dmpi_encode_buf(iam,kmx,mmk,kh,ia,tfsbuf(idi)%p(soff:), &
                                       it,ij,idi,afac,a,b,c)
                elseif (afac == 4) then
                  call dmpi_encode_buf(iam,kmx,mmk,kh,ia,tfsbuf(idi)%p(soff:), &
                                       it,ij,idi,afac,a,b,c,d)
                elseif (afac == 5) then
                  call dmpi_encode_buf(iam,kmx,mmk,kh,ia,tfsbuf(idi)%p(soff:), &
                                       it,ij,idi,afac,a,b,c,d,e)
                elseif (afac == 6) then
                  call dmpi_encode_buf(iam,kmx,mmk,kh,ia,tfsbuf(idi)%p(soff:), &
                                       it,ij,idi,afac,a,b,c,d,e,f)
                endif
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! SEND:
    if (FLAG == 4 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                ! send buffer
!$OMP MASTER
                soff  = (isbuf-1)*b_size + 1
                nreqs = nreqs + 1
                sreq  = (itagoff-1)*max_nb + nreqs
                call dmpi_isend( tfsbuf(idi)%p(soff:), afac*b_size, it-1,      &
                                 itag, ireqsh(sreq))
!$OMP END MASTER
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! WAIT for RECV requests and DECODE:
    if (FLAG == 2 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! reconstruct nreqr:
        if (dd(ia)%recv_buf_srf(id) > 0) nreqr = nreqr + 1
      enddo
      if (nreqr > 0) then 
!$OMP MASTER
        rreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqr, ireqrh(rreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (kmx > 1) then
            b_size = dd(ia)%recv_bufsize(id)
          else
            b_size = dd(ia)%recv_buf_srf(id)
          endif
          if (b_size > 0) then
            roff = (irbuf-1)*b_size + 1
            ! decode buffer
            if (afac == 1) then
              call dmpi_decode_buf(iam,id,kmx,mmk,kh,ia,tfrbuf(id)%p(roff:),   &
                                   afac,a)
            elseif (afac == 2) then
              call dmpi_decode_buf(iam,id,kmx,mmk,kh,ia,tfrbuf(id)%p(roff:),   &
                                   afac,a,b)
            elseif (afac == 3) then
              call dmpi_decode_buf(iam,id,kmx,mmk,kh,ia,tfrbuf(id)%p(roff:),   &
                                   afac,a,b,c)
            elseif (afac == 4) then
              call dmpi_decode_buf(iam,id,kmx,mmk,kh,ia,tfrbuf(id)%p(roff:),   &
                                   afac,a,b,c,d)
            elseif (afac == 5) then
              call dmpi_decode_buf(iam,id,kmx,mmk,kh,ia,tfrbuf(id)%p(roff:),   &
                                   afac,a,b,c,d,e)
            elseif (afac == 6) then
              call dmpi_decode_buf(iam,id,kmx,mmk,kh,ia,tfrbuf(id)%p(roff:),   &
                                   afac,a,b,c,d,e,f)
            endif
!$OMP BARRIER
          endif
        enddo
      endif
    endif

    ! WAIT for SEND requests:
    if (FLAG == 5 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it .and. dd(ia)%send_buf_srf(idi) > 0) nreqs = nreqs + 1
          endif
        enddo
      enddo   ! it
      if (nreqs > 0) then 
!$OMP MASTER
        sreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqs, ireqsh(sreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
      endif
    endif

  end subroutine dmpi_distribute_halo_nb

!===============================================================================

  subroutine dmpi_distribute_halo_nb_col(kmx, msrf, mcol, kh, ia,              &
                                         irbuf, isbuf, itagoff, FLAG,          &
                                         a, b, c, d, e, f)
    implicit none

    integer(4), intent(in)              :: ia, kmx
    integer(4), intent(in)              :: msrf(0:,0:), mcol(0:), kh(0:)
    integer(4), intent(in)              :: irbuf, isbuf, itagoff, FLAG
    real(8),    intent(inout)           :: a(0:)
    real(8),    intent(inout), optional :: b(0:), c(0:), d(0:), e(0:), f(0:)


    integer(4) :: itag, it, nb, b_size, id, iam, nbi, idi, ij, afac
    integer(4) :: roff, soff, nreqr, rreq, nreqs, sreq

    iam  = mpi_rank + 1
    itag = NextFreeItag + (itagoff-1)

    afac = 1
    if (present(b)) then
      afac = 2
      if (present(c)) then
        afac = 3
        if (present(d)) then
          afac = 4
          if (present(e)) then
            afac = 5
            if (present(f)) afac = 6
          endif
        endif
      endif
    endif

    ! RECV:
    if (FLAG == 1 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! nb is the sender, iam is the receiver
        nb = mpi_tt(ia,iam)%hnb(id)
        if (kmx > 1) then
          b_size = dd(ia)%recv_bufsize(id)
        else
          b_size = dd(ia)%recv_buf_srf(id)
        endif
        if (b_size > 0) then
          ! recv buffer
!$OMP MASTER
          roff  = (irbuf-1)*b_size + 1
          nreqr = nreqr + 1
          rreq  = (itagoff-1)*max_nb + nreqr
          call dmpi_irecv( tfrbuf(id)%p(roff:), afac*b_size, nb-1, itag,       &
                           ireqrh(rreq))
!$OMP END MASTER
!$OMP BARRIER
        endif
      enddo
    endif

    ! ENCODE:
    if (FLAG == 3 .or. FLAG == 6) then
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                soff = (isbuf-1)*b_size + 1
                ! fill buffer
                if (afac == 1) then
                  call dmpi_encode_buf(iam,kmx,msrf,mcol,kh,ia,                &
                                       tfsbuf(idi)%p(soff:),                   &
                                       it,ij,idi,afac,a)
                elseif (afac == 2) then
                  call dmpi_encode_buf(iam,kmx,msrf,mcol,kh,ia,                &
                                       tfsbuf(idi)%p(soff:),                   &
                                       it,ij,idi,afac,a,b)
                elseif (afac == 3) then
                  call dmpi_encode_buf(iam,kmx,msrf,mcol,kh,ia,                &
                                       tfsbuf(idi)%p(soff:),                   &
                                       it,ij,idi,afac,a,b,c)
                elseif (afac == 4) then
                  call dmpi_encode_buf(iam,kmx,msrf,mcol,kh,ia,                &
                                       tfsbuf(idi)%p(soff:),                   &
                                       it,ij,idi,afac,a,b,c,d)
                elseif (afac == 5) then
                  call dmpi_encode_buf(iam,kmx,msrf,mcol,kh,ia,                &
                                       tfsbuf(idi)%p(soff:),                   &
                                       it,ij,idi,afac,a,b,c,d,e)
                elseif (afac == 6) then
                  call dmpi_encode_buf(iam,kmx,msrf,mcol,kh,ia,                &
                                       tfsbuf(idi)%p(soff:),                   &
                                       it,ij,idi,afac,a,b,c,d,e,f)
                endif
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! SEND:
    if (FLAG == 4 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                ! send buffer
!$OMP MASTER
                soff  = (isbuf-1)*b_size + 1
                nreqs = nreqs + 1
                sreq  = (itagoff-1)*max_nb + nreqs
                call dmpi_isend( tfsbuf(idi)%p(soff:), afac*b_size, it-1,      &
                                 itag, ireqsh(sreq))
!$OMP END MASTER
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! WAIT for RECV requests and DECODE:
    if (FLAG == 2 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! reconstruct nreqr:
        if (dd(ia)%recv_buf_srf(id) > 0) nreqr = nreqr + 1
      enddo
      if (nreqr > 0) then 
!$OMP MASTER
        rreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqr, ireqrh(rreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (kmx > 1) then
            b_size = dd(ia)%recv_bufsize(id)
          else
            b_size = dd(ia)%recv_buf_srf(id)
          endif
          if (b_size > 0) then
            roff = (irbuf-1)*b_size + 1
            ! decode buffer
            if (afac == 1) then
              call dmpi_decode_buf(iam,id,kmx,msrf,mcol,kh,ia,                 &
                                   tfrbuf(id)%p(roff:),                        &
                                   afac,a)
            elseif (afac == 2) then
              call dmpi_decode_buf(iam,id,kmx,msrf,mcol,kh,ia,                 &
                                   tfrbuf(id)%p(roff:),                        &
                                   afac,a,b)
            elseif (afac == 3) then
              call dmpi_decode_buf(iam,id,kmx,msrf,mcol,kh,ia,                 &
                                   tfrbuf(id)%p(roff:),                        &
                                   afac,a,b,c)
            elseif (afac == 4) then
              call dmpi_decode_buf(iam,id,kmx,msrf,mcol,kh,ia,                 &
                                   tfrbuf(id)%p(roff:),                        &
                                   afac,a,b,c,d)
            elseif (afac == 5) then
              call dmpi_decode_buf(iam,id,kmx,msrf,mcol,kh,ia,                 &
                                   tfrbuf(id)%p(roff:),                        &
                                   afac,a,b,c,d,e)
            elseif (afac == 6) then
              call dmpi_decode_buf(iam,id,kmx,msrf,mcol,kh,ia,                 &
                                   tfrbuf(id)%p(roff:),                        &
                                   afac,a,b,c,d,e,f)
            endif
!$OMP BARRIER
          endif
        enddo
      endif
    endif

    ! WAIT for SEND requests:
    if (FLAG == 5 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it .and. dd(ia)%send_buf_srf(idi) > 0) nreqs = nreqs + 1
          endif
        enddo
      enddo   ! it
      if (nreqs > 0) then 
!$OMP MASTER
        sreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqs, ireqsh(sreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
      endif
    endif

  end subroutine dmpi_distribute_halo_nb_col

!===============================================================================

  subroutine dmpi_distribute_halo_TF_nb(kmx, mmk, kh, ia,                      &
                                        irbuf, isbuf, itagoff, FLAG,           &
                                        a, b, c)
    implicit none

    integer(4), intent(in)              :: ia, kmx, mmk(:,0:,0:), kh(0:)
    integer(4), intent(in)              :: irbuf, isbuf, itagoff, FLAG
    real(8),    intent(inout)           :: a(:,0:)
    real(8),    intent(inout), optional :: b(:,0:), c(:,0:)


    integer(4) :: itag, it, nb, b_size, id, iam, nbi, idi, ij, nc, afac
    integer(4) :: roff, soff, nreqr, rreq, nreqs, sreq

    iam  = mpi_rank + 1
    nc   = mpi_nc
    itag = NextFreeItag + (itagoff-1)

    afac = 1
    if (present(b)) then
      afac = 2
      if (present(c)) afac = 3
    endif

    ! RECV:
    if (FLAG == 1 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! nb is the sender, iam is the receiver
        nb = mpi_tt(ia,iam)%hnb(id)
        if (kmx > 1) then
          b_size = dd(ia)%recv_bufsize(id)
        else
          b_size = dd(ia)%recv_buf_srf(id)
        endif
        if (b_size > 0) then
          ! recv buffer
!$OMP MASTER
          roff  = (irbuf-1)*nc*b_size + 1
          nreqr = nreqr + 1
          rreq  = (itagoff-1)*max_nb + nreqr
          call dmpi_irecv( tfrbuf(id)%p(roff:), afac*nc*b_size, nb-1, itag,    &
                           ireqrh(rreq))
!$OMP END MASTER
!$OMP BARRIER
        endif
      enddo
    endif

    ! ENCODE:
    if (FLAG == 3 .or. FLAG == 6) then
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                soff = (isbuf-1)*nc*b_size + 1
                ! fill buffer
                if (afac == 1) then
                  call dmpi_encode_buf_nc(iam,kmx,mmk,kh,ia,a,                 &
                                          tfsbuf(idi)%p(soff:),it,ij,idi,nc)
                elseif (afac == 2) then
                  call dmpi_encode_buf_nc(iam,kmx,mmk,kh,ia,a,b,               &
                                          tfsbuf(idi)%p(soff:),it,ij,idi,nc)
                elseif (afac == 3) then
                  call dmpi_encode_buf_nc(iam,kmx,mmk,kh,ia,a,b,c,             &
                                          tfsbuf(idi)%p(soff:),it,ij,idi,nc)
                endif
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! SEND:
    if (FLAG == 4 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                ! send buffer
!$OMP MASTER
                soff  = (isbuf-1)*nc*b_size + 1
                nreqs = nreqs + 1
                sreq  = (itagoff-1)*max_nb + nreqs
                call dmpi_isend( tfsbuf(idi)%p(soff:), afac*nc*b_size, it-1,   &
                                 itag, ireqsh(sreq))
!$OMP END MASTER
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! WAIT for RECV requests and DECODE:
    if (FLAG == 2 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! reconstruct nreqr:
        if (dd(ia)%recv_buf_srf(id) > 0) nreqr = nreqr + 1
      enddo
      if (nreqr > 0) then 
!$OMP MASTER
        rreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqr, ireqrh(rreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (kmx > 1) then
            b_size = dd(ia)%recv_bufsize(id)
          else
            b_size = dd(ia)%recv_buf_srf(id)
          endif
          if (b_size > 0) then
            roff = (irbuf-1)*nc*b_size + 1
            ! decode buffer
            if (afac == 1) then
              call dmpi_decode_buf_nc(iam,id,kmx,mmk,kh,ia,a,                  &
                                      tfrbuf(id)%p(roff:),nc)
            elseif (afac == 2) then
              call dmpi_decode_buf_nc(iam,id,kmx,mmk,kh,ia,a,b,                &
                                      tfrbuf(id)%p(roff:),nc)
            elseif (afac == 3) then
              call dmpi_decode_buf_nc(iam,id,kmx,mmk,kh,ia,a,b,c,              &
                                      tfrbuf(id)%p(roff:),nc)
            endif
!$OMP BARRIER
          endif
        enddo
      endif
    endif

    ! WAIT for SEND requests:
    if (FLAG == 5 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it .and. dd(ia)%send_buf_srf(idi) > 0) nreqs = nreqs + 1
          endif
        enddo
      enddo   ! it
      if (nreqs > 0) then 
!$OMP MASTER
        sreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqs, ireqsh(sreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
      endif
    endif

  end subroutine dmpi_distribute_halo_TF_nb

!===============================================================================

  subroutine dmpi_distribute_halo_TF_nb2(kmx, mmk, kh, ia,                     &
                                         irbuf, isbuf, itagoff, FLAG,          &
                                         a, nci)
    implicit none

    integer(4), intent(in)    :: ia, kmx, mmk(:,0:,0:), kh(0:)
    integer(4), intent(in)    :: irbuf, isbuf, itagoff, FLAG, nci
    real(8),    intent(inout) :: a(:,0:)


    integer(4) :: itag, it, nb, b_size, id, iam, nbi, idi, ij, nc
    integer(4) :: roff, soff, nreqr, rreq, nreqs, sreq

    iam  = mpi_rank + 1
    nc   = nci
    itag = NextFreeItag + (itagoff-1)

    ! RECV:
    if (FLAG == 1 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! nb is the sender, iam is the receiver
        nb = mpi_tt(ia,iam)%hnb(id)
        if (kmx > 1) then
          b_size = dd(ia)%recv_bufsize(id)
        else
          b_size = dd(ia)%recv_buf_srf(id)
        endif
        if (b_size > 0) then
          ! recv buffer
!$OMP MASTER
          roff  = (irbuf-1)*nc*b_size + 1
          nreqr = nreqr + 1
          rreq  = (itagoff-1)*max_nb + nreqr
          call dmpi_irecv( tfrbuf(id)%p(roff:), nc*b_size, nb-1, itag,         &
                           ireqrh(rreq))
!$OMP END MASTER
!$OMP BARRIER
        endif
      enddo
    endif

    ! ENCODE:
    if (FLAG == 3 .or. FLAG == 6) then
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                soff = (isbuf-1)*nc*b_size + 1
                ! fill buffer
                call dmpi_encode_buf_nc(iam,kmx,mmk,kh,ia,a,                   &
                                        tfsbuf(idi)%p(soff:),it,ij,idi,nc)
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! SEND:
    if (FLAG == 4 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #it and I must send to him

            ! Find my neighbour which corresponds to task #it and
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it) then
              if (kmx > 1) then
                b_size = dd(ia)%send_bufsize(idi)
              else
                b_size = dd(ia)%send_buf_srf(idi)
              endif
              if (b_size > 0) then
                call dmpi_roff (idi, ia, iam, it, ij)
                ! send buffer
!$OMP MASTER
                soff  = (isbuf-1)*nc*b_size + 1
                nreqs = nreqs + 1
                sreq  = (itagoff-1)*max_nb + nreqs
                call dmpi_isend( tfsbuf(idi)%p(soff:), nc*b_size, it-1,        &
                                 itag, ireqsh(sreq))
!$OMP END MASTER
!$OMP BARRIER
              endif
            endif  ! nbi == it
          endif
        enddo
      enddo   ! it
    endif

    ! WAIT for RECV requests and DECODE:
    if (FLAG == 2 .or. FLAG == 6) then
      ! I should not send to myself, but receive from up to 8 of my neighbours:
      nreqr = 0
      do id=1,max_nb
        ! reconstruct nreqr:
        if (dd(ia)%recv_buf_srf(id) > 0) nreqr = nreqr + 1
      enddo
      if (nreqr > 0) then 
!$OMP MASTER
        rreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqr, ireqrh(rreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (kmx > 1) then
            b_size = dd(ia)%recv_bufsize(id)
          else
            b_size = dd(ia)%recv_buf_srf(id)
          endif
          if (b_size > 0) then
            roff = (irbuf-1)*nc*b_size + 1
            ! decode buffer
            call dmpi_decode_buf_nc(iam,id,kmx,mmk,kh,ia,a,                    &
                                    tfrbuf(id)%p(roff:),nc)
!$OMP BARRIER
          endif
        enddo
      endif
    endif

    ! WAIT for SEND requests:
    if (FLAG == 5 .or. FLAG == 6) then
      nreqs = 0
      do it=1,mpi_size
        if (it == iam) cycle
        do id=1,max_nb
          nb = mpi_tt(ia,it)%hnb(id)
          if (nb == iam) then
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == it .and. dd(ia)%send_buf_srf(idi) > 0) nreqs = nreqs + 1
          endif
        enddo
      enddo   ! it
      if (nreqs > 0) then 
!$OMP MASTER
        sreq = (itagoff-1)*max_nb + 1
        call MPI_Waitall(nreqs, ireqsh(sreq:), lstatush, ierr)
!$OMP END MASTER
!$OMP BARRIER
      endif
    endif

  end subroutine dmpi_distribute_halo_TF_nb2

!===============================================================================

  subroutine dmpi_distribute_halo_log(mmk, ia, a)
    implicit none

    integer(4), intent(in)    :: ia, mmk(:,0:,0:)
    logical,    intent(inout) :: a(0:)

    integer(4) :: itag, it_recv, itp1, nb, b_size, id, iam, nbi, idi
    integer(4) :: ij

    iam = mpi_rank + 1

    ! FIXME: Not necessary to run through ALL mpi-tasks, is it? Guess, we 
    !        should only run through those tasks which could possibly have
    !        anything to communicate with actual task #iam. Make a look-up
    !        table or something ...
    !        Also, it should be possible to do the receive fork (itp1 == iam) 
    !        outside the it_recv loop. 

    do it_recv=0,mpi_size-1
      itp1 = it_recv + 1

      ! I should not send to myself, but receive from up to 8 of my neighbours:
      if (itp1 == iam) then
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb     = mpi_tt(ia,iam)%hnb(id)
          b_size = dd(ia)%recv_buf_srf(id)
          if (b_size > 0) then 
            itag = 1
!$OMP MASTER
            ! recv buffer
            call dmpi_recv( lbuf, b_size, nb-1, itag )
!$OMP END MASTER
!$OMP BARRIER
            ! decode buffer
            call dmpi_decode_buf_log(iam, id, mmk, ia, a, lbuf)
!$OMP BARRIER
          endif
        enddo

      ! I should not receive from myself, but send to up to 8 of my neighbours:
      else
        do id=1,max_nb
          nb = mpi_tt(ia,itp1)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #itp1 and I must send to him

            ! Find my neighbour which corresponds to task #itp1 and 
            ! send to it
            idi    = n2n(id)
            b_size = 0
            nbi    = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == itp1) then
              b_size = dd(ia)%send_buf_srf(idi)
              if (b_size > 0) then 
                call dmpi_roff (idi, ia, iam, itp1, ij)
                itag = 1
                ! fill buffer
                call dmpi_encode_buf_log(iam, mmk, ia, a, lbuf, itp1, ij, idi)
!$OMP BARRIER
!$OMP MASTER
                ! send buffer
                call dmpi_send( lbuf, b_size, it_recv, itag )
!$OMP END MASTER
!$OMP BARRIER
              endif
            endif  ! nbi == itp1
          endif
        enddo

      endif
    enddo   ! it_recv, itp1

  end subroutine dmpi_distribute_halo_log

!===============================================================================

  subroutine dmpi_gather_global_srf (ia,ind_l,mmk_l,mmk,a,trecv)
    implicit none

    integer(4), intent(in)    :: ia, trecv
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), ind_l(:,:) 
    real(8),    intent(inout) :: a(0:)

    integer(4) :: iw, itag, ii, it, bufsize, itp1, i, j

    itag = 1

    if ((trecv >= mpi_size) .or. (trecv <0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    if (mpi_rank /= trecv) then

      if (dd(ia)%up_ws /= 0) then
        bufsize = dd(ia)%up_ws-dd(ia)%low_ws+1
        ii=0
        do iw=1,bufsize
          i = ind_l(1,iw)
          j = ind_l(2,iw)
          if (mmk_l(1,i,j) == 0) cycle
          ii       = ii + 1
          gbuf(ii) = a(mmk(1,i,j)) 
        enddo
        call dmpi_send( gbuf, bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        if (it /= trecv) then 
          itp1 = it + 1
          if (mpi_tt(ia,itp1)%up_ws /= 0) then
            bufsize = mpi_tt(ia,itp1)%up_ws-mpi_tt(ia,itp1)%low_ws+1
            call dmpi_recv( gbuf, bufsize, it, itag )
            ii=0
            do j=mpi_tt(ia,itp1)%low_j,mpi_tt(ia,itp1)%up_j
              do i=mpi_tt(ia,itp1)%low_i,mpi_tt(ia,itp1)%up_i
                if (mmk(1,i,j) == 0) cycle
                ii            = ii + 1
                a(mmk(1,i,j)) = gbuf(ii)
              enddo
            enddo
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_global_srf

!===============================================================================

  subroutine dmpi_gather_global_all(ia,kmx,ind_l,mmk_l,mmk,a,trecv)
    implicit none
    integer(4), intent(in)    :: ia, kmx, trecv
    integer(4), intent(in)    :: mmk(:,0:,0:), mmk_l(:,0:,0:), ind_l(:,:)
    real(8),    intent(inout) :: a(0:)

    integer(4) :: iw, itag, ii, it, itp1, bufsize, i, j, k

    if ((trecv >= mpi_size) .or. (trecv < 0)) then
      call exitme(1,'Called dmpi_gather with invalid trecv parameter')
    endif

    itag = 1

    if (mpi_rank /= trecv) then

      if (dd(ia)%nwet3 > 0) then
        bufsize = dd(ia)%nwet3
        ii = 0
        do iw=1,dd(ia)%up_ws - dd(ia)%low_ws + 1
          i = ind_l(1,iw)
          j = ind_l(2,iw)
          if (mmk_l(1,i,j) == 0) cycle
          do k=1,kmx
            if (mmk_l(k,i,j) == 0) exit
            ii       = ii + 1
            gbuf(ii) = a(mmk(k,i,j)) 
          enddo
        enddo
        call dmpi_send( gbuf, bufsize, trecv, itag )
      endif

    else

      do it=0,mpi_size-1
        if (it /= trecv) then 
          itp1 = it + 1
          if (mpi_tt(ia,itp1)%nwet3 > 0) then
            bufsize = mpi_tt(ia,itp1)%nwet3
            call dmpi_recv( gbuf, bufsize, it, itag )
            ii = 0
            do j=mpi_tt(ia,itp1)%low_j,mpi_tt(ia,itp1)%up_j
              do i=mpi_tt(ia,itp1)%low_i,mpi_tt(ia,itp1)%up_i
                if (mmk(1,i,j) == 0) cycle
                do k=1,kmx
                  if (mmk(k,i,j) == 0) exit
                  ii            = ii + 1
                  a(mmk(k,i,j)) = gbuf(ii)
                enddo
              enddo
            enddo
          endif
        endif
      enddo

    endif

  end subroutine dmpi_gather_global_all

!===============================================================================

  subroutine dmpi_gather_bnd_data_out(ia, nc, kmx, nz1, nz2, krz, a)
    implicit none
    integer(4), intent(in)           :: ia, nc, kmx, nz1, nz2, krz(:,0:)
    real(8),    intent(inout)        :: a(:,:,0:)           ! (1:nc,1:kmx,0:nz)

    integer(4) :: it, iam, itag, b_size, iz, ik, ic, ibnd, jbnd, nzz, i, j

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  nothing to do here? ----------------------------------------------------
    if (nz2 > 0) then
      nzz = nz2 - nz1 + 1
    else
      nzz = 0
    endif
    if (nzz < 1) return

    !-  misc vars --------------------------------------------------------------
    itag = 1

    !-  recv on I/O task from all tasks having any data ------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam) then
          !  find buffer size:
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            b_size = b_size + kmx*nc
          enddo

          !  receive and decode:
          if (b_size > 0) call dmpi_recv(gbuf, b_size, it-1, itag)
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            do ik=1,kmx
              do ic=1,nc
                b_size      = b_size + 1
                a(ic,ik,iz) = gbuf(b_size)
              enddo
            enddo
          enddo
        endif
      enddo

    !-  if I have data, send it to I/O task ------------------------------------
    else
      !  encode and send:
      b_size = 0
      do iz=nz1,nz2
        ! skip if this point is not on this task:
        i = krz(1,iz)
        j = krz(2,iz)
        if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                         &
            j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
        do ik=1,kmx
          do ic=1,nc
            b_size       = b_size + 1
            gbuf(b_size) = a(ic,ik,iz)
          enddo
        enddo
      enddo
      if (b_size > 0) call dmpi_send(gbuf, b_size, mpi_io_rank, itag)

    endif

  end subroutine dmpi_gather_bnd_data_out

!===============================================================================

  subroutine dmpi_gather_bnd_data(ia, nc, kmx, nz1, nz2, krz, a, iga, jga,     &
                                  ii, FLAG)

    !  FALG = 1 : pre-post RECV
    !       = 2 : encode and pre-post SEND
    !       = 3 : WAIT for RECVs and decode
    !       = 4 : WAIT for SENDs
    !       = 5 : all the above in order 1, 3, 2, 4

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: ia, nc, kmx, nz1, nz2, krz(:,0:)
    integer(4), intent(in)    :: iga(:), jga(:), ii, FLAG
    real(8),    intent(inout) :: a(:,:,0:)           ! (1:nc,1:kmx,0:nz)

    integer(4) :: it, iam, itag, iz, ik, ibnd, jbnd, nzz, i, j, ierr
    integer(4) :: nzl, nzu, btmp, tnum, itr, ireqi, ireql, irequ, nnreq, imask
    integer(4) :: itl, itu, irrl, irru
    ! valid on MASTER only:
    integer(4) :: b_size

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1
    call domp_get_thread_no( tnum )

    !-  nothing to do here? ----------------------------------------------------
    if (nz2 > 0) then
      nzz = nz2 - nz1 + 1
    else
      nzz = 0
    endif
    if (nzz < 1) return

    !-  misc vars --------------------------------------------------------------
    if (do_comm) itag = nestitags7(ia)%p(ii)

    !-  Prepost recv on I/O task from all tasks having any data ----------------
    if (mpi_io_serial .and. (FLAG == 1 .or. FLAG == 5) .and. do_comm) then
      irequ = itagreq7(itag,1,2)
      if (irequ > 0) then
        ireqi = itagreq7(itag,1,1)
        call domp_get_domain(1, mpi_size, itl, itu, .true.)
        call domp_get_domain(nz1, nz2, nzl, nzu)
        do it=1,mpi_size
          if (it == iam) cycle
          bsiz2(tnum,it,1,1) = 0
          ! receive but not from myself:
          imask = rsmask7(ia,it)%p(ii)
          if (imask == 0 .or. imask == 10) cycle
          !  find buffer size:
          btmp = 0
          do iz=nzl,nzu
            ! skip if this point is not on the relevant task:
            ibnd = krz(1,iz)
            jbnd = krz(2,iz)
            i = iga(1)+int((ibnd-1)/iga(2),4)
            j = jga(1)+int((jbnd-1)/jga(2),4)
            if     (krz(3,iz) == 1) then
              j = j-1
            elseif (krz(3,iz) == 2) then
              i = i-1
            elseif (krz(3,iz) == 3) then
              j = j+1
            elseif (krz(3,iz) == 4) then
              i = i+1
            endif
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            btmp = btmp + kmx*nc
          enddo
          bsiz2(tnum,it,1,1) = btmp
!$OMP BARRIER
!$OMP MASTER
          b_size = bsiz2(1,it,1,1)
          bsiz2(1,it,2,1) = 0
          do itr=2,mpi_nt
            b_size = b_size + bsiz2(itr,it,1,1)
            bsiz2(itr,it,2,1) = bsiz2(itr-1,it,1,1) + bsiz2(itr-1,it,2,1)
          enddo
          thr_b_size = b_size
!$OMP END MASTER
!$OMP BARRIER
          if (itl <= it .and. it <= itu) then
            call dmpi_irecv(rbufnb7(ireqi)%p,thr_b_size,it-1,itag,ireq7(ireqi))
          endif
          ireqi = ireqi + 1
        enddo
      endif
    endif

    !-  Wait for recv on I/O task and decode -----------------------------------
    if (mpi_io_serial .and. (FLAG == 3 .or. FLAG == 5) .and. do_comm) then
      irequ = itagreq7(itag,1,2)
      if (irequ > 0) then
        ireql = itagreq7(itag,1,1)
        call domp_get_domain(nz1, nz2, nzl, nzu)
        call domp_get_domain(ireql, irequ, irrl, irru)

        nnreq = irru - irrl + 1
        if (nnreq > 0) then
          call MPI_Waitall(nnreq,ireq7(irrl:irru),lstat7(:,irrl:irru),ierr)
        endif

!$OMP BARRIER
        ireqi = ireql
        do it=1,mpi_size
          ! recieve but not from myself:
          imask = rsmask7(ia,it)%p(ii)
          if (imask == 0 .or. imask == 10) cycle
          !  decode:
          if (bsiz2(tnum,it,1,1) > 0) then
            btmp = bsiz2(tnum,it,2,1)
            do iz=nzl,nzu
              ! skip if this point is not on the relevant task:
              ibnd = krz(1,iz)
              jbnd = krz(2,iz)
              i = iga(1)+int((ibnd-1)/iga(2),4)
              j = jga(1)+int((jbnd-1)/jga(2),4)
              if     (krz(3,iz) == 1) then
                j = j-1
              elseif (krz(3,iz) == 2) then
                i = i-1
              elseif (krz(3,iz) == 3) then
                j = j+1
              elseif (krz(3,iz) == 4) then
                i = i+1
              endif
              if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.     &
                  j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
              do ik=1,kmx
                a(1:nc,ik,iz) = rbufnb7(ireqi)%p(btmp+1:btmp+nc)
                btmp = btmp + nc
              enddo
            enddo
          endif
          ireqi = ireqi + 1
        enddo
      endif
    endif

    !-  if I have data, send it to I/O task ------------------------------------
    if ((.not.mpi_io_serial) .and. (FLAG == 2 .or. FLAG == 5) .and. do_comm)then
      irequ = itagreq7(itag,2,2)
      if (irequ > 0) then
        call domp_get_domain(nz1, nz2, nzl, nzu)
        call domp_get_domain(1, mpi_size, itl, itu, .true.)
        it    = mpi_io_rank + 1
        ireqi = itagreq7(itag,2,1)
        imask = rsmask7(ia,iam)%p(ii)
        if (imask /= 0 .and. imask /= 1) then

          !  first, find thread offset:
          offset(tnum,1:2) = 0
          btmp = 0
          do iz=nzl,nzu
            ! skip if this point is not on this task:
            ibnd = krz(1,iz)
            jbnd = krz(2,iz)
            i = iga(1)+int((ibnd-1)/iga(2),4)
            j = jga(1)+int((jbnd-1)/jga(2),4)
            if     (krz(3,iz) == 1) then
              j = j-1
            elseif (krz(3,iz) == 2) then
              i = i-1
            elseif (krz(3,iz) == 3) then
              j = j+1
            elseif (krz(3,iz) == 4) then
              i = i+1
            endif
            if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                     &
                j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
            btmp = btmp + kmx*nc
          enddo
          offset(tnum,1) = btmp
!$OMP BARRIER
!$OMP MASTER
          b_size = offset(1,1)
          offset(1,2) = 0
          do itr=2,mpi_nt
            b_size = b_size + offset(itr,1)
            offset(itr,2) = offset(itr-1,1) + offset(itr-1,2)
          enddo
          thr_b_size = b_size
!$OMP END MASTER
!$OMP BARRIER

          !  then, encode and send:
          if (offset(tnum,1) > 0) then
            btmp = offset(tnum,2)
            do iz=nzl,nzu
              ! skip if this point is not on this task:
              ibnd = krz(1,iz)
              jbnd = krz(2,iz)
              i = iga(1)+int((ibnd-1)/iga(2),4)
              j = jga(1)+int((jbnd-1)/jga(2),4)
              if     (krz(3,iz) == 1) then
                j = j-1
              elseif (krz(3,iz) == 2) then
                i = i-1
              elseif (krz(3,iz) == 3) then
                j = j+1
              elseif (krz(3,iz) == 4) then
                i = i+1
              endif
              if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                   &
                  j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
              do ik=1,kmx
                sbufnb7(ireqi)%p(btmp+1:btmp+nc) = a(1:nc,ik,iz)
                btmp = btmp + nc
              enddo
            enddo
          endif
!$OMP BARRIER
          if (thr_b_size > 0 .and. itl <= iam .and. iam <= itu) then
            call dmpi_isend(sbufnb7(ireqi)%p,thr_b_size,it-1,itag,irqs7(ireqi))
          endif
          ireqi = ireqi + 1
        endif
      endif
    endif

    !-  complete the send from this task ---------------------------------------
    if ((.not.mpi_io_serial) .and. (FLAG == 4 .or. FLAG == 5) .and. do_comm)then
      !  wait for the send to finish
      irequ = itagreq7(itag,2,2)
      if (irequ > 0) then
!$OMP MASTER
        ireql = itagreq7(itag,2,1)
        nnreq = irequ - ireql + 1
        call MPI_Waitall(nnreq, irqs7(ireql:irequ), lstat7(:,ireql:irequ), ierr)
!$OMP END MASTER
      endif
    endif

!$OMP BARRIER

  end subroutine dmpi_gather_bnd_data

!===============================================================================

  subroutine dmpi_gather_bnd_data_2(ia, na, nz1, nz2, krz, a, iga, jga)
    implicit none
    integer(4), intent(in)           :: ia, na, nz1, nz2, krz(:,0:)
    integer(4), intent(in), optional :: iga(:), jga(:)
    real(8),    intent(inout)        :: a(:,0:)   ! (1:na,0:nz)

    integer(4) :: it, iam, itag, b_size, iz, ic, ibnd, jbnd, nzz, i, j
    logical    :: bc

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  nothing to do here? ----------------------------------------------------
    if (nz2 > 0) then
      nzz = nz2 - nz1 + 1
    else
      nzz = 0
    endif
    if (nzz < 1) return

    !-  misc vars --------------------------------------------------------------
    itag = 1
    bc   = (present(iga))

    !-  recv on I/O task from all tasks having any data ------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam) then
          !  find buffer size:
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            ibnd = krz(1,iz)
            jbnd = krz(2,iz)
            if (bc) then
              i = iga(1)+int((ibnd-1)/iga(2),4)
              j = jga(1)+int((jbnd-1)/jga(2),4)
              if     (krz(3,iz) == 1) then
                j = j-1
              elseif (krz(3,iz) == 2) then
                i = i-1
              elseif (krz(3,iz) == 3) then
                j = j+1
              elseif (krz(3,iz) == 4) then
                i = i+1
              endif
            else
              i = ibnd
              j = jbnd
            endif
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            b_size = b_size + na
          enddo

          !  receive and decode:
          if (b_size > 0) call dmpi_recv(gbuf, b_size, it-1, itag)
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            ibnd = krz(1,iz)
            jbnd = krz(2,iz)
            if (bc) then
              i = iga(1)+int((ibnd-1)/iga(2),4)
              j = jga(1)+int((jbnd-1)/jga(2),4)
              if     (krz(3,iz) == 1) then
                j = j-1
              elseif (krz(3,iz) == 2) then
                i = i-1
              elseif (krz(3,iz) == 3) then
                j = j+1
              elseif (krz(3,iz) == 4) then
                i = i+1
              endif
            else
              i = ibnd
              j = jbnd
            endif
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            do ic=1,na
              b_size   = b_size + 1
              a(ic,iz) = gbuf(b_size)
            enddo
          enddo
        endif
      enddo

    !-  if I have data, send it to I/O task ------------------------------------
    else
      !  encode and send:
      b_size = 0
      do iz=nz1,nz2
        ! skip if this point is not on this task:
        ibnd = krz(1,iz)
        jbnd = krz(2,iz)
        if (bc) then
          i = iga(1)+int((ibnd-1)/iga(2),4)
          j = jga(1)+int((jbnd-1)/jga(2),4)
          if     (krz(3,iz) == 1) then
            j = j-1
          elseif (krz(3,iz) == 2) then
            i = i-1
          elseif (krz(3,iz) == 3) then
            j = j+1
          elseif (krz(3,iz) == 4) then
            i = i+1
          endif
        else
          i = ibnd
          j = jbnd
        endif
        if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                         &
            j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
        do ic=1,na
          b_size       = b_size + 1
          gbuf(b_size) = a(ic,iz)
        enddo
      enddo
      if (b_size > 0) call dmpi_send(gbuf, b_size, mpi_io_rank, itag)

    endif

  end subroutine dmpi_gather_bnd_data_2

!===============================================================================

  subroutine dmpi_scatter_bnd_data(ia, iia, ii, nc, kmx, nz1, nz2, krz, a)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    integer(4), intent(in)    :: ia, iia, ii, nc, kmx, nz1, nz2, krz(:,0:)
    real(8),    intent(inout) :: a(:,:,0:)           ! (1:nc,1:kmx,0:nz)

    integer(4) :: it, iam, itag, iz, ik, nzz, i, j, tnum, btmp, itr, nzl, nzu
    integer(4) :: b_size
!Fixme: we must find the max size needed and allocate bbbb only once
!       the present code is only proof-of-concept
    real(8), allocatable :: bbbb(:)

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1
    call domp_get_thread_no( tnum )

    !-  nothing to do here? ----------------------------------------------------
    if (nz2 > 0) then
      nzz = nz2 - nz1 + 1
    else
      nzz = 0
    endif
    if (nzz < 1) return

    !-  misc vars --------------------------------------------------------------
    itag = NextFreeItag + tnum - 1 

    !-  send from I/O task to all tasks wanting data ---------------------------
    if (mpi_io_serial) then
      call domp_get_domain(nz1, nz2, nzl, nzu)

      do it=1,mpi_size
        if (it == iam) cycle
        !  encode and send:
        b_size = 0
        do iz=nzl,nzu
          ! skip if this point is not on the relevant task:
          i = krz(1,iz)
          j = krz(2,iz)
          if (i < mpi_tt(iia,it)%low_i .or. i > mpi_tt(iia,it)%up_i .or.       &
              j < mpi_tt(iia,it)%low_j .or. j > mpi_tt(iia,it)%up_j) cycle
          b_size = b_size + kmx*nc
        enddo

        !  then, encode and send:
        if (b_size > 0) then
          allocate( bbbb(b_size) )
          btmp = 0
          do iz=nzl,nzu
            ! skip if this point is not on the relevant task:
            i = krz(1,iz)
            j = krz(2,iz)
            if (i < mpi_tt(iia,it)%low_i .or. i > mpi_tt(iia,it)%up_i .or.     &
                j < mpi_tt(iia,it)%low_j .or. j > mpi_tt(iia,it)%up_j) cycle
            do ik=1,kmx
              bbbb(btmp+1:btmp+nc) = a(1:nc,ik,iz)
              btmp = btmp + nc
            enddo
          enddo
          call dmpi_send(bbbb, b_size, it-1, itag)
          deallocate( bbbb )
        endif
      enddo

    !-  if I want data, recv it from I/O task ----------------------------------
    else
      call domp_get_domain(nz1, nz2, nzl, nzu)

      it = mpi_io_rank + 1
      !  find buffer size:
      bsiz2(tnum,it,1,2) = 0
      b_size = 0
      do iz=nzl,nzu
        ! skip if this point is not on the relevant task:
        i = krz(1,iz)
        j = krz(2,iz)
        if (i < dd(iia)%low_i .or. i > dd(iia)%up_i .or.                       &
            j < dd(iia)%low_j .or. j > dd(iia)%up_j       ) cycle
        b_size = b_size + kmx*nc
      enddo

      !  receive and decode:
      if (b_size > 0) then
        allocate( bbbb(b_size) )
        call dmpi_recv(bbbb, b_size, it-1, itag)
        btmp = 0
        do iz=nzl,nzu
          ! skip if this point is not on the relevant task:
          i = krz(1,iz)
          j = krz(2,iz)
          if (i < dd(iia)%low_i .or. i > dd(iia)%up_i .or.                     &
              j < dd(iia)%low_j .or. j > dd(iia)%up_j       ) cycle
          do ik=1,kmx
            a(1:nc,ik,iz) = bbbb(btmp+1:btmp+nc)
            btmp = btmp + nc
          enddo
        enddo
        deallocate( bbbb )
      endif

    endif

  end subroutine dmpi_scatter_bnd_data

!===============================================================================

  subroutine dmpi_scatter_bnd_data_2(ia, na, nz1, nz2, krz, a, iga, jga)
    implicit none
    integer(4), intent(in)           :: ia, na,nz1, nz2, krz(:,0:)
    integer(4), intent(in), optional :: iga(:), jga(:)
    real(8),    intent(inout)        :: a(:,0:)    ! (1:na,0:nz)

    integer(4) :: it, iam, itag, b_size, iz, ic, ibnd, jbnd, nzz, i, j
    logical    :: bc

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  nothing to do here? ----------------------------------------------------
    if (nz2 > 0) then
      nzz = nz2 - nz1 + 1
    else
      nzz = 0
    endif
    if (nzz < 1) return

    !-  misc vars --------------------------------------------------------------
    itag = 1
    bc   = present(iga)

    !-  send from I/O task to all tasks wanting data ---------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam) then
          !  encode and send:
          b_size = 0
          do iz=nz1,nz2
            ! skip if this point is not on the relevant task:
            ibnd = krz(1,iz)
            jbnd = krz(2,iz)
            if (bc) then
              i = iga(1)+int((ibnd-1)/iga(2),4)
              j = jga(1)+int((jbnd-1)/jga(2),4)
              if     (krz(3,iz) == 1) then
                j = j-1
              elseif (krz(3,iz) == 2) then
                i = i-1
              elseif (krz(3,iz) == 3) then
                j = j+1
              elseif (krz(3,iz) == 4) then
                i = i+1
              endif
            else
              i = ibnd
              j = jbnd
            endif
            if (i < mpi_tt(ia,it)%low_i .or. i > mpi_tt(ia,it)%up_i .or.       &
                j < mpi_tt(ia,it)%low_j .or. j > mpi_tt(ia,it)%up_j) cycle
            do ic=1,na
              b_size       = b_size + 1
              gbuf(b_size) = a(ic,iz)
            enddo
          enddo
          if (b_size > 0) call dmpi_send(gbuf, b_size, it-1, itag)

        endif
      enddo

    !-  if I want data, recv it from I/O task ----------------------------------
    else
      !  find buffer size:
      b_size = 0
      do iz=nz1,nz2
        ! skip if this point is not on the relevant task:
        ibnd = krz(1,iz)
        jbnd = krz(2,iz)
        if (bc) then
          i = iga(1)+int((ibnd-1)/iga(2),4)
          j = jga(1)+int((jbnd-1)/jga(2),4)
          if     (krz(3,iz) == 1) then
            j = j-1
          elseif (krz(3,iz) == 2) then
            i = i-1
          elseif (krz(3,iz) == 3) then
            j = j+1
          elseif (krz(3,iz) == 4) then
            i = i+1
          endif
        else
          i = ibnd
          j = jbnd
        endif
        if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                         &
            j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
        b_size = b_size + na
      enddo

      !  receive and decode:
      if (b_size > 0) call dmpi_recv(gbuf, b_size, mpi_io_rank, itag)
      b_size = 0
      do iz=nz1,nz2
        ! skip if this point is not on the relevant task:
        ibnd = krz(1,iz)
        jbnd = krz(2,iz)
        if (bc) then
          i = iga(1)+int((ibnd-1)/iga(2),4)
          j = jga(1)+int((jbnd-1)/jga(2),4)
          if     (krz(3,iz) == 1) then
            j = j-1
          elseif (krz(3,iz) == 2) then
            i = i-1
          elseif (krz(3,iz) == 3) then
            j = j+1
          elseif (krz(3,iz) == 4) then
            i = i+1
          endif
        else
          i = ibnd
          j = jbnd
        endif
        if (i < dd(ia)%low_i .or. i > dd(ia)%up_i .or.                         &
            j < dd(ia)%low_j .or. j > dd(ia)%up_j       ) cycle
        do ic=1,na
          b_size   = b_size + 1
          a(ic,iz) = gbuf(b_size)
        enddo
      enddo

    endif

  end subroutine dmpi_scatter_bnd_data_2

!===============================================================================

  subroutine dmpi_barrier(comm)
    implicit none
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm
    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif
    call MPI_Barrier(mycomm, ierr)
  end subroutine dmpi_barrier

!===============================================================================

  subroutine dmpi_comm_error( sr, sr_error )
    implicit none

    integer(4), intent(in) :: sr, sr_error

    if (sr == 1) then
      write (iu06,'(a)')    ' MPI_Send failed.'
      write (iu06,'(a,i4)') ' MPI_Send Error =  ', sr_error
    elseif (sr == 2) then
      write (iu06,'(a)')    ' MPI_Recv failed.'
      write (iu06,'(a,i4)') ' MPI_Recv Error =  ', sr_error
    else
      write (iu06,'(a)')    ' MPI_Reduce failed.'
      write (iu06,'(a,i4)') ' MPI_Reduce Error =  ', sr_error
    endif
    call exitme(1,'Bailing out due to MPI errors')
  end subroutine dmpi_comm_error

!===============================================================================

  subroutine dmpi_send_real8( sbuf, bsize, trecv, itag, comm )
    implicit none

    integer(4), intent(in) :: bsize, trecv, itag
    real(8),    intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif

    call MPI_Send(sbuf, bsize, MPI_REAL8, trecv, itag, mycomm, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 1, ierr )
  end subroutine dmpi_send_real8

  subroutine dmpi_isend_real8( sbuf, bsize, trecv, itag, ireq, comm )
    implicit none
    integer(4), intent(in) :: bsize, trecv, itag
    integer(4), intent(out) :: ireq
    real(8),    intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm

    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif

    call MPI_ISEND(sbuf, bsize, MPI_REAL8, trecv, itag, mycomm, ireq, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 1, ierr )
  end subroutine dmpi_isend_real8


  subroutine dmpi_send_integer4( sbuf, bsize, trecv, itag, comm )
    implicit none

    integer(4), intent(in) :: bsize, trecv, itag
    integer(4), intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif
    call MPI_Send(sbuf, bsize, MPI_integer4, trecv, itag, mycomm, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 1, ierr )
  end subroutine dmpi_send_integer4

  subroutine dmpi_isend_integer4( sbuf, bsize, trecv, itag, ireq, comm )
    implicit none
    integer(4), intent(in)  :: bsize, trecv, itag
    integer(4), intent(out) :: ireq
    integer(4), intent(in)  :: sbuf(:)
    integer(4), intent(in), optional :: comm

    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif

    call MPI_ISEND(sbuf, bsize, MPI_INTEGER4, trecv, itag, mycomm, ireq, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 1, ierr )
  end subroutine dmpi_isend_integer4


  subroutine dmpi_send_logical( sbuf, bsize, trecv, itag )
    implicit none

    integer(4), intent(in) :: bsize, trecv, itag
    logical,    intent(in) :: sbuf(:)

    call MPI_Send(sbuf, bsize, MPI_LOGICAL, trecv, itag, mpi_comm_model, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 1, ierr )
  end subroutine dmpi_send_logical

!===============================================================================

  subroutine dmpi_irecv_real8( rbuf, bsize, tsend, itag, ireq, comm )
    implicit none

    integer(4), intent(in)  :: bsize, tsend, itag
    integer(4), intent(out) :: ireq
    real(8),    intent(out) :: rbuf(:)
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif

    call MPI_IRecv(rbuf, bsize, MPI_REAL8, tsend, itag, mycomm, ireq, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_irecv_real8

!===============================================================================

  subroutine dmpi_recv_real8( rbuf, bsize, tsend, itag, comm )
    implicit none

    integer(4), intent(in)  :: bsize, tsend, itag
    real(8),    intent(out) :: rbuf(:)
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif

    call MPI_Recv(rbuf, bsize, MPI_REAL8, tsend, itag, mycomm, istatus, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_recv_real8

  subroutine dmpi_recv_integer4( rbuf, bsize, tsend, itag )
    implicit none

    integer(4), intent(in)  :: bsize, tsend, itag
    integer(4), intent(out) :: rbuf(:)

    call MPI_Recv(rbuf, bsize, MPI_integer4, tsend, itag, mpi_comm_model,      &
                  istatus, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_recv_integer4

  subroutine dmpi_recv_logical( rbuf, bsize, tsend, itag )
    implicit none
    integer(4), intent(in)  :: bsize, tsend, itag
    logical,    intent(out) :: rbuf(:)

    call MPI_Recv(rbuf, bsize, MPI_LOGICAL, tsend, itag, mpi_comm_model,       &
                  istatus, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_recv_logical

!===============================================================================

  subroutine dmpi_bcast_real8( bbuf, bsize, iroot, comm )
    implicit none

    integer(4), intent(in)    :: bsize, iroot
    real(8),    intent(inout) :: bbuf(:)
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif
    call MPI_BCAST(bbuf, bsize, MPI_REAL8, iroot, mycomm, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_bcast_real8

  subroutine dmpi_bcast_integer4( bbuf, bsize, iroot, comm )
    implicit none

    integer(4), intent(in)    :: bsize, iroot
    integer(4), intent(inout) :: bbuf(:)
    integer(4), intent(in), optional :: comm
    integer(4) :: mycomm

    if (present(comm)) then
      mycomm=comm
    else
      mycomm=mpi_comm_model
    endif
    call MPI_BCAST(bbuf, bsize, MPI_integer4, iroot, mycomm, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_bcast_integer4

  subroutine dmpi_bcast_logical( bbuf, bsize, iroot )
    implicit none
    integer(4), intent(in)    :: bsize, iroot
    logical,    intent(inout) :: bbuf(:)

    call MPI_BCAST(bbuf, bsize, MPI_LOGICAL, iroot, mpi_comm_model, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_bcast_logical

  subroutine dmpi_bcast_char( bbuf, bsize, iroot )
    implicit none
    integer(4),   intent(in)    :: bsize, iroot
    character(*), intent(inout) :: bbuf

    call MPI_BCAST(bbuf, bsize, MPI_CHARACTER, iroot, mpi_comm_model, ierr)
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 2, ierr )
  end subroutine dmpi_bcast_char

!===============================================================================

  subroutine dmpi_broadcast_met_info(cm, tm, it, md)
    implicit none
    character(*), intent(inout)           :: cm
    real(8),      intent(inout)           :: tm
    integer(4),   intent(inout), optional :: it
    character(*), intent(inout), optional :: md

    ! dummy buffers to survive MPI_BCAST calls:
    real(8)    :: rbuf(1)
    integer(4) :: ibuf(1)

    call dmpi_bcast( cm, len_trim(cm), mpi_io_rank )
    rbuf(1) = tm
    call dmpi_bcast( rbuf, 1, mpi_io_rank )
    tm = rbuf(1)
    if (present(it)) then
      ibuf(1) = it
      call dmpi_bcast( ibuf, 1, mpi_io_rank )
      it = ibuf(1)
    endif
    if (present(md)) call dmpi_bcast( md, len_trim(md), mpi_io_rank ) 
  end subroutine dmpi_broadcast_met_info

  subroutine dmpi_broadcast_met_data(iw2, pl, wu, wv, at, hu, cl, pr, pe)
    implicit none
    integer(4), intent(in)                   :: iw2
    logical,    intent(in)                   :: pe
    real(8),    intent(inout), dimension(0:) :: pl, wu, wv, at, hu, cl, pr

    call dmpi_bcast( pl, iw2+1, mpi_io_rank )
    call dmpi_bcast( wu, iw2+1, mpi_io_rank )
    call dmpi_bcast( wv, iw2+1, mpi_io_rank )
    call dmpi_bcast( at, iw2+1, mpi_io_rank )
    call dmpi_bcast( hu, iw2+1, mpi_io_rank )
    call dmpi_bcast( cl, iw2+1, mpi_io_rank )
    if (pe) call dmpi_bcast( pr, iw2+1, mpi_io_rank )
  end subroutine dmpi_broadcast_met_data

!===============================================================================

  subroutine dmpi_broadcast_logical( l )
    implicit none
    logical, intent(inout) :: l

    logical :: lbuf(1)

    lbuf(1) = l
    call dmpi_bcast( lbuf, 1, mpi_io_rank )
    l = lbuf(1)
  end subroutine dmpi_broadcast_logical

!===============================================================================

  subroutine dmpi_broadcast_bnd_data(nc, kmx, nz, a)
    implicit none
    integer(4), intent(in)    :: nc, kmx, nz
    real(8),    intent(inout) :: a(:,:,0:)      ! (1:nc,1:kmx,0:nz)

    integer(4) :: it, iam, itag, b_size, i, iz, ik, ic

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  nothing to do here? ----------------------------------------------------
    if (nz < 1 .or. nz /= nbnd(iam)) return

    !-  misc vars --------------------------------------------------------------
    itag   = 1
    b_size = nz*kmx*nc

    !-  send from I/O task to all tasks needing data ---------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam .and. nbnd(it) > 0) then
          !  encode and send:
          i = 0
          do iz=1,nz
            do ik=1,kmx
              do ic=1,nc
                i = i + 1
                gbuf(i) = a(ic,ik,iz)
              enddo
            enddo
          enddo
          call dmpi_send(gbuf, b_size, it-1, itag)
        endif
      enddo

    !-  if I need data, receive from I/O task ----------------------------------
    else
      !  receive and decode:
      call dmpi_recv(gbuf, b_size, mpi_io_rank, itag)
      i = 0
      do iz=1,nz
        do ik=1,kmx
          do ic=1,nc
            i = i + 1
            a(ic,ik,iz) = gbuf(i)
          enddo
        enddo
      enddo

    endif

  end subroutine dmpi_broadcast_bnd_data

!===============================================================================

  subroutine dmpi_broadcast_bnd_data_2(nc, nz, a)
    implicit none
    integer(4), intent(in)    :: nc, nz
    real(8),    intent(inout) :: a(:,0:)      ! (1:nc,0:nz)

    integer(4) :: it, iam, itag, b_size, i, iz, ic

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  nothing to do here? ----------------------------------------------------
    if (nz < 1 .or. nz /= nbnd(iam)) return

    !-  misc vars --------------------------------------------------------------
    itag   = 1
    b_size = nz*nc

    !-  send from I/O task to all tasks needing data ---------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
        if (it /= iam .and. nbnd(it) > 0) then
          !  encode and send:
          i = 0
          do iz=1,nz
            do ic=1,nc
              i = i + 1
              gbuf(i) = a(ic,iz)
            enddo
          enddo
          call dmpi_send(gbuf, b_size, it-1, itag)
        endif
      enddo

    !-  if I need data, receive from I/O task ----------------------------------
    else
      !  receive and decode:
      call dmpi_recv(gbuf, b_size, mpi_io_rank, itag)
      i = 0
      do iz=1,nz
        do ic=1,nc
          i = i + 1
          a(ic,iz) = gbuf(i)
        enddo
      enddo

    endif

  end subroutine dmpi_broadcast_bnd_data_2

!===============================================================================

  subroutine dmpi_broadcast_bnd_1d(nz, a, nb1, nb2, b)
    implicit none
    integer(4), intent(in)              :: nz
    real(8),    intent(inout)           :: a(0:)  ! (0:nz)
    integer(4), intent(inout), optional :: nb1, nb2
    logical,    intent(in),    optional :: b

    integer(4) :: it, iam, itag, b_size
    integer(4) :: ibuf(1) 
    logical    :: bb

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1

    !-  figure out what to do --------------------------------------------------
    if (present(b)) then
      bb = b   ! allow for a special case selected by b=F
    else
      bb = .true.
    endif

    !-  nothing to do here? ----------------------------------------------------
    if (nz < 1 .or. nz /= nbnd(iam)) return

    !-  misc vars --------------------------------------------------------------
    itag   = 1
    b_size = nz + 1

    !-  send from I/O task to all tasks needing data ---------------------------
    if (mpi_io_serial) then
      do it=1,mpi_size
! FIXME: This one doesn't always work, but is needed in cases without 
!        any open boundary points. ...
!       if (it /= iam) then
        if (it /= iam .and. nbnd(it) > 0) then
          if (bb) call dmpi_send(a, b_size, it-1, itag)
          if (present(nb1)) then
            ibuf(1) = nb1
            call dmpi_send(ibuf, (1), it-1, itag)
          endif
          if (present(nb2)) then
            ibuf(1) = nb2
            call dmpi_send(ibuf, (1), it-1, itag)
          endif
        endif
      enddo

    !-  if I need data, receive from I/O task ----------------------------------
    else
      if (bb) call dmpi_recv(a, b_size, mpi_io_rank, itag)
      if (present(nb1)) then
        call dmpi_recv(ibuf, (1), mpi_io_rank, itag)
        nb1 = ibuf(1) 
      endif
      if (present(nb2)) then
        call dmpi_recv(ibuf, (1), mpi_io_rank, itag)
        nb2 = ibuf(1) 
      endif

    endif

  end subroutine dmpi_broadcast_bnd_1d

!===============================================================================

  subroutine dmpi_allreduce_real8_scalar( r, op )
    implicit none
    real(8),      intent(inout) :: r
    character(3), intent(in)    :: op

    real(8)    :: recv(1), send(1)
    integer(4) :: mycomm

    mycomm = mpi_comm_model

    send(1) = r

    if (op == 'max' .or. op == 'MAX') then
      call MPI_Allreduce(send, recv, 1, MPI_REAL8, MPI_MAX, mycomm, ierr)
    else
      ! sorry, only max is implmeneted at present.
    endif
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 3, ierr )

    r = recv(1)
  end subroutine dmpi_allreduce_real8_scalar

  subroutine dmpi_allreduce_real8( r, bsize, op )
    implicit none
    real(8),      intent(inout) :: r(1:)
    integer(4),   intent(in)    :: bsize
    character(3), intent(in)    :: op

    real(8)    :: recv(1:bsize)
    integer(4) :: mycomm

    mycomm = mpi_comm_model

    if (op == 'max' .or. op == 'MAX') then
      call MPI_Allreduce(r, recv, bsize, MPI_REAL8, MPI_MAX, mycomm, ierr)
    else
      ! sorry, only max is implmeneted at present.
    endif
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 3, ierr )

    r(1:bsize) = recv(1:bsize)
  end subroutine dmpi_allreduce_real8

  subroutine dmpi_allreduce_integer4_scalar( i, op )
    implicit none
    integer(4),   intent(inout) :: i
    character(3), intent(in)    :: op

    integer(4) :: recv(1), send(1), mycomm

    mycomm = mpi_comm_model

    send(1) = i

    if (op == 'max' .or. op == 'MAX') then
      call MPI_Allreduce(send, recv, 1, MPI_INTEGER4, MPI_MAX, mycomm, ierr)
    else
      ! sorry, only max is implmeneted at present.
    endif
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 3, ierr )

    i = recv(1)
  end subroutine dmpi_allreduce_integer4_scalar

  subroutine dmpi_allreduce_integer4( i, bsize, op )
    implicit none
    integer(4),   intent(inout) :: i(1:)
    integer(4),   intent(in)    :: bsize
    character(3), intent(in)    :: op

    integer(4) :: recv(1:bsize), mycomm

    mycomm = mpi_comm_model

    if (op == 'max' .or. op == 'MAX') then
      call MPI_Allreduce(i, recv, bsize, MPI_INTEGER4, MPI_MAX, mycomm, ierr)
    else
      ! sorry, only max is implmeneted at present.
    endif
    if (ierr /= MPI_SUCCESS) call dmpi_comm_error( 3, ierr )

    i(1:bsize) = recv(1:bsize)
  end subroutine dmpi_allreduce_integer4

!===============================================================================

#else
  !-----------------------------------------------------------------------------
  !
  ! Please note: The subroutines below are just interfaces that enables 
  !              compilation of the source code even without MPI.  We will allow
  !              unused variables to be present in our code for these intent(in)
  !              arguments here and 
  !                                  ONLY HERE
  !
  !-----------------------------------------------------------------------------

  subroutine dmpi_gather_global_srf (ia,ind_l,mmk_l,mmk,a,trecv)
    integer(4), intent(in) :: ia, trecv
    integer(4), intent(in) :: mmk(:,0:,0:), mmk_l(:,0:,0:), ind_l(:,:)
    real(8),    intent(in) :: a(0:)
    return
  end subroutine dmpi_gather_global_srf

  subroutine dmpi_gather_global_all (ia,kmx,ind_l,mmk_l,mmk,a,trecv)
    integer(4), intent(in) :: ia, kmx, trecv
    integer(4), intent(in) :: mmk(:,0:,0:), mmk_l(:,0:,0:), ind_l(:,:)
    real(8),    intent(in) :: a(0:)
  end subroutine dmpi_gather_global_all

  subroutine dmpi_gather_bnd_data_out(ia, nc, kmx, nz1, nz2, krz, a)
    implicit none
    integer(4), intent(in) :: ia, nc, kmx, nz1, nz2, krz(:,0:)
    real(8),    intent(in) :: a(:,:,0:)      ! (1:nc,1:kmx,0:nz)
    return
  end subroutine dmpi_gather_bnd_data_out

  subroutine dmpi_gather_bnd_data(ia, nc, kmx, nz1, nz2, krz, a, iga, jga, i, F)
    implicit none
    integer(4), intent(in) :: ia, nc, kmx, nz1, nz2, krz(:,0:)
    integer(4), intent(in) :: iga(:), jga(:), i, F
    real(8),    intent(in) :: a(:,:,0:)      ! (1:nc,1:kmx,0:nz)
    return
  end subroutine dmpi_gather_bnd_data

  subroutine dmpi_gather_bnd_data_2(ia, na, nz1, nz2, krz, a, iga, jga)
    implicit none
    integer(4), intent(in)           :: ia, na, nz1, nz2, krz(:,0:)
    integer(4), intent(in), optional :: iga(:), jga(:)
    real(8),    intent(in)           :: a(:,0:)   ! (1:nc,0:nz)
    return
  end subroutine dmpi_gather_bnd_data_2

  subroutine dmpi_barrier()
    implicit none
    return
  end subroutine dmpi_barrier

  subroutine dmpi_distribute_halo_nb(kmx, mmk, kh, ia,                         &
                                     irbuf, isbuf, itagoff, FLAG,              &
                                     a, b, c, d, e, f)
    implicit none
    integer(4), intent(in)           :: ia, kmx, mmk(:,0:,0:), kh(0:)
    integer(4), intent(in)           :: irbuf, isbuf, itagoff, FLAG
    real(8),    intent(in)           :: a(0:)
    real(8),    intent(in), optional :: b(0:)
    real(8),    intent(in), optional :: c(0:)
    real(8),    intent(in), optional :: d(0:)
    real(8),    intent(in), optional :: e(0:)
    real(8),    intent(in), optional :: f(0:)
    return
  end subroutine dmpi_distribute_halo_nb

  subroutine dmpi_distribute_halo_nb_col (kmx, msrf, mcol, kh, ia,             &
                                          irbuf, isbuf, itagoff, FLAG,         &
                                          a, b, c, d, e, f)
    implicit none
    integer(4), intent(in)           :: ia, kmx, msrf(0:,0:), mcol(0:), kh(0:)
    integer(4), intent(in)           :: irbuf, isbuf, itagoff, FLAG
    real(8),    intent(in)           :: a(0:)
    real(8),    intent(in), optional :: b(0:)
    real(8),    intent(in), optional :: c(0:)
    real(8),    intent(in), optional :: d(0:)
    real(8),    intent(in), optional :: e(0:)
    real(8),    intent(in), optional :: f(0:)
    return
  end subroutine dmpi_distribute_halo_nb_col

  subroutine dmpi_distribute_halo_TF_nb2(kmx, mmk, kh, ia,                     &
                                         irbuf, isbuf, itagoff, FLAG,          &
                                         a, nci)
    implicit none
    integer(4), intent(in) :: ia, kmx, irbuf, isbuf, itagoff, FLAG , nci
    real(8),    intent(in) :: a(:,0:)
    integer(4), intent(in) :: mmk(:,0:,0:), kh(0:)
    return
  end subroutine dmpi_distribute_halo_TF_nb2

  subroutine dmpi_distribute_halo_TF_nb(kmx, mmk, kh, ia,                      &
                                        irbuf, isbuf, itagoff, FLAG,           &
                                        a, b, c)
    implicit none
    integer(4), intent(in)           :: ia, kmx, irbuf, isbuf, itagoff, FLAG
    real(8),    intent(in)           :: a(:,0:)
    real(8),    intent(in), optional :: b(:,0:), c(:,0:)
    integer(4), intent(in)           :: mmk(:,0:,0:), kh(0:)
    return
  end subroutine dmpi_distribute_halo_TF_nb

  subroutine dmpi_distribute_halo_log(mmk, ia, a)
    implicit none
    integer(4), intent(in) :: ia
    logical,    intent(in) :: a(0:)
    integer(4), intent(in) :: mmk(:,0:,0:)
    return
  end subroutine dmpi_distribute_halo_log

  subroutine dmpi_broadcast_met_info(cm, tm, it, md)
    implicit none
    character(*), intent(in)           :: cm
    real(8),      intent(in)           :: tm
    integer(4),   intent(in), optional :: it
    character(*), intent(in), optional :: md
    return
  end subroutine dmpi_broadcast_met_info

  subroutine dmpi_broadcast_met_data(iw2, pl, wu, wv, at, hu, cl, pr, pe)
    implicit none
    integer(4), intent(in)                :: iw2
    logical,    intent(in)                :: pe
    real(8),    intent(in), dimension(0:) :: pl, wu, wv, at, hu, cl, pr
    return
  end subroutine dmpi_broadcast_met_data

  subroutine dmpi_broadcast_logical( l )
    implicit none
    logical, intent(inout) :: l
    return
  end subroutine dmpi_broadcast_logical

  subroutine dmpi_broadcast_bnd_data(nc, kmx, nz, a)
    implicit none
    integer(4), intent(in) :: nc, kmx, nz
    real(8),    intent(in) :: a(:,:,0:)      ! (1:nc,1:kmx,0:nz)
    return
  end subroutine dmpi_broadcast_bnd_data

  subroutine dmpi_broadcast_bnd_data_2(nc, nz, a)
    implicit none
    integer(4), intent(in) :: nc, nz
    real(8),    intent(in) :: a(:,0:)      ! (1:nc,0:nz)
    return
  end subroutine dmpi_broadcast_bnd_data_2

  subroutine dmpi_broadcast_bnd_1d(nz, a, nb1, nb2, b)
    implicit none
    integer(4), intent(in)           :: nz
    real(8),    intent(in)           :: a(0:)  ! (0:nz)
    integer(4), intent(in), optional :: nb1, nb2
    logical,    intent(in), optional :: b
    return
  end subroutine dmpi_broadcast_bnd_1d

  subroutine dmpi_scatter_bnd_data(ia, iia, ii, nc, kmx, nz1, nz2, krz, a)
    implicit none
    integer(4), intent(in) :: ia, iia, ii, nc, kmx, nz1, nz2, krz(:,0:)
    real(8),    intent(in) :: a(:,:,0:)           ! (1:nc,1:kmx,0:nz)
    return
  end subroutine dmpi_scatter_bnd_data

  subroutine dmpi_scatter_bnd_data_2(ia, na, nz1, nz2, krz, a, iga, jga)
    implicit none
    integer(4), intent(in)           :: ia, na, nz1, nz2, krz(:,0:)
    integer(4), intent(in), optional :: iga(:), jga(:)
    real(8),    intent(in)           :: a(:,0:)    ! (1:nc,0:nz)
    return
  end subroutine dmpi_scatter_bnd_data_2

  subroutine dmpi_send_real8( sbuf, bsize, trecv, itag, comm )
    implicit none
    integer(4), intent(in) :: bsize, trecv, itag
    real(8),    intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_send_real8

  subroutine dmpi_send_integer4( sbuf, bsize, trecv, itag, comm )
    implicit none
    integer(4), intent(in) :: bsize, trecv, itag
    integer(4), intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_send_integer4

  subroutine dmpi_send_logical( sbuf, bsize, trecv, itag )
    implicit none
    integer(4), intent(in) :: bsize, trecv, itag
    logical,    intent(in) :: sbuf(:)
  end subroutine dmpi_send_logical

  subroutine dmpi_recv_real8( rbuf, bsize, tsend, itag, comm )
    implicit none
    integer(4), intent(in) :: bsize, tsend, itag
    real(8),    intent(in) :: rbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_recv_real8

  subroutine dmpi_recv_integer4( rbuf, bsize, tsend, itag )
    implicit none
    integer(4), intent(in) :: bsize, tsend, itag
    integer(4), intent(in) :: rbuf(:)
  end subroutine dmpi_recv_integer4

  subroutine dmpi_recv_logical( rbuf, bsize, tsend, itag )
    implicit none
    integer(4), intent(in) :: bsize, tsend, itag
    logical,    intent(in) :: rbuf(:)
  end subroutine dmpi_recv_logical

  subroutine dmpi_irecv_real8( rbuf, bsize, tsend, itag, ireq, comm )
    implicit none
    integer(4), intent(in) :: bsize, tsend, itag
    integer(4), intent(in) :: ireq
    real(8),    intent(in) :: rbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_irecv_real8

  subroutine dmpi_isend_real8( sbuf, bsize, trecv, itag, ireq, comm )
    implicit none
    integer(4), intent(in) :: bsize, trecv, itag
    integer(4), intent(in) :: ireq
    real(8),    intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_isend_real8

  subroutine dmpi_isend_integer4( sbuf, bsize, trecv, itag, ireq, comm )
    implicit none
    integer(4), intent(in) :: bsize, trecv, itag
    integer(4), intent(in) :: ireq
    integer(4), intent(in) :: sbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_isend_integer4

  subroutine MPI_Waitall(nreq, ireq, lstatus, ierr)
    implicit none
    integer(4), intent(in) :: nreq, ireq(:), ierr, lstatus(:,:)
  end subroutine MPI_Waitall

  subroutine dmpi_bcast_real8( bbuf, bsize, iroot, comm )
    implicit none
    integer(4), intent(in) :: bsize, iroot
    real(8),    intent(in) :: bbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_bcast_real8

  subroutine dmpi_bcast_integer4( bbuf, bsize, iroot, comm )
    implicit none
    integer(4), intent(in) :: bsize, iroot, bbuf(:)
    integer(4), intent(in), optional :: comm
  end subroutine dmpi_bcast_integer4

  subroutine dmpi_bcast_logical( bbuf, bsize, iroot )
    implicit none
    integer(4), intent(in) :: bsize, iroot
    logical,    intent(in) :: bbuf(:)
  end subroutine dmpi_bcast_logical

  subroutine dmpi_bcast_char( bbuf, bsize, iroot )
    implicit none
    integer(4),   intent(in) :: bsize, iroot
    character(*), intent(in) :: bbuf
  end subroutine dmpi_bcast_char

  subroutine dmpi_allreduce_real8_scalar( r, op )
    implicit none
    real(8),      intent(in) :: r
    character(3), intent(in) :: op
  end subroutine dmpi_allreduce_real8_scalar

  subroutine dmpi_allreduce_real8( r, bsize, op )
    implicit none
    real(8),      intent(in) :: r(1:)
    integer(4),   intent(in) :: bsize
    character(3), intent(in) :: op
  end subroutine dmpi_allreduce_real8

  subroutine dmpi_allreduce_integer4_scalar( i, op )
    implicit none
    integer(4),   intent(in) :: i
    character(3), intent(in) :: op
  end subroutine dmpi_allreduce_integer4_scalar

  subroutine dmpi_allreduce_integer4( i, bsize, op )
    implicit none
    integer(4),   intent(in) :: i(1:)
    integer(4),   intent(in) :: bsize
    character(3), intent(in) :: op
  end subroutine dmpi_allreduce_integer4

#endif

!===============================================================================

  subroutine dmpi_gather_uv_bdr_nb(ia, iia, n1, n2, kr, ijpm, recv, send,      &
                                   kmx, mmk_l, mmk, mmk_c, kh_l, kh, iga, jga, &
                                   uvdir, ii, FLAG,                            &
                                   a_l, a, b_l, b, c_l, c, d_l, d)

    use dmi_omp, only : domp_get_domain, domp_get_thread_no

    implicit none

    logical,      intent(in)    :: recv, send
    integer(4),   intent(in)    :: ia, iia, n1, n2, kr(:,0:), ijpm, kmx
    integer(4),   intent(in)    :: uvdir, ii, FLAG
    integer(4),   intent(in)    :: mmk_l(:,0:,0:), mmk(:,0:,0:), iga(:), jga(:)
    integer(4),   intent(in)    :: mmk_c(:,0:,0:), kh_l(0:), kh(0:)
    real(8),      intent(in)    :: a_l(0:)
    real(8),      intent(inout) :: a(0:)
    real(8),      intent(in),    optional :: b_l(0:), c_l(0:), d_l(0:)
    real(8),      intent(inout), optional :: b(0:), c(0:), d(0:)

    integer(4) :: it, iam, itag, b_size, iz, ik, nuv, ig, jg, il, iu, jl, ju, kb
    integer(4) :: ig1, jg1, ig2, jg2, iii, jjj, iff, jff, midx, i, j, badd
    integer(4) :: btmp, nzl, nzu, mm2, ml2, mmdx, tnum, itr
    integer(4) :: imask, nnreq, ireql, irequ, ireqi

    !-  nothing to do here? ----------------------------------------------------
    if (n2 > 0) then
      nuv = n2 - n1 + 1
    else
      nuv = 0
    endif
    if (nuv < 1) return

    !-  who am I? --------------------------------------------------------------
    iam = mpi_rank + 1
    call domp_get_thread_no( tnum )

    !-  misc vars --------------------------------------------------------------
    if (do_comm) itag = nestitags2(ia)%p(ii,uvdir)
    ig1  = iga(1)
    ig2  = iga(2)
    jg1  = jga(1)
    jg2  = jga(2)
    badd = 1
    if (present(b) .and. present(b_l)) then
      badd = 2
      if (present(c) .and. present(c_l)) then
        badd = 3
        if (present(d) .and. present(d_l)) badd = 4
      endif
    endif
    call domp_get_domain(n1, n2, nzl, nzu, .true.)


    !-  Set up partial receive buffer sizes and   ------------------------------
    !   prepost non blocking receive call to iam:
    if (recv .and. (FLAG == 1 .or. FLAG == 4) .and. do_comm) then
      irequ = itagreq2(itag,1,2,uvdir)
      if (irequ > 0) then
        ireqi = itagreq2(itag,1,1,uvdir)

        do it=1,mpi_size
          bsiz2(tnum,it,1:2,uvdir) = 0

          ! receive but not from myself:
          imask = rsmask2(ia,it)%p(ii,uvdir)
          if (imask /= 0 .and. imask /= 10) then

            !  find buffer size:
            btmp = 0
            do iz=nzl,nzu
              ig = kr(1,iz)
              jg = kr(2,iz)
              if (kr(3,iz) == 3) then
                jg = jg+1
              elseif (kr(3,iz) == 4) then
                ig = ig+1
              endif
              if (.not. ((dd(ia)%low_i <= ig .and. ig <= dd(ia)%up_i) .and.    &
                         (dd(ia)%low_j <= jg .and. jg <= dd(ia)%up_j))  ) cycle
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=-ijpm,jg2-1+ijpm
                jff = j + jjj
                if (jff < mpi_tt(iia,it)%low_j .or.                            &
                    jff > mpi_tt(iia,it)%up_j       ) cycle
                do iii=-ijpm,ig2-1+ijpm
                  iff = i + iii
                  if (iff < mpi_tt(iia,it)%low_i .or.                          &
                      iff > mpi_tt(iia,it)%up_i       ) cycle
                  midx = mmk(1,iff,jff)
                  if (midx > 0) btmp = btmp + badd*min(kmx,kh(midx))
                enddo
              enddo
            enddo
            bsiz2(tnum,it,1,uvdir) = btmp
!$OMP BARRIER
!$OMP MASTER
            b_size = bsiz2(1,it,1,uvdir)
            bsiz2(1,it,2,uvdir) = 0
            do itr=2,mpi_nt
              b_size = b_size + bsiz2(itr,it,1,uvdir)
              bsiz2(itr,it,2,uvdir) = bsiz2(itr-1,it,1,uvdir)                  &
                                    + bsiz2(itr-1,it,2,uvdir)
            enddo

            !  recv a dummy in case of degenerated buffer, FIXME
            if (b_size == 0) b_size = 1

            !  prepost the receive:
            call dmpi_irecv(rbufnb2(ireqi)%p, b_size, it-1, itag, ireq2(ireqi))
            ireqi = ireqi + 1
!$OMP END MASTER
          endif
        enddo
      endif
    endif


    !-  if I have data, send it to the relevant task ---------------------------
    if (send .and. (FLAG == 2 .or. FLAG == 4) .and. do_comm) then
      irequ = itagreq2(itag,2,2,uvdir)
      if (irequ > 0) then
        ireqi = itagreq2(itag,2,1,uvdir)

        do it=1,mpi_size
          ! send but not to myself:
          imask = rsmask2(ia,it)%p(ii,uvdir)
          if (imask /= 0 .and. imask /= 1) then

            !  coarse grid dims on task #it
            il = mpi_tt(ia,it)%low_i
            iu = mpi_tt(ia,it)%up_i
            jl = mpi_tt(ia,it)%low_j
            ju = mpi_tt(ia,it)%up_j

            !  find thread offset:
            offset(tnum,1:2) = 0
            btmp = 0
            do iz=nzl,nzu
              ig = kr(1,iz)
              jg = kr(2,iz)
              if (kr(3,iz) == 3) then
                jg = jg+1
              elseif (kr(3,iz) == 4) then
                ig = ig+1
              endif
              if (ig < il .or. iu < ig) cycle
              if (jg < jl .or. ju < jg) cycle
              if (mmk_c(1,ig,jg) <= 0) cycle
              i = ig2*(ig-ig1) + 1
              j = jg2*(jg-jg1) + 1
              do jjj=-ijpm,jg2-1+ijpm
                jff = j + jjj
                if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                do iii=-ijpm,ig2-1+ijpm
                  iff = i + iii
                  if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                  midx = mmk_l(1,iff,jff)
                  if (midx > 0) btmp = btmp + badd*min(kmx,kh_l(midx))
                enddo
              enddo
            enddo
            offset(tnum,1) = btmp
!$OMP BARRIER
!$OMP MASTER
            b_size = offset(1,1)
            offset(1,2) = 0
            do itr=2,mpi_nt
              b_size = b_size + offset(itr,1)
              offset(itr,2) = offset(itr-1,1) + offset(itr-1,2)
            enddo
!$OMP END MASTER
!$OMP BARRIER

            if (offset(tnum,1) > 0) then
              ! encode and send:
              btmp = offset(tnum,2)
              do iz=nzl,nzu
                ig = kr(1,iz)
                jg = kr(2,iz)
                if (kr(3,iz) == 3) then
                  jg = jg+1
                elseif (kr(3,iz) == 4) then
                  ig = ig+1
                endif
                if (ig < il .or. iu < ig) cycle
                if (jg < jl .or. ju < jg) cycle
                if (mmk_c(1,ig,jg) <= 0) cycle
                i = ig2*(ig-ig1) + 1
                j = jg2*(jg-jg1) + 1
                do jjj=-ijpm,jg2-1+ijpm
                  jff = j + jjj
                  if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
                  do iii=-ijpm,ig2-1+ijpm
                    iff = i + iii
                    if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
                    midx = mmk_l(1,iff,jff)
                    if (midx <= 0) cycle
                    btmp       = btmp + 1
                    sbufnb2(ireqi)%p(btmp) = a_l(midx)
                    if (badd >= 2) then
                      btmp       = btmp + 1
                      sbufnb2(ireqi)%p(btmp) = b_l(midx)
                      if (badd >= 3) then
                        btmp       = btmp + 1
                        sbufnb2(ireqi)%p(btmp) = c_l(midx)
                        if (badd == 4) then
                          btmp       = btmp + 1
                          sbufnb2(ireqi)%p(btmp) = d_l(midx)
                        endif
                      endif
                    endif
                    kb = min(kmx,kh_l(midx))
                    if (kb >= 2) then
                      ml2 = mmk_l(2,iff,jff) - 2
                      sbufnb2(ireqi)%p(btmp+1:btmp+kb-1) = a_l(ml2+2:ml2+kb)
                      btmp = btmp + kb - 1
                      if (badd >= 2) then
                        sbufnb2(ireqi)%p(btmp+1:btmp+kb-1) = b_l(ml2+2:ml2+kb)
                        btmp = btmp + kb - 1
                        if (badd >= 3) then
                          sbufnb2(ireqi)%p(btmp+1:btmp+kb-1) = c_l(ml2+2:ml2+kb)
                          btmp = btmp + kb - 1
                          if (badd == 4) then
                            sbufnb2(ireqi)%p(btmp+1:btmp+kb-1) =               &
                                                               d_l(ml2+2:ml2+kb)
                            btmp = btmp + kb - 1
                          endif
                        endif
                      endif
                    endif
                  enddo
                enddo
              enddo
            endif
!$OMP BARRIER
!$OMP MASTER
            !  send a dummy in case of degenerated buffer, FIXME
            if (b_size == 0) then
              b_size = 1
              sbufnb2(ireqi)%p(1) = dummy_value
            endif
            call dmpi_isend(sbufnb2(ireqi)%p, b_size, it-1, itag, irqs2(ireqi))
!$OMP END MASTER
            ireqi = ireqi + 1
          endif
        enddo
      endif
    endif


    !-  recv on this task from all tasks having any data -----------------------
    if (recv .and. (FLAG == 3 .or. FLAG == 4)) then
      ! copy from local to global array
      if (send .and. recv) then
        it = iam
        do iz=nzl,nzu
          ig = kr(1,iz)
          jg = kr(2,iz)
          if (kr(3,iz) == 3) then
            jg = jg+1
          elseif (kr(3,iz) == 4) then
            ig = ig+1
          endif
          if (ig < dd(ia)%low_i .or. dd(ia)%up_i < ig) cycle
          if (jg < dd(ia)%low_j .or. dd(ia)%up_j < jg) cycle
          if (mmk_c(1,ig,jg) <= 0) cycle
          i = ig2*(ig-ig1) + 1
          j = jg2*(jg-jg1) + 1
          do jjj=-ijpm,jg2-1+ijpm
            jff = j + jjj
            if (jff < dd(iia)%low_j .or. jff > dd(iia)%up_j) cycle
            do iii=-ijpm,ig2-1+ijpm
              iff = i + iii
              if (iff < dd(iia)%low_i .or. iff > dd(iia)%up_i) cycle
              midx = mmk_l(1,iff,jff)
              if (midx <= 0) cycle
              mmdx = mmk(1,iff,jff)
              if (mmdx <= 0) cycle
              a(mmdx) = a_l(midx)
              if (badd >= 2) then 
                b(mmdx) = b_l(midx)
                if (badd >= 3) then
                  c(mmdx) = c_l(midx)
                  if (badd == 4) d(mmdx) = d_l(midx)
                endif
              endif
              kb = min(kmx,kh_l(midx))
              if (kb >= 2) then
                ml2 = mmk_l(2,iff,jff) - 2
                mm2 = mmk(2,iff,jff)   - 2
                a(mm2+2:mm2+kb) = a_l(ml2+2:ml2+kb)
                if (badd >= 2) then
                  b(mm2+2:mm2+kb) = b_l(ml2+2:ml2+kb)
                  if (badd >= 3) then
                    c(mm2+2:mm2+kb) = c_l(ml2+2:ml2+kb)
                    if (badd == 4) d(mm2+2:mm2+kb) = d_l(ml2+2:ml2+kb)
                  endif
                endif
              endif
            enddo
          enddo
        enddo
      endif

      !  wait for the buffers to get filled and then decode
      if (do_comm) then
        irequ = itagreq2(itag,1,2,uvdir)
        if (irequ > 0) then
          ireql = itagreq2(itag,1,1,uvdir)
          nnreq = irequ - ireql + 1
!$OMP MASTER
          call MPI_Waitall(nnreq,ireq2(ireql:irequ),lstat2(:,ireql:irequ),ierr)
!$OMP END MASTER
!$OMP BARRIER

          ireqi = ireql
          do it=1,mpi_size
            ! recieve but not from myself:
            imask = rsmask2(ia,it)%p(ii,uvdir)
            if (imask /= 0 .and. imask /= 10) then
              if (bsiz2(tnum,it,1,uvdir) > 0) then
                ! decode buffer:
                btmp = bsiz2(tnum,it,2,uvdir)
                do iz=nzl,nzu
                  ig = kr(1,iz)
                  jg = kr(2,iz)
                  if (kr(3,iz) == 3) then
                    jg = jg+1
                  elseif (kr(3,iz) == 4) then
                    ig = ig+1
                  endif
                  if (.not.((dd(ia)%low_i <= ig .and. ig <= dd(ia)%up_i) .and. &
                            (dd(ia)%low_j <= jg .and. jg <= dd(ia)%up_j))) cycle
                  if (mmk_c(1,ig,jg) <= 0) cycle
                  i = ig2*(ig-ig1) + 1
                  j = jg2*(jg-jg1) + 1
                  do jjj=-ijpm,jg2-1+ijpm
                    jff = j + jjj
                    if (jff < mpi_tt(iia,it)%low_j .or.                        &
                        jff > mpi_tt(iia,it)%up_j       ) cycle
                    do iii=-ijpm,ig2-1+ijpm
                      iff = i + iii
                      if (iff < mpi_tt(iia,it)%low_i .or.                      &
                        iff > mpi_tt(iia,it)%up_i       ) cycle
                      midx = mmk(1,iff,jff)
                      if (midx <= 0) cycle
                      btmp  = btmp + 1
                      a(midx) = rbufnb2(ireqi)%p(btmp)
                      if (badd >= 2) then
                        btmp  = btmp + 1
                        b(midx) = rbufnb2(ireqi)%p(btmp)
                        if (badd >= 3) then
                          btmp  = btmp + 1
                          c(midx) = rbufnb2(ireqi)%p(btmp)
                          if (badd == 4) then
                            btmp  = btmp + 1
                            d(midx) = rbufnb2(ireqi)%p(btmp)
                          endif
                        endif
                      endif
                      kb = min(kmx,kh(midx))
                      if (kb >= 2) then
                        mm2 = mmk(2,iff,jff) - 2
                        a(mm2+2:mm2+kb) = rbufnb2(ireqi)%p(btmp+1:btmp+kb-1)
                        btmp = btmp + kb - 1
                        if (badd >= 2) then
                          b(mm2+2:mm2+kb) = rbufnb2(ireqi)%p(btmp+1:btmp+kb-1)
                          btmp = btmp + kb - 1
                          if (badd >= 3) then
                            c(mm2+2:mm2+kb) = rbufnb2(ireqi)%p(btmp+1:btmp+kb-1)
                            btmp = btmp + kb - 1
                            if (badd == 4) then
                              d(mm2+2:mm2+kb) =                                &
                                              rbufnb2(ireqi)%p(btmp+1:btmp+kb-1)
                              btmp = btmp + kb - 1
                            endif
                          endif
                        endif
                      endif
                    enddo
                  enddo
                enddo
              endif
              ireqi = ireqi + 1
            endif
          enddo
        endif
      endif
    endif


    !-  complete the send from this task ---------------------------------------
    if (send .and. (FLAG == 3 .or. FLAG == 4) .and. do_comm) then
      !  wait for the send to finish
      irequ = itagreq2(itag,2,2,uvdir)
      if (irequ > 0) then
!$OMP MASTER
        ireql = itagreq2(itag,2,1,uvdir)
        nnreq = irequ - ireql + 1
        call MPI_Waitall(nnreq, irqs2(ireql:irequ), lstat2(:,ireql:irequ), ierr)
!$OMP END MASTER
      endif
    endif

  end subroutine dmpi_gather_uv_bdr_nb

#if defined (MPI)
  subroutine dmpi_halosize_srf(kmx, mm, kh, il, iu, jl, ju, ihl, ihu, jhl, jhu,&
                               nh3, nh2)

    implicit none

    integer(4), intent(in)  :: kmx, mm(0:,0:), kh(:)
    integer(4), intent(in)  :: il, iu, jl, ju, ihl, ihu, jhl, jhu
    integer(4), intent(out) :: nh3, nh2

    integer(4) :: i, j, k

    nh3 = 0
    nh2 = 0
    !  w, n/w, s/w
    do j=jhl,jl-1
      do i=ihl,ihu
        if (mm(i,j) == 0) cycle
        nh2 = nh2 + 1
        nh3=nh3+kh(mm(i,j))
      enddo
    enddo
    !  e, n/e, s/e
    do j=ju+1,jhu
      do i=ihl,ihu
        if (mm(i,j) == 0) cycle
        nh2 = nh2 + 1
        nh3=nh3+kh(mm(i,j))
      enddo
    enddo
    !  n
    do j=jl,ju
      do i=ihl,il-1
        if (mm(i,j) == 0) cycle
        nh2 = nh2 + 1
        nh3=nh3+kh(mm(i,j))
      enddo
    enddo
    !  s
    i = ihu
    do j=jl,ju
      do i=iu+1,ihu
        if (mm(i,j) == 0) cycle
        nh2 = nh2 + 1
        nh3=nh3+kh(mm(i,j))
      enddo
    enddo

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_halosize_srf

  subroutine dmpi_bufsize_srf(kmx, mm, kh, ia, recv, send, rsrf, ssrf )

    implicit none

    integer(4), intent(in)  :: kmx, mm(0:,0:), ia, kh(1:)
    integer(4), intent(out) :: recv(:), send(:), rsrf(:), ssrf(:)

    integer(4) :: wrecv, wsend, wrsrf, wssrf
    integer(4) :: iam, it_recv, itp1, id, nb, idi, nbi
    integer(4) :: i, j, k, ij
    integer(4) :: irl, iru, jrl, jru, isl, isu, jsl, jsu
    integer(4) :: il, iu, jl, ju, ihl, ihu, jhl, jhu

    iam = mpi_rank + 1

    il  = mpi_tt(ia,iam)%low_i
    iu  = mpi_tt(ia,iam)%up_i
    jl  = mpi_tt(ia,iam)%low_j
    ju  = mpi_tt(ia,iam)%up_j
    ihl = mpi_tt(ia,iam)%low_hi
    ihu = mpi_tt(ia,iam)%up_hi
    jhl = mpi_tt(ia,iam)%low_hj
    jhu = mpi_tt(ia,iam)%up_hj

    do it_recv=0,mpi_size-1
      itp1 = it_recv + 1

      ! I should not send to myself, but receive from up to 8 of my neighbours:
      if (itp1 == iam) then
        do id=1,max_nb
          ! nb is the sender, iam is the receiver
          nb = mpi_tt(ia,iam)%hnb(id)
          if (nb > 0) then
            ! recv-buffer coordinates:
            select case (id)
              case (1) ! west
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhl
                jru = jhl
                ij  = max(irl,ihl) - 1
              case (2) ! north
                irl = ihl
                iru = ihl
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
              case (3) ! east
                irl = mpi_tt(ia,nb)%low_i
                iru = mpi_tt(ia,nb)%up_i
                jrl = jhu
                jru = jhu
                ij  = max(irl,ihl) - 1
              case (4) ! south
                irl = ihu
                iru = ihu
                jrl = mpi_tt(ia,nb)%low_j
                jru = mpi_tt(ia,nb)%up_j
                ij  = max(jrl,jhl) - 1
              case (5) ! north-west
                irl = ihl
                iru = ihl
                jrl = jhl
                jru = jhl
                ij  =  0
              case (6) ! north-east
                irl = ihl
                iru = ihl
                jrl = jhu
                jru = jhu
                ij  =  0
              case (7) ! south-east
                irl = ihu
                iru = ihu
                jrl = jhu
                jru = jhu
                ij  =  0
              case (8) ! south-west
                irl = ihu
                iru = ihu
                jrl = jhl
                jru = jhl
                ij  =  0
            end select

            ! calculate No. of wet points to recv/decode:
            wrsrf = 0
            wrecv = 0
            do j=max(jrl,jhl),min(jru,jhu)
              do i=max(irl,ihl),min(iru,ihu)
                ij = ij + 1
                if (decode(ia,id)%mask(ij) == 0) cycle
                wrsrf = wrsrf + 1
                if (mm(i,j) == 0) cycle
                wrecv = wrecv + 1
                wrecv = wrecv + kh(mm(i,j))
              enddo
            enddo
            if (wrecv > 0) then
              rsrf(id) = wrsrf
              recv(id) = wrecv
            endif

          endif
        enddo  ! id

      ! I should not receive from myself, but send to up to 8 of my neighbours:
      else

        do id=1,max_nb
          nb = mpi_tt(ia,itp1)%hnb(id)
          if (nb == iam) then
            ! I am a neighbour to task #itp1 and I must send to him

            ! Find my neighbour which corresponds to task #itp1
            idi = n2n(id)
            nbi = mpi_tt(ia,iam)%hnb(idi)
            if (nbi == itp1) then
              ! send-buffer coordinates:
              select case (idi)
                case (1) ! west
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = jl
                  jsu = jl
                  ij  = max(isl,il) - 1
                case (2) ! north
                  isl = il
                  isu = il
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                case (3) ! east
                  isl = mpi_tt(ia,itp1)%low_hi
                  isu = mpi_tt(ia,itp1)%up_hi
                  jsl = ju
                  jsu = ju
                  ij  = max(isl,il) - 1
                case (4) ! south
                  isl = iu
                  isu = iu
                  jsl = mpi_tt(ia,itp1)%low_hj
                  jsu = mpi_tt(ia,itp1)%up_hj
                  ij  = max(jsl,jl) - 1
                case (5) ! north-west
                  isl = il
                  isu = il
                  jsl = jl
                  jsu = jl
                  ij  =  0
                case (6) ! north-east
                  isl = il
                  isu = il
                  jsl = ju
                  jsu = ju
                  ij  =  0
                case (7) ! south-east
                  isl = iu
                  isu = iu
                  jsl = ju
                  jsu = ju
                  ij  =  0
                case (8) ! south-west
                  isl = iu
                  isu = iu
                  jsl = jl
                  jsu = jl
                  ij  =  0
              end select

              ! calculate No. of wet points to send/encode:
              wssrf = 0
              wsend = 0
              do j=max(jsl,jl),min(jsu,ju)
                do i=max(isl,il),min(isu,iu)
                  ij = ij + 1
                  if (encode(ia,idi)%mask(ij) == 0) cycle
                  wssrf = wssrf + 1
                  if (mm(i,j) == 0) cycle
                  wsend = wsend + 1
                  wsend = wsend + kh(mm(i,j))
                enddo
              enddo
              if (wsend > 0) then
                ssrf(idi) = wssrf
                send(idi) = wsend
              endif

            endif  ! nbi == itp1
          endif    ! nb == iam
        enddo      ! id

      endif
    enddo   ! it_recv, itp1

    ! That's it! ---------------------------------------------------------------

  end subroutine dmpi_bufsize_srf

#endif
!===============================================================================

end module dmi_mpi

