module probin_module

  use bl_types

  implicit none

  private

  real (kind=dp_t), save, public :: rfire = 0.01d0
  integer, save, public :: verbose = 0
  real (kind=dp_t), save, public :: stop_time = -1.d0
  integer, save, public :: max_step = 1
  real (kind=dp_t), save, public :: cflfac = 0.5d0
  integer, save, public :: cfl_int = 10
  real (kind=dp_t), save, public :: init_shrink = 1.d0
  real (kind=dp_t), save, public :: small_dt = 1.d-30
  real (kind=dp_t), save, public :: max_dt_growth = 1.1d0
  real (kind=dp_t), save, public :: max_dt = 1.0d33
  real (kind=dp_t), save, public :: fixed_dt = -1.0d0
  real (kind=dp_t), save, public :: prob_lo_x = 0.d0
  real (kind=dp_t), save, public :: prob_lo_y = 0.d0
  real (kind=dp_t), save, public :: prob_lo_z = 0.d0
  real (kind=dp_t), save, public :: prob_hi_x = 1.d0
  real (kind=dp_t), save, public :: prob_hi_y = 1.d0
  real (kind=dp_t), save, public :: prob_hi_z = 1.d0
  integer, save, public :: max_grid_size = 64
  integer, save, public :: n_cellx = -1
  integer, save, public :: n_celly = -1
  integer, save, public :: n_cellz = -1
  integer, save, public :: tb_split_dim = 2
  logical, save, public :: tb_collapse_boxes = .false.
  integer, save, public :: tb_idim_more = 2
  integer, save, public :: tb_idim_less = 1
  integer, save, public :: tb_blocksize_x = -1
  integer, save, public :: tb_blocksize_y = 16
  integer, save, public :: tb_blocksize_z = 16

  ! These will be allocated and defined below
  logical,    allocatable, save, public :: pmask(:)
  real(dp_t), allocatable, save, public :: prob_lo(:)
  real(dp_t), allocatable, save, public :: prob_hi(:)

end module probin_module


module runtime_init_module

  use bl_types
  use probin_module

  implicit none

  namelist /probin/ rfire
  namelist /probin/ verbose
  namelist /probin/ stop_time
  namelist /probin/ max_step
  namelist /probin/ cflfac
  namelist /probin/ cfl_int
  namelist /probin/ init_shrink
  namelist /probin/ small_dt
  namelist /probin/ max_dt_growth
  namelist /probin/ max_dt
  namelist /probin/ fixed_dt
  namelist /probin/ prob_lo_x
  namelist /probin/ prob_lo_y
  namelist /probin/ prob_lo_z
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
  namelist /probin/ max_grid_size
  namelist /probin/ n_cellx
  namelist /probin/ n_celly
  namelist /probin/ n_cellz
  namelist /probin/ tb_split_dim
  namelist /probin/ tb_collapse_boxes
  namelist /probin/ tb_idim_more
  namelist /probin/ tb_idim_less
  namelist /probin/ tb_blocksize_x
  namelist /probin/ tb_blocksize_y
  namelist /probin/ tb_blocksize_z

  private

  public :: probin

  public :: runtime_init, runtime_close

contains

  subroutine runtime_init()

    use parallel
    use bl_IO_module
    use bl_error_module
    
    logical    :: lexist, need_inputs
    integer    :: natonce, myproc, nprocs, nsets, myset, iset, ibuff(1)
    integer    :: wakeuppid, waitforpid, tag, un
    real(dp_t) :: pistart, piend, pitotal, pistartall, piendall, pitotalall
    real(dp_t) :: pitotal_max, pitotalall_max

    need_inputs = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !
    ! Don't have more than 64 processes trying to read from disk at once.
    !
    natonce = min(64,parallel_nprocs())
    myproc  = parallel_myproc()
    nprocs  = parallel_nprocs()
    nsets   = ((nprocs + (natonce - 1)) / natonce)
    myset   = (myproc / natonce)
    pistart = 0.d0
    piend   = 0.d0

    pistartall = parallel_wtime()

    ! loop over the processor groups (sets) and read the inputs file.
    ! We look first for an environment variable "PROBIN", then at the
    ! first argument to the executable, and finally for a file named
    ! inputs_SMC
    do iset = 0, nsets-1

       if (myset .eq. iset) then

          pistart = parallel_wtime()
          
          inquire(file = 'inputs_SMC', exist = lexist)
          if ( need_inputs .AND. lexist ) then
             un = unit_new()
             open(unit=un, file = 'inputs_SMC', status = 'old', action = 'read')
             read(unit=un, nml = probin)
             close(unit=un)
             need_inputs = .false.
          end if

          piend = parallel_wtime()

          ibuff(1)  = 0
          wakeuppid = myproc + natonce
          tag       = mod(myproc,natonce)
          
          if (wakeuppid < nprocs) call parallel_send(ibuff, wakeuppid, tag)

       end if

      if (myset .eq. (iset + 1)) then

         tag        = mod(myproc,natonce)
         waitforpid = myproc - natonce

         call parallel_recv(ibuff, waitforpid, tag)
      endif

    end do

    if (need_inputs) then
       call bl_error("Cannot find the 'inputs_SMC' file in current directory")
    end if

    piendall   = parallel_wtime()
    pitotal    = piend - pistart
    pitotalall = piendall - pistartall

    call parallel_reduce(pitotal_max,    pitotal,    MPI_MAX, &
                         proc = parallel_IOProcessorNode())
    call parallel_reduce(pitotalall_max, pitotalall, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    if (parallel_IOProcessor()) then
      print*, "PROBINIT max time   = ", pitotal_max
      print*, "PROBINIT total time = ", pitotalall_max
    endif

    !-------------------------------------------------------------------------
    ! some sanity checks and set some defaults
    !-------------------------------------------------------------------------

    ! initialize pmask
    allocate(pmask(3))
    pmask = .TRUE.

    ! initialize prob_lo and prob_hi
    allocate(prob_lo(3))
    prob_lo(1) = prob_lo_x
    prob_lo(2) = prob_lo_y
    prob_lo(3) = prob_lo_z
    allocate(prob_hi(3))
    prob_hi(1) = prob_hi_x
    prob_hi(2) = prob_hi_y
    prob_hi(3) = prob_hi_z

  end subroutine runtime_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine runtime_close()

    use probin_module

    deallocate(pmask)
    deallocate(prob_lo)
    deallocate(prob_hi)

  end subroutine runtime_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module runtime_init_module
