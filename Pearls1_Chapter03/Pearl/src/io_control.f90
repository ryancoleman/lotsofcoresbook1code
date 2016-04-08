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

module io_control
  !-----------------------------------------------------------------------------
  ! Module for controlling the asynchronous IO
  !
  ! TODO: * Need to condider isend-wait if more than one fh is choosen
  !       * Need to write sponge files
  !       * IO-server should take care of the Validation
  !       * Headers in restart and output files
  !       * Read restart (and sponge in the restart files)
  !       * Metread? 
  !       * wlfiles?
  !       * iota* variables
  !       * boundary/rivers
  !-----------------------------------------------------------------------------
  use cmod_mem, only : cmr1
  implicit none

  !- Public modules variables and arrays
  logical,    public, save :: io_server, in_hbm, io_debug, io_readrestart
  integer(4), public, save :: io_ntask, io_maxfh, io_maxfh_coupler

  !- Private modules variables and arrays
  integer(4), parameter, private :: maxiotasks = 5 
  
  !- codes for io_signals
  integer(4), parameter, private :: io_finish         = 1 ! prepare to die
  integer(4), parameter, private :: io_open           = 2 ! open file
  integer(4), parameter, private :: io_write          = 3 ! write output
  integer(4), parameter, private :: io_read           = 4 ! read input
  integer(4), parameter, private :: io_close          = 5 ! close file
  integer(4), parameter, private :: io_writerestart   = 6 ! write restart
  integer(4), parameter, private :: io_validate       = 7 ! write validation 
  integer(4), parameter, private :: io_open_c         = 8 ! open coupler file
  integer(4), parameter, private :: io_write_c        = 9  ! open coupler file
  integer(4), parameter, private :: io_close_c        = 10 ! open coupler file
  integer(4), parameter, private :: io_writerestart_c = 11 ! write coupler rest
  !- unused messages passed  
  real(8),    parameter, private :: unused = -1.0_8

  !- arrays for ouput information
  integer(4), allocatable, private, save :: output_fields(:,:), output_freq(:),&
                                            output_fh(:,:), output_dims(:,:),  &
                                            restart_fields(:), restart_dims(:),&
                                            restart_fh(:)
  character(len=20), allocatable, private, save :: output_fname(:)

  !- arrays for coupler ouput information
  integer(4), allocatable, private, save :: output_fields_coupler(:,:),        &
                                            output_freq_coupler(:),            &
                                            output_fh_coupler(:,:),            &
                                            output_dims_coupler(:,:),          &
                                            restart_fields_coupler(:),         &
                                            restart_dims_coupler(:),           &
                                            restart_fh_coupler(:)
  character(len=20), allocatable, private, save :: output_fname_coupler(:)

  !- bookkepping arrays for IO-server
  integer(4), allocatable, private, save :: io_rank(:)
  integer(4), allocatable, private, save :: bufsize2d(:,:), bufsize3d(:,:)
  integer(4), allocatable, private, save :: sbufsize2d(:), sbufsize3d(:)
  integer(4),              private, save :: hbm_lb, hbm_ub
  real(8),                 private, save :: msg_io(4)
  logical,                 private, save :: wait_msgsend       = .false.
  logical,    allocatable, private, save :: wait_outscatter(:)
  logical,                 private, save :: wait_restscatter   = .false.
  logical,                 private, save :: wait_valscatter    = .false.
  logical,    allocatable, private, save :: wait_outscatter_c(:)
  logical,                 private, save :: wait_restscatter_c = .false.

  !- Public routines 
  public  :: ios_init, ios_usage
             
contains

!===============================================================================

!===============================================================================

  subroutine ios_init()
    !---------------------------------------------------------------------------
    ! setting up MPI-communication
    !---------------------------------------------------------------------------
    use io_mpi,         only : set_comm_models, inhbm, my_iotask, my_rank
    use dmi_mpi,        only : mpi_size
    use dmi_mpi_global, only : iu06
    implicit none
    !- Let's split the tasks into two groups HBM and IO
    call set_comm_models(io_server, io_ntask)
  end subroutine ios_init

!===============================================================================

  subroutine ios_usage(iu06,coupler_active, coupler_3d, coupler_2d, ldefault)
    ! modules used
    use io_subs,     only : flush_unit, io_new_unit
    use io_mpi,      only : color, first_model_task
    use io_miscsubs, only : f_tempdat, f_coupler,                              &
                            max_io_fields_def, set_predefined, set_restart
    use dmi_mpi,     only : mpi_rank, mpi_size
    use exits,       only : exitme
    use constants,   only : coldst

    implicit none

    !- Arguments
    integer(4), intent(in)           :: iu06
    logical,    intent(in)           :: coupler_active
    integer(4), intent(in)           :: coupler_3d, coupler_2d
    logical,    intent(in), optional :: ldefault
    
    !- Local variables and params
    character(5), parameter :: w = '(a57)'
    integer(4),  parameter :: uninitialized=-99999
    integer(4)    :: i, ios, j, f, funit
    integer(4)    :: io_rank_in(maxiotasks)
    integer(4)    :: fields_in(max_io_fields_def), freq_in
    character(10) :: filename
   
    !- MAybe all these should go in module header instead
    namelist /iospec_general/     io_server, io_ntask, io_maxfh, io_debug,     &
                                  io_readrestart, io_maxfh_coupler
    namelist /iospec_rank/        io_rank_in
    namelist /iospec_out/         filename, fields_in,freq_in
    namelist /iospec_out_coupler/ filename, fields_in,freq_in

    !-----Defaults -------------------------------------------------------------
    io_server        = .false.
    io_ntask         = 0
    io_debug         = .false.
    io_readrestart   = .true.
    io_maxfh         = 1
    io_maxfh_coupler = 0
    do i = 1, maxiotasks 
      io_rank_in(i) = uninitialized 
    enddo
 
    if (present(ldefault)) then
      return
    endif

  end subroutine ios_usage

!===============================================================================

end module io_control
