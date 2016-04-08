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

module io_mpi
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !- modules 
#if defined (MPI)
  use mpi
#endif
  use cmod_mem, only : cmr1, cmi1, cmi3

  implicit none

  !- Some definitions
  integer(4), parameter, public :: msgtag         = 100
  integer(4), parameter, public :: scattag        = 200
  integer(4), parameter, public :: valtag         = 500

  !- Public MPI related variables
  logical,    public, save      :: inhbm
  integer(4), public, save      :: my_rank, my_size, color, comm_model, comm_io
  integer(4), public, save      :: comm_inter, first_model_task
  integer(4), public, save      :: ireq_msg, ireq_val, my_iotask
  type(cmi1), allocatable, public, save :: ireq_scatter(:,:)

  !- Private module variables and parameters
  integer(4), private, allocatable, save :: iw2io(:), iw3io(:)
  !- Extra communicators
  integer(4), private, save :: comm_all, comm_dupworld

  !- more MPI variables
  integer(4), allocatable, private, save :: filetype(:,:), restart_filetype(:)
  integer(4), allocatable, private, save :: filetype_coupler(:,:)
  integer(4), allocatable, private, save :: restart_filetype_coupler(:)
  real(8),    allocatable, private, save :: rbuf(:), wbuf(:)
  type(cmi3), allocatable, private, save :: mo(:), mop(:)
  type(cmr1), allocatable, private, save :: hbmbuf(:,:,:)
  integer(4), allocatable, private, save :: n2dt(:), n3dlb(:), n3ds(:),        &
                                            n3dlen(:), n3dub(:)
  integer(4), allocatable, private, save :: wrk1(:), wrk2(:)
  integer(4), allocatable, private, save :: commtags(:,:,:) 

  !- Public routines
  public  :: set_comm_models
contains

#if defined (MPI)
!===============================================================================
  

  subroutine set_comm_models(io_server, io_ntask)
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    !- modules
    use dmi_mpi, only : mpi_rank, mpi_size, dmpi_reinit, dmpi_reset_size

    !- directives
    implicit none

    !- arguments
    integer(4), intent(in) :: io_ntask
    logical,    intent(in) :: io_server

    !- locals
    integer(4) :: i, ierr

    inhbm = .true.

    call dmpi_reinit(comm_model)
 
  end subroutine set_comm_models

!===============================================================================

#else
!-------------------------------------------------------------------------------
! The subroutines below are interfaces which are used when running without MPI.
! Obviously, the ioserver is not activated in this case.
!-------------------------------------------------------------------------------
  subroutine set_comm_models(io_server, io_ntask)
    use dmi_mpi, only : dmpi_reinit
    implicit none
    integer(4), intent(in) :: io_ntask
    logical,    intent(in) :: io_server
    integer(4) :: dummy
    inhbm = .true.
    call dmpi_reinit(dummy)
  end subroutine set_comm_models
  !---------
#endif

end module io_mpi
