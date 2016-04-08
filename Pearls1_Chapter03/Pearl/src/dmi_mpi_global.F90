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

module dmi_mpi_global

#if defined (MPI)
  use mpi
#endif

  !- implicit directives -------------------------------------------------------
  implicit none
  private

  !- public vars ---------------------------------------------------------------
  logical,            save, public :: debug_print
  integer(4),         save, public :: iu06, mpi_comm_model
  integer(4),         save, public :: mpi_io_rank = 0
  ! A "Initialisation expression for MPI_DECOMP_FILE padded with blanks"-warning
  ! is accepted here:
  character(LEN=256), save, public :: mpi_decomp_file = 'mpi_decompo.txt'

  !- public methods ------------------------------------------------------------
  public :: dmpi_close, dmpi_abort

contains

!===============================================================================

  subroutine dmpi_close ( flag )
    implicit none

    logical, intent(in) :: flag

#if defined (MPI)
    integer(4) :: ierr

    if (flag) then
      write(*,'(a19)') 'Time to die for MPI'
    else
      write(iu06,'(a19)') 'Time to die for MPI'
    endif
#endif

    if (.not.flag) close (iu06)

#if defined (MPI)
    call mpi_finalize(ierr)
#endif

  end subroutine dmpi_close

  subroutine dmpi_abort ( flag )
    implicit none

    logical, intent(in) :: flag

#if defined (MPI)
    integer(4) :: ierr

    if (flag) then
      write(*,'(a19)') 'Time to abort MPI'
    else
      write(iu06,'(a19)') 'Time to abort MPI'
    endif
#endif

    if (.not.flag) close (iu06)

#if defined (MPI)
    call mpi_abort(MPI_COMM_WORLD, 0, ierr)
#endif

  end subroutine dmpi_abort


!===============================================================================

end module dmi_mpi_global

