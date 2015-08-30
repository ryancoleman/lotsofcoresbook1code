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

module io_subs
  !-----------------------------------------------------------------------------
  !
  ! Module containing methods for handling I/O
  !
  !   io_new_unit()
  !   flush_unit()
  !
  ! Per Berg, DMI.
  !
  !-----------------------------------------------------------------------------

  !- modules -------------------------------------------------------------------
  !  none so far ...

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- interfaces  ---------------------------------------------------------------
  interface flush_unit
    module procedure flush_unit_orig
    module procedure flush_unit_debug
  end interface 

  !- public routines -----------------------------------------------------------
  public :: io_new_unit, flush_unit

  !- private routines ----------------------------------------------------------
  !  none so far ...

  !- public vars ---------------------------------------------------------------
  !  none so far ...

  !- private vars --------------------------------------------------------------
  !  none so far ...

contains

!- Public routines -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  integer(4) function io_new_unit ()

    use exits, only : exitme

    integer(4) :: i
    logical    :: file_is_open
    integer(4), parameter :: io_min=7, io_max=130

    file_is_open = .true.
    i = io_min-1
    do while (file_is_open)
      i = i+1
      if (i > io_max) then
        call exitme('ERROR: Could not find new I/O unit number')
      endif
      inquire(unit=i,opened=file_is_open)
    enddo

    io_new_unit = i

  end function io_new_unit

  !-----------------------------------------------------------------------------

  subroutine flush_unit_orig ( lun )
    use port_flush, only : portable_flush
    use exits, only : exitme
    integer(4), intent(in) :: lun
    logical :: file_is_open

    inquire(unit=lun,opened=file_is_open)
    if (.not.file_is_open) call exitme('ERROR: Cannot flush non-opened file')
    call portable_flush( lun ) 
  end subroutine flush_unit_orig

  !-----------------------------------------------------------------------------

  subroutine flush_unit_debug ( lun, debug )
    use port_flush, only : portable_flush
    use exits, only : exitme
    integer(4), intent(in) :: lun
    logical,    intent(in) :: debug

    logical :: file_is_open

    if (debug) then
      inquire(unit=lun,opened=file_is_open)
      if (.not.file_is_open) call exitme('ERROR: Cannot flush non-opened file')
      call portable_flush( lun )
    endif
  end subroutine flush_unit_debug
end module io_subs
