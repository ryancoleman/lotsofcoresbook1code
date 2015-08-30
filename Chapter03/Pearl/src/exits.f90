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

module exits

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- interfaces ----------------------------------------------------------------
  interface exitme
    module procedure exitme_txt
    module procedure exitme_val
    module procedure exitme_val_txt
    module procedure exitme_stdout
  end interface !f90 std exitme

  !- public vars & methods -----------------------------------------------------
  public  :: exitme

  !- private vars & methods ----------------------------------------------------
  private :: exitme_val, exitme_txt, exitme_val_txt, exitme_stdout, exit_hack
  private :: get_fmt
  integer(4),   private :: oval
  character(8), private :: cfmt

contains

  ! ----------------------------------------------------------------------------

  subroutine exitme_txt(txt)
    use dmi_mpi_global, only : iu06

    character (len=*), intent(in) :: txt

    write(iu06,'(a)') txt
    write(iu06,'(a14)') 'Program exited'
    call exit_hack(1_4)
  end subroutine exitme_txt

  ! ----------------------------------------------------------------------------

  subroutine exitme_val(val)
    use dmi_mpi_global, only : iu06

    integer(4), intent(in) :: val

    call get_fmt( val )

    write(iu06,cfmt) 'Program exited, return value=', oval
    call exit_hack(val)
  end subroutine exitme_val

  ! ----------------------------------------------------------------------------

  subroutine exitme_val_txt(val,txt)
    use dmi_mpi_global, only : iu06

    integer(4), intent(in)        :: val
    character (len=*), intent(in) :: txt

    call get_fmt( val )

    write(iu06,'(a)') txt
    write(iu06,cfmt) 'Program exited, return value=', oval
    call exit_hack(val)
  end subroutine exitme_val_txt

  ! ----------------------------------------------------------------------------

  subroutine exitme_stdout(val,txt,flag)
    use dmi_mpi_global, only : iu06

    integer(4), intent(in)        :: val
    character (len=*), intent(in) :: txt
    logical, intent(in)           :: flag

    call get_fmt( val )

    if (flag) then
      write(*,'(a)') txt
      write(*,cfmt) 'Program exited, return value=', oval
    else
      write(iu06,'(a)') txt
      write(iu06,cfmt) 'Program exited, return value=', oval
    endif
    call exit_hack(val,flag)
  end subroutine exitme_stdout

  ! ----------------------------------------------------------------------------

  subroutine exit_hack( returncode, flag )
    use port_exit, only : portable_exit
    use dmi_mpi_global, only : dmpi_close, dmpi_abort

    ! This is a hack to have only one single ansi f90/f95 warning about exit()
    integer(4), intent(in), optional :: returncode
    logical,    intent(in), optional :: flag
    integer(4) :: exit_status

    if (present(returncode)) then
      exit_status = returncode
    else
      exit_status = 0
    endif

    if (present(flag)) then
      if (exit_status == 0) then
        call dmpi_close( flag )
      else
        call dmpi_abort( flag )
      endif
    else
      if (exit_status == 0) then
        call dmpi_close( .false. )
      else
        call dmpi_abort( .false. )
      endif
    endif

    call portable_exit( exit_status ) 

  end subroutine exit_hack

  ! ----------------------------------------------------------------------------

  subroutine get_fmt( val )
    use dmi_mpi_global, only : iu06

    implicit none
    integer(4), intent(in) :: val

    oval = val
    if (val < 10) then
      cfmt = '(a29,'//'i1'//')'
    elseif (val < 100) then
      cfmt = '(a29,'//'i2'//')'
    elseif (val < 1000) then
      cfmt = '(a29,'//'i3'//')'
    elseif (val < 10000) then
      cfmt = '(a29,'//'i4'//')'
    else 
      cfmt = '(a29,'//'i4'//')'
      write(iu06,'(a39)') 'exit-value too large, set to 9999'
      oval = 9999
    endif

  end subroutine get_fmt

end module exits

