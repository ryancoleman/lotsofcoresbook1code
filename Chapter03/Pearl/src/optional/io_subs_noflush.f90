module port_flush
  implicit none
  private
  public :: portable_flush
contains
  subroutine portable_flush ( lun )
    ! A "Unused dummy variable LUN"-warning is accepted here:
    integer(4), intent(in) :: lun
    return
  end subroutine portable_flush
end module port_flush

