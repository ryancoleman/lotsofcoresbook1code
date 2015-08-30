module port_flush
  implicit none
  private
  public :: portable_flush
contains
  subroutine portable_flush ( lun )
    integer(4), intent(in) :: lun
    ! wrap the F2003 version of 'FLUSH'
    flush (lun)
  end subroutine portable_flush
end module port_flush
