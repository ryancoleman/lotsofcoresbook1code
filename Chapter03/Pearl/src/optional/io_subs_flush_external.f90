module port_flush
  implicit none
  private
  public :: portable_flush
contains
  subroutine portable_flush ( lun )
    integer(4), intent(in) :: lun
    ! wrap the "usual" external flush
    call flush(lun)
  end subroutine portable_flush
end module port_flush
