module port_exit
  implicit none
  private
  public  :: portable_exit
contains
  subroutine portable_exit( exit_status )
    integer(4),intent(in) :: exit_status
    call exit (exit_status)
  end subroutine portable_exit
end module port_exit
