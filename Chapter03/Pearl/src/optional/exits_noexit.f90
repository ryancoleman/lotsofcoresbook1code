module port_exit
  implicit none
  private
  public  :: portable_exit
contains
  subroutine portable_exit( exit_status )
    integer(4),intent(in) :: exit_status
    write (*,*) 'I am going to kill myself - my exit status is ', exit_status
    stop 'Missing exit so using stop instead'
  end subroutine portable_exit
end module port_exit
