program main

  use BoxLib

  implicit none

  call boxlib_initialize()

  call smc()

  call boxlib_finalize()

end program main
