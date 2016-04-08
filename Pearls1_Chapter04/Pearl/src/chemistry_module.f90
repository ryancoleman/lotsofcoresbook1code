module chemistry_module

  use bl_types

  implicit none

  integer, parameter :: nspecies = 9   ! number of species

  logical, save :: chemistry_initialized = .false.

  double precision, save :: Ru, Ruc, Patm
  double precision, save :: molecular_weight(nspecies), inv_mwt(nspecies)

contains

  subroutine chemistry_init()
    integer :: iwrk
    double precision :: rwrk

    call ckrp(iwrk, rwrk, Ru, Ruc, Patm)

    call ckwt(iwrk, rwrk, molecular_weight)
    inv_mwt = 1.d0 / molecular_weight

    chemistry_initialized = .true.
  end subroutine chemistry_init

  subroutine chemistry_close()
    return
  end subroutine chemistry_close

end module chemistry_module
