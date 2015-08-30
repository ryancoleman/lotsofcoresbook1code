
module sfc_module
  use box_module
  !
  ! This is more than a bit of a hack.  We're going to store a few values here
  ! while inside the sfc_i() routine, so that they're accessible via
  ! sfc_greater_i().
  !
  integer            :: dm, mpower
  type(box), pointer :: pbxs(:)
end module sfc_module

module knapsack_module

  use bl_types

  implicit none

  logical, private :: knapsack_verbose = .false.

  real(kind=dp_t), private :: knapsack_threshold = 0.9_dp_t

  private :: sfc_greater_i

contains

  subroutine knapsack_set_verbose(yesorno)
    logical, intent(in) :: yesorno
    knapsack_verbose = yesorno
  end subroutine knapsack_set_verbose

  subroutine sfc_i(prc, ibxs, bxs, np, verbose)

    use sfc_module
    use parallel
    use box_module
    use sort_i_module
    use bl_error_module

    integer,   intent(out)           :: prc(:)
    integer,   intent(in )           :: ibxs(:)
    type(box), intent(in ), target   :: bxs(:)
    integer,   intent(in )           :: np
    logical,   intent(in ), optional :: verbose

    logical         :: lverb
    integer         ::  i
    real(kind=dp_t) :: t1, t2, efficiency

    integer, allocatable :: iorder(:), whichcpu(:)

    if ( np < 1 ) call bl_error('sfc_i(): np < 1')

    if ( size(ibxs) < 1 ) call bl_error('sfc_i(): size(ibxs) < 1')

    call cpu_time(t1)

    lverb = knapsack_verbose ; if ( present(verbose) ) lverb = verbose
    !
    ! Set dm, mpower & pbxs in sfc_module so they're accessible to sfc_greater_i().
    !
    dm     =  get_dim(bxs(1))
    mpower =  maxpower()
    pbxs   => bxs

    allocate(iorder(size(ibxs)), whichcpu(size(ibxs)))
    !
    ! Set to "bad" value that we can check for later to ensure array filled correctly.
    !
    whichcpu = -1

    do i = 1, size(ibxs)
       iorder(i) = i
    end do

    call sort(iorder, sfc_greater_i)
    !
    ! "iorder" now indexes the boxes in morton space-filling-curve order.
    !
    call distribute()
    !
    ! "whichcpu(i)" now contains the CPU on which to place the iorder(i)'th box.
    !
    if ( minval(whichcpu) < 0 ) call bl_error('sfc_i(): improper CPU number')

    do i = 1, size(ibxs)
       prc(iorder(i)) = whichcpu(i)
    end do

    call cpu_time(t2)

    if ( lverb .and. np > 1 .and. parallel_ioprocessor() ) then
       print *, 'SFC effi = ', efficiency
       print *, 'SFC time = ', t2-t1
    end if

    contains

      function maxpower() result(r)

        integer :: r, maxijk, i, d

        maxijk = lwb(bxs(1),1)

        do i = 1,size(bxs)
           do d = 1, dm
              maxijk = max(maxijk,lwb(bxs(i),d))
           end do
        end do

        r = 0
        do while ( ishft(1,r) <= maxijk )
           r = r + 1
        end do

      end function maxpower

      subroutine distribute()

        integer         :: k, cnt, sz
        real(kind=dp_t) :: totalvol, volpercpu, vol, maxvol

        maxvol = -Huge(1_dp_t)

        volpercpu = 0_dp_t
        do i = 1, size(ibxs)
           volpercpu = volpercpu + ibxs(i)
        end do
        volpercpu = volpercpu / np

        k        = 1
        sz       = size(ibxs)
        totalvol = 0_dp_t

        do i = 1, np

           cnt = 0
           vol = 0_dp_t

           do while ( (k <= sz) .and. ((i == np) .or. (vol < volpercpu)) )
              whichcpu(k) = i
              vol = vol + ibxs(iorder(k))
              k   = k   + 1
              cnt = cnt + 1
           end do

           totalvol = totalvol + vol

           if ( (totalvol/i) > volpercpu .and. (cnt > 1) .and. (k <= sz) ) then
              k        = k - 1
              vol      = vol - ibxs(iorder(k))
              totalvol = totalvol - ibxs(iorder(k))
           endif

           maxvol = max(maxvol,vol)

        end do
        !
        ! Force "whichcpu" values to be zero-based instead of one-based.
        !
        whichcpu = whichcpu - 1

        efficiency = volpercpu / maxvol

      end subroutine distribute

  end subroutine sfc_i

  function sfc_greater_i(ilhs,irhs) result(r)

    use sfc_module

    logical             :: r
    integer, intent(in) :: ilhs,irhs

    integer m,d,NNN,llwb,rlwb

    do m = mpower,0,-1

       NNN = ishft(1,m)

       do d = 1, dm

          llwb = lwb(pbxs(ilhs),d)/NNN
          rlwb = lwb(pbxs(irhs),d)/NNN

          if ( llwb < rlwb ) then
             r = .true.
             return
          else if ( llwb > rlwb ) then
             r = .false.
             return
          end if
       end do
    end do

    r = .false.

  end function sfc_greater_i

end module knapsack_module
