module sort_box_module
  ! Adapted from Meissner, Adams, et.al., and NR.

  use box_module

  implicit none

  interface box_sort
     module procedure heapsort_box
  end interface box_sort

contains

  subroutine heapsort_box( array, cmp)
    type(box), dimension(:), intent(in out) :: array
    interface
       function cmp(a,b) result(r)
         use box_module
         implicit none
         logical :: r
         type(box), intent(in) :: a, b
       end function cmp
    end interface
    optional cmp

    if ( present(cmp) ) then
       call he_sort_cmp(array,cmp)
    else
       call he_sort_no_cmp(array)
    end if

  contains

    subroutine box_swap(a,b)
      type(box), intent(inout) :: a, b
      type(box) :: t
      t = a; a = b; b = t
    end subroutine box_swap

    subroutine he_sort_cmp(array,cmp)
      type(box), dimension(:), intent(in out) :: array
      interface
         function cmp(a,b) result(r)
           use box_module
           implicit none
           logical :: r
           type(box), intent(in) :: a, b
         end function cmp
      end interface
      integer :: n, i

      n = size (array)
      do i = n / 2, 1, -1
         call peck_cmp(array,i,n,cmp)
      end do
      do i = n, 2, -1
         call box_swap(array(1),array(i))
         call peck_cmp(array,1, i-1,cmp)
      end do

    end subroutine he_sort_cmp

    subroutine he_sort_no_cmp(array)
      type(box), dimension(:), intent(in out) :: array
      integer :: n, i

      n = size (array)
      do i = n / 2, 1, -1
         call peck_no_cmp(array,i,n)
      end do
      do i = n, 2, -1
         call box_swap(array(1),array(i))
         call peck_no_cmp(array,1, i-1)
      end do

    end subroutine he_sort_no_cmp

    subroutine peck_cmp(array,l,r,cmp)
      type(box), dimension(:), intent(inout) :: array
      integer, intent(in) :: r
      integer, intent(in) :: l
      interface
         function cmp(a,b) result(r)
           use box_module
           implicit none
           logical :: r
           type(box), intent(in) :: a, b
         end function cmp
      end interface
      integer :: i
      type(box) :: next

      next = array(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( cmp(array(i), array(i+1)) ) i = i+1
         end if
         if ( .not. cmp( next, array(i)) ) exit
         array(i/2) = array(i)
         i = 2*i
      end do
      array(i/2) = next

    end subroutine peck_cmp

    subroutine peck_no_cmp(array,l,r)
      type(box), dimension(:), intent(inout) :: array
      integer, intent(in) :: r
      integer, intent(in) :: l
      integer :: i

      type(box) :: next
      next = array(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( box_less(array(i), array(i+1)) ) i = i+1
         end if
         if ( .not. box_less( next, array(i)) ) exit
         array(i/2) = array(i)
         i = 2*i
      end do
      array(i/2) = next

    end subroutine peck_no_cmp

  end subroutine heapsort_box

end module sort_box_module

