module list_box_module

  use box_module  

  implicit none

  type list_box_node
     private
     type(box) :: v ! = 
     type(list_box_node), pointer :: prev => NULL()
     type(list_box_node), pointer :: next => NULL()
  end type list_box_node

  type list_box
     private
     integer :: size = 0
     type(list_box_node), pointer:: head => NULL()
     type(list_box_node), pointer:: tail => NULL()
  end type list_box

  interface value
     module procedure list_node_value_box
  end interface value

  interface set
     module procedure list_node_set_box
  end interface set

  interface next
     module procedure list_node_next_box
  end interface next

  interface build
     module procedure list_build_v_box
  end interface build

  interface destroy
     module procedure list_destroy_box
  end interface destroy

  interface size
     module procedure list_size_box
  end interface size

  interface empty
     module procedure list_empty_box
  end interface empty

  interface begin
     module procedure list_begin_box
  end interface begin

  interface push_back
     module procedure list_push_back_box
  end interface

  interface clear
     module procedure list_clear_box
  end interface clear

  interface erase
     module procedure list_erase_box
     module procedure list_erase_range_box
  end interface erase

  interface splice
     module procedure list_splice_box
  end interface splice

contains

  ! Node accessors.

  pure function list_node_value_box(n) result(r)
    type(list_box_node), intent(in) :: n
    type(box) :: r

    r = n%v

  end function list_node_value_box

  subroutine list_node_set_box(n, v)
    type(list_box_node), intent(inout) :: n
    type(box), intent(in) ::  v

    n%v = v

  end subroutine list_node_set_box

  function list_node_next_box(n) result(r)
    type(list_box_node), pointer :: n, r

    r => n%next

  end function list_node_next_box

  subroutine list_build_v_box(r, d)
    type(list_box), intent(out) :: r
    type(box), dimension(:), intent(in) :: d
    integer :: i

    do i = 1, size(d)
       call push_back(r, d(i))
    end do
    r%size = size(d)

  end subroutine list_build_v_box

  subroutine list_destroy_box(r)
    type(list_box), intent(inout) :: r

    call clear(r)
    r%size = 0

  end subroutine list_destroy_box

  pure function list_size_box(l) result (r)
    type(list_box), intent(in) :: l
    integer :: r

    r = l%size

  end function list_size_box

  function list_begin_box(l) result(r)
    type(list_box), intent(in) :: l
    type(list_box_node), pointer :: r

    r => l%head

  end function list_begin_box

  pure function list_empty_box(l) result (r)
    type(list_box), intent(in) :: l
    logical:: r

    r = .not.(associated(l%head) .AND. associated(l%tail))

  end function list_empty_box

  subroutine list_push_back_box(l,v)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(box), intent(in) :: v
    type(list_box_node), pointer :: n
    integer :: ierr
    allocate(n, stat = ierr)
    if ( ierr /= 0 ) then
       call bl_error("list_push_back_box: failed to allocate memory")
    end if
    n%v = v
    n%prev => l%tail
    l%tail => n
    if ( .not. associated (l%head ) ) then
       l%head => l%tail
    else
       l%tail%prev%next => n
    end if
    l%size = l%size + 1

  end subroutine list_push_back_box

  function list_erase_box(l, beg) result (r)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: beg, r
    integer :: ierr

    if ( associated(beg, l%tail) ) then
       r => Null()
    else
       r => beg%next
    end if
    if ( associated(beg, l%head) .AND. associated(beg, l%tail) ) then
       l%head => Null()
       l%tail => Null()
    else if ( associated(beg, l%head) ) then
       l%head => beg%next
       l%head%prev => Null()
    else if ( associated(beg, l%tail) ) then
       l%tail => beg%prev
       l%tail%next => Null()
    else
       beg%next%prev => beg%prev
       beg%prev%next => beg%next
    end if
    deallocate(beg, stat = ierr)
    if ( ierr /= 0 ) then
       call bl_error("list_erase_box: failed to deallocate")
    end if
    l%size = l%size - 1

  end function list_erase_box

  function list_erase_range_box(l, beg, end) result(rr)
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: beg, end
    type(list_box_node), pointer :: r, e, rr

    r => beg
    e => end
    if ( .not. associated(e) ) then
       do while ( associated(r) )
          r => erase(l,r)
       end do
    else
       do while ( .not. associated(r, e) )
          r => erase(l,r)
       end do
    end if
    rr => r

  end function list_erase_range_box

  subroutine list_clear_box(l)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: r

    r => erase(l, l%head, Null(r))
    if ( .not. empty(l) ) call bl_error('list_clear_box(): clear is broken')
    if ( l%size /= 0 ) call bl_error('list_clear_box(): clear is broken')

  end subroutine list_clear_box

  subroutine list_splice_box(l1, l2)
    type(list_box), intent(inout) :: l1, l2

    if ( empty(l2) ) return
    if ( empty(l1) ) then
       l1%head => l2%head
       l1%tail => l2%tail
    else
       l1%tail%next => l2%head
       l2%head%prev => l1%tail
       l1%tail => l2%tail
    end if
    nullify(l2%head, l2%tail)
    l1%size = l2%size + l1%size
    l2%size = 0

  end subroutine list_splice_box

  function boxlist_box_diff(bx1, b2) result(r)
    type(box), intent(in) :: bx1, b2
    type(list_box) :: r
    type(box) :: b1, bn
    integer, dimension(bx1%dim) :: b2lo, b2hi, b1lo, b1hi
    integer :: i, dm

    dm = bx1%dim
    b1 = bx1
    if ( .not. contains(b2,b1) ) then
       if ( .not. intersects(b1,b2) ) then
          call push_back(r, b1)
       else
          b2lo = lwb(b2); b2hi = upb(b2)
          do i = 1, dm
             b1lo = lwb(b1); b1hi = upb(b1)
             if ( b1lo(i) < b2lo(i) .AND. b2lo(i) <= b1hi(i) ) then
                bn = b1
                call set_lwb(bn, i, b1lo(i))
                call set_upb(bn, i, b2lo(i)-1)
                call push_back(r, bn)
                call set_lwb(b1, i, b2lo(i))
             end if
             if ( b1lo(i) <= b2hi(i) .AND. b2hi(i) < b1hi(i) ) then
                bn = b1
                call set_lwb(bn, i, b2hi(i)+1)
                call set_upb(bn, i, b1hi(i))
                call push_back(r, bn)
                call set_upb(b1, i, b2hi(i))
             end if
          end do
       end if
    end if

  end function boxlist_box_diff

  ! r = bx - bxl
  function boxlist_boxlist_diff(bx, bxl) result(r)

    type(box), intent(in) :: bx
    type(list_box), intent(in) :: bxl
    type(list_box) :: r, bl
    type(list_box_node), pointer :: blp, bp

    call push_back(r, bx)
    blp => begin(bxl)
    do while ( associated(blp) .AND. .NOT. empty(r) )
       bp => begin(r)
       do while ( associated(bp) )
          if ( intersects(value(bp), value(blp)) ) then
             bl = boxlist_box_diff(value(bp), value(blp))
             call splice(r, bl)
             bp => erase(r, bp)
          else
             bp => next(bp)
          end if
       end do
       blp => next(blp)
    end do

  end function boxlist_boxlist_diff

end module list_box_module
