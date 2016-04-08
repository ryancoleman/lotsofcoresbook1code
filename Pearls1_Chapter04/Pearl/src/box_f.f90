!! Support for box calculus
!!
!! _Boxs_ are define rectangular domains in space.  They can also be
!! thought of as providing a definition of a Cell-Centered index region.
module box_module

  use bl_types

  implicit none

  integer, parameter :: MAX_SPACEDIM = 3  

  !! Boxes consist of
  !!   - A dimensionality, and
  !!   - A low and high end of cells
  !!
  !! Initially, each box is completely empty, with a low end at positive
  !! infinity and a low end at negative infinity.  Almost all box operations
  !! are ill-defined on default boxes.
  type box
     integer :: dim  = 0
     integer :: lo(MAX_SPACEDIM) =  Huge(1)
     integer :: hi(MAX_SPACEDIM) = -Huge(1)
  end type box

  !! Returns a default box
  interface nobox
     module procedure box_nobox
  end interface

  !! Returns the dimensionality of the box.
  interface get_dim
     module procedure box_dim
  end interface

  !! Returns the volume of the box as an integer; aborts on overflow.
  interface volume
     module procedure box_volume
  end interface

  !! Returns the upper bound of a box, or of a box in
  !! a given direction
  interface upb
     module procedure box_upb
     module procedure box_upb_d
  end interface
  !! Sets the upper bound of a box, or sets the upper
  !! bound of a box in a given direction
  interface set_upb
     module procedure box_set_upb_d
  end interface set_upb

  !! Returns the lower bound of a box, or of a box in
  !! a given direction
  interface lwb
     module procedure box_lwb
     module procedure box_lwb_d
  end interface
  !! Sets the lower bound of a box, or sets the lower
  !! bound of a box in a given direction
  interface set_lwb
     module procedure box_set_lwb_d
  end interface set_lwb

  !! returns a coarsened version of a box
  interface coarsen
     module procedure box_coarsen_i
  end interface

  !! grows a box
  interface grow
     module procedure box_grow_n
     module procedure box_grow_n_d_f
  end interface

  !! shifts a box
  interface shift
     module procedure box_shift_d
     module procedure box_shift_v
  end interface shift

  !! Returns true if the lwb(bx1) < lwb(bx2)
  interface operator(.LT.)
     module procedure box_less
  end interface

  !! Returns the intersection of two boxes
  interface intersection
     module procedure box_intersection
  end interface
  interface operator(.INTERSECT.)
     module procedure box_intersection
  end interface

  !! Returns true if two boxes intersect
  interface intersects
     module procedure box_intersects
  end interface

  !! returns the bounding box of two boxes; that is
  !! the smallest box the contains them both.
  interface bbox
     module procedure box_bbox
  end interface

  !! Returns the size, or extent, of a box or of a box 
  !! in a given direction.
  interface extent
     module procedure box_extent
     module procedure box_extent_d
  end interface

  !! Prints a box.
  interface print
     module procedure box_print
  end interface

  interface contains
     module procedure box_contains
     module procedure box_contains_iv
  end interface

  interface empty
     module procedure box_empty
  end interface

contains

  elemental function int_coarsen(v, i) result(r)
    integer, intent(in) :: v
    integer, intent(in) :: i
    integer :: r
    if ( v < 0 ) then
       r = -abs(v+1)/i - 1
    else
       r = v/i
    end if
  end function int_coarsen

  !! Makes a box given a low and high end.
  !! If given a single point (index), it returs a unit box located
  !! at the that point.
  function make_box(lo, hi) result(r)
    integer, intent(in), dimension(:) :: lo
    integer, intent(in), dimension(:), optional :: hi
    type(box) :: r
    if ( present(hi) ) then
       call box_build_2(r, lo, hi)
    else
       call box_build_2(r, lo, lo)
    end if
  end function make_box

  !! Returns a completely empty box [+Inf,-Inf]; this box is also
  !! the default box of dimension _dim_.
  pure function box_nobox(dim) result(r)
    integer, intent(in) :: dim
    type(box) :: r
    r%dim = dim
    r%lo(1:dim) =  Huge(1)
    r%hi(1:dim) = -Huge(1)
  end function box_nobox

  !! Returns .TRUE. if the box _bx_ contains no cells, or if it is of dimension
  !! zero.
  pure function box_empty(bx) result(r)
    logical :: r
    type(box), intent(in) :: bx
    if ( bx%dim == 3 ) then
        r = bx%lo(1) > bx%hi(1) .or. bx%lo(2) > bx%hi(2) .or. bx%lo(3) > bx%hi(3)
    else if ( bx%dim == 2 ) then
        r = bx%lo(1) > bx%hi(1) .or. bx%lo(2) > bx%hi(2)
    else if ( bx%dim == 1 ) then
        r = bx%lo(1) > bx%hi(1)
    else
        r = .true.
    end if
  end function box_empty

  !! Builds a box given a low and high end of index space.  It is an error
  !! for the sizes of _lo_ and _hi_ to be unequal or for the size of _lo_
  !! to exceed MAX_SPACEDIM.
  subroutine box_build_2(bx, lo, hi)
    use bl_error_module
    type(box), intent(out) :: bx
    integer, intent(in) :: lo(:), hi(:)
    if ( size(lo) /= size(hi)    ) call bl_error("BOX_BUILD_2: lo /= hi")
    if ( size(lo) > MAX_SPACEDIM ) call bl_error("BOX_BUILD_2: size(lo,hi) > MAX_SPACEDIM")
    bx%dim = size(lo)
    bx%lo(1:bx%dim) = lo
    bx%hi(1:bx%dim) = hi
  end subroutine box_build_2

  !! Returns the dimensionality of the box _bx_.
  pure function box_dim(bx) result(r)
    type(box), intent(in) :: bx
    integer :: r
    r = bx%dim
  end function box_dim

  !! Returns a coarsened version of the box _bx_, coarsened by the integer
  !! _ci_.
  pure function box_coarsen_i(bx, ci) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ci
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = int_coarsen(bx%lo(1:bx%dim), ci)
    r%hi(1:r%dim) = int_coarsen(bx%hi(1:bx%dim), ci)
  end function box_coarsen_i

  !! Box_nodalize, perhaps inaptly named, grows a box
  !! on the hiside by one if the correspending directions of
  !! nodal are true, or else not at all if nodal is not present.
  !! this last because by default fabs, multifabs, and etc, are cell-
  !! centered.
  function box_nodalize(bx, nodal) result(r)
    type(box) :: r
    type(box), intent(in) :: bx
    logical, intent(in), optional :: nodal(:)
    integer :: i
    r = bx
    if ( .not. present(nodal) ) return
    do i = 1, bx%dim
       if ( nodal(i) ) r = grow(r, 1, i, +1)
    end do
  end function box_nodalize

  function box_grow_n_d_f(bx, ri, dim, face) result(r)
    use bl_error_module
    type(box), intent(in) :: bx
    integer, intent(in) :: ri, dim
    integer, intent(in) :: face
    type(box) :: r
    r = bx
    if ( face == -1 ) then
       r%lo(dim) = bx%lo(dim) - ri
    else if ( face == 1 ) then
       r%hi(dim) = bx%hi(dim) + ri
    else 
       call bl_error("BOX_GROW_N_D_F: unexpected face: ", face)
    end if
  end function box_grow_n_d_f

  pure function box_grow_n(bx, ri) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: ri
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = bx%lo(1:bx%dim) - ri
    r%hi(1:r%dim) = bx%hi(1:bx%dim) + ri
  end function box_grow_n

  !! lexigraphic ordering in the lower-bound
  pure function box_less(bx1, bx2) result(r)
    logical :: r
    type(box), intent(in) :: bx1, bx2
    integer :: i
    do i = 1, bx1%dim
       if ( bx1%lo(i) == bx2%lo(i) ) cycle
       if ( bx1%lo(i) < bx2%lo(i) ) r = .TRUE.
       if ( bx1%lo(i) > bx2%lo(i) ) r = .FALSE.
       return
    end do
    r = .FALSE.
  end function box_less

  function box_intersection(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    type(box) :: r
    r%dim = min(bx1%dim, bx2%dim)
    if ( r%dim == 3 ) then
       r%lo(1) = max(bx1%lo(1), bx2%lo(1))
       r%hi(1) = min(bx1%hi(1), bx2%hi(1))
       r%lo(2) = max(bx1%lo(2), bx2%lo(2))
       r%hi(2) = min(bx1%hi(2), bx2%hi(2))
       r%lo(3) = max(bx1%lo(3), bx2%lo(3))
       r%hi(3) = min(bx1%hi(3), bx2%hi(3))
    else if ( r%dim == 2 ) then
       r%lo(1) = max(bx1%lo(1), bx2%lo(1))
       r%hi(1) = min(bx1%hi(1), bx2%hi(1))
       r%lo(2) = max(bx1%lo(2), bx2%lo(2))
       r%hi(2) = min(bx1%hi(2), bx2%hi(2))
    else if ( r%dim == 1 ) then
       r%lo(1) = max(bx1%lo(1), bx2%lo(1))
       r%hi(1) = min(bx1%hi(1), bx2%hi(1))
    end if
  end function box_intersection

  pure function box_intersects(bx1, bx2) result(r)
    type(box), intent(in) :: bx1, bx2
    logical :: r
    integer :: dm
    dm = min(bx1%dim, bx2%dim)
    r = all( (min(bx1%hi(1:dm),bx2%hi(1:dm)) &
         - max(bx1%lo(1:dm),bx2%lo(1:dm))) >= 0 )
  end function box_intersects

  pure function box_contains(bx1, bx2, strict) result(r)
    logical :: r
    type(box), intent(in) :: bx1, bx2
    logical, intent(in), optional :: strict
    logical :: lstrict
    lstrict = .false.; if ( present(strict) ) lstrict = strict
    if ( lstrict ) then
       r = all(bx1%lo(1:bx1%dim) < bx2%lo(1:bx2%dim)) .and. &
            all(bx1%lo(1:bx1%dim) < bx2%hi(1:bx2%dim)) .and. &
            all(bx1%hi(1:bx1%dim) > bx2%lo(1:bx2%dim)) .and. &
            all(bx1%hi(1:bx1%dim) > bx2%hi(1:bx2%dim))
    else
       r = all(bx1%lo(1:bx1%dim) <= bx2%lo(1:bx2%dim)) .and. &
            all(bx1%lo(1:bx1%dim) <= bx2%hi(1:bx2%dim)) .and. &
            all(bx1%hi(1:bx1%dim) >= bx2%lo(1:bx2%dim)) .and. &
            all(bx1%hi(1:bx1%dim) >= bx2%hi(1:bx2%dim))
    end if
  end function box_contains

  pure function box_contains_iv(bx1, iv, strict) result(r)
    logical :: r
    type(box), intent(in) :: bx1
    integer, intent(in) :: iv(:)
    logical, intent(in), optional :: strict
    logical :: lstrict
    lstrict = .false.; if ( present(strict) ) lstrict = strict
    if ( lstrict ) then
       r = all(bx1%lo(1:bx1%dim) < iv) .and. all(bx1%hi(1:bx1%dim) > iv)
    else
       r = all(bx1%lo(1:bx1%dim) <= iv) .and. all(bx1%hi(1:bx1%dim) >= iv)
    end if
  end function box_contains_iv

  pure function box_shift_v(bx, iv) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: iv(:)
    type(box) :: r
    r%dim = bx%dim
    r%lo(1:r%dim) = bx%lo(1:r%dim) + iv(1:r%dim)
    r%hi(1:r%dim) = bx%hi(1:r%dim) + iv(1:r%dim)
  end function box_shift_v
  pure function box_shift_d(bx, i, dim) result(r)
    type(box), intent(in) :: bx
    integer, intent(in) :: i
    integer, intent(in) :: dim
    type(box) :: r
    r = bx
    r%lo(dim) = bx%lo(dim) + i
    r%hi(dim) = bx%hi(dim) + i
  end function box_shift_d

  pure function box_extent(bx) result(r)
    type (box), intent(in) :: bx
    integer, dimension(bx%dim) :: r
    r = bx%hi(1:bx%dim) - bx%lo(1:bx%dim) + 1
  end function box_extent

  pure function box_extent_d(bx, dim) result(r)
    type (box), intent(in) :: bx
    integer, intent(in) :: dim
    integer :: r
    r = bx%hi(dim)-bx%lo(dim)+1
  end function box_extent_d

  pure function box_lwb(bx) result(r)
    type (box), intent(in) :: bx
    integer, dimension(bx%dim) :: r
    r = Huge(1)
    r(1:bx%dim) = bx%lo(1:bx%dim)
  end function box_lwb
  pure function box_lwb_d(bx, dim) result(r)
    type (box), intent(in) :: bx
    integer, intent(in) :: dim
    integer :: r
    r = bx%lo(dim)
  end function box_lwb_d
  subroutine box_set_lwb_d(bx, dim, v)
    integer, intent(in) :: v
    integer, intent(in) :: dim
    type (box), intent(inout) :: bx
    bx%lo(dim) = v
  end subroutine box_set_lwb_d

  pure function box_upb(bx) result(r)
    type (box), intent(in) :: bx
    integer, dimension(bx%dim) :: r
    r = -Huge(1)
    r(1:bx%dim) = bx%hi(1:bx%dim)
  end function box_upb
  pure function box_upb_d(bx, dim) result(r)
    type (box), intent(in) :: bx
    integer, intent(in) :: dim
    integer :: r
    r = bx%hi(dim)
  end function box_upb_d
  subroutine box_set_upb_d(bx, dim, v)
    integer, intent(in) :: v
    integer, intent(in) :: dim
    type (box), intent(inout) :: bx
    bx%hi(dim) = v
  end subroutine box_set_upb_d

  !! Box_chop chops a box into two boxes in the coordinate direction
  !! dir at the cell ichop.  The resulting boxes are  bxl and bxr where
  !! bxl is to the left of the chop line and bxr is to the right.  If
  !! the chop line is outside the extent of the box; the coresponding 
  !! is an empty, default box.
  subroutine box_chop(bx, bxl, bxr, dim, ichop)
    type(box), intent(in) :: bx
    type(box), intent(out) :: bxl, bxr
    integer, intent(in) :: dim, ichop
    if ( ichop <= lwb(bx,dim) ) then
       bxr = bx
    else if ( ichop > upb(bx,dim) ) then
       bxl = bx
    else
       bxl = bx
       bxr = bx
       bxl%hi(dim) = ichop-1
       bxr%lo(dim) = ichop
    end if
  end subroutine box_chop

  function box_volume(bx) result(r)
    use bl_error_module
    type(box), intent(in) :: bx
    integer :: r, i, l
    r = 1
    do i = 1, bx%dim
       l = (bx%hi(i)-bx%lo(i)+1)
       if ( r .le. Huge(r) / l ) then
          r = r*l
       else
          call print(bx, 'BAD BOX')
          call bl_error('box_volume(): overflow')
       end if
    end do
  end function box_volume

  pure function box_bbox(b1, b2) result(r)
    type(box), intent(in) :: b1, b2
    type(box) :: r
    r%dim = b1%dim
    r%lo(1:r%dim) = min(b1%lo(1:b1%dim),b2%lo(1:b2%dim))
    r%hi(1:r%dim) = max(b1%hi(1:b1%dim),b2%hi(1:b2%dim))
  end function box_bbox

  subroutine box_print(bx, str, unit, advance, legacy, nodal, skip)
    use bl_IO_module
    use bl_string_module
    type(box), intent(in) :: bx
    character(len=*), intent(in), optional :: str
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: skip
    logical, intent(in), optional :: legacy
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: nodal(:)
    integer :: un
    character(len=3) :: adv
    integer :: i
    logical :: llegacy
    un = unit_stdout(unit)
    adv = unit_advance(advance)
    llegacy = .FALSE.; if ( present(legacy) ) llegacy = legacy
    if ( llegacy ) then
       call write_a_legacy_box(nodal)
    else
       call write_a_box()
    end if
    if ( eq_i(adv, 'YES') ) write(unit=un, fmt='()')
  contains
    subroutine write_a_box
      call unit_skip(un, skip)
      write(unit=un, fmt='("BOX[")', advance='no')
      if ( present(str) ) then
         write(unit=un, fmt='("(* ", A, " *)")', advance = 'no') str
      end if
      if ( bx%dim > 0 ) then
         write(unit=un, fmt='("{", 3(I0,:,", "))', advance = 'no') bx%lo(1:bx%dim)
         write(unit=un, fmt='("}, {", 3(I0,:,", "))', advance = 'no') bx%hi(1:bx%dim)
         write(unit=un, fmt='("}]")', advance = 'no' )
      else
         write(unit=un, fmt='("]")', advance = 'no' )
      end if
    end subroutine write_a_box
    subroutine write_a_legacy_box(nodal)
      logical, intent(in), optional :: nodal(:)
      integer :: tp(bx%dim), thi(bx%dim)
      tp = 0
      if ( present(nodal) ) then
         do i = 1, bx%dim
            if ( nodal(i) ) tp(i) = 1
         end do
      end if
      if ( bx%dim == 0 ) then
         write(unit=un, fmt='("()")', advance = 'no')
      else
         thi = bx%hi(1:bx%dim)
         if ( present(nodal) ) then
            where ( nodal ) thi = thi + 1
         end if
         write(unit=un, fmt='("(")', advance='no')

         call write_a_legacy_intvect(bx%lo); write(unit=un, fmt='(" ")', advance='no')
         call write_a_legacy_intvect(thi); write(unit=un, fmt='(" ")', advance='no')
         call write_a_legacy_intvect(tp)
         write(unit=un, fmt='(")")', advance='no')
      end if
    end subroutine write_a_legacy_box

    subroutine write_a_legacy_intvect(vc)
      integer, intent(in) :: vc(:)
      write(unit=un, fmt='("(")', advance='no')
      do i = 1, bx%dim-1
         write(unit=un, fmt='(I0, ",")', advance = 'no' ) vc(i)
      end do
      write(unit=un, fmt='(I0, ")")', advance = 'no' ) vc(bx%dim)
    end subroutine write_a_legacy_intvect
  end subroutine box_print

  !
  ! Used by layout.
  !
  subroutine box_periodic_shift(dmn,b,nodal,pmask,ng,shft,cnt,bxs,sync_shift)

    type(box), intent(in)            :: dmn,b
    logical,   intent(in)            :: nodal(:),pmask(:)
    integer,   intent(in)            :: ng
    integer,   intent(out)           :: shft(:,:),cnt
    type(box), intent(out), optional :: bxs(:)
    logical,   intent(in), optional  :: sync_shift

    type(box) :: dom, bx, src
    integer   :: nbeg(3),nend(3),ldom(3),r(3),ri,rj,rk,l(3),i
    logical   :: lsync_shift

    lsync_shift = .false.; if ( present(sync_shift) ) lsync_shift = sync_shift

    if ( .not. lsync_shift ) cnt = 0

    if ( all(pmask .eqv. .false.) ) return

    dom = box_nodalize(dmn,nodal)

    if ( lsync_shift ) then
       bx = b
    else
       bx = grow(box_nodalize(b,nodal),ng)
    end if

    if ( contains(dom,bx,strict = any(nodal)) ) return

    l = 0

    do i = 1, dom%dim
       if (nodal(i)) l(i) = 1
    end do

    nbeg           = 0
    nend           = 0
    nbeg(1:bx%dim) = -1
    nend(1:bx%dim) = +1
    ldom           = (/ extent(dom,1)-l(1), extent(dom,2)-l(2), extent(dom,3)-l(3) /)

    do ri = nbeg(1), nend(1)
       if (ri /= 0 .and. (.not. is_periodic(1))) cycle
       if (ri /= 0 .and. is_periodic(1)) bx = shift(bx,ri*ldom(1),1)

       do rj = nbeg(2), nend(2)
          if (rj /= 0 .and. (.not. is_periodic(2))) cycle
          if (rj /= 0 .and. is_periodic(2)) bx = shift(bx,rj*ldom(2),2)

          do rk = nbeg(3), nend(3)
             if (rk /= 0 .and. (.not. is_periodic(3))) cycle
             if (rk /= 0 .and. is_periodic(3)) bx = shift(bx,rk*ldom(3),3)

             if (ri == 0 .and. rj == 0 .and. rk == 0) cycle

             src = intersection(dom,bx)

             if ( .not. empty(src) ) then
                cnt = cnt + 1
                if ( present(bxs) ) bxs(cnt) = src
                r = (/ri*ldom(1),rj*ldom(2),rk*ldom(3)/)
                shft(cnt,1:bx%dim) = r(1:bx%dim)
             end if

             if (rk /= 0 .and. is_periodic(3)) bx = shift(bx,-rk*ldom(3),3)
          end do

          if (rj /= 0 .and. is_periodic(2)) bx = shift(bx,-rj*ldom(2),2)
       end do

       if (ri /= 0 .and. is_periodic(1)) bx = shift(bx,-ri*ldom(1),1)
    end do

  contains

    function is_periodic(i) result(r)
      integer, intent(in) :: i
      logical             :: r
      r = .false.
      if (i >= 1 .and. i <= bx%dim) r = pmask(i) .eqv. .true.
    end function is_periodic

  end subroutine box_periodic_shift

end module box_module
