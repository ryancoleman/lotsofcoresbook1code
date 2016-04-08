!!
!! A _BoxArray_ is an array of boxes.
!!
module boxarray_module

  use bl_types
  use box_module
  use list_box_module

  implicit none

  type boxarray
     private
     integer :: dim = 0
     integer :: nboxes = 0
     type(box), pointer :: bxs(:) => Null()
  end type boxarray

  interface dataptr
     module procedure boxarray_dataptr
  end interface

  interface get_dim
     module procedure boxarray_dim
  end interface

  interface built_q
     module procedure boxarray_built_q
  end interface

  interface copy
     module procedure boxarray_build_copy
  end interface

  interface build
     module procedure boxarray_build_v
     module procedure boxarray_build_l
     module procedure boxarray_build_bx
  end interface

  interface destroy
     module procedure boxarray_destroy
  end interface

  interface nboxes
     module procedure boxarray_nboxes
  end interface

  interface volume
     module procedure boxarray_volume
  end interface

  interface get_box
     module procedure boxarray_get_box
  end interface

  interface boxarray_maxsize
     module procedure boxarray_maxsize_i
     module procedure boxarray_maxsize_v
  end interface

  interface boxarray_grow
     module procedure boxarray_grow_n
  end interface

  interface bbox
     module procedure boxarray_bbox
  end interface

  private :: boxarray_maxsize_l

contains

  function boxarray_dataptr(ba) result(r)
    type(boxarray), intent(in) :: ba
    type(box), pointer :: r(:)
    r => ba%bxs
  end function boxarray_dataptr
 
  pure function boxarray_built_q(ba) result(r)
    logical :: r
    type(boxarray), intent(in) :: ba
    r = ba%dim /= 0
  end function boxarray_built_q

  pure function boxarray_dim(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%dim
  end function boxarray_dim

  pure function boxarray_get_box(ba, i) result(r)
    type(boxarray), intent(in) :: ba
    integer, intent(in) :: i
    type(box) :: r
    r = ba%bxs(i)
  end function boxarray_get_box

  ! kepp
  subroutine boxarray_build_copy(ba, ba1)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(boxarray), intent(in) :: ba1
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_COPY: already built")
    if ( .not. built_q(ba1) ) return
    ba%nboxes = size(ba1%bxs)
    allocate(ba%bxs(size(ba1%bxs)))
    ba%bxs = ba1%bxs
    ba%dim = ba1%dim
    call boxarray_verify_dim(ba)
  end subroutine boxarray_build_copy

  subroutine boxarray_build_v(ba, bxs, sort)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(box), intent(in), dimension(:) :: bxs
    logical, intent(in), optional :: sort
    logical :: lsort
    
    lsort = .false. ; if (present(sort)) lsort = sort
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_V: already built")
    ba%nboxes = size(bxs)
    allocate(ba%bxs(size(bxs)))
    ba%bxs = bxs
    if ( ba%nboxes > 0 ) then
       ba%dim = ba%bxs(1)%dim
    end if
    call boxarray_verify_dim(ba)
    if (lsort) call boxarray_sort(ba) !! make sure all grids are sorted
  end subroutine boxarray_build_v

  subroutine boxarray_build_bx(ba, bx)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(box), intent(in) :: bx
    
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_BX: already built")
    ba%nboxes = 1
    allocate(ba%bxs(1))
    ba%bxs(1) = bx
    ba%dim = bx%dim
    call boxarray_verify_dim(ba)
  end subroutine boxarray_build_bx

  subroutine boxarray_build_l(ba, bl, sort)
    use bl_error_module
    type(boxarray), intent(inout) :: ba
    type(list_box), intent(in) :: bl
    logical, intent(in), optional :: sort
    type(list_box_node), pointer :: bln
    logical :: lsort
    integer :: i
    !
    ! Default is to sort.
    !
    lsort = .true. ; if ( present(sort) ) lsort = sort
    if ( built_q(ba) ) call bl_error("BOXARRAY_BUILD_L: already built")
    ba%nboxes = size(bl)
    allocate(ba%bxs(ba%nboxes))
    bln => begin(bl)
    i = 1
    do while (associated(bln))
       ba%bxs(i) = value(bln)
       i = i + 1
       bln=>next(bln)
    end do
    if ( ba%nboxes > 0 ) then
       ba%dim = ba%bxs(1)%dim
    end if
    call boxarray_verify_dim(ba)
    if ( lsort ) call boxarray_sort(ba)
  end subroutine boxarray_build_l

  subroutine boxarray_destroy(ba)
    type(boxarray), intent(inout) :: ba
    if ( associated(ba%bxs) ) then
       deallocate(ba%bxs) 
       ba%bxs => Null()
    end if
    ba%dim = 0
    ba%nboxes = 0
  end subroutine boxarray_destroy

  subroutine boxarray_sort(ba)
    use sort_box_module
    type(boxarray), intent(inout) :: ba
    call box_sort(ba%bxs)
  end subroutine boxarray_sort

  subroutine boxarray_verify_dim(ba, stat)
    use bl_error_module
    type(boxarray), intent(in) :: ba
    integer, intent(out), optional :: stat
    integer :: i, dm
    if ( present(stat) ) stat = 0
    if ( ba%nboxes < 1 ) return
    dm = ba%dim
    if ( dm == 0 ) then
       dm = ba%bxs(1)%dim
    end if
    if ( dm == 0 ) then
       call bl_error("BOXARRAY_VERIFY_DIM: dim is zero!")
    end if
    do i = 1, ba%nboxes
       if ( ba%dim /= ba%bxs(i)%dim ) then
          if ( present(stat) ) then
             stat = 1
             return
          else
             call bl_error("BOXARRAY_VERIFY_DIM: " // &
                  "ba%dim not equal to some boxes dim: ", ba%dim)
          end if
       end if
    end do
  end subroutine boxarray_verify_dim

  subroutine boxarray_grow_n(ba, n)
    type(boxarray), intent(inout) :: ba
    integer, intent(in) :: n
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = grow(ba%bxs(i), n)
    end do
  end subroutine boxarray_grow_n

  subroutine boxarray_nodalize(ba, nodal)
    type(boxarray), intent(inout) :: ba
    logical, intent(in), optional :: nodal(:)
    integer :: i
    do i = 1, ba%nboxes
       ba%bxs(i) = box_nodalize(ba%bxs(i), nodal)
    end do
  end subroutine boxarray_nodalize

  subroutine boxarray_box_boundary_n(bao, bx, n)
    type(boxarray), intent(out) :: bao
    type(box), intent(in)  :: bx
    integer, intent(in) :: n
    type(boxarray) :: baa
    call boxarray_build_bx(baa, bx)
    call boxarray_boundary_n(bao, baa, n)
    call boxarray_destroy(baa)
  end subroutine boxarray_box_boundary_n

  subroutine boxarray_boundary_n(bao, ba, n)
    type(boxarray), intent(out) :: bao
    type(boxarray), intent(in)  :: ba
    integer,        intent(in)  :: n
    call boxarray_build_copy(bao, ba)
    call boxarray_grow(bao, n)
    call boxarray_diff(bao, ba)
  end subroutine boxarray_boundary_n

  pure function boxarray_nboxes(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer :: r
    r = ba%nboxes
  end function boxarray_nboxes

  function boxarray_volume(ba) result(r)
    type(boxarray), intent(in) :: ba
    integer(kind=ll_t) :: r
    integer :: i
    r = 0_ll_t
    do i = 1, ba%nboxes
       r = r + box_volume(ba%bxs(i))
    end do
  end function boxarray_volume

  pure function boxarray_bbox(ba) result(r)
    type(boxarray), intent(in) :: ba
    type(box) :: r
    integer :: i
    r = nobox(ba%dim)
    do i = 1, ba%nboxes
       r = bbox(r, ba%bxs(i))
    end do
  end function boxarray_bbox

  subroutine boxarray_diff(bao, ba)
    type(boxarray), intent(inout) :: bao
    type(boxarray), intent(in) :: ba
    type(list_box) :: bl, bl1, bl2
    integer :: i
    call build(bl1, ba%bxs)
    do i = 1, bao%nboxes
       bl2 = boxlist_boxlist_diff(bao%bxs(i), bl1)
       call splice(bl, bl2)
    end do
    call boxarray_destroy(bao)
    call boxarray_build_l(bao, bl)
    call destroy(bl)
    call destroy(bl1)
  end subroutine boxarray_diff

  subroutine boxarray_maxsize_i(bxa, chunk) 
    type(boxarray), intent(inout) :: bxa
    integer, intent(in) :: chunk
    integer :: vchunk(bxa%dim) 
    vchunk = chunk
    call boxarray_maxsize_v(bxa, vchunk)
  end subroutine boxarray_maxsize_i

  subroutine boxarray_maxsize_v(bxa, chunk) 
    type(boxarray), intent(inout) :: bxa
    integer, intent(in), dimension(:) :: chunk
    type(list_box) :: bl
    bl = boxarray_maxsize_l(bxa, chunk)
    call boxarray_destroy(bxa)
    call boxarray_build_l(bxa, bl)
    call destroy(bl)
  end subroutine boxarray_maxsize_v

  function boxarray_maxsize_l(bxa, chunk) result(r)
    type(list_box) :: r
    type(boxarray), intent(in) ::  bxa
    integer, intent(in), dimension(:) :: chunk
    integer :: i,k
    type(list_box_node), pointer :: li
    integer :: len(bxa%dim)
    integer :: nl, bs, rt, nblk, sz, ex, ks, ps
    type(box) :: bxr, bxl

    do i = 1, bxa%nboxes
       call push_back(r, bxa%bxs(i))
    end do
    li => begin(r)
    do while ( associated(li) )
       len = extent(value(li))
       do i = 1, bxa%dim
          if ( len(i) > chunk(i) ) then
             rt = 1
             bs = chunk(i)
             nl = len(i)
             do while ( mod(bs,2) == 0 .AND. mod(nl,2) == 0)
                rt = rt * 2
                bs = bs/2
                nl = nl/2
             end do
             nblk = nl/bs
             if ( mod(nl,bs) /= 0 ) nblk = nblk + 1
             sz   = nl/nblk
             ex   = mod(nl,nblk)
             do k = 0, nblk-2
                if ( k < ex ) then
                   ks = (sz+1)*rt
                else
                   ks = sz*rt
                end if
                ps = upb(value(li), i) - ks + 1
                call box_chop(value(li), bxr, bxl, i, ps)
                call set(li, bxr)
                call push_back(r, bxl)
             end do
          end if
       end do
       li => next(li)
    end do

  end function boxarray_maxsize_l

end module boxarray_module
