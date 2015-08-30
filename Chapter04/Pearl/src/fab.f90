module fab_module

  use bl_types
  use box_module

  implicit none

  !! Each Fab type, integer, double, and logical, consists of
  !! DIM:    Dimension
  !  BX:     The Box in index space for which this FAB is defined
  !! IBX:    The index range of the valid data for the FAB
  !! PBX:    The physical box for the FAB
  !! NC:     Number of components
  !! NG:     Number of ghost cells
  
  !! When a FAB is created IBX = BX, unless it is nodal, in which case
  !! IBX = grow(BX, FACE=hi, 1 (in the nodal directions).
  !! PBX = grow(IBX, NG)
  
  !! For parallel systems the BX, IBX, PBX, etc are all defined, but the
  !! underlying pointer will not be allocated.

  !! All FABS are 'Four' dimensional, conventially, (NX,NY,NZ,NC) in size.
  !! NY = 1, NZ = 1, when DIM =1, NZ = 1, when DIM = 2.

  type fab
     private
     integer   :: dim = 0
     type(box) :: bx
     type(box) :: pbx
     type(box) :: ibx
     integer   :: nc = 1
     real(kind = dp_t), pointer, dimension(:,:,:,:) :: p => Null()
  end type fab

  interface build
     module procedure fab_build
  end interface

  interface destroy
     module procedure fab_destroy
  end interface

  interface dataptr
     module procedure fab_dataptr
     module procedure fab_dataptr_c
     module procedure fab_dataptr_bx_c
  end interface

  interface setval
     module procedure fab_setval
     module procedure fab_setval_c
     module procedure fab_setval_bx_c
  end interface

  interface get_box
     module procedure fab_get_box
  end interface

  interface get_pbox
     module procedure fab_get_pbox
  end interface

  interface get_ibox
     module procedure fab_get_ibox
  end interface

  interface volume
     module procedure fab_volume
  end interface

  logical,       private, save ::      do_init_fabs = .false.
  real(dp_t),    private, save ::  fab_default_init = -Huge(1.0_dp_t)

contains

  function fab_volume(fb, all) result(r)
    integer(kind=ll_t) :: r
    type(fab), intent(in) :: fb
    logical, intent(in), optional :: all
    if ( all ) then
       r = volume(get_pbox(fb))
    else
       r = volume(get_box(fb))
    end if
    r = r * fb%nc
  end function fab_volume

  pure function fab_get_box(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%bx
  end function fab_get_box

  pure function fab_get_pbox(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%pbx
  end function fab_get_pbox

  pure function fab_get_ibox(fb) result(r)
    type(fab), intent(in) :: fb
    type(box) :: r
    r = fb%ibx
  end function fab_get_ibox
  
  subroutine fab_build(fb, bx, nc, ng, nodal, alloc, stencil)
    type(fab), intent(out) :: fb
    type(box), intent(in)  :: bx
    integer, intent(in), optional :: ng, nc
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: alloc
    logical, intent(in), optional :: stencil
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: lnc, lng
    logical :: lal, lst
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    lal = .true. ; if ( present(alloc)   ) lal = alloc
    lst = .false.; if ( present(stencil) ) lst = stencil
    lo = 1
    hi = 1
    fb%dim = bx%dim
    fb%bx = bx
    fb%nc = lnc
    fb%ibx = box_nodalize(bx, nodal)
    fb%pbx = grow(fb%ibx, lng)
    lo(1:fb%dim) = fb%pbx%lo(1:fb%dim)
    hi(1:fb%dim) = fb%pbx%hi(1:fb%dim)
    if ( lal ) then
       if ( lst) then
          allocate(fb%p(1:lnc,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
       else
          allocate(fb%p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1:lnc))
       end if
       if ( do_init_fabs ) call setval(fb, fab_default_init)
    end if
  end subroutine fab_build

  subroutine fab_destroy(fb)
    type(fab), intent(inout) :: fb
    if ( associated(fb%p) ) then
       deallocate(fb%p)
    end if
    fb%bx  = nobox(fb%dim)
    fb%dim = 0
    fb%nc  = 0
  end subroutine fab_destroy

  function fab_dataptr(fb) result(r)
    type(fab), intent(in) :: fb
    real(kind=dp_t), pointer :: r(:,:,:,:)
    r => fb%p
  end function fab_dataptr

  function fab_dataptr_c(fb, c, nc) result(r)
    use bl_error_module
    type(fab), intent(in) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('fab_dataptr_c: not enough components')
    r => fb%p(:,:,:,c:c+lnc-1)
  end function fab_dataptr_c

  function fab_dataptr_bx_c(fb, bx, c, nc) result(r)
    use bl_error_module
    type(fab), intent(in) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), pointer :: r(:,:,:,:)
    integer :: lnc
    lnc = 1; if ( present(nc) ) lnc = nc
    if ( (c+lnc-1) > fb%nc ) call bl_error('fab_dataptr_bx_c: not enough components')
    if ( .not. contains(fb%pbx, bx) ) call bl_error('fab_dataptr_bx_c: bx is too large')
    select case (fb%dim)
    case (1)
       r => fb%p(bx%lo(1):bx%hi(1),:,:,c:c+lnc-1)
    case (2)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),:,c:c+lnc-1)
    case (3)
       r => fb%p(bx%lo(1):bx%hi(1),bx%lo(2):bx%hi(2),bx%lo(3):bx%hi(3),c:c+lnc-1)
    end select
  end function fab_dataptr_bx_c

  subroutine fab_setval(fb, val)
    use bl_error_module
    type(fab), intent(inout) :: fb
    real(kind=dp_t), intent(in) :: val
    if ( .not. associated(fb%p) ) call bl_error("FAB_SETVAL: not associated")
    fb%p = val
  end subroutine fab_setval

  subroutine fab_setval_c(fb, val, c, nc)
    use bl_error_module
    type(fab), intent(inout) :: fb
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_c(fb, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_c

  subroutine fab_setval_bx_c(fb, val, bx, c, nc)
    use bl_error_module
    type(fab), intent(inout) :: fb
    type(box), intent(in) :: bx
    integer, intent(in) :: c
    integer, intent(in), optional :: nc
    real(kind=dp_t), intent(in) :: val
    real(kind=dp_t), pointer :: p(:,:,:,:)
    p => fab_dataptr_bx_c(fb, bx, c, nc)
    if ( .not. associated(p) ) call bl_error("FAB_SETVAL: not associated")
    p = val
  end subroutine fab_setval_bx_c

end module fab_module
