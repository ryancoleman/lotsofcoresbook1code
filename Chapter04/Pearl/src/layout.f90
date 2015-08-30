module layout_module

  use parallel
  use boxarray_module

  implicit none

  integer, private, parameter :: LA_UNDF = 0
  integer, private, parameter :: LA_BASE = 1
  integer, private, parameter :: LA_CRSN = 2
  integer, private, parameter :: LA_PCHD = 3
  integer, private, parameter :: LA_PPRT = 4

  integer, parameter :: LA_KNAPSACK   = 101
  integer, parameter :: LA_ROUNDROBIN = 102
  integer, parameter :: LA_LOCAL      = 103
  integer, parameter :: LA_EXPLICIT   = 104

  integer, private :: verbose = 0

  integer, private :: def_mapping = LA_KNAPSACK

  type comm_dsc
     integer   :: nd = 0                 ! dst box number
     integer   :: ns = 0                 ! src box number
     type(box) :: dbx                    ! dst sub-box
     type(box) :: sbx                    ! src sub-box
     integer   :: pv = 0                 ! number of points in buf prior to this
     integer   :: sh(MAX_SPACEDIM+1) = 0 ! shape for data from rcvbuf
     integer   :: pr                     ! Processors number of src or dest
  end type comm_dsc

  type trns_dsc
     integer :: sz = 0              ! Size of chunk 
     integer :: pv = -1             ! Number points in buf prior to this
     integer :: pr = MPI_ANY_SOURCE ! src or destination processor
  end type trns_dsc

  type remote_conn
     integer                 :: svol = 0          ! numpts in snd volume
     integer                 :: rvol = 0          ! numpts in rcv volume
     integer                 :: nsnd = 0          ! Number of snd chunks
     integer                 :: nrcv = 0          ! Number of rcv chunks
     type(comm_dsc), pointer :: snd(:) => Null()
     type(comm_dsc), pointer :: rcv(:) => Null()
     integer                 :: nrp  = 0          ! Number of processes receiving from
     integer                 :: nsp  = 0          ! Number of processes sending to
     type(trns_dsc), pointer :: str(:) => Null()
     type(trns_dsc), pointer :: rtr(:) => Null()
  end type remote_conn

  type local_copy_desc
     integer   :: ns = 0    ! Source box in layout
     integer   :: nd = 0    ! Destination box in layout
     type(box) :: sbx       ! Sub-box for this copy
     type(box) :: dbx       ! Sub-box for this copy
  end type local_copy_desc

  type local_conn
     integer                        :: ncpy   ! Number of cpy chunks
     type(local_copy_desc), pointer :: cpy(:) => Null()
  end type local_conn

  type boxassoc
     integer                 :: dim      = 0       ! spatial dimension 1, 2, or 3
     integer                 :: nboxes   = 0       ! number of boxes
     integer                 :: grwth    = 0       ! growth factor
     integer                 :: idim     = 0       ! if 1, 2, or 3, grow only in that direction
     logical, pointer        :: nodal(:) => Null() ! nodal flag
     logical                 :: cross    = .false. ! cross/full stencil?
     type(local_conn)        :: l_con
     type(remote_conn)       :: r_con
     type(boxassoc), pointer :: next     => Null()
  end type boxassoc

  type copyassoc
     integer                   :: dim        = 0       ! spatial dimension 1, 2, or 3
     integer                   :: reused     = 0
     integer                   :: hash
     logical, pointer          :: nd_dst(:)  => Null() ! dst nodal flag
     logical, pointer          :: nd_src(:)  => Null() ! src nodal flag
     type(local_conn)          :: l_con
     type(remote_conn)         :: r_con
     type(copyassoc),  pointer :: next       => Null()
     type(boxarray)            :: ba_src
     type(boxarray)            :: ba_dst
     integer, pointer          :: prc_src(:) => Null()
     integer, pointer          :: prc_dst(:) => Null()
  end type copyassoc
  !
  ! Used by layout_get_box_intersector().
  !
  type box_intersector
     integer   :: i
     type(box) :: bx
  end type box_intersector

  type box_hash_bin
     integer, pointer :: iv(:) => Null()
  end type box_hash_bin
  !
  ! Global list of copyassoc's used by multifab copy routines.
  !
  type(copyassoc), pointer, save, private :: the_copyassoc_head => Null()

  integer, save, private :: the_copyassoc_cnt = 0  ! Count of copyassocs on list.
  integer, save, private :: the_copyassoc_max = 25 ! Maximum # copyassocs allowed on list.

  type layout
     integer                   :: la_type =  LA_UNDF
     type(layout_rep), pointer :: lap     => Null()
  end type layout

  !! Defines the box distribution and box connectivity of a boxarray.
  type layout_rep
     integer                         :: dim    = 0            ! Spatial dimension: 1, 2, or 3.
     integer                         :: id     = 0
     integer                         :: nboxes = 0            ! Total count of boxes.
     integer                         :: nlocal = 0            ! Number of boxes we own.
     type(box)                       :: pd                    ! Problem Domain.
     logical, pointer                :: pmask(:)    => Null() ! Periodic mask.
     integer, pointer                :: prc(:)      => Null() ! Which processor owns each box.
     integer, pointer                :: idx(:)      => Null() ! Global indices of boxes we own.
     type(boxarray)                  :: bxa
     type(boxassoc), pointer         :: bxasc       => Null()
     type(coarsened_layout), pointer :: crse_la     => Null()
     type(pn_layout), pointer        :: pn_children => Null()
     ! Box Hashing
     integer                         :: crsn                = -1
     integer                         :: plo(MAX_SPACEDIM)   = 0
     integer                         :: phi(MAX_SPACEDIM)   = 0
     integer                         :: vshft(MAX_SPACEDIM) = 0
     type(box_hash_bin), pointer     :: bins(:,:,:)         => Null()
  end type layout_rep

  !! A layout that is derived by coarsening an existing layout,
  !! The processor distribution and the number of boxes will be
  !! the same as for the parent layout.  The intent is to be used
  !! in multigrid solvers that keep coarsened grids on the same
  !! processor as their parent in the hierarchy.
  type coarsened_layout
     integer                         :: dim = 0
     integer, pointer                :: crse(:) => Null()
     type(layout)                    :: la
     type(coarsened_layout), pointer :: next => Null()
  end type coarsened_layout

  type pn_layout
     integer                  :: dim = 0
     integer, pointer         :: refr(:) => Null()
     type(layout)             :: la
     type(pn_layout), pointer :: next => Null()
  end type pn_layout

  integer, private :: g_layout_next_id = 0

  interface built_q
     module procedure boxassoc_built_q
  end interface

  interface build
     module procedure layout_build_ba
  end interface

  interface destroy
     module procedure layout_destroy
  end interface

  interface local
     module procedure layout_local
  end interface

  interface remote
     module procedure layout_remote
  end interface

  interface nboxes
     module procedure layout_nboxes
  end interface

  interface nlocal
     module procedure layout_nlocal
  end interface

  interface local_index
     module procedure layout_local_index
  end interface

  interface global_index
     module procedure layout_global_index
  end interface

  interface get_box
     module procedure layout_get_box
  end interface

  interface get_boxarray
     module procedure layout_boxarray
  end interface

  interface get_proc
     module procedure layout_get_proc
  end interface

  interface get_dim
     module procedure layout_dim
  end interface

  private :: greater_i, layout_next_id, layout_rep_build, layout_rep_destroy

contains

  function layout_next_id() result(r)
    integer :: r
    g_layout_next_id = g_layout_next_id + 1
    r = g_layout_next_id
  end function layout_next_id

  pure function layout_dim(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%dim
  end function layout_dim

  pure function layout_nboxes(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%nboxes
  end function layout_nboxes

  pure function layout_nlocal(la) result(r)
    integer :: r
    type(layout), intent(in) :: la
    r = la%lap%nlocal
  end function layout_nlocal

  !
  ! Given a global index "i" into the boxarray on which
  ! the layout is built, returns the corresponding local
  ! index "r" at which that index is stored in la&lap%idx(:).
  ! It's an error if the layout does not "own" index "r"
  ! or if the index "i" is not in the range of boxes for the
  ! layout. nboxes is the number of boxes in the associated layout.
  !
  function layout_local_index(la,i) result(r)

    use bl_error_module

    integer,      intent(in) :: i
    type(layout), intent(in) :: la
    integer                  :: r

    call bl_assert(i >= 1 .and. i <= la%lap%nboxes, "layout_local_index: invalid global index")

    if (parallel_nprocs() == 1) then
       call bl_assert(la%lap%idx(i) == i, "layout_local_index: how did this happen?")
       r = i
    else
       r = bsearch(la%lap%idx,i)
       call bl_assert(r >= 1, "layout_local_index: no corresponding local index")
    endif

  contains
      !
      ! Returns -1 if "val" is not found in array "arr".
      !
      ! "arr" is assumed to be sorted from smallest to largest.
      !
      pure function bsearch (arr,val) result (r)
        integer, intent(in) :: arr(:), val
        integer             :: r, lo, hi, mid
        r  = -1
        lo = lbound(arr,1)
        hi = ubound(arr,1)
        do while (lo <= hi)
           mid = (lo + hi) / 2
           if (arr(mid) == val) then
              r = mid
              exit
           else if (arr(mid) > val) then
              hi = mid - 1
           else
              lo = mid + 1
           end if
        end do
      end function bsearch

  end function layout_local_index

  pure function layout_global_index(la,i) result(r)
    integer,      intent(in) :: i
    type(layout), intent(in) :: la
    integer                  :: r
    r = la%lap%idx(i)
  end function layout_global_index
  
  function layout_boxarray(la) result(r)
    type(layout), intent(in) :: la
    type(boxarray) :: r
    r = la%lap%bxa
  end function layout_boxarray

  subroutine layout_rep_build(lap, ba, pd, pmask, mapping, explicit_mapping)
    use bl_error_module
    type(layout_rep), intent(out) :: lap
    type(boxarray), intent(in) :: ba
    type(box), intent(in) :: pd
    logical, intent(in) :: pmask(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)
    integer :: lmapping, i, j

    lmapping = def_mapping; if ( present(mapping) ) lmapping = mapping
    if ( present(explicit_mapping) ) then
       if ( present(mapping) ) then
          if ( mapping /= LA_EXPLICIT ) then
             call bl_error("layout_rep_build():explicit_mapping doesn't match mapping")
          end if
       end if
       lmapping = LA_EXPLICIT
    end if

    call boxarray_build_copy(lap%bxa, ba)

    lap%dim    = get_dim(lap%bxa)
    lap%nboxes = nboxes(lap%bxa)
    lap%id     = layout_next_id()
    lap%pd     = pd
    allocate(lap%pmask(lap%dim))
    lap%pmask = pmask

    allocate(lap%prc(lap%nboxes))

    call layout_knapsack(lap%prc, ba)

    lap%nlocal = 0
    do i = 1, lap%nboxes
       if (parallel_myproc() == lap%prc(i)) lap%nlocal = lap%nlocal + 1
    end do

    allocate(lap%idx(lap%nlocal))

    j = 1
    do i = 1, lap%nboxes
       if (parallel_myproc() == lap%prc(i)) then
          lap%idx(j) = i
          j = j + 1
       end if
    end do
  end subroutine layout_rep_build

  recursive subroutine layout_rep_destroy(lap, la_type)
    type(layout_rep), pointer :: lap
    integer, intent(in) :: la_type
    type(coarsened_layout), pointer :: clp, oclp
    type(pn_layout), pointer :: pnp, opnp
    type(boxassoc),  pointer :: bxa, obxa
    integer :: i, j, k
    if ( la_type /= LA_CRSN ) then
       deallocate(lap%prc)
       deallocate(lap%idx)
    end if
    call destroy(lap%bxa)
    if ( (la_type == LA_BASE) .or. (la_type == LA_PCHD) ) then
       deallocate(lap%pmask)
    end if
    clp => lap%crse_la
    do while ( associated(clp) )
       oclp => clp%next
       deallocate(clp%crse)
       call layout_rep_destroy(clp%la%lap, LA_CRSN)
       deallocate(clp)
       clp => oclp
    end do
    pnp => lap%pn_children
    do while ( associated(pnp) )
       opnp => pnp%next
       deallocate(pnp%refr)
       call layout_rep_destroy(pnp%la%lap, LA_PCHD)
       deallocate(pnp)
       pnp  => opnp
    end do
    !
    ! Get rid of boxassocs.
    !
    bxa => lap%bxasc
    do while ( associated(bxa) )
       obxa => bxa%next
       call boxassoc_destroy(bxa)
       deallocate(bxa)
       bxa => obxa
    end do
    !
    ! Remove any boxarray hash.
    !
    if ( associated(lap%bins) ) then
       do k = lbound(lap%bins,3), ubound(lap%bins,3)
          do j = lbound(lap%bins,2), ubound(lap%bins,2)
             do i = lbound(lap%bins,1), ubound(lap%bins,1)
                deallocate(lap%bins(i,j,k)%iv)
             end do
          end do
       end do
       deallocate(lap%bins)
    end if

    deallocate(lap)
  end subroutine layout_rep_destroy

  subroutine layout_build_ba(la, ba, pd, pmask, mapping, explicit_mapping)
    type(layout)  , intent(  out) :: la
    type(boxarray), intent(in   ) :: ba
    type(box)     , intent(in   ) :: pd
    logical, intent(in), optional :: pmask(:)
    integer, intent(in), optional :: mapping
    integer, intent(in), optional :: explicit_mapping(:)

    logical :: lpmask(get_dim(ba))
    lpmask = .false.; if ( present(pmask) ) lpmask = pmask
    allocate(la%lap)
    la%la_type = LA_BASE
    call layout_rep_build(la%lap, ba, pd, lpmask, mapping, explicit_mapping)
  end subroutine layout_build_ba

  subroutine layout_destroy(la)
    use bl_error_module
    type(layout), intent(inout) :: la
    if ( la%la_type /= LA_BASE ) call bl_error("layout_destroy(): confused")
    call layout_rep_destroy(la%lap, LA_BASE)
  end subroutine layout_destroy

  pure function layout_remote(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    logical :: r
    r = la%lap%prc(i) /= parallel_myproc()
  end function layout_remote

  pure function layout_local(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    logical :: r
    r = la%lap%prc(i) == parallel_myproc()
  end function layout_local

  pure function layout_get_box(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(la%lap%bxa, i)
  end function layout_get_box

  pure function layout_get_proc(la, i) result(r)
    type(layout), intent(in) :: la
    integer, intent(in) :: i
    integer :: r
    r = la%lap%prc(i)
  end function layout_get_proc

  function greater_i(a,b) result(r)
    logical :: r
    integer, intent(in) :: a, b
    r = a > b
  end function greater_i

  subroutine layout_knapsack(prc, bxs)
    use fab_module
    use bl_error_module
    use knapsack_module
    integer,        intent(out), dimension(:) :: prc
    type(boxarray), intent(in )               :: bxs

    integer              :: i
    integer, allocatable :: ibxs(:)
    type(box), pointer   :: pbxs(:)

    if ( nboxes(bxs) /= size(prc) ) call bl_error('layout_knapsack: how did this happen?')

    allocate(ibxs(nboxes(bxs)))

    do i = 1, size(ibxs,1)
       ibxs(i) = volume(get_box(bxs,i))
    end do

    !
    ! Always use Morton space-filling-curve distribution in this stripped down version
    !
    pbxs => dataptr(bxs)
    call sfc_i(prc, ibxs, pbxs, parallel_nprocs())
    deallocate(ibxs)
  end subroutine layout_knapsack

  pure function boxassoc_check(bxa, ng, nodal, cross, idim) result(r)
    type(boxassoc), intent(in) :: bxa
    integer,        intent(in) :: ng
    logical,        intent(in) :: nodal(:)
    logical,        intent(in) :: cross
    integer, intent(in), optional :: idim
    logical                    :: r
    r = (bxa%grwth == ng) .and. all(bxa%nodal .eqv. nodal) .and. (bxa%cross .eqv. cross)
    if (present(idim)) then
       r = r .and. (bxa%idim == idim)
    end if
  end function boxassoc_check

  function layout_boxassoc(la, ng, nodal, cross, idim) result(r)
    type(boxassoc)               :: r
    type(layout) , intent(inout) :: la
    integer, intent(in)          :: ng
    logical, intent(in)          :: nodal(:)
    logical, intent(in)          :: cross
    integer, intent(in),optional :: idim

    type(boxassoc), pointer :: bp

    bp => la%lap%bxasc
    do while ( associated(bp) )
       if ( boxassoc_check(bp, ng, nodal, cross, idim) ) then
          r = bp
          return
       end if
       bp => bp%next
    end do
    !
    ! Have to build one.
    !
    allocate (bp)
    call boxassoc_build(bp, la%lap, ng, nodal, cross, idim=idim)
    bp%next      => la%lap%bxasc
    la%lap%bxasc => bp
    r = bp
  end function layout_boxassoc

  pure function boxassoc_built_q(bxasc) result(r)
    logical :: r
    type(boxassoc), intent(in) :: bxasc
    r = bxasc%dim /= 0
  end function boxassoc_built_q

  subroutine boxarray_bndry_periodic(bxai, dmn, b, nodal, pmask, ng, shfts, cross, idim)
    type(boxarray), intent(out) :: bxai
    type(box),      intent(in)  :: dmn, b
    logical,        intent(in)  :: nodal(:), pmask(:)
    integer,        intent(in)  :: ng
    integer,        intent(out) :: shfts(:,:)
    logical,        intent(in)  :: cross
    integer,        intent(in)  :: idim

    integer               :: i, cnt
    type(box)             :: bxs(3**get_dim(b)), gbx
    type(box),allocatable :: bv(:)
    integer               :: shft(3**get_dim(b),get_dim(b))
    integer               :: upbx(1:get_dim(b)), lwbx(1:get_dim(b))
    type(boxarray)        :: tba
    type(list_box)        :: bl

    if ( cross ) then
       gbx = box_nodalize(b,nodal)
       do i = 1, gbx%dim
          if (idim .ne. 0 .and. idim .ne. i) cycle
          !
          ! lo face
          !
          upbx    = upb(gbx)
          lwbx    = lwb(gbx)
          upbx(i) = lwbx(i) - 1
          lwbx(i) = lwbx(i) - ng
          call push_back(bl, make_box(lwbx,upbx))
          !
          ! hi face
          !
          upbx    = upb(gbx)
          lwbx    = lwb(gbx)
          lwbx(i) = upbx(i) + 1
          upbx(i) = upbx(i) + ng
          call push_back(bl, make_box(lwbx,upbx))
       end do
       call build(tba, bl, sort = .false.)
       call destroy(bl)
    else
       call boxarray_box_boundary_n(tba, box_nodalize(b,nodal), ng)       
    end if

    shfts = 0

    call box_periodic_shift(dmn, b, nodal, pmask, ng, shft, cnt, bxs)

    if ( cnt > 0 ) then
       allocate(bv(nboxes(tba)+cnt))
       do i = 1, nboxes(tba)
          bv(i) = get_box(tba,i)
       end do
       bv(nboxes(tba)+1:nboxes(tba)+cnt) = bxs(1:cnt)
       shfts(nboxes(tba)+1:nboxes(tba)+cnt,:) = shft(1:cnt,:)
       call destroy(tba)
       call boxarray_build_v(tba, bv, sort = .false.)
    end if

    bxai = tba

  end subroutine boxarray_bndry_periodic

  subroutine sumassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, lcnt_r, cnt_s, cnt_r, pvol, parr)
    type(boxassoc),   intent(inout) :: bxasc
    type(layout),     intent(inout) :: la
    type(layout_rep), intent(in)    :: lap
    type(layout),     intent(inout) :: latmp
    integer,          intent(in)    :: ng
    logical,          intent(in)    :: anynodal
    integer,          intent(inout) :: lcnt_r, cnt_s, cnt_r
    integer,          intent(inout) :: pvol(0:,:), parr(0:,:)

    integer                         :: i, j, ii, jj, i_r, i_s, li_r
    type(boxarray)                  :: bxa, bxai
    type(box)                       :: sbx, dbx
    type(local_copy_desc), pointer  :: n_cpy(:) => Null()
    type(comm_dsc), pointer         :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(box_intersector), pointer  :: bi(:)
    integer                         :: shft(2*3**la%lap%dim,la%lap%dim), sh(MAX_SPACEDIM+1)
    integer, parameter              :: ChunkSize = 50

    bxa = get_boxarray(la)

    li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, nboxes(bxa)
       call boxarray_bndry_periodic(bxai, lap%pd, get_box(bxa,i), bxasc%nodal, lap%pmask, ng, shft, .false., bxasc%idim)
       do ii = 1, nboxes(bxai)
          if ( anynodal ) then
             bi => layout_get_box_intersector(latmp, get_box(bxai,ii))
          else
             bi => layout_get_box_intersector(la,    get_box(bxai,ii))
          end if
          do jj = 1, size(bi)
             j   = bi(jj)%i
             if ( remote(la,i) .and. remote(la,j) ) cycle
             dbx = bi(jj)%bx
             sbx = shift(dbx,-shft(ii,:))
             if ( local(la,i) .and. local(la, j) ) then
                if ( li_r > size(bxasc%l_con%cpy) ) then
                   allocate(n_cpy(size(bxasc%l_con%cpy) + ChunkSize))
                   n_cpy(1:li_r-1) = bxasc%l_con%cpy(1:li_r-1)
                   deallocate(bxasc%l_con%cpy)
                   bxasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                bxasc%l_con%cpy(li_r)%nd  = j
                bxasc%l_con%cpy(li_r)%ns  = i
                bxasc%l_con%cpy(li_r)%sbx = sbx
                bxasc%l_con%cpy(li_r)%dbx = dbx
                li_r                      = li_r + 1
             else if ( local(la, j) ) then
                if ( i_r > size(bxasc%r_con%rcv) ) then
                   allocate(n_rcv(size(bxasc%r_con%rcv) + ChunkSize))
                   n_rcv(1:i_r-1) = bxasc%r_con%rcv(1:i_r-1)
                   deallocate(bxasc%r_con%rcv)
                   bxasc%r_con%rcv => n_rcv
                end if
                cnt_r                    = cnt_r + 1
                parr(lap%prc(i), 1)      = parr(lap%prc(i), 1) + 1
                pvol(lap%prc(i), 1)      = pvol(lap%prc(i), 1) + volume(dbx)
                bxasc%r_con%rcv(i_r)%nd  = j
                bxasc%r_con%rcv(i_r)%ns  = i
                bxasc%r_con%rcv(i_r)%sbx = sbx
                bxasc%r_con%rcv(i_r)%dbx = dbx
                bxasc%r_con%rcv(i_r)%pr  = get_proc(la, i)
                sh                       = 1
                sh(1:bxasc%dim)          = extent(dbx)
                bxasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             else if ( local(la, i) ) then
                if ( i_s > size(bxasc%r_con%snd) ) then
                   allocate(n_snd(size(bxasc%r_con%snd) + ChunkSize))
                   n_snd(1:i_s-1) = bxasc%r_con%snd(1:i_s-1)
                   deallocate(bxasc%r_con%snd)
                   bxasc%r_con%snd => n_snd
                end if
                cnt_s                    = cnt_s + 1
                parr(lap%prc(j), 2)      = parr(lap%prc(j), 2) + 1
                pvol(lap%prc(j), 2)      = pvol(lap%prc(j), 2) + volume(dbx)
                bxasc%r_con%snd(i_s)%nd  = j
                bxasc%r_con%snd(i_s)%ns  = i
                bxasc%r_con%snd(i_s)%sbx = sbx
                bxasc%r_con%snd(i_s)%dbx = dbx
                bxasc%r_con%snd(i_s)%pr  = get_proc(la, j)
                i_s                      = i_s + 1
             end if
          end do
          deallocate(bi)
       end do
       call boxarray_destroy(bxai)
    end do
  end subroutine sumassoc_build_innards

  subroutine boxassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, cross, lcnt_r, cnt_s, cnt_r, pvol, parr)
    type(boxassoc),   intent(inout) :: bxasc
    type(layout),     intent(inout) :: la
    type(layout_rep), intent(in)    :: lap
    type(layout),     intent(inout) :: latmp
    integer,          intent(in)    :: ng
    logical,          intent(in)    :: anynodal
    logical,          intent(in)    :: cross
    integer,          intent(inout) :: lcnt_r, cnt_s, cnt_r
    integer,          intent(inout) :: pvol(0:,:), parr(0:,:)

    integer                         :: i, j, ii, jj, i_r, i_s, li_r
    type(boxarray)                  :: bxa, bxai
    type(box)                       :: sbx, dbx
    type(local_copy_desc), pointer  :: n_cpy(:) => Null()
    type(comm_dsc), pointer         :: n_snd(:) => Null(), n_rcv(:) => Null()
    type(box_intersector), pointer  :: bi(:)
    integer                         :: shft(2*3**la%lap%dim,la%lap%dim), sh(MAX_SPACEDIM+1)
    integer, parameter              :: ChunkSize = 50

    bxa = get_boxarray(la)

    li_r = 1; i_r = 1; i_s = 1
    !
    ! Consider all copies I <- J.
    !
    do i = 1, nboxes(bxa)
       call boxarray_bndry_periodic(bxai, lap%pd, get_box(bxa,i), bxasc%nodal, lap%pmask, ng, shft, cross, bxasc%idim)
       do ii = 1, nboxes(bxai)
          if ( anynodal ) then
             bi => layout_get_box_intersector(latmp, get_box(bxai,ii))
          else
             bi => layout_get_box_intersector(la,    get_box(bxai,ii))
          end if
          do jj = 1, size(bi)
             j   = bi(jj)%i
             if ( remote(la,i) .and. remote(la,j) ) cycle
             sbx = bi(jj)%bx
             dbx = shift(sbx,-shft(ii,:))
             if ( local(la,i) .and. local(la, j) ) then
                if ( li_r > size(bxasc%l_con%cpy) ) then
                   allocate(n_cpy(size(bxasc%l_con%cpy) + ChunkSize))
                   n_cpy(1:li_r-1) = bxasc%l_con%cpy(1:li_r-1)
                   deallocate(bxasc%l_con%cpy)
                   bxasc%l_con%cpy => n_cpy
                end if
                lcnt_r                    = lcnt_r + 1
                bxasc%l_con%cpy(li_r)%nd  = i
                bxasc%l_con%cpy(li_r)%ns  = j
                bxasc%l_con%cpy(li_r)%sbx = sbx
                bxasc%l_con%cpy(li_r)%dbx = dbx
                li_r                      = li_r + 1
             else if ( local(la, j) ) then
                if ( i_s > size(bxasc%r_con%snd) ) then
                   allocate(n_snd(size(bxasc%r_con%snd) + ChunkSize))
                   n_snd(1:i_s-1) = bxasc%r_con%snd(1:i_s-1)
                   deallocate(bxasc%r_con%snd)
                   bxasc%r_con%snd => n_snd
                end if
                cnt_s                    = cnt_s + 1
                parr(lap%prc(i), 2)      = parr(lap%prc(i), 2) + 1
                pvol(lap%prc(i), 2)      = pvol(lap%prc(i), 2) + volume(sbx)
                bxasc%r_con%snd(i_s)%nd  = i
                bxasc%r_con%snd(i_s)%ns  = j
                bxasc%r_con%snd(i_s)%sbx = sbx
                bxasc%r_con%snd(i_s)%dbx = dbx
                bxasc%r_con%snd(i_s)%pr  = get_proc(la, i)
                i_s                      = i_s + 1
             else if ( local(la, i) ) then
                if ( i_r > size(bxasc%r_con%rcv) ) then
                   allocate(n_rcv(size(bxasc%r_con%rcv) + ChunkSize))
                   n_rcv(1:i_r-1) = bxasc%r_con%rcv(1:i_r-1)
                   deallocate(bxasc%r_con%rcv)
                   bxasc%r_con%rcv => n_rcv
                end if
                cnt_r                    = cnt_r + 1
                parr(lap%prc(j), 1)      = parr(lap%prc(j), 1) + 1
                pvol(lap%prc(j), 1)      = pvol(lap%prc(j), 1) + volume(sbx)
                bxasc%r_con%rcv(i_r)%nd  = i
                bxasc%r_con%rcv(i_r)%ns  = j
                bxasc%r_con%rcv(i_r)%sbx = sbx
                bxasc%r_con%rcv(i_r)%dbx = dbx
                bxasc%r_con%rcv(i_r)%pr  = get_proc(la, j)
                sh                       = 1
                sh(1:bxasc%dim)          = extent(sbx)
                bxasc%r_con%rcv(i_r)%sh  = sh
                i_r                      = i_r + 1
             end if
          end do
          deallocate(bi)
       end do
       call boxarray_destroy(bxai)
    end do
  end subroutine boxassoc_build_innards

  subroutine boxassoc_build(bxasc, lap, ng, nodal, cross, do_sum_boundary, idim)
    use bl_error_module

    integer,          intent(in)           :: ng
    logical,          intent(in)           :: nodal(:)
    type(layout_rep), intent(in), target   :: lap
    type(boxassoc),   intent(inout)        :: bxasc
    logical,          intent(in)           :: cross
    logical,          intent(in), optional :: do_sum_boundary
    integer,          intent(in), optional :: idim

    integer                         :: pv, rpv, spv, pi_r, pi_s, pcnt_r, pcnt_s
    type(boxarray)                  :: bxa, batmp
    type(layout)                    :: la, latmp
    integer                         :: lcnt_r, cnt_r, cnt_s, i_r, i_s, np, i
    integer, parameter              :: ChunkSize = 50
    integer, allocatable            :: pvol(:,:), ppvol(:,:), parr(:,:)
    type(local_copy_desc), pointer  :: n_cpy(:) => Null()
    type(comm_dsc), pointer         :: n_snd(:) => Null(), n_rcv(:) => Null()
    logical                         :: anynodal, ldsm

    if ( built_q(bxasc) ) call bl_error("boxassoc_build(): already built")

    la%lap       => lap
    bxa          =  get_boxarray(la)
    bxasc%dim    =  get_dim(bxa)
    bxasc%grwth  =  ng
    bxasc%nboxes =  nboxes(bxa)
    bxasc%cross  =  cross
    np           =  parallel_nprocs()
    anynodal     =  any( nodal .eqv. .true. )

    ldsm = .false.; if ( present(do_sum_boundary) ) ldsm = do_sum_boundary

    bxasc%idim = 0; if ( present(idim) ) bxasc%idim = idim
    if (bxasc%idim < 0 .or. bxasc%idim > bxasc%dim) then
       call bl_error("BOXASSOC_BUILD:", bxasc%idim)
    end if
    
    if (bxasc%idim .ne. 0) then
       bxasc%cross = .true.   ! idim > 0 is a special case of cross
    end if

    allocate(bxasc%nodal(bxasc%dim))
    allocate(parr(0:np-1,2))
    allocate(pvol(0:np-1,2))
    allocate(ppvol(0:np-1,2))
    allocate(bxasc%l_con%cpy(ChunkSize))
    allocate(bxasc%r_con%snd(ChunkSize))
    allocate(bxasc%r_con%rcv(ChunkSize))

    if ( anynodal ) then
       !
       ! Build a temporary layout to be used in intersection tests below.
       !
       call copy(batmp, bxa)
       call boxarray_nodalize(batmp, nodal)
       call build(latmp, batmp, boxarray_bbox(batmp), mapping = LA_LOCAL)  ! LA_LOCAL ==> bypass processor distribution calculation.
       call destroy(batmp)
    end if

    bxasc%nodal = nodal

    parr = 0; pvol = 0; lcnt_r = 0; cnt_r = 0; cnt_s = 0

    if ( ldsm ) then
       if (bxasc%idim .ne. 0) then
          call bl_error("BOXASSOC_BUILD: sumassoc_build_innards not supported for nonzero idim")
       end if
       call sumassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, &
            lcnt_r, cnt_s, cnt_r, pvol, parr)
    else
       call boxassoc_build_innards(bxasc, la, lap, latmp, ng, anynodal, bxasc%cross, &
            lcnt_r, cnt_s, cnt_r, pvol, parr)
    end if

    if ( anynodal ) call destroy(latmp)

    bxasc%l_con%ncpy = lcnt_r
    bxasc%r_con%nsnd = cnt_s
    bxasc%r_con%nrcv = cnt_r
    !
    ! Trim off unused space.
    !
    allocate(n_cpy(lcnt_r))
    n_cpy(1:lcnt_r) = bxasc%l_con%cpy(1:lcnt_r)
    deallocate(bxasc%l_con%cpy)
    bxasc%l_con%cpy => n_cpy

    allocate(n_snd(cnt_s))
    n_snd(1:cnt_s) = bxasc%r_con%snd(1:cnt_s)
    deallocate(bxasc%r_con%snd)
    bxasc%r_con%snd => n_snd

    allocate(n_rcv(cnt_r))
    n_rcv(1:cnt_r) = bxasc%r_con%rcv(1:cnt_r)
    deallocate(bxasc%r_con%rcv)
    bxasc%r_con%rcv => n_rcv
    !
    ! This region packs the src/recv boxes into processor order
    !
    do i = 0, np-1
       ppvol(i,1) = sum(pvol(0:i-1,1))
       ppvol(i,2) = sum(pvol(0:i-1,2))
    end do
    !
    ! Pack Receives maintaining original ordering
    !
    do i_r = 1, cnt_r
       i = bxasc%r_con%rcv(i_r)%pr
       bxasc%r_con%rcv(i_r)%pv = ppvol(i,1)
       pv = volume(bxasc%r_con%rcv(i_r)%dbx)
       ppvol(i,1) = ppvol(i,1) + pv
    end do
    !
    ! Pack Sends maintaining original ordering
    !
    do i_s = 1, cnt_s
       i = bxasc%r_con%snd(i_s)%pr
       bxasc%r_con%snd(i_s)%pv = ppvol(i,2)
       pv = volume(bxasc%r_con%snd(i_s)%dbx)
       ppvol(i,2) = ppvol(i,2) + pv
    end do
    !
    ! Now compute the volume of data that each processor expects.
    !
    pcnt_r = count(parr(:,1) /= 0 )
    pcnt_s = count(parr(:,2) /= 0 )
    bxasc%r_con%nrp  = pcnt_r
    bxasc%r_con%nsp  = pcnt_s
    bxasc%r_con%rvol = sum(pvol(:,1))
    bxasc%r_con%svol = sum(pvol(:,2))
    allocate(bxasc%r_con%str(pcnt_s))
    allocate(bxasc%r_con%rtr(pcnt_r))
    pi_r = 1; pi_s = 1; rpv  = 0; spv  = 0
    do i = 0, size(pvol,dim=1)-1
       if ( pvol(i,1) /= 0 ) then
          bxasc%r_con%rtr(pi_r)%sz = pvol(i,1)
          bxasc%r_con%rtr(pi_r)%pr = i
          bxasc%r_con%rtr(pi_r)%pv = rpv
          rpv  = rpv + pvol(i,1)
          pi_r = pi_r + 1
       end if
       if ( pvol(i,2) /= 0 ) then
          bxasc%r_con%str(pi_s)%sz = pvol(i,2)
          bxasc%r_con%str(pi_s)%pr = i
          bxasc%r_con%str(pi_s)%pv = spv
          spv  = spv + pvol(i,2)
          pi_s = pi_s + 1
       end if
    end do

  end subroutine boxassoc_build

  subroutine boxassoc_destroy(bxasc)
    use bl_error_module
    type(boxassoc), intent(inout) :: bxasc
    if ( .not. built_q(bxasc) ) call bl_error("boxassoc_destroy(): not built")
    deallocate(bxasc%nodal)
    deallocate(bxasc%l_con%cpy)
    deallocate(bxasc%r_con%snd)
    deallocate(bxasc%r_con%rcv)
    deallocate(bxasc%r_con%str)
    deallocate(bxasc%r_con%rtr)
  end subroutine boxassoc_destroy

  subroutine init_box_hash_bin(la, crsn)
    use bl_error_module
    type(layout), intent(inout) :: la
    integer, intent(in), optional :: crsn
    type(boxarray) :: ba
    integer, dimension(MAX_SPACEDIM) :: ext, vsz
    integer :: dm, i, j, k, n
    type(box) :: bx, cbx
    integer :: lcrsn
    integer :: sz
    type(box_hash_bin), pointer :: bins(:,:,:)
    integer, pointer :: ipv(:)

    dm = la%lap%dim
    ba = get_boxarray(la)
    vsz = 0; vsz(1:dm) = -Huge(1)
    do n = 1, nboxes(ba)
       vsz(1:dm) = max(vsz(1:dm),extent(get_box(ba,n)))
    end do
    if ( present(crsn) ) then
       lcrsn = crsn
    else
       lcrsn = maxval(vsz)
    end if
    la%lap%crsn = lcrsn
    bx = boxarray_bbox(ba)
    cbx = coarsen(bx, lcrsn)
    la%lap%plo = 0; la%lap%plo(1:dm) = lwb(cbx)
    la%lap%phi = 0; la%lap%phi(1:dm) = upb(cbx)
    la%lap%vshft = int_coarsen(vsz, lcrsn+1)
    allocate(la%lap%bins(la%lap%plo(1):la%lap%phi(1),la%lap%plo(2):la%lap%phi(2),la%lap%plo(3):la%lap%phi(3)))
    bins => la%lap%bins
    do k = la%lap%plo(3), la%lap%phi(3)
       do j = la%lap%plo(2), la%lap%phi(2)
          do i = la%lap%plo(1), la%lap%phi(1)
             allocate(bins(i,j,k)%iv(0))
          end do
       end do
    end do
    do n = 1, nboxes(ba)
       ext = 0; ext(1:dm) = int_coarsen(lwb(get_box(ba,n)), lcrsn)
       if ( .not. contains(cbx, ext(1:dm)) ) then
          call bl_error("init_box_hash_bin(): not contained!")
       end if
       sz = size(bins(ext(1),ext(2),ext(3))%iv)
       allocate(ipv(sz+1))
       ipv(1:sz) = bins(ext(1),ext(2),ext(3))%iv(1:sz)
       ipv(sz+1) = n
       deallocate(bins(ext(1),ext(2),ext(3))%iv)
       bins(ext(1),ext(2),ext(3))%iv => ipv
    end do
  end subroutine init_box_hash_bin

  function layout_get_box_intersector(la, bx) result(bi)
    use bl_error_module
    type(box_intersector), pointer :: bi(:)
    type(layout), intent(inout) :: la
    type(box), intent(in) :: bx
    type(box_hash_bin), pointer :: bins(:,:,:)
    integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
    integer :: dm
    type(box) :: bx1
    integer :: i, j, k, n
    type(boxarray) :: ba
    integer, parameter :: ChunkSize = 50
    integer :: cnt
    type(box_intersector), pointer :: tbi(:)  => Null()

    if (.not. associated(la%lap%bins)) call init_box_hash_bin(la)

    allocate(bi(ChunkSize))

    dm = la%lap%dim
    ba = get_boxarray(la)
    bins => la%lap%bins
    bx1 = coarsen(bx, la%lap%crsn)
    lo = 0; lo(1:dm) = lwb(bx1)
    hi = 0; hi(1:dm) = upb(bx1)
    cnt = 0
    select case ( dm ) 
    case (3)
       do k = max(lo(3)-la%lap%vshft(3)-1,la%lap%plo(3)), min(hi(3)+la%lap%vshft(3), la%lap%phi(3))
          do j = max(lo(2)-la%lap%vshft(2)-1,la%lap%plo(2)), min(hi(2)+la%lap%vshft(2), la%lap%phi(2))
             do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))
                do n = 1, size(bins(i,j,k)%iv)
                   bx1 = intersection(bx, get_box(ba,bins(i,j,k)%iv(n)))
                   if ( empty(bx1) ) cycle
                   cnt = cnt + 1

                   if (cnt > size(bi)) then
                      allocate(tbi(size(bi) + ChunkSize))
                      tbi(1:cnt-1) = bi(1:cnt-1)
                      deallocate(bi)
                      bi => tbi
                   end if

                   bi(cnt)%i  = bins(i,j,k)%iv(n)
                   bi(cnt)%bx = bx1
                end do
             end do
          end do
       end do
    case (2)
       k = 0
       do j = max(lo(2)-la%lap%vshft(2)-1,la%lap%plo(2)), min(hi(2)+la%lap%vshft(2), la%lap%phi(2))
          do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))
             do n = 1, size(bins(i,j,k)%iv)
                bx1 = intersection(bx, get_box(ba,bins(i,j,k)%iv(n)))
                if ( empty(bx1) ) cycle
                cnt = cnt + 1

                if (cnt > size(bi)) then
                   allocate(tbi(size(bi) + ChunkSize))
                   tbi(1:cnt-1) = bi(1:cnt-1)
                   deallocate(bi)
                   bi => tbi
                end if

                bi(cnt)%i  = bins(i,j,k)%iv(n)
                bi(cnt)%bx = bx1
             end do
          end do
       end do
    case (1)
       k = 0
       j = 0
       do i = max(lo(1)-la%lap%vshft(1)-1,la%lap%plo(1)), min(hi(1)+la%lap%vshft(1), la%lap%phi(1))
          do n = 1, size(bins(i,j,k)%iv)
             bx1 = intersection(bx, get_box(ba,bins(i,j,k)%iv(n)))
             if ( empty(bx1) ) cycle
             cnt = cnt + 1

             if (cnt > size(bi)) then
                allocate(tbi(size(bi) + ChunkSize))
                tbi(1:cnt-1) = bi(1:cnt-1)
                deallocate(bi)
                bi => tbi
             end if

             bi(cnt)%i  = bins(i,j,k)%iv(n)
             bi(cnt)%bx = bx1
          end do
       end do
    end select

    allocate(tbi(cnt))
    tbi(1:cnt) = bi(1:cnt)
    deallocate(bi)
    bi => tbi

  end function layout_get_box_intersector

  end module layout_module
