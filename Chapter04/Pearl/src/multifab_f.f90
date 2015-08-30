module multifab_module

  use bl_error_module
  use layout_module
  use fab_module

  implicit none

  type multifab
     !private
     integer      :: dim   = 0
     integer      :: nc    = 1
     integer      :: ng    = 0
     type(layout) :: la
     logical,   pointer :: nodal(:) => Null()
     type(fab), pointer :: fbs(:)   => Null()
  end type multifab

  interface built_q
     module procedure multifab_built_q
  end interface

  interface build
     module procedure multifab_build
  end interface

  interface local_index
     module procedure  multifab_local_index
  end interface

  interface global_index
     module procedure  multifab_global_index
  end interface

  interface destroy
     module procedure multifab_destroy
  end interface

  interface dataptr
     module procedure multifab_dataptr
     module procedure multifab_dataptr_bx_c
  end interface

  interface get_boxarray
     module procedure multifab_get_boxarray
  end interface

  interface get_ibox
     module procedure multifab_get_ibox
  end interface
  
  interface get_box
     module procedure multifab_get_box
  end interface

  interface nfabs
     module procedure multifab_nfabs
  end interface

  interface nghost
     module procedure multifab_nghost
  end interface

  interface get_layout
     module procedure multifab_get_layout
  end interface

  interface fill_boundary
     module procedure multifab_fill_boundary
     module procedure multifab_fill_boundary_c
  end interface

  interface ncomp
     module procedure multifab_ncomp
  end interface

  private :: cpy_d
  private :: reshape_d_4_1, reshape_d_1_4
  private :: mf_fb_fancy_double, mf_fb_easy_double

contains

  function multifab_local_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(multifab),  intent(in) :: mf
    integer                     :: r
    r = layout_local_index(mf%la,i)
  end function multifab_local_index

  pure function multifab_global_index(mf,i) result(r)
    integer,         intent(in) :: i
    type(multifab),  intent(in) :: mf
    integer                     :: r
    r = layout_global_index(mf%la,i)
  end function multifab_global_index

  pure function multifab_ncomp(mf) result(r)
    integer :: r
    type(multifab), intent(in) :: mf
    r = mf%nc
  end function multifab_ncomp

  pure function multifab_built_q(mf) result(r)
    logical :: r
    type(multifab), intent(in) :: mf
    r = mf%dim /= 0
  end function multifab_built_q

  function multifab_get_layout(mf) result(r)
    type(layout) :: r
    type(multifab), intent(in) :: mf
    r = mf%la
  end function multifab_get_layout

  subroutine multifab_build(mf, la, nc, ng, nodal, stencil)
    type(multifab), intent(out)   :: mf
    type(layout),   intent(in )   :: la
    integer, intent(in), optional :: nc, ng
    logical, intent(in), optional :: nodal(:), stencil
    integer :: i, lnc, lng
    if ( built_q(mf) ) call bl_error("MULTIFAB_BUILD: already built")
    lng = 0; if ( present(ng) ) lng = ng
    lnc = 1; if ( present(nc) ) lnc = nc
    mf%dim = get_dim(la)
    mf%la  = la
    mf%nc  = lnc
    mf%ng  = lng
    allocate(mf%nodal(mf%dim))
    mf%nodal = .False.; if ( present(nodal) ) mf%nodal = nodal(1:mf%dim)
    allocate(mf%fbs(nlocal(mf%la)))

    do i = 1, nlocal(mf%la)
      call build( &
           mf%fbs(i), get_box(mf%la, global_index(mf%la,i)), &
           mf%nc, mf%ng, mf%nodal,  &
           alloc = .true., stencil = stencil)
    end do
  end subroutine multifab_build

  subroutine multifab_destroy(mf)
    type(multifab), intent(inout) :: mf
    integer :: i
    do i = 1, nlocal(mf%la)
       call destroy(mf%fbs(i))
    end do
    deallocate(mf%fbs)
    deallocate(mf%nodal)
    mf%dim = 0
    mf%nc  = 0
    mf%ng  = 0
  end subroutine multifab_destroy

  pure function multifab_nfabs(mf) result(r)
    type(multifab), intent(in) :: mf
    integer :: r
    r = nlocal(mf%la)
  end function multifab_nfabs

  pure function multifab_nghost(mf) result(r)
    type(multifab), intent(in) :: mf
    integer :: r
    r = mf%ng
  end function multifab_nghost

  function multifab_get_boxarray(mf) result(r)
    type(boxarray) :: r
    type(multifab), intent(in) :: mf
    r = get_boxarray(mf%la)
  end function multifab_get_boxarray

  pure function multifab_get_box(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_box(mf%fbs(i))
  end function multifab_get_box

  pure function multifab_get_ibox(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    type(box) :: r
    r = get_ibox(mf%fbs(i))
  end function multifab_get_ibox

  function multifab_dataptr(mf, i) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i
    real(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i))
  end function multifab_dataptr

  function multifab_dataptr_bx_c(mf, i, bx, c, nc) result(r)
    type(multifab), intent(in) :: mf
    integer, intent(in) :: i, c
    integer, intent(in), optional :: nc
    type(box), intent(in) :: bx
    real(dp_t), pointer :: r(:,:,:,:)
    r => dataptr(mf%fbs(i), bx, c, nc)
  end function multifab_dataptr_bx_c

  subroutine cpy_d(out, in, filter)
     real(dp_t), intent(inout) :: out(:,:,:,:)
     real(dp_t), intent(in   ) ::  in(:,:,:,:)
     integer                   :: i, j, k, n
     integer :: s(4)
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( present(filter) ) then
       call filter(out,in)
    else
       s = shape(out)
       !$OMP PARALLEL DO PRIVATE(i,j,k,n) COLLAPSE(3)
       do n = 1, s(4)
          do k = 1, s(3)
             do j = 1, s(2)
                do i = 1, s(1)
                   out(i,j,k,n) = in(i,j,k,n)
                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine cpy_d

  subroutine reshape_d_4_1(dst,ic,src)
    real(dp_t),intent(in)    :: src(:,:,:,:)
    real(dp_t),intent(inout) :: dst(:)
    integer,intent(in)       :: ic
    integer                  :: i,j,k,n,c,s(4),s12,s123
    s = shape(src)
    s12 = s(1)*s(2)
    s123 = s12*s(3)
    !$OMP PARALLEL DO PRIVATE(i,j,k,n,c) COLLAPSE(3)
    do n = 1, s(4)
       do k = 1, s(3)
          do j = 1, s(2)
             do i = 1, s(1) 
                c = ic + (i-1) + (j-1)*s(1) + (k-1)*s12 + (n-1)*s123
                dst(c) = src(i,j,k,n)
             end do
          end do
       end do
    end do
    !$OMP END PARALLEL DO
  end subroutine reshape_d_4_1

  subroutine reshape_d_1_4(dst,src,ic,sh,filter)
    real(dp_t),intent(in)    :: src(:)
    real(dp_t),intent(inout) :: dst(:,:,:,:)
    integer,intent(in)       :: ic
    integer,intent(in)       :: sh(:)
    integer                  :: i,j,k,n,c,s(4),s12,s123
    real(dp_t), allocatable  :: ptmp(:,:,:,:)
    interface
       subroutine filter(out, in)
         use bl_types
         real(dp_t), intent(inout) :: out(:,:,:,:)
         real(dp_t), intent(in   ) ::  in(:,:,:,:)
       end subroutine filter
    end interface
    optional filter
    if ( size(sh) /= 4 ) call bl_error("reshape_d_1_4: how did this happen?")
    s = shape(dst)
    s12 = s(1)*s(2)
    s123 = s12*s(3)
    if ( present(filter) ) then
       allocate(ptmp(sh(1),sh(2),sh(3),sh(4)))
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,c) COLLAPSE(3)
       do n = 1, s(4)
          do k = 1, s(3)
             do j = 1, s(2)
                do i = 1, s(1)
                   c = ic + (i-1) + (j-1)*s(1) + (k-1)*s12 + (n-1)*s123
                   ptmp(i,j,k,n) = src(c)
                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO
       call filter(dst, ptmp)
       deallocate(ptmp)
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,n,c) COLLAPSE(3)
       do n = 1, s(4)
          do k = 1, s(3)
             do j = 1, s(2)
                do i = 1, s(1)
                   c = ic + (i-1) + (j-1)*s(1) + (k-1)*s12 + (n-1)*s123
                   dst(i,j,k,n) = src(c)
                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if
  end subroutine reshape_d_1_4

  subroutine mf_fb_fancy_double(mf, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer     :: p(:,:,:,:), p1(:,:,:,:), p2(:,:,:,:)
    integer,    allocatable :: rst(:)
    integer,    parameter   :: tag = 1102
    integer                 :: i, ii, jj, np, sh(MAX_SPACEDIM+1)
    type(boxassoc)          :: bxasc
    real(dp_t), allocatable :: g_snd_d(:), g_rcv_d(:)

    bxasc = layout_boxassoc(mf%la, ng, mf%nodal, lcross, idim)

    do i = 1, bxasc%l_con%ncpy
       ii  =  local_index(mf,bxasc%l_con%cpy(i)%nd)
       jj  =  local_index(mf,bxasc%l_con%cpy(i)%ns)
       p1  => dataptr(mf%fbs(ii), bxasc%l_con%cpy(i)%dbx, c, nc)
       p2  => dataptr(mf%fbs(jj), bxasc%l_con%cpy(i)%sbx, c, nc)
       call cpy_d(p1,p2)
    end do

    np = parallel_nprocs()

    if (np == 1) return

    allocate(g_snd_d(nc*bxasc%r_con%svol))
    allocate(g_rcv_d(nc*bxasc%r_con%rvol))

    do i = 1, bxasc%r_con%nsnd
       p => dataptr(mf, local_index(mf,bxasc%r_con%snd(i)%ns), bxasc%r_con%snd(i)%sbx, c, nc)
       call reshape_d_4_1(g_snd_d, 1 + nc*bxasc%r_con%snd(i)%pv, p)
    end do

    allocate(rst(bxasc%r_con%nrp))
    do i = 1, bxasc%r_con%nrp
       rst(i) = parallel_irecv_dv(g_rcv_d(1+nc*bxasc%r_con%rtr(i)%pv:), &
            nc*bxasc%r_con%rtr(i)%sz, bxasc%r_con%rtr(i)%pr, tag)
    end do
    do i = 1, bxasc%r_con%nsp
       call parallel_send_dv(g_snd_d(1+nc*bxasc%r_con%str(i)%pv), &
            nc*bxasc%r_con%str(i)%sz, bxasc%r_con%str(i)%pr, tag)
    end do
    call parallel_wait(rst)

    do i = 1, bxasc%r_con%nrcv
       sh = bxasc%r_con%rcv(i)%sh
       sh(4) = nc
       p => dataptr(mf, local_index(mf,bxasc%r_con%rcv(i)%nd), bxasc%r_con%rcv(i)%dbx, c, nc)
       call reshape_d_1_4(p, g_rcv_d, 1 + nc*bxasc%r_con%rcv(i)%pv, sh)
    end do

  end subroutine mf_fb_fancy_double


  subroutine mf_fb_easy_double(mf, c, nc, ng, lcross, idim)
    type(multifab), intent(inout) :: mf
    integer,        intent(in)    :: c, nc, ng
    logical,        intent(in)    :: lcross
    integer, intent(in), optional :: idim

    real(dp_t), pointer :: pdst(:,:,:,:), psrc(:,:,:,:)
    type(boxarray)      :: bxai
    type(box)           :: abx
    integer             :: i, j, ii, proc, sz, il, jl
    integer             :: shft(2*3**mf%dim,mf%dim), sh(MAX_SPACEDIM+1)
    integer, parameter  :: tag = 1101
    real(dp_t), allocatable :: buf(:)

    allocate(buf(64))

    do i = 1, nboxes(mf%la)
       if (.not.remote(mf%la,i)) il = local_index(mf,i)
       call boxarray_bndry_periodic(bxai, mf%la%lap%pd, get_box(mf%la,i), mf%nodal, &
            mf%la%lap%pmask, ng, shft, lcross, idim)
       do j = 1, nboxes(mf%la)
          if (.not.remote(mf%la,j)) jl = local_index(mf,j)
          if ( remote(mf%la,i) .and. remote(mf%la,j) ) cycle
          do ii = 1, nboxes(bxai)
             abx = intersection(get_box(mf%la,j), get_box(bxai,ii))
             if ( .not. empty(abx) ) then
                if ( local(mf%la,i) .and. local(mf%la,j) ) then
                   psrc => dataptr(mf, jl, abx, c, nc)
                   pdst => dataptr(mf, il, shift(abx,-shft(ii,:)), c, nc)
                   pdst = psrc
                else if ( local(mf%la,j) ) then ! must send
                   psrc => dataptr(mf, jl, abx, c, nc)
                   sz = size(psrc)
                   if (sz .gt. size(buf)) then
                      deallocate(buf)
                      allocate(buf(sz))
                   end if
                   call reshape_d_4_1(buf,1,psrc)
                   proc = get_proc(mf%la, i)
                   call parallel_send_dv(buf, sz, proc, tag)
                else if ( local(mf%la,i) ) then  ! must recv
                   pdst => dataptr(mf, il, shift(abx,-shft(ii,:)), c, nc)
                   sz = size(pdst)
                   if (sz .gt. size(buf)) then
                      deallocate(buf)
                      allocate(buf(sz))
                   end if
                   proc = get_proc(mf%la,j)
                   call parallel_recv_dv(buf, sz, proc, tag)
                   call reshape_d_1_4(pdst,buf,1,sh)
                end if
             end if
          end do
       end do
       call destroy(bxai)
    end do
  end subroutine mf_fb_easy_double


  subroutine multifab_fill_boundary_c(mf, c, nc, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    integer, intent(in)           :: c, nc
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    integer :: lng
    logical :: lcross
    logical, parameter :: aggregate = .true.
    lcross  = .false.; if ( present(cross)  ) lcross  = cross
    lng     = mf%ng;   if ( present(ng)     ) lng     = ng
    if ( lng > mf%ng      ) call bl_error("MULTIFAB_FILL_BOUNDARY_C: ng too large", lng)
    if ( mf%nc < (c+nc-1) ) call bl_error('MULTIFAB_FILL_BOUNDARY_C: nc too large', nc)
    if ( present(idim) ) then
       if (idim > 0) lcross = .true. 
    end if
   
    if (aggregate) then
       call mf_fb_fancy_double(mf, c, nc, lng, lcross, idim)
    else
       call mf_fb_easy_double(mf, c, nc, lng, lcross, idim)
    end if
  end subroutine multifab_fill_boundary_c

  subroutine multifab_fill_boundary(mf, ng, cross, idim)
    type(multifab), intent(inout) :: mf
    integer, intent(in), optional :: ng, idim
    logical, intent(in), optional :: cross
    call multifab_fill_boundary_c(mf, 1, mf%nc, ng, cross, idim)
  end subroutine multifab_fill_boundary
  
end module multifab_module
