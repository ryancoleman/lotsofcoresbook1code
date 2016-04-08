module threadbox_module

  use bl_error_module
  use parallel
  use omp_module
  use layout_module
  use multifab_module

  implicit none

  integer, save :: numthreads, ng, nb, ndim
  integer, save :: numboxgroups, boxgroupsize, nthreadsperbox
  integer, save :: nthreads_d(3)

  logical, allocatable, save :: worktodo(:)
  integer, allocatable, save :: tb_lo(:,:), tb_hi(:,:) 
  integer, allocatable, save :: tb_glo(:,:), tb_ghi(:,:) 

  type blocks  ! in_a_threadbox
     integer :: nblocks
     integer, pointer :: lo(:,:) => Null()
     integer, pointer :: hi(:,:) => Null()
  end type blocks

  type(blocks), allocatable, save :: allblocks(:)

  !$omp threadprivate(worktodo,tb_lo,tb_hi,tb_glo,tb_ghi,allblocks)

  private

  public ::build_threadbox, destroy_threadbox,  &
       tb_get_valid_lo, tb_get_valid_hi, tb_get_grown_lo, tb_get_grown_hi, &
       tb_get_block_lo, tb_get_block_hi, tb_get_nblocks, &
       tb_multifab_setval, tb_worktodo

contains

  subroutine destroy_threadbox()
    integer :: ibox
    !$omp parallel private(ibox)
    if (allocated(allblocks)) then
       do ibox=1,nb
          if (associated(allblocks(ibox)%lo)) then
             deallocate(allblocks(ibox)%lo)
          end if
          if (associated(allblocks(ibox)%hi)) then
             deallocate(allblocks(ibox)%hi)
          end if
       end do
       deallocate(allblocks)
    end if
    if (allocated(tb_lo)) deallocate(tb_lo)
    if (allocated(tb_hi)) deallocate(tb_hi)
    if (allocated(tb_glo)) deallocate(tb_glo)
    if (allocated(tb_ghi)) deallocate(tb_ghi)
    if (allocated(worktodo)) deallocate(worktodo)
    !$omp end parallel
  end subroutine destroy_threadbox

  subroutine build_threadbox(la, ng_in)
    use probin_module, only : tb_split_dim, tb_collapse_boxes, tb_idim_more, tb_idim_less, &
         tb_blocksize_x, tb_blocksize_y, tb_blocksize_z
    implicit none
    type(layout), intent(in) :: la
    integer, intent(in) :: ng_in

    call destroy_threadbox()

    if (nboxes(la) < parallel_nprocs()) then
       call bl_error('ERROR: more MPI ranks than boxes')
    end if

    ndim = la%lap%dim
    if (ndim .eq. 2) then
       call bl_error('2D not supported')
    end if
    ng = ng_in
    nb = nlocal(la) ! number of local boxes
    numthreads = omp_get_max_threads()

    if (tb_idim_more .eq. tb_idim_less) then
       call bl_error("threadbox_module: tb_idim_more .eq. tb_idim_less")
    end if
    if (tb_idim_more < 1 .or. tb_idim_more > 3) then
       call bl_error("threadbox_module: invalid tb_idim_more")
    end if
    if (tb_idim_less < 1 .or. tb_idim_less > 3) then
       call bl_error("threadbox_module: invalid tb_idim_less")
    end if

    !$omp parallel

    allocate(tb_lo (3,nb))
    allocate(tb_hi (3,nb))
    allocate(tb_glo(3,nb))
    allocate(tb_ghi(3,nb))
    allocate(worktodo(nb))
    allocate(allblocks(nb))

    !$omp single
    call setup_boxgroups(tb_collapse_boxes) 

    call init_thread_topology(nthreadsperbox, nthreads_d, tb_split_dim, &
         tb_idim_more, tb_idim_less)
    !$omp end single

    ! use the masteter thread because it uses bl_error, which in turn uses MPI
    !$omp master
    call check_boxsize(la)
    !$omp end master
    !$omp barrier

    call init_worktodo() 

    call init_threadbox(la)

    call init_allblocks(tb_blocksize_x,tb_blocksize_y,tb_blocksize_z)

    !$omp end parallel

  end subroutine build_threadbox


  ! Boxes are divided into groups.
  ! Then boxes in a group are divided by the number of threads.
  subroutine setup_boxgroups(collapse)
    implicit none
    logical, intent(in) :: collapse 
    if (collapse) then
       boxgroupsize = greatest_common_factor(nb, numthreads)
    else
       boxgroupsize = 1
    end if
    numboxgroups = nb / boxgroupsize
    nthreadsperbox = numthreads / boxgroupsize
  contains
    recursive function greatest_common_factor(a, b) result(c)
      implicit none
      integer :: c
      integer, intent(in) :: a, b
      if (a.eq.b) then
         c = a
      else if (a.gt.b) then
         c = greatest_common_factor(a-b,b)
      else
         c = greatest_common_factor(a,b-a)
      end if
    end function greatest_common_factor
  end subroutine setup_boxgroups


  subroutine init_thread_topology(n, n3d, d, imore, iless)
    implicit none
    integer, intent(in) :: n, d, imore, iless
    integer, intent(out) :: n3d(3)

    integer, allocatable :: facs(:)
    integer :: lfacs, f, nfac, rem, i, j, nmin, nmax
    integer :: n2d(2)

    if (d.eq.1) then
       n3d = 1
       n3d(imore) = n
       return
    end if

    lfacs = int(log(dble(n))/log(2.d0))
    allocate(facs(lfacs))
    
    nfac = 0
    f = 2
    rem = n
    do while (rem .ne. 1)
       if (mod(rem, f) .eq. 0) then
          rem = rem/f
          nfac = nfac+1
          facs(nfac) = f
       else
          f = f + 1
       end if
    end do

    if (d .eq. 3) then
       n3d = 1
       do i = nfac, 1, -1
          j = minloc(n3d,1)
          n3d(j) = n3d(j) * facs(i)
       end do
    else
       n2d = 1
       do i = nfac, 1, -1
          j = minloc(n2d,1)
          n2d(j) = n2d(j) * facs(i)
       end do
       n3d(1) = 1
       n3d(2:3) = n2d
    end if

    nmin = minval(n3d)
    nmax = maxval(n3d)

    n3d = n / (nmin*nmax)
    n3d(iless) = nmin
    n3d(imore) = nmax

    deallocate(facs)

  end subroutine init_thread_topology


  subroutine check_boxsize(la)
    implicit none
    type(layout),intent(in) :: la

    integer :: ilocal, iglobal, idim
    integer :: box_size_1(la%lap%dim), box_size_i(la%lap%dim)
    type(box) :: bx

    ilocal = 1
    iglobal = global_index(la, ilocal)
    bx = get_box(la, iglobal)
    box_size_1 = box_extent(bx)

    do ilocal = 2, nb

       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       box_size_i = box_extent(bx)

       do idim = 1, ndim
          if (box_size_1(idim) .ne. box_size_i(idim)) then
             call bl_error('All boxes in the same MPI task must have the same size')
             ! make sure all boxes have the same size
             ! otherwise, the load is unbalanced
          end if
       end do

    end do

    do ilocal = 1, nb

       iglobal = global_index(la, ilocal)
       bx = get_box(la, iglobal)
       box_size_i = box_extent(bx)

       do idim = 1, ndim
          if (box_size_i(idim) < nthreads_d(idim)) then
             print *, 'Box #', iglobal, 'idim =', idim, ' box size = ', box_size_i(idim), &
                  '  threads in this direction', nthreads_d(idim)
             call bl_error("threadbox_module: Too many threads for such small box") 
          end if
       end do

    end do

  end subroutine check_boxsize


  ! All threads will work on a group, but not necessarily on a box.
  subroutine init_worktodo()
    implicit none
    integer :: tid, ibox, ib2, ib3
    tid = omp_get_thread_num()
    ib3 = (tid/nthreadsperbox) + 1 ! box index in boxgroup for this thread
    do ibox=1,nb
       ib2 = mod(ibox-1,boxgroupsize) + 1 ! box index in a group
       if (ib2 .eq. ib3) then
          worktodo(ibox) = .true.
       else
          worktodo(ibox) = .false.
       end if
    end do
  end subroutine init_worktodo


  subroutine init_threadbox(la)
    implicit none
    type(layout), intent(in) :: la
    
    integer :: tid
    integer, allocatable :: xsize(:), ysize(:), zsize(:)
    integer, allocatable :: xstart(:), ystart(:), zstart(:)
    integer :: my_box_size(3), my_box_lo(3), zero_lo(3), zero_hi(3)
    integer :: i,j,k, itbox
    integer :: igroup, ibg, ibox, istart
    double precision :: stride
    type(box) :: bx

    tid = omp_get_thread_num()

    allocate(xsize (nthreads_d(1)))
    allocate(ysize (nthreads_d(2)))
    allocate(zsize (nthreads_d(3)))
    allocate(xstart(nthreads_d(1)))
    allocate(ystart(nthreads_d(2)))
    allocate(zstart(nthreads_d(3)))

    stride = dble(nthreadsperbox) / dble(numboxgroups)

    ibox=0
    do igroup=1,numboxgroups

       istart = nint((igroup-1)*stride)

       ! For the group of threads working on the same box, 
       ! this is its new ID in the group, ranging from 0 to nthreadsperbox-1.
       itbox = mod(tid+istart, nthreadsperbox)

       ! convert the 1D index into 3D indices
       ! Here the 1D index itbox starts from 0; the 3D indices start from 1
       i = mod(itbox,nthreads_d(1)) + 1
       j = mod(int(itbox/nthreads_d(1)), nthreads_d(2)) + 1
       k = int(itbox / (nthreads_d(1)*nthreads_d(2))) + 1

       do ibg=1,boxgroupsize

          ibox = ibox+1

          if (.not. worktodo(ibox)) then
             tb_lo(:,ibox) = HUGE(0)
             tb_hi(:,ibox) = -(HUGE(0)-1)
             tb_glo(:,ibox) = HUGE(0)
             tb_ghi(:,ibox) = -(HUGE(0)-1)
          else

             bx = get_box(la, global_index(la, ibox))
             my_box_size = box_extent(bx)
             my_box_lo = box_lwb(bx)

             !
             ! valid box
             !
             call split_domain(my_box_size(1), nthreads_d(1), xsize, xstart)
             call split_domain(my_box_size(2), nthreads_d(2), ysize, ystart)
             call split_domain(my_box_size(3), nthreads_d(3), zsize, zstart)
             
             zero_lo(1) = xstart(i)
             zero_lo(2) = ystart(j)
             zero_lo(3) = zstart(k)
             zero_hi(1) = xstart(i) + xsize(i) - 1
             zero_hi(2) = ystart(j) + ysize(j) - 1
             zero_hi(3) = zstart(k) + zsize(k) - 1
             
             tb_lo(:,ibox) = zero_lo + my_box_lo
             tb_hi(:,ibox) = zero_hi + my_box_lo

             !
             ! grown box
             !
             call split_domain(my_box_size(1)+2*ng, nthreads_d(1), xsize, xstart)
             call split_domain(my_box_size(2)+2*ng, nthreads_d(2), ysize, ystart)
             call split_domain(my_box_size(3)+2*ng, nthreads_d(3), zsize, zstart)

             zero_lo(1) = xstart(i)
             zero_lo(2) = ystart(j)
             zero_lo(3) = zstart(k)
             zero_hi(1) = xstart(i) + xsize(i) - 1
             zero_hi(2) = ystart(j) + ysize(j) - 1
             zero_hi(3) = zstart(k) + zsize(k) - 1
             
             tb_glo(:,ibox) = zero_lo + my_box_lo - ng
             tb_ghi(:,ibox) = zero_hi + my_box_lo - ng

          end if
          
       end do
    end do

    deallocate(xsize, ysize, zsize, xstart, ystart, zstart)

  end subroutine init_threadbox

  subroutine split_domain(totalsize, n, isize, start)
    implicit none
    integer, intent(in) :: totalsize, n
    integer, intent(out) :: isize(n), start(n)
    
    integer :: i, iavg, ileft
    
    iavg = int(totalsize/n)
    ileft = totalsize - iavg*n
    
    start(1) = 0
    isize(1) = iavg
    do i=2, n
       start(i) = start(i-1) + isize(i-1)
       if (ileft > 0) then
          isize(i) = iavg + 1
          ileft = ileft - 1
       else
          isize(i) = iavg
       end if
    end do
    
  end subroutine split_domain


  subroutine init_allblocks(blksizex,blksizey,blksizez)
    implicit none
    integer,intent(in) :: blksizex,blksizey,blksizez
    integer :: idim, ibox, nbk(3), tbsize(3)
    integer :: i,j,k,ibk,nblk
    integer :: bksize(3), zero_lo(3), zero_hi(3)
    integer, allocatable :: xsize(:), ysize(:), zsize(:)
    integer, allocatable :: xstart(:), ystart(:), zstart(:)

    bksize(1) = blksizex
    bksize(2) = blksizey
    bksize(3) = blksizez

    do ibox=1,nb

       if (.not. worktodo(ibox)) then

          allblocks(ibox)%nblocks = 0

       else

          tbsize = tb_hi(:,ibox) - tb_lo(:,ibox) + 1

          do idim=1,3
             if (bksize(idim) <= 0) then
                nbk(idim) = 1
             else
                nbk(idim) = max(int(tbsize(idim)/bksize(idim)),1)
             end if
          end do

          nblk = nbk(1)*nbk(2)*nbk(3)
          allblocks(ibox)%nblocks = nblk
          allocate(allblocks(ibox)%lo(3,nblk))
          allocate(allblocks(ibox)%hi(3,nblk))

          allocate(xsize (nbk(1)))
          allocate(ysize (nbk(2)))
          allocate(zsize (nbk(3)))
          allocate(xstart(nbk(1)))
          allocate(ystart(nbk(2)))
          allocate(zstart(nbk(3)))

          call split_domain(tbsize(1), nbk(1), xsize, xstart)
          call split_domain(tbsize(2), nbk(2), ysize, ystart)
          call split_domain(tbsize(3), nbk(3), zsize, zstart)

          ibk = 1
          do k = 1, nbk(3)
             do j = 1, nbk(2)
                do i = 1, nbk(1)

                   zero_lo(1) = xstart(i)
                   zero_lo(2) = ystart(j)
                   zero_lo(3) = zstart(k)
                   zero_hi(1) = xstart(i) + xsize(i) - 1
                   zero_hi(2) = ystart(j) + ysize(j) - 1
                   zero_hi(3) = zstart(k) + zsize(k) - 1

                   allblocks(ibox)%lo(:,ibk) = zero_lo + tb_lo(:,ibox)
                   allblocks(ibox)%hi(:,ibk) = zero_hi + tb_lo(:,ibox)

                   ibk = ibk + 1
                end do
             end do
          end do

          deallocate(xsize,ysize,zsize,xstart,ystart,zstart)

       end if
    end do

  end subroutine init_allblocks


  function tb_get_valid_lo(ilocal) result (lo)
    implicit none
    integer, intent(in) :: ilocal
    integer, dimension(3) :: lo
    lo = tb_lo(:,ilocal)
  end function tb_get_valid_lo

  function tb_get_valid_hi(ilocal) result (hi)
    implicit none
    integer, intent(in) :: ilocal
    integer, dimension(3) :: hi
    hi = tb_hi(:,ilocal)
  end function tb_get_valid_hi


  function tb_get_grown_lo(ilocal) result (lo)
    implicit none
    integer, intent(in) :: ilocal
    integer, dimension(3) :: lo
    lo = tb_glo(:,ilocal)
  end function tb_get_grown_lo

  function tb_get_grown_hi(ilocal) result (hi)
    implicit none
    integer, intent(in) :: ilocal
    integer, dimension(3) :: hi
    hi = tb_ghi(:,ilocal)
  end function tb_get_grown_hi


  function tb_get_nblocks(ilocal) result (nblk)
    implicit none
    integer, intent(in) :: ilocal
    integer :: nblk
    nblk = allblocks(ilocal)%nblocks
  end function tb_get_nblocks

  function tb_get_block_lo(iblock, ilocal) result(lo)
    implicit none
    integer, intent(in) :: iblock, ilocal
    integer, dimension(3) :: lo
    lo = allblocks(ilocal)%lo(:,iblock)
  end function tb_get_block_lo

  function tb_get_block_hi(iblock, ilocal) result(hi)
    implicit none
    integer, intent(in) :: iblock, ilocal
    integer, dimension(3) :: hi
    hi = allblocks(ilocal)%hi(:,iblock)
  end function tb_get_block_hi


  subroutine tb_multifab_setval(mf, val, all)
    type(multifab), intent(inout) :: mf
    real(dp_t), intent(in) :: val
    logical, intent(in), optional :: all

    logical :: lall
    integer :: ib, wlo(3), whi(3), ngmf
    double precision, pointer :: p(:,:,:,:)

    if (present(all)) then
       lall = all
    else
       lall = .false.
    end if

    ngmf = nghost(mf)

    if (lall) then
       if (ngmf .ne. ng .and. ngmf .ne. 0) then
          print *, 'ng, ngmf =', ng, ngmf
          call bl_error("threadbox: do not know how to do setval in this case.")
       end if
    end if

    if (lall .and. ngmf.eq.ng) then
       !$omp parallel private(ib, wlo, whi, p)
       do ib = 1, nfabs(mf)
          if (.not.tb_worktodo(ib)) cycle
          p => dataptr(mf, ib)
          wlo = tb_get_grown_lo(ib)
          whi = tb_get_grown_hi(ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    else
       !$omp parallel private(ib, wlo, whi, p)
       do ib = 1, nfabs(mf)
          if (.not.tb_worktodo(ib)) cycle
          p => dataptr(mf, ib)
          wlo = tb_get_valid_lo(ib)
          whi = tb_get_valid_hi(ib)
          p(wlo(1):whi(1),wlo(2):whi(2),wlo(3):whi(3),:) = val
       end do
       !$omp end parallel
    end if

  end subroutine tb_multifab_setval


  function tb_worktodo(ilocal) result(r)
    implicit none
    logical :: r
    integer,intent(in) :: ilocal
    r = worktodo(ilocal)
  end function tb_worktodo


end module threadbox_module
