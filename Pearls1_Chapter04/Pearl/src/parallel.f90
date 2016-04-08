
      module mpi
      implicit none
      include 'mpif.h'
      end module mpi

!! MPI wrappers
!! These wrappers are used so that a non-MPI version can coexist
!! with the MPI version.
!! Not all of MPI is supported.  Only sufficient to execute
!! fBoxLib calls.

module parallel

  ! Assumption:
  ! 1) The user has not replaced the default error handler which
  !    is MPI_ERRORS_ARE_FATAL

  use bl_types; use mpi

  implicit none

  integer, parameter :: parallel_root = 0
  integer, parameter, private :: io_processor_node = parallel_root
  integer, private :: m_nprocs = -1
  integer, private :: m_myproc = -1
  integer, private :: m_comm   = -1

  interface parallel_reduce
     module procedure parallel_reduce_d
     module procedure parallel_reduce_i
  end interface parallel_reduce

  interface parallel_send
     module procedure parallel_send_i1
  end interface parallel_send

  interface parallel_recv
     module procedure parallel_recv_i1
  end interface parallel_recv

  interface parallel_wait
     module procedure parallel_wait_vec_vec
  end interface parallel_wait

  interface parallel_bcast
     module procedure parallel_bcast_l
  end interface parallel_bcast

  interface parallel_gather
     !
     ! Gather fixed size blocks to specified processor
     !
     module procedure parallel_gather_i
  end interface parallel_gather

  interface parallel_allgather
     module procedure parallel_allgather_iv
  end interface parallel_allgather

contains

  subroutine parallel_initialize(comm)
    integer, intent(in), optional :: comm
    integer ierr
    logical flag
    external MPI_Init, MPI_Comm_Dup, MPI_Comm_Size, MPI_Comm_Rank
    call MPI_Initialized(flag, ierr)
    if ( .not. flag ) call MPI_Init(ierr)
    if ( present(comm) ) then
       call MPI_Comm_Dup(comm, m_comm, ierr)
    else
       call MPI_Comm_Dup(MPI_COMM_WORLD, m_comm, ierr)
    endif
    call MPI_Comm_Size(m_comm, m_nprocs, ierr)
    call MPI_Comm_Rank(m_comm, m_myproc, ierr)
    call parallel_barrier()
  end subroutine parallel_initialize
  subroutine parallel_finalize(do_finalize_MPI)
    logical, intent(in), optional :: do_finalize_MPI
    integer ierr
    external MPI_Comm_Free, MPI_Finalize
    !call MPI_Comm_Free(m_comm, ierr)  !Note: This is *supposed* to be the right way to do this, but it crashes on Linux.  comment out leads to small mem leak
    m_comm = MPI_COMM_WORLD
    if (present(do_finalize_MPI) ) then
       if (do_finalize_MPI) call MPI_Finalize(ierr)
    else
       call MPI_Finalize(ierr)
    endif
    
  end subroutine parallel_finalize

  subroutine parallel_abort(str)
    character(len=*), optional :: str
    external MPI_Abort
    integer :: ierr
    if ( parallel_IOProcessor() ) then
       if ( present(str) ) then
          print*, 'parallel_abort(): ', str
       else
          print*, 'parallel_abort() !!!'
       end if
    end if
    call MPI_Abort(m_comm, -1, ierr)
  end subroutine parallel_abort

  pure function parallel_nprocs() result(r)
    integer r
    r = m_nprocs
  end function parallel_nprocs
  pure function parallel_myproc() result(r)
    integer r
    r = m_myproc
  end function parallel_myproc
  pure function parallel_IOProcessor() result(r)
    logical :: r
    r = m_myproc == io_processor_node
  end function parallel_IOProcessor
  pure function parallel_IOProcessorNode() result(r)
    integer :: r
    r = io_processor_node
  end function parallel_IOProcessorNode
  function parallel_wtime() result(r)
    real(kind=dp_t) :: r
    r = MPI_Wtime()
  end function parallel_wtime

  ! REDUCE:
  subroutine parallel_reduce_d(r, a, op, proc, comm)
    real(kind=dp_t), intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    real(kind=dp_t), intent(out) :: r
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, 1, MPI_DOUBLE_PRECISION, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, 1, MPI_DOUBLE_PRECISION, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_d
  subroutine parallel_reduce_i(r, a, op, proc, comm)
    integer, intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc
    integer, intent(in), optional :: comm
    integer, intent(out) :: r
    integer ierr
    external MPI_Reduce, MPI_AllReduce
    integer :: l_comm
    l_comm = m_comm; if ( present(comm) ) l_comm = comm
    if ( present(proc) ) then
       CALL MPI_Reduce(a, r, 1, MPI_INTEGER, op, proc, l_comm, ierr) 
    else
       call MPI_AllReduce(a, r, 1, MPI_INTEGER, op, l_comm, ierr)
    end if
  end subroutine parallel_reduce_i

  ! RECV:
  ! blocking receive.
  subroutine parallel_recv_dv(a, n, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_dv

  ! IRECV:
  ! non-blocking receive.
  function parallel_irecv_dv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_IRecv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_IRecv(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, r, ierr)
  end function parallel_irecv_dv


  ! WAIT:
  ! completes the isend/iwait calls
  subroutine parallel_wait_vec_vec(req, status)
    integer, intent(inout)  :: req(:)
    integer, intent(out), optional :: status(:,:)
    integer :: ierr, lstatus(MPI_STATUS_SIZE,size(req))
    external MPI_WaitAll
    call MPI_WaitAll(size(req), req, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_vec_vec


  ! SEND:
  ! Blocking Send.
  subroutine parallel_send_dv(a, n, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_DOUBLE_PRECISION, proc, tag, l_comm, ierr)
  end subroutine parallel_send_dv
  subroutine parallel_send_i1(a, proc, tag, comm)
    integer, intent(in) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    call parallel_send_iv(a, size(a), proc, tag, comm)
  end subroutine parallel_send_i1
  subroutine parallel_send_iv(a, n, proc, tag, comm)
    integer, intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm
    external MPI_Send
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Send(a, n, MPI_INTEGER, proc, tag, l_comm, ierr)
  end subroutine parallel_send_iv

  ! RECV:
  ! Blocking Receive.
  subroutine parallel_recv_i1(a, proc, tag, comm, status)
    integer, intent(out) :: a(:)
    integer, intent(in) :: proc, tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    call parallel_recv_iv(a, size(a), proc, tag, comm, status)
  end subroutine parallel_recv_i1
  subroutine parallel_recv_iv(a, n, proc, tag, comm, status)
    integer, intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer ierr, l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    external MPI_Recv
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Recv(a, n, MPI_INTEGER, proc, tag, l_comm, lstatus, ierr)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_iv

  ! Broadcast
  subroutine parallel_bcast_l(a, root, comm)
    logical, intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer ierr, l_comm, l_root
    external MPI_Bcast
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    CALL MPI_Bcast(a, 1, MPI_LOGICAL, l_root, l_comm, ierr)
  end subroutine parallel_bcast_l


  ! Barrier:
  subroutine parallel_barrier(comm)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Barrier
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Barrier(l_comm, ierr)
  end subroutine parallel_barrier


  !
  ! Gather fixed size blocks to specified processor
  !
  subroutine parallel_gather_i(snd, rcv, n, root, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer :: ierr, l_root, l_comm
    external MPI_Gather
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Gather(snd, n, MPI_INTEGER, &
         rcv, n, MPI_INTEGER, &
         l_root, l_comm, ierr)
  end subroutine parallel_gather_i


  ! Allgather:
  subroutine parallel_allgather_iv(snd, rcv, n, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: ierr, l_comm
    external MPI_Allgather
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call MPI_Allgather(snd, n, MPI_INTEGER, rcv, n, MPI_INTEGER, l_comm, ierr)
  end subroutine parallel_allgather_iv

end module parallel
