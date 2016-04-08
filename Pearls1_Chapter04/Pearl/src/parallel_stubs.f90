! Parallel wrappers
!! These wrappers are used so that a non-MPI version can coexist
!! with the MPI version.

module parallel

  use bl_types

  implicit none


  ! Some selected values based on MPICH/1.2.5 unix implementation

  integer, parameter :: MPI_UNDEFINED   = -32766
  integer, parameter :: MPI_STATUS_SIZE = 4
  integer, parameter :: MPI_SOURCE      = 2
  integer, parameter :: MPI_TAG         = 3
  integer, parameter :: MPI_ERROR       = 4

  integer, parameter :: MPI_MAX       = 100
  integer, parameter :: MPI_MIN       = 101
  integer, parameter :: MPI_SUM       = 102
  integer, parameter :: MPI_PROD      = 103
  integer, parameter :: MPI_LAND      = 104
  integer, parameter :: MPI_BAND      = 105
  integer, parameter :: MPI_LOR       = 106
  integer, parameter :: MPI_BOR       = 107
  integer, parameter :: MPI_LXOR      = 108
  integer, parameter :: MPI_BXOR      = 109
  integer, parameter :: MPI_MINLOC    = 110
  integer, parameter :: MPI_MAXLOC    = 111

  integer, parameter :: MPI_ANY_SOURCE = -2
  integer, parameter :: MPI_ANY_TAG    = -1

  integer, parameter, private :: io_processor_node = 0
  integer, private :: m_nprocs = 1
  integer, private :: m_myproc = 0
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

  interface
     subroutine sys_abort()
     end subroutine sys_abort
  end interface

  logical, private :: g_init = .False.

contains

  subroutine parallel_initialize(comm)
    integer, intent(in), optional :: comm
    g_init = .True.
  end subroutine parallel_initialize

  subroutine parallel_finalize(do_finalize_MPI)
    logical, intent(in), optional :: do_finalize_MPI
    g_init = .False.
  end subroutine parallel_finalize

  subroutine parallel_abort(str)
    character(len=*), optional :: str
    if ( present(str) ) then
       print*, 'parallel_abort(): ', str
    else
       print*, 'parallel_abort() !!!'
    end if
    call sys_abort()
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
    r = parallel_myproc() == io_processor_node
  end function parallel_IOProcessor
  pure function parallel_IOProcessorNode() result(r)
    integer :: r
    r = io_processor_node
  end function parallel_IOProcessorNode
  function parallel_wtime() result(r)
    real(kind=dp_t) :: r
    interface
       subroutine wall_second(s)
         use bl_types
         real(kind=dp_t), intent(out) :: s
       end subroutine wall_second
    end interface
    call wall_second(r)
  end function parallel_wtime

  subroutine parallel_reduce_d(r, a, op, proc, comm)
    real(kind=dp_t), intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    real(kind=dp_t), intent(out) :: r
    r = a
  end subroutine parallel_reduce_d
  subroutine parallel_reduce_i(r, a, op, proc, comm)
    integer, intent(in) :: a
    integer, intent(in) :: op
    integer, intent(in), optional :: proc, comm
    integer, intent(out) :: r
    r = a
  end subroutine parallel_reduce_i


  subroutine parallel_recv_dv(a, n, proc, tag, comm, status)
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = HUGE(a)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_dv

  function parallel_irecv_dv(a, n, proc, tag, comm) result(r)
    integer :: r
    real(kind=dp_t), intent(out) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_ISEND_RV')
    a(1:n) = HUGE(a)
    r = -1
  end function parallel_irecv_dv

  subroutine parallel_wait_vec_vec(req, status)
    integer, intent(in)  :: req(:)
    integer, intent(out), optional :: status(:,:)
    integer :: lstatus(MPI_STATUS_SIZE,size(req))
    lstatus = 0
    if ( size(req) > 0 ) call parallel_abort('PARALLEL_WAIT')
    if ( present(status) ) status = lstatus
  end subroutine parallel_wait_vec_vec


  subroutine parallel_send_dv(a, n, proc, tag, comm)
    real(kind=dp_t), intent(in) :: a(*)
    integer, intent(in) :: proc, tag, n
    integer, intent(in), optional :: comm
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
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
    integer l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_SEND')
  end subroutine parallel_send_iv

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
    integer l_comm, lstatus(MPI_STATUS_SIZE)
    integer, intent(out), optional, dimension(MPI_STATUS_SIZE) :: status
    lstatus = 0
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    call parallel_abort('PARALLEL_RECV')
    a(1:n) = HUGE(a)
    if ( present(status) ) status = lstatus
  end subroutine parallel_recv_iv


  subroutine parallel_bcast_l(a, root, comm)
    logical, intent(in) :: a
    integer, intent(in), optional :: root
    integer, intent(in), optional :: comm
    integer l_comm, l_root
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
  end subroutine parallel_bcast_l


  subroutine parallel_barrier(comm)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
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
    integer :: l_root, l_comm
    l_root = io_processor_node
    if ( present(root) ) l_root = root
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_gather_i


  subroutine parallel_allgather_iv(snd, rcv, n, comm)
    integer, intent(in) :: n
    integer, intent(in) :: snd(*)
    integer, intent(out) :: rcv(*)
    integer, intent(in), optional :: comm
    integer :: l_comm
    l_comm = m_comm
    if ( present(comm) ) l_comm = comm
    rcv(1:n) = snd(1:n)
  end subroutine parallel_allgather_iv

end module parallel
