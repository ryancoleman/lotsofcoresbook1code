!Copyright (c) 2014, Per Berg and Jacob Weismann Poulsen, DMI
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met: 
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer. 
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution. 
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!The views and conclusions contained in the software and documentation are those
!of the authors and should not be interpreted as representing official policies, 
!either expressed or implied, of the FreeBSD Project.

module e_timer

  !=============================================================================
  !
  ! Module containing methods for doing timer experiments. 
  ! Please, use it during testing and developement, and turn it off again 
  ! before actual production runs.
  ! 
  !
  ! Example of usage:
  ! 
  !   The first thing to do is to turn the timer stuff on and reset it. You do 
  !   this by calling timer_init() with first argument set to .TRUE. and 
  !   follow it by a call to timer_rest(). These to calls must not be from 
  !   inside a parallel section. See, the beginning of the main program,
  !   main.f90, for an example.
  ! 
  !   Then you choose the block of statements that you want to time and wrap it
  !   with calls to timer(). In the call preceeding the block, the first
  !   argument to timer() is 1. You identify and distinguish your chosen block
  !   from other timer-blocks by giving it a unique name (here 'myblock'). In 
  !   the call proceeding the block, the first argument to timer() is 2. You may
  !   print the result of the timing by calling timer_print(). Something like
  !   this:
  !
  ! 
  !     call timer(1, 'myblock')
  ! 
  !     ! ... here goes the chosen block of statements 
  !
  !     call timer(2, 'myblock')
  !     call timer_print('My chosen block took: ', 'myblock')
  ! 
  ! 
  ! Enjoy!
  !
  !=============================================================================


  !- use exitme() and task specific unit throughout ----------------------------
  use dmi_mpi_global, only : iu06
  use exits,          only : exitme
  use constants,      only : zero

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- module variables ----------------------------------------------------------
  integer(4),   private, parameter :: idlen=10
  character(2), private, parameter :: clen='10'       ! idlen as a char
  integer(4),   private, parameter :: ntimers=50
  integer(4),   private, parameter :: rc = 1          ! error return code
  integer(4),   private, save      :: nthreads
  logical,      private, save      :: timer_active
  integer(4),   private, save      :: ntim
  integer(4),   private, save      :: nclass
  character(len=idlen), allocatable, save, private :: cc(:,:) 
  real(8),              allocatable, save, private :: t1(:,:), t2(:,:)
  real(8), parameter, private      :: tmin0 = 1000000.0_8
  real(8), parameter, private      :: days = 24.0_8*3600.0_8
  real(8), parameter, private      :: hrs  = 3600.0_8, mins = 60.0_8
  real(8), parameter, private      :: secs = 1.0_8, msec = 0.001_8
!$OMP THREADPRIVATE(ntim)

  !- public interface ----------------------------------------------------------
  public :: timer, timer_init, timer_destroy, timer_print

  !- private interface ---------------------------------------------------------
  public :: timer_reset
contains

!-------------------------------------------------------------------------------

  subroutine timer_init (lflag,nt,ncl)

#if defined (_OPENMP)
    use omp_lib, only : omp_in_parallel
#endif
    implicit none
    logical, intent(in)              :: lflag
    integer(4), intent(in), optional :: nt
    integer(4), intent(in), optional :: ncl
    integer(4)                       :: s

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      call exitme(rc,'Calling timer_init() in parrallel region is not allowed')
    endif
#endif

    timer_active = lflag
    if (lflag) then
      if (.not.present(nt)) then
        nthreads=1
      else
        nthreads=nt
      endif
      if (.not.present(ncl)) then
        nclass=1
      else
        nclass=ncl
      endif
      allocate (t1(ntimers,nthreads), t2(ntimers,nthreads),                    &
                cc(ntimers,nthreads), STAT=s)
      if (s /= 0) then
        call exitme(rc,'Allocation error')
      endif
      write(iu06,*) 'Setting up timer with ', nthreads, ' threads and ',       &
                 'reporting for classes <= ', nclass
    endif
    call timer_reset()

  end subroutine timer_init

!-------------------------------------------------------------------------------

  subroutine timer_reset (nt)

#if defined (_OPENMP)
    use omp_lib, only : omp_in_parallel
#endif

    implicit none
    integer(4), intent(in), optional :: nt
    integer(4)                       :: i,j

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      call exitme(rc,'Calling timer_reset() in parrallel region is not allowed')
    endif
#endif

    ! initialize the timing
    ntim = 0
    do j=1,nthreads
      do i=1,ntimers
        t1(i,j) = zero
        t2(i,j) = zero
      enddo
    enddo
    if (present(nt)) then
      write(iu06,*) 'Resetting timer to handle ', nthreads, ' threads'
      nthreads=nt
    endif

  end subroutine timer_reset

!-------------------------------------------------------------------------------

  subroutine timer_destroy ()
    implicit none
    integer(4) :: s
   
    if (allocated (t1)) then   ! admitted, the test is sloppy - so sue me
      deallocate (t1, t2, cc, STAT=s)
      if (s /= 0) then
        call exitme(rc,'Deallocation error')
      endif
    endif
  end subroutine timer_destroy

!-------------------------------------------------------------------------------

  subroutine timer_print (string,id,class)

#if defined (_OPENMP)
    use omp_lib, only : omp_in_parallel
#endif

    implicit none
    character(len=*), intent(in), optional :: string
    character(len=*), intent(in), optional :: id
    integer(4),       intent(in), optional :: class

    real(8), parameter :: tfmtmax = 100000.0_8
    integer(4) :: i
#if defined (_OPENMP)
    integer(4) :: j
    real(8) :: t_max(ntim), t_min(ntim),t_ave(ntim)
#endif

    if (present(class)) then
      if (class > nclass) return

#if defined (_OPENMP)
!$OMP BARRIER
!$OMP MASTER
      if (omp_in_parallel()) then
        do i = 1, ntim
          t_max(i) = zero
          t_min(i) = tmin0
          t_ave(i) = zero
        enddo
        do i = 1, ntim
          do j = 1, nthreads
            if ( t2(i,j) > t_max(i) ) then
              t_max(i)=t2(i,j)
            endif
            if ( t2(i,j) < t_min(i) ) then
              t_min(i)=t2(i,j)
            endif
            t_ave(i)=t_ave(i)+t2(i,j)
          enddo
          t_ave(i)=t_ave(i)/nthreads
        enddo
        do i = 1, ntim
          write (iu06,*) 'Thread timing information (min, max, average) for ', &
                      cc(i,1), ' : ', t_min(i), t_max(i), t_ave(i)
        enddo
      else
        if (present(string) .and. present(id)) then
          do i = 1, ntim
            if (cc(i,1) == id) then
              if (t2(i,1) < tfmtmax) then
                write (iu06,'(a,a1,f13.7,a8)') string, ' ', t2(i,1), ' seconds'
              else
                write (iu06,*) string, ' ', t2(i,1), ' seconds'
              endif
              t1(i,1)=zero
              t2(i,1)=zero
            endif
          enddo
        else
          if (present(id)) then
            do i = 1, ntim
              if (cc(i,1) == id) then
                write(iu06,*) 'Timing information for ',id,' : ',cc(i,1),t2(i,1)
              endif
            enddo
          else
            do i = 1, ntim
              write(iu06,*) 'Timing information for ', cc(i,1), ' : ', t2(i,1)
            enddo
          endif
        endif
      endif
!$OMP END MASTER
#else
      if (present(string) .and. present(id)) then
        do i = 1, ntim
          if (cc(i,1) == id) then
            if (t2(i,1) < tfmtmax) then
              write (iu06,'(a,a1,f13.7,a8)') string, ' ', t2(i,1), ' seconds'
            else
              write (iu06,*) string, ' ', t2(i,1), ' seconds'
            endif
            t1(i,1)=zero
            t2(i,1)=zero
          endif
        enddo
      else
        if (present(id)) then
          do i = 1, ntim
            if (cc(i,1) == id) then
              write(iu06,*) 'Timing information for ', id, ' : ',cc(i,1),t2(i,1)
            endif
          enddo
        else
          do i = 1, ntim
            write(iu06,*) 'Timing information for ', cc(i,1), ' : ', t2(i,1)
          enddo
        endif
      endif
#endif
    endif

  end subroutine timer_print

!-------------------------------------------------------------------------------

  subroutine timer (mode,id,lrestart)

#if defined (_OPENMP)
    use omp_lib, only : omp_get_thread_num
#endif

    implicit none

    ! input/output parameters

    integer(4),       intent(in)           :: mode
    character(len=*), intent(in)           :: id
    logical,          intent(in), optional :: lrestart

    ! local parameters

    integer(4) :: ta(8)
    integer(4) :: i,tid
    real(8)    :: mtime
    logical    :: llrestart

    if (.not. timer_active) return

    if (len(id) > idlen) then
      call exitme(rc,'size of '//trim(id)//' exceeds the allowed size '//      &
                     clen//' - bailing out')
    endif

    if (.not.present(lrestart)) then
      llrestart=.true. ! default is to restart the timer
    else
      llrestart=lrestart 
    endif

#if defined (_OPENMP)
    tid=omp_get_thread_num()+1
#else
    tid=1
#endif

    select case (mode)

    !---------------------------------------------------------------------------
    !    mode = 1: start timing of 'id'
    !           2: stop timing  of 'id'
    !---------------------------------------------------------------------------

    case (1)
      ! start timing of 'id'
      call date_and_time(values=ta)
      mtime = ta(3)*(days) + & ! number of days
              ta(5)*(hrs)  + & ! number of hours
              ta(6)*(mins) + & ! number of minutes
              ta(7)*(secs) + & ! number of seconds
              ta(8)*(msec)     ! number of milliseconds
      do i = 1, ntim
        if (cc(i,tid) == id) then
          if (llrestart) then
            ! restarting timer
            t1(i,tid) = mtime
          endif
          return
        end if
      end do

      if (ntim+1 > ntimers) then
        call exitme(rc,'Trying to add more timers than we have reserved '//    &
                       'space for, please adjust ntimers')
      endif

      ntim         = ntim + 1
      t1(ntim,tid) = mtime
      t2(ntim,tid) = zero
      cc(ntim,tid) = id

    case (2)
      ! stop timing  of 'id'
      call date_and_time(values=ta)
      mtime = ta(3)*(days) + & ! number of days
              ta(5)*(hrs)  + & ! number of hours
              ta(6)*(mins) + & ! number of minutes
              ta(7)*(secs) + & ! number of seconds
              ta(8)*(msec)     ! number of milliseconds
      do i = 1, ntim
        if (cc(i,tid) == id) then
          t2(i,tid) = t2(i,tid) + mtime - t1(i,tid)
          return
        end if
      end do
      if (i == ntim+1) then 
        call exitme(rc,'Trying to stop timer '//trim(id)//' which has not '//  &
                       'been started - bailing out')
      endif

   end select

   end subroutine timer

end module e_timer


