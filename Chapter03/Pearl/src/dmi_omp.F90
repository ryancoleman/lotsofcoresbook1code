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

module dmi_omp
  use constants, only : zero, one

  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- interfaces ----------------------------------------------------------------
  interface domp_get_domain
    module procedure domp_get_domain_rlu
    module procedure domp_get_domain_wet
    module procedure domp_get_domain_wet_idx 
    module procedure domp_get_domain_wet_idx_stream
    module procedure domp_get_domain_wet_halo_idx 
  end interface !F90 std domp_get_domain

  !- private vars & methods ----------------------------------------------------
  ! Please note, this constant will create a compiler info for a constant
  ! expression in IF statements:
  logical(4), private, save :: omp_debug=.true.
  integer(4), private :: domp_iam, domp_nt
  integer(4), public :: domp_ns
  integer(4), allocatable, save, private :: tsize2d(:), tsize3d(:)
  integer(4), allocatable, save, private :: hsize2d(:), hsize3d(:)
  integer(4), allocatable, save, private :: lubounds(:,:)
  private :: domp_get_domain_rlu, domp_get_domain_wet, domp_get_domain_wet_idx
  private :: domp_get_domain_wet_halo_idx, domp_get_domain_wet_idx_stream
  private :: domp_get_domain_exact_idx, reconstruct_solution

  !- public vars & methods -----------------------------------------------------
  public  :: domp_init, domp_get_domain, domp_get_thread_no,                   &
             domp_estimate_weighted_load, domp_load_decompo, domp_debug_info

#if defined (_OPENMP)
  real(8), private :: rdomp_iam, rdomp_nt
!$OMP THREADPRIVATE(domp_iam,domp_nt,rdomp_iam,rdomp_nt) 
#endif

contains
    
  ! ----------------------------------------------------------------------------

  subroutine domp_init(nt_out)

    use dmi_mpi_global, only : iu06

#if defined (_OPENMP)
    use omp_lib, only : omp_get_thread_num, omp_get_num_threads
#endif

    !- argument(s) -------------------------------------------------------------
    integer(4), intent(out) :: nt_out

!$OMP PARALLEL DEFAULT(none)
#if defined (_OPENMP)
    domp_iam  = omp_get_thread_num()
    rdomp_iam = real(domp_iam,8)
    domp_nt   = omp_get_num_threads()
    rdomp_nt  = real(domp_nt,8)
#else
    domp_iam  = 0
    domp_nt   = 1
    domp_ns   = 0
#endif
!$OMP END PARALLEL

#if defined (_OPENACC)
    domp_ns=16 ! FIXME this should be a tunable parameter that we read in 
             ! from a namelist
    write(iu06,'(a27)') 'Build with openACC support'
!#elif defined (_OPENMP)
!    write(iu06,'(a26)') 'Build with openMP support'
!#else
!    write(iu06,'(a41)') 'Build without openMP and openACC support'
#endif

#if defined (_OPENACC)
    write(iu06,'(a20,i5,a8)') 'Running openACC with ', domp_ns, ' streams'
#endif

    !- echo #threads:
    if (domp_nt > 1) then
      write(iu06,'(a20,i5,a8)') 'Running openMP with ', domp_nt, ' threads'
    else
#if defined (_OPENMP)
      write(iu06,'(a35)') 'Running openMP with a single thread'
#else
      write(iu06,'(a22)') 'Running without openMP'
#endif
    endif

    !- return value of #threads:
    nt_out = domp_nt
    if (omp_debug) then
      allocate(tsize2d(1:domp_nt),tsize3d(1:domp_nt), lubounds(1:domp_nt,1:2), &
               hsize2d(1:domp_nt),hsize3d(1:domp_nt))
    endif

  end subroutine domp_init
 
  ! ----------------------------------------------------------------------------

  subroutine domp_load_decompo(idx,ia,narea)
    use io_subs,        only : io_new_unit, flush_unit
    use dmi_mpi_global, only : iu06
    use exits,          only : exitme
    implicit none
    integer(4), intent(in)   :: ia,narea
    integer(4), intent(out)  :: idx(1:,1:)
    integer(4)               :: lun, ios, nt, it, iit, itl, itu, iia, ial
    character(LEN=18), parameter  :: domp_decomp_file = 'openmp_decompo.txt'
    idx(:,:) = -1
    lun = io_new_unit()
    open (unit=lun, file=trim(domp_decomp_file), status='old', iostat=ios)
    if (ios /= 0) call exitme(1,'Cannot open '//trim(domp_decomp_file))
    read(lun, '(i10)', iostat=ios) nt
    if (ios /= 0) then
      call exitme(1,'Cannot read number of threads from '                      &
                //trim(domp_decomp_file))
    endif
    if (nt /= domp_nt) call exitme(1,'Invalid number of threads specified')
    naloop: do ial=1,narea
      do it=1,nt
        read(lun, '(4i10)', iostat=ios) iia, iit, itl, itu
        write (*,*) it, ia, iia, iit, itl, itu
        if (ios /= 0) then
          call exitme(1,'Cannot read grid extent data from '                   &
                               //trim(domp_decomp_file))
        endif
        if (ia == iia) then
          idx(1,iit)=itl
          idx(2,iit)=itu
        endif
      enddo 
      if (iia==ia) exit naloop
    enddo naloop
    close(lun)
    if (omp_debug) then
      write(iu06,'(a14,i3)') 'Subdomain ', ia
      do it=1,nt
        write(iu06,'(a14,i3,a24,i10,i10)') 'openMP thread ', it,               &
         ' handles surface range: ', idx(1,it), idx(2,it)
      enddo
      call flush_unit(iu06)
    endif
  end subroutine domp_load_decompo

  subroutine domp_debug_info(idx,kh)
    use dmi_mpi_global, only : iu06
    use io_subs,        only : flush_unit
    implicit none
    integer(4), intent(in)   :: idx(1:,1:)
    integer(4), intent(in)   :: kh(0:)
    integer(4)               :: csize,it,n
    do it=1,domp_nt
      csize=0
      do n=idx(1,it), idx(2,it)
        csize = csize + kh(n)
      enddo
      write(iu06,'(a14,i3,a24,i10,i10,i10)') 'openMP thread ', it,             &
        ' handles surface range: ', idx(1,it), idx(2,it), csize
    enddo
    call flush_unit(iu06)
  end subroutine domp_debug_info

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_rlu(lower,upper,d_lower,d_upper,reverse)

#if defined (_OPENMP)
    use omp_lib,   only : omp_in_parallel
    use constants, only : half
#endif

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)  :: lower,upper
    integer(4), intent(out) :: d_lower,d_upper
    logical,    intent(in), optional :: reverse

#if defined (_OPENMP)
    !  local variables ---------------------------------------------------------
    real(8)    :: dlen
    integer(4) :: lr, ur
#endif

    ! proper action in "null" cases:
    if (upper <= 0 .or. upper < lower) then
      d_lower = 0
      d_upper = -1
      return
    endif

    ! proper action in serial sections
    d_lower = lower
    d_upper = upper

#if defined (_OPENMP)
    if (omp_in_parallel()) then
      dlen    = real(upper-lower+1, 8)
      d_lower = lower    + floor((rdomp_iam*dlen+half)/rdomp_nt)   !F90 std, 4)
      d_upper = lower -1 + floor((rdomp_iam*dlen+dlen+half)/rdomp_nt)
      if (present(reverse)) then
        if (reverse) then
          ! reverse order:
          lr = upper + lower - d_upper
          ur = upper + lower - d_lower
          d_lower = lr
          d_upper = ur
        endif
      endif
    endif
#endif

  end subroutine domp_get_domain_rlu

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_wet(Nwet, kh, lower, upper, d_lower, d_upper)
    ! attempt to use an equal number of wet points for each thread

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)  :: Nwet,lower,upper
    integer(4), intent(in)  :: kh(0:)
    integer(4), intent(out) :: d_lower,d_upper

    !  local variables ---------------------------------------------------------
    integer(4) :: chunk, csize, tnum, n


    !  check if there is nothing to do:
    if (Nwet <= 0) then
      d_lower = 0
      d_upper = -1
      return
    endif

    !  OK, we should if lower/upper loop limits:
    chunk   = Nwet/domp_nt
    d_lower = lower
    csize   = 0
    tnum    = 1

    n_loop: do n=lower,upper
      csize = csize + kh(n)
      if (csize >= chunk .or. n == upper) then
        ! this chunk is large enough; set up limits

        d_upper = min(n,upper)

        if (tnum == domp_iam + 1) then
          ! got what I came for
          exit n_loop
        elseif (n < upper) then
          ! try next thread number
          tnum    = tnum + 1
          d_lower = n+1
        else
          ! oh no, I'm trapped !!! gotta do something
          d_lower = d_upper + 1
          exit n_loop
        endif

        ! start next chunk
        csize = 0

      elseif (tnum == domp_nt) then
        ! no need to proceed; no more threads anyway so we 
        ! put the rest into the last thread
        d_upper = upper
        exit n_loop
      endif
    enddo n_loop

  end subroutine domp_get_domain_wet

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_wet_idx(kh, lower, upper, d_lower, d_upper, idx)
    ! attempt to use an equal number of wet points for each thread 
    ! storing lower/upper index pair for easy look-up during subsequent calls

    !- modules -----------------------------------------------------------------
    use constants,      only : ompsplit_exact, ompsplit_op
#if defined (_OPENMP)
    use dmi_mpi_global, only : iu06
    use io_subs,        only : flush_unit
#endif

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: lower,upper
    integer(4), intent(in)    :: kh(0:)
    integer(4), intent(out)   :: d_lower,d_upper
    integer(4), intent(inout) :: idx(1:,1:)

    !  local variables ---------------------------------------------------------
    integer(4) :: chunk, csize, tnum, n, nwet, csizenext, penalty_scal, csize_n
    logical    :: got_trapped, better_than_next
    real(8)    :: rounding_offset, op
#if defined (_OPENMP)
    integer(4) :: icount
    real(8)    :: wmean
#endif

    ! check if idx has already been set and terminate --------------------------
    tnum = domp_iam + 1
    if (idx(1,tnum) >= 0) then
      d_lower = idx(1,tnum)
      d_upper = idx(2,tnum)
      return
    endif

    if (ompsplit_op) then
      penalty_scal = 1
    else
      penalty_scal = 0
    endif


    ! calc No of wet points ----------------------------------------------------
    nwet = 0
    do n=lower,upper
      nwet = nwet + kh(n)
    enddo
    chunk = nwet/domp_nt

    ! Treat trivial case -------------------------------------------------------
    if (nwet <= 0) then
      ! ooppsss ..., nothing to do but ok
      d_lower = 0
      d_upper = -1
      ! store d_lower and d_upper for subsequent entries:
    elseif (ompsplit_exact) then
      !- split domain with exhausive search (exact solution)
      call domp_get_domain_exact_idx(kh, lower, upper, d_lower, d_upper, csize)
    else
      ! get domains of nearly equal chunk sizes --------------------------------
      rounding_offset = zero
      d_lower = lower
      csize   = 0
      tnum    = 1

      got_trapped = .false.
      n_loop: do n=lower,upper
        csize = csize + kh(n)
        if (n < upper .and. csize < chunk) then
          csize_n   = kh(n+1)
          csizenext = csize + csize_n
         !if (.not.(csizenext > chunk .and. csizenext-chunk > chunk-csize)) then

          !- if we made a too small domain last time, op should make iteasier to
          !  make a larger split this time. We do this so we don't end up with 
          !  a very small or very large last domain
          op               = rounding_offset*penalty_scal
          better_than_next = csizenext-chunk > chunk-csize+op
          if (.not.(csizenext > chunk .and. better_than_next)) then
            csizenext = 0
            csize_n   = 0
          endif
        else
          csizenext = 0
          csize_n   = 0
        endif
        if (csize >= chunk .or. n == upper .or. csizenext > 0) then
        ! this chunk is large enough; set up limits

          d_upper = min(n,upper)

          if (tnum == domp_iam + 1) then
            ! got what I came for
            if (tnum == domp_nt) d_upper = upper
            exit n_loop
          elseif (n < upper) then
            ! try next thread number
            rounding_offset = rounding_offset + chunk - csize - csize_n
            tnum    = tnum + 1
            d_lower = n+1
          else
            ! oh no, I'm trapped !!! gotta do something
            got_trapped = .true.
            d_lower = d_upper + 1
            exit n_loop
          endif

          ! start next chunk
          csize = 0
  
        elseif (tnum == domp_nt) then
          ! no need to proceed; no more threads anyway so we 
          ! put the rest into the last thread
          d_upper = upper
          exit n_loop
        endif
      enddo n_loop
    endif ! trivail nwet < 0 or domp_exact true or false

    ! store d_lower and d_upper to avoid doing the above again ---------------
    tnum = domp_iam + 1
    idx(1,tnum) = d_lower
    idx(2,tnum) = d_upper
    if (omp_debug) then
      if (got_trapped) then
        tsize3d(tnum) = 0
      else
        tsize3d(tnum) = csize
      endif
      tsize2d(tnum)    = d_upper-d_lower
      lubounds(tnum,1) = d_lower
      lubounds(tnum,2) = d_upper
    endif
#if defined (_OPENMP)
    if (omp_debug) then
!$OMP BARRIER
!$OMP MASTER
      ! compute the total number of 3D wetpoints for the last thread
      csize = 0
      do n=lubounds(domp_nt,1), lubounds(domp_nt,2)
        csize = csize + kh(n)
      enddo
      tsize3d(domp_nt)=csize
      icount = count(tsize3d > zero)
      if (icount > 0) then
        wmean =  sum(tsize3d, tsize3d > zero) / real(icount,8)
      else
        wmean = zero
      endif
      write(iu06,'(a64)')                                                      &
        'OpenMP 3D decomposition statistics (ideal mean, mean, min, max):'
      write(iu06,'(i10,f10.0,i13,a2,i3,a2,i13,a2,i3,a2)')                      &
        chunk, wmean,                                                          &
        minval(tsize3d), '(', minloc(tsize3d), ')',                            &
        maxval(tsize3d), '(', maxloc(tsize3d), ')' 
      icount = count(tsize2d > zero)
      if (icount > 0) then
        wmean =  sum(tsize2d, tsize2d > zero) / real(icount,8)
      else
        wmean = zero
      endif
      write(iu06,'(a64)')                                                      &
        'OpenMP 2D decomposition statistics (ideal mean, mean, min, max):'
      write(iu06,'(2f10.0,i13,a2,i3,a2,i13,a2,i3,a2)')                         &
        real((upper-lower)/domp_nt,8), wmean,                                  &
        minval(tsize2d), '(', minloc(tsize2d), ')',                            &
        maxval(tsize2d), '(', maxloc(tsize2d), ')'
      do tnum=1,domp_nt
        write(iu06,'(a14,i3,a24,i10,i10)') 'openMP thread ', tnum,             &
         ' handles surface range: ', lubounds(tnum,1), lubounds(tnum,2)
      enddo
      call flush_unit(iu06)
!$OMP END MASTER
!$OMP BARRIER
    endif
#endif
  end subroutine domp_get_domain_wet_idx

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_exact_idx(kh, lower, upper, d_lower, d_upper,     &
                                       my_cost)
    implicit none

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: lower, upper
    integer(4), intent(in)    :: kh(0:)
    integer(4), intent(out)   :: d_lower, d_upper, my_cost

    !- locals -----------------------------------------------------------------
    integer(4) :: n2d, iw, threads, x, s, tnum
    integer(4), allocatable :: psum(:), cost(:,:), solution(:,:), d_allupper(:)

    n2d = upper - lower + 1

    allocate( psum(n2d), cost(n2d,domp_nt), solution(n2d, domp_nt),            &
              d_allupper(domp_nt) )

    ! compute prefix sums
    psum(1) = kh(1)
    do iw=2,n2d
      psum(iw) = psum(iw-1) + kh(iw)
    enddo

    !=======
    ! BC
    !=======

    ! boundary condition 1
    do iw=1,n2d
      !- the smallest cost of the first partition is just the total iw3 
      cost(iw,1) = psum(iw)
    enddo

    !- boundary condition 2  
    !  the smallest cost of splitting i=1 is always idist(1) since we cannot
    !  divide it. 
    cost(1,1:domp_nt) = kh(1)

    !=======
    ! Do the exhausive search
    !=======
    do iw=2,n2d
      do threads=2,domp_nt
        x = 1

        !- unroll this from loop
        cost(iw,threads)     = max( cost(x,threads-1), psum(iw)-psum(x) )
        solution(iw,threads) = x

        do x=2,iw-1
          
          s = max( cost(x,threads-1), psum(iw)-psum(x) )

          !- NOTE: if we chage '>=' to '>' we might get empty tasks 'cause
          !- then it stops splitting the the tasks (i+1...N) if the best
          !- possible solution has been found with largest subdomain at i. 
          if (cost(iw,threads) >= s ) then
            cost(iw,threads)     = s
            solution(iw,threads) = x
          endif
        enddo
      enddo
    enddo
 
    !=======
    ! Reconstruct partitions from the solution array
    !
    ! Explanation:
    ! D(mmx,ni) = where sould the (ni-1)'s divider be placed if we have mmx
    ! ilines which should be divided into ni subdomains
    !
    ! We solve this problem with help of a recursive routine:
    ! D(mmx,ni) = z  says where the (ni-1)'s divider should be placed. Use this
    ! info to get the placement of the (ni-2)'s divider since we then now the 
    ! size of this problem [D(z,ni-1)]
    
    !
    ! I.E.: to get the placement of (ni-3)'s divider D( D( D(...D(iw,1)...) ) )
    ! 
    !=======
    d_allupper(:) = 0
    call reconstruct_solution(solution, n2d, domp_nt, d_allupper)

    ! calculate upper and lower for this thread

    !- NOTE if the thread is dry, upper_d = 0 and lower_d = 1 
    tnum = domp_iam + 1
    d_upper = d_allupper(tnum)
    
    if (tnum > 1) then
      d_lower = d_allupper(tnum-1)+1 ! calculated from prev threads' upper lim
    else
      d_lower = 1
    endif

    !- calculate cost for this thread
    my_cost = 0
    if (omp_debug) then
      do iw = d_lower, d_upper
        my_cost = my_cost + kh(iw)
      enddo
    endif

    !- clean up
    deallocate( psum, cost, solution, d_allupper )
  end subroutine domp_get_domain_exact_idx

  ! ----------------------------------------------------------------------------

  recursive subroutine reconstruct_solution(solution, n, k, upper)
    implicit none

    !- arguments
    integer(4), intent(in)    :: solution(:,:), n, k
    integer(4), intent(inout) :: upper(:)

    if (k==1 .or. n == 1) then
      upper(k) = n
    else
      call reconstruct_solution(solution, solution(n,k), k-1, upper)
      upper(k) = n
    endif

  end subroutine reconstruct_solution

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_wet_idx_stream(kh, lower, upper, d_lower, d_upper,&
                                            idx, iam)
#if defined (_OPENACC) 
    use dmi_mpi_global, only : iu06
#endif
    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: lower,upper,iam
    integer(4), intent(in)    :: kh(0:)
    integer(4), intent(out)   :: d_lower,d_upper
    integer(4), intent(inout) :: idx(1:,1:)

#if defined (_OPENACC) 
    !  local variables ---------------------------------------------------------
    integer(4) :: chunk, csize, tnum, n, nwet, csizenext
    logical    :: got_trapped

    ! check if idx has already been set and terminate --------------------------
    if (idx(1,iam) >= 0) then
      d_lower = idx(1,iam)
      d_upper = idx(2,iam)
      return
    endif

    ! calc No of wet points ----------------------------------------------------
    nwet = 0
    do n=lower,upper
      nwet = nwet + kh(n)
    enddo

    ! Treat trivial case -------------------------------------------------------
    if (nwet <= 0) then
      ! ooppsss ..., nothing to do but ok
      d_lower = 0
      d_upper = -1
      ! store d_lower and d_upper for subsequent entries:
      idx(1,iam) = d_lower
      idx(2,iam) = d_upper
    else
      ! get domains of nearly equal chunk sizes --------------------------------
      chunk   = nwet/domp_ns
      d_lower = lower
      csize   = 0
      tnum    = 1

      got_trapped = .false.
      n_loop: do n=lower,upper
        csize = csize + kh(n)
        if (n < upper .and. csize < chunk) then
          csizenext = csize + kh(n+1)
          if (.not.(csizenext > chunk .and. csizenext-chunk > chunk-csize)) then
            csizenext = 0
          endif
        else
          csizenext = 0
        endif
        if (csize >= chunk .or. n == upper .or. csizenext > 0) then
        ! this chunk is large enough; set up limits

          d_upper = min(n,upper)

          if (tnum == iam) then
            ! got what I came for
            if (tnum == domp_ns) d_upper = upper
            exit n_loop
          elseif (n < upper) then
            ! try next thread number
            tnum    = tnum + 1
            d_lower = n+1
          else
            ! oh no, I'm trapped !!! gotta do something
            got_trapped = .true.
            d_lower = d_upper + 1
            exit n_loop
          endif

          ! start next chunk
          csize = 0
  
        elseif (tnum == domp_ns) then
          ! no need to proceed; no more threads anyway so we 
          ! put the rest into the last thread
          d_upper = upper
          exit n_loop
        endif
      enddo n_loop

      ! store d_lower and d_upper to avoid doing the above again ---------------
      idx(1,iam) = d_lower
      idx(2,iam) = d_upper
    endif

    write(iu06,'(a14,i3,a24,i10,i10)') 'openACC stream ', iam,                 &
       ' handles surface range: ', idx(1,iam), idx(2,iam)
#else
    call domp_get_domain(kh, lower, upper, d_lower, d_upper, idx)
#endif
  end subroutine domp_get_domain_wet_idx_stream

  ! ----------------------------------------------------------------------------

  subroutine domp_get_domain_wet_halo_idx(kh, lower, upper, d_lower, d_upper,  &
                                          mmk, ind, a0, a1, a2, a3, a4, idx)
    ! Attempt to use an equal number of a weighted sum of wet points and halo 
    ! points for each thread
    ! storing lower/upper index pair for easy look-up during subsequent calls
    !
    ! First, you should call domp_get_domain_wet_idx() to get idx for the even
    ! split of wetpoints, then call this domp_get_domain_wet_halo_idx() a
    ! number of times to get each iteration 
    !
    ! If coef a1 is zero, we can reuse the previous idx and no new decompo is
    ! calculated.

    !- modules -----------------------------------------------------------------
#if defined (_OPENMP)
    use dmi_mpi_global, only : iu06
    use io_subs,        only : flush_unit
#endif

    !- arguments ---------------------------------------------------------------
    integer(4), intent(in)    :: lower,upper
    integer(4), intent(in)    :: kh(0:), mmk(1:,0:,0:), ind(1:,1:)
    integer(4), intent(out)   :: d_lower,d_upper
    integer(4), intent(inout) :: idx(1:,1:)
    real(8),    intent(in)    :: a0, a1, a2, a3, a4

    !  local variables ---------------------------------------------------------
    integer(4) :: tnum, n, nwet, hsize, thrdsize, halosize, i, j, it
    integer(4) :: nh, ih, jh, hsize2, hsiz2
    real(8)    :: csize, csizenext, chunk, d1, d2
    logical    :: got_trapped
#if defined (_OPENMP)
    integer(4) :: icount
    real(8)    :: wmean
#endif

    ! check if idx has already been set and terminate --------------------------
    tnum = domp_iam + 1
    if (idx(1,tnum) >= 0 .and. a1 == zero) then
      d_lower = idx(1,tnum)
      d_upper = idx(2,tnum)
      return
    endif

    ! calc No of wet points ----------------------------------------------------
    nwet  = 0
    do n=lower,upper
      nwet = nwet + kh(n)
    enddo

    ! calc hsize from previous decompo -----------------------------------------
    hsize = 0
    hsiz2 = 0
    do it=1,domp_nt
      halosize = 0
      hsize2   = 0
      do n=idx(1,it),idx(2,it)
        i = ind(1,n)
        j = ind(2,n)
        if (mmk(1,i+1,j) > idx(2,it)) then
          halosize = halosize + kh(mmk(1,i+1,j))
          hsize2   = hsize2 + 1
        endif
        if (mmk(1,i,j+1) > idx(2,it)) then
          halosize = halosize + kh(mmk(1,i,j+1))
          hsize2   = hsize2 + 1
        endif
        if (mmk(1,i-1,j) < idx(1,it)) then
          halosize = halosize + kh(mmk(1,i-1,j))
          hsize2   = hsize2 + 1
        endif
        if (mmk(1,i,j-1) < idx(1,it)) then
          halosize = halosize + kh(mmk(1,i,j-1))
          hsize2   = hsize2 + 1
        endif
      enddo
      hsize = hsize + halosize
      hsiz2 = hsiz2 + hsize2
    enddo

    ! Print sizes --------------------------------------------------------------
#if defined (_OPENMP)
    if (omp_debug) then
!$OMP BARRIER
!$OMP MASTER
      write(iu06,'(a23,i10)') '3D wet points:         ', nwet
      write(iu06,'(a23,i10)') '2D halo size at entry: ', hsiz2
      write(iu06,'(a23,i10)') '3D halo size at entry: ', hsize
!$OMP END MASTER
!$OMP BARRIER
    endif
#endif

    ! Treat trivial case -------------------------------------------------------
    if (nwet <= 0) then
      ! ooppsss ..., nothing to do but ok
      d_lower = 0
      d_upper = -1
      ! store d_lower and d_upper for subsequent entries:
      tnum = domp_iam + 1
      idx(1,tnum) = d_lower
      idx(2,tnum) = d_upper
    else
      ! get domains of nearly equal chunk sizes --------------------------------
      chunk    = a0 + (  a1*real(nwet,8)  + a2*real(upper-lower+1,8)           &
                       + a3*real(hsiz2,8) + a4*real(hsize,8) )/real(domp_nt,8)
      d_lower  = lower
      thrdsize = 0
      tnum     = 1

      got_trapped = .false.
      n_loop: do n=lower,upper
        i = ind(1,n)
        j = ind(2,n)

        ! running size of thread:
        thrdsize = thrdsize + kh(n)

        ! running size of halo:
        halosize = 0
        hsize2   = 0
        do nh=d_lower,n
          ih = ind(1,nh)
          jh = ind(2,nh)
          if (mmk(1,ih+1,jh) > n)       then
            halosize = halosize + kh(mmk(1,ih+1,jh))
            hsize2   = hsize2 + 1
          endif
          if (mmk(1,ih,jh+1) > n)       then
            halosize = halosize + kh(mmk(1,ih,jh+1))
            hsize2   = hsize2 + 1
          endif
          if (mmk(1,ih-1,jh) < d_lower) then
            halosize = halosize + kh(mmk(1,ih-1,jh))
            hsize2   = hsize2 + 1
          endif
          if (mmk(1,ih,jh-1) < d_lower) then
            halosize = halosize + kh(mmk(1,ih,jh-1))
            hsize2   = hsize2 + 1
          endif
        enddo

        ! running size of chunk and discrepnacy:
        csize = a0 + a1*real(thrdsize,8) + a2*real(n-d_lower+1,8)              &
                   + a3*real(hsize2,8)   + a4*real(halosize,8)
        d1    = abs(csize - chunk)

        ! running size of chunk, including next point:
        if (n < upper) then
          csizenext = csize + a1*real(kh(n+1),8) + a2
          ih = ind(1,n+1)
          jh = ind(2,n+1)
          if (mmk(1,ih+1,jh) > n+1)     then
            csizenext = csizenext + a3 + a4*real(kh(mmk(1,ih+1,jh)),8)
          endif
          if (mmk(1,ih,jh+1) > n+1)     then
            csizenext = csizenext + a3 + a4*real(kh(mmk(1,ih,jh+1)),8)
          endif
          if (mmk(1,ih-1,jh) < d_lower) then
            csizenext = csizenext + a3 + a4*real(kh(mmk(1,ih-1,jh)),8)
          endif
          if (mmk(1,ih,jh-1) < d_lower) then
            csizenext = csizenext + a3 + a4*real(kh(mmk(1,ih,jh-1)),8)
          endif
          d2 = abs(csizenext - chunk)
!         if (domp_nt > 1) then
!           d2 = d2*((0.08_8/real(tnum,8)+0.92_8)*real(domp_nt,8)-one)         &
!                  /(real(domp_nt,8)-one)
!         endif
          d2 = d2*0.92_8
        else
          d2 = d1
        endif

        if (d2 > d1 .or. n == upper) then
        ! this chunk is large enough; set up limits

          d_upper = min(n,upper)

          if (tnum == domp_iam + 1) then
            ! got what I came for
            if (tnum == domp_nt) d_upper = upper
            exit n_loop
          elseif (n < upper) then
            ! try next thread number
            tnum    = tnum + 1
            d_lower = n+1
          else
            ! oh no, I'm trapped !!! gotta do something
            got_trapped = .true.
            d_lower = d_upper + 1
            exit n_loop
          endif

          ! start next chunk
          thrdsize = 0
  
        elseif (tnum == domp_nt) then
          ! no need to proceed; no more threads anyway so we 
          ! put the rest into the last thread
          d_upper = upper
          exit n_loop
        endif
      enddo n_loop

      ! store d_lower and d_upper to avoid doing the above again ---------------
      tnum = domp_iam + 1
      idx(1,tnum) = d_lower
      idx(2,tnum) = d_upper
      if (omp_debug) then
        if (got_trapped) then
          tsize3d(tnum) = 0
          hsize2d(tnum) = 0
          hsize3d(tnum) = 0
        else
          tsize3d(tnum) = thrdsize
          hsize2d(tnum) = hsize2
          hsize3d(tnum) = halosize
        endif
        tsize2d(tnum)    = d_upper-d_lower+1
        lubounds(tnum,1) = d_lower
        lubounds(tnum,2) = d_upper
      endif
    endif
#if defined (_OPENMP)
    if (omp_debug) then
!$OMP BARRIER
!$OMP MASTER
      ! compute the sizes for the last thread:
      thrdsize = 0
      halosize = 0
      hsize2   = 0
      do n=lubounds(domp_nt,1), lubounds(domp_nt,2)
        i = ind(1,n)
        j = ind(2,n)
        thrdsize = thrdsize + kh(n)
        if ( mmk(1,i+1,j) > lubounds(domp_nt,2) ) then
          halosize = halosize + kh(mmk(1,i+1,j))
          hsize2   = hsize2 + 1
        endif
        if ( mmk(1,i,j+1) > lubounds(domp_nt,2) ) then
          halosize = halosize + kh(mmk(1,i,j+1))
          hsize2   = hsize2 + 1
        endif
        if ( mmk(1,i-1,j) < lubounds(domp_nt,1) ) then
          halosize = halosize + kh(mmk(1,i-1,j))
          hsize2   = hsize2 + 1
        endif
        if ( mmk(1,i,j-1) < lubounds(domp_nt,1) ) then
          halosize = halosize + kh(mmk(1,i,j-1))
          hsize2   = hsize2 + 1
        endif
      enddo
      tsize3d(domp_nt) = thrdsize
      hsize2d(domp_nt) = hsize2
      hsize3d(domp_nt) = halosize
 
      ! print target and attained chunk sizes:
      write(iu06,'(a23,f10.4)') 'target load size:      ', chunk
      write(iu06,'(a23,i10)')  '2D halo size on exit:  ', sum(hsize2d)
      write(iu06,'(a23,i10)')  '3D halo size on exit:  ', sum(hsize3d)
      do tnum=1,domp_nt
        write(iu06,'(a14,i3,a16,f10.4)') 'openMP thread ', tnum,               &
         ' has load size: ', (a0 + a1*tsize3d(tnum) + a2*tsize2d(tnum)         &
                                 + a3*hsize2d(tnum) + a4*hsize3d(tnum))
      enddo
      
      ! print 3D stats
      icount = count(tsize3d > zero)
      if (icount > 0) then
        wmean =  sum(tsize3d, tsize3d > zero) / real(icount,8)
      else
        wmean = zero
      endif
      write(iu06,'(a52)')                                                      &
        'OpenMP 3D decomposition statistics (mean, min, max):'
      write(iu06,'(f10.0,i10,a2,i3,a2,i10,a2,i3,a2)')                          &
        wmean,                                                                 &
        minval(tsize3d), ' (', minloc(tsize3d), ') ',                          &
        maxval(tsize3d), ' (', maxloc(tsize3d), ') ' 

      ! print 3D halo stats:
      icount = count(hsize3d > zero)
      if (icount > 0) then
        wmean =  sum(hsize3d, hsize3d > zero) / real(icount,8)
      else
        wmean = zero
      endif
      write(iu06,'(a43)') 'OpenMP 3D halo statistics (mean, min, max):'
      write(iu06,'(f10.0,i10,a2,i3,a2,i10,a2,i3,a2)')                          &
        wmean,                                                                 &
        minval(hsize3d), ' (', minloc(hsize3d), ') ',                          &
        maxval(hsize3d), ' (', maxloc(hsize3d), ') ' 
      
      ! print 2D stats
      icount = count(tsize2d > zero)
      if (icount > 0) then
        wmean =  sum(tsize2d, tsize2d > zero) / real(icount,8)
      else
        wmean = zero
      endif
      write(iu06,'(a52)')                                                      &
        'OpenMP 2D decomposition statistics (mean, min, max):'
      write(iu06,'(f10.0,i10,a2,i3,a2,i10,a2,i3,a2)')                          &
        wmean,                                                                 &
        minval(tsize2d), ' (', minloc(tsize2d), ') ',                          &
        maxval(tsize2d), ' (', maxloc(tsize2d), ') '

      ! print 2D halo stats:
      icount = count(hsize2d > zero)
      if (icount > 0) then
        wmean =  sum(hsize2d, hsize2d > zero) / real(icount,8)
      else
        wmean = zero
      endif
      write(iu06,'(a43)') 'OpenMP 2D halo statistics (mean, min, max):'
      write(iu06,'(f10.0,i10,a2,i3,a2,i10,a2,i3,a2)')                          &
        wmean,                                                                 &
        minval(hsize2d), ' (', minloc(hsize2d), ') ',                          &
        maxval(hsize2d), ' (', maxloc(hsize2d), ') ' 

      ! print 2D ranges:
      do it=1,domp_nt
        write(iu06,'(a14,i3,a24,i10,i10)') 'openMP thread ', it,               &
         ' handles surface range: ', lubounds(it,1), lubounds(it,2)
      enddo
      call flush_unit(iu06)
!$OMP END MASTER
!$OMP BARRIER
    endif
#endif
  end subroutine domp_get_domain_wet_halo_idx

  ! ----------------------------------------------------------------------------

  subroutine domp_estimate_weighted_load (iu06, ia, iwet2, mmk, ind, kh, idx)
    use io_subs, only : flush_unit
    implicit none

    integer(4), intent(in)  :: iu06, ia, iwet2
    integer(4), intent(in)  :: mmk(1:,0:,0:), ind(1:,1:), kh(0:)
    integer(4), intent(out) :: idx(1:,1:)

    integer(4)            :: ib, it, nl, nu
    integer(4), parameter :: itmax = 2
    ! either do the experiment with b(:) or with a_i coefs:
    logical,    parameter :: do_a = .true.
    ! number of load sets to be estimated, 1 or 11:
    integer(4), parameter :: nb = 1 
    ! b(:) if do_a=F
    real(8)               :: b(1:nb)
    real(8),    parameter :: b0 = 0.05_8, db = 0.1_8
    ! coefs:                 k   N3  N2  H2  H3
    real(8)               :: a0, a1, a2, a3, a4

    ! first, make sure to allocate and turn on printouts:
    if (.not.omp_debug) then
      omp_debug = .true.
      if (.not.allocated(tsize3d))                                             &
        allocate(tsize2d(1:domp_nt),tsize3d(1:domp_nt),lubounds(1:domp_nt,1:2),&
                 hsize2d(1:domp_nt),hsize3d(1:domp_nt))
    endif

    ! then, define the weights (hardcoded experiment):
    if (do_a) then   ! NOTE: hardcoded IF construct is constant.
      b(:) = b0
      ! From Larry's f:
      !a0 = 118.88585_8
      !a1 = -0.00234967_8
      !a2 = -0.01324462_8
      !a3 = -0.00159484_8
      !a4 =  0.000052243881_8
      ! From Larry's 1/f, using 5 coefs:
      a0 = -0.079020449_8
      a1 =  2.2999119e-06_8
      a2 =  3.9998704e-05_8
      a3 =  5.8200283e-07_8
      a4 =  7.6242388e-08_8
      ! From Larry's 1/f, using only 3 coefs:
      !a0 = -0.077383095_8
      !a1 =  2.5610786e-06_8
      !a2 =  3.4530795e-05_8
      !a3 =  zero
      !a4 =  zero
    else
      b(1) = b0
      do ib=2,nb ! NOTE: loop has zero iteration counts when nb is hardwired =1
        b(ib) = (ib-1)*db
      enddo
      a0 = zero
      a1 = one
      a2 = zero
      a3 = zero
    endif

    write(iu06,*) '============================================================'
    write(iu06,*) 'Loadbalancing experiment, area:', ia
    do ib=1,nb
      write(iu06,*) '----------------------------------------------------------'
      if (do_a) then  ! NOTE: hardcoded IF construct is constant.
        write(iu06,'(a9,f10.5)') 'coef. a0=', a0
        write(iu06,'(a9,f10.7)') 'coef. a1=', a1
        write(iu06,'(a9,f10.7)') 'coef. a2=', a2
        write(iu06,'(a9,f10.7)') 'coef. a3=', a3
        write(iu06,'(a9,f10.7)') 'coef. a4=', a4
      else
        write(iu06,'(a9,f8.4)') 'coef. b= ', b(ib)
        a4 = b(ib)
      endif

      ! make sure to redefine idx before each initial guess:
      idx(:,:) = -1

      ! initial guess stored as usual in idx:
      write(iu06,*) 'Initial guess: -------------------------------------------'
!$OMP PARALLEL DEFAULT(shared) PRIVATE(nl, nu)
      call domp_get_domain(kh, 1, iwet2, nl, nu, idx)
!$OMP END PARALLEL
      call flush_unit(iu06)

      write(iu06,*) 'Iterations: ----------------------------------------------'
      it = 1
      do while (it <= itmax) ! .or. some_measure > some_threshold)
        write(iu06,*) 'it=', it
        ! obtain refined guesses (idx is inout):
!$OMP PARALLEL DEFAULT(shared) PRIVATE(nl, nu)
        call domp_get_domain(kh,1,iwet2,nl,nu,mmk,ind,a0,a1,a2,a3,a4,idx)
!$OMP END PARALLEL
        it = it + 1
        if (it <= itmax)  write(iu06,*)                                        &
                    '----------------------------------------------------------'
      enddo
    enddo
    write(iu06,*) '============================================================'
    call flush_unit(iu06)

    do it=1,domp_nt
      write(iu06,*) 'nt, n3d(nt):', it, tsize3d(it)
    enddo
    do it=1,domp_nt
      write(iu06,*) 'nt, n2d(nt):', it, tsize2d(it)
    enddo
    do it=1,domp_nt
      write(iu06,*) 'nt, h2d(nt):', it, hsize2d(it)
    enddo
    do it=1,domp_nt
      write(iu06,*) 'nt, h3d(nt):', it, hsize3d(it)
    enddo
    write(iu06,*) '============================================================'

  end subroutine domp_estimate_weighted_load

  ! ----------------------------------------------------------------------------

  subroutine domp_get_thread_no (tnum)
    implicit none

    integer(4), intent(out) :: tnum

    tnum = domp_iam + 1

  end subroutine domp_get_thread_no

  ! ----------------------------------------------------------------------------

end module dmi_omp


