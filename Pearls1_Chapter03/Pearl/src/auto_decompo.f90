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

module auto_decompo

  !- modules: ------------------------------------------------------------------
  !  none

  !- implicit directives: ------------------------------------------------------
  implicit none
  private
  
  real(8), parameter :: zero = 0.0_8

  !- arrays needeed when calculating the exact solution
  integer(4), private, allocatable, save :: solution(:,:)
  real(8),    private, allocatable, save :: cost(:,:)

  !- public stuff --------------------------------------------------------------
  public :: autodecompose
  
  !- private stuff -------------------------------------------------------------
  private :: decompose_optimal_Mx1, decompose_write_MxN, decompose_optimal_1xN
contains

!===============================================================================

  subroutine autodecompose(ntaskmax, narea, mmx, nmx, iw3, msrf, kh, iu06)

    use constants, only : only_islices, autodecomp_op, decomp_version
    use cmod_mem,  only : cmi1,cmi2
    use exits,     only : exitme

    implicit none

    integer(4),  intent(in) :: ntaskmax, narea, iu06
    integer(4),  intent(in) :: mmx(:), nmx(:), iw3(:)
    type (cmi1), intent(in) :: kh(:)
    type (cmi2), intent(in) :: msrf(:)

    integer(4)              :: nproci, nprocj, nact, maxit, maxtt, ia, it
    integer(4)              :: ntasksj, max_mmx
    integer(4), allocatable :: jtask(:,:), jlows(:,:), jhghs(:,:)
    integer(4), allocatable :: ilows(:,:), ihghs(:,:)
    real(8),    allocatable :: cmax(:), cmin(:), hmax(:), hmin(:)

    !- set No. of tasks in J-direction:
    if (only_islices) then
      ntasksj = 1
    else
      ntasksj = ntaskmax
    endif

    !- do some allocations
    allocate(jtask(ntasksj,narea),jlows(ntasksj,narea),jhghs(ntasksj,narea),   &
             ilows(ntaskmax,narea),ihghs(ntaskmax,narea),                      &
             cmax(narea), cmin(narea), hmax(narea), hmin(narea))
    cmax(:) = zero
    cmin(:) = zero

    if (decomp_version == 3) then
      max_mmx = maxval(mmx(:))
      allocate(cost(max_mmx, ntaskmax), solution(max_mmx,ntaskmax))
    endif

    !- run through number of i-slices
    do nproci=1,ntaskmax

      !- autodecompose with attempt to keep the tasks well-balanced ------------
      do nprocj=1,ntasksj

  
        !- find j-slices:
        nact = 0
        jtask(1:nprocj,:) = 0
        jlows(1:nprocj,:) = 0
        jhghs(1:nprocj,:) = 0

        do ia=1,narea
          call decompose_optimal_1xN(ia, nprocj, mmx(ia), nmx(ia),             &
                                     iw3(ia), msrf(ia)%p, kh(ia)%p, nact,      &
                                     jtask, jlows, jhghs)

          ! find sliced task with max No. of wet points:
          maxit = 0
          maxtt = 0
          do it=1,nprocj
            if (maxtt < jtask(it,ia)) then
              maxtt = jtask(it,ia)
              maxit = it
            endif
          enddo

          ! sub-divide the largest slice:
          ilows(1:nproci,ia) = 0
          ihghs(1:nproci,ia) = 0
          if (decomp_version == 1) then
            call decompose_optimal_Mx1(nproci,mmx(ia),msrf(ia)%p,              &
                             kh(ia)%p,jlows(maxit,ia), jhghs(maxit,ia), maxtt, &
                                       ilows(:,ia), ihghs(:,ia),autodecomp_op)
          else
            call exitme(1,"decomp_verision in options.nml can only be 1,2 or 3")
          endif                         
        enddo        

        ! sub-divide each j-slice according to the largest one, and 
        ! write this nproci by nprocj decomposition: 
        call decompose_write_MxN( nproci, nprocj, narea, msrf, kh, iw3(1),     &
                                 ilows, ihghs, jlows, jhghs, ntaskmax, iu06,   &
                                 cmax(:), cmin(:), hmax(:), hmin(:) ) 
      enddo
    enddo

    !- clean up
    deallocate( jtask, jlows, jhghs, ilows, ihghs, cmax, cmin, hmax, hmin )
    if (decomp_version == 3) deallocate(cost, solution)

  end subroutine autodecompose

!===============================================================================

  subroutine decompose_optimal_1xN(ia, nj, mmx, nmx, iw3, msrf, kh,            &
                                   nact, jtask, jlows, jhghs, jnsrf)
    implicit none

    integer(4), intent(in)    :: ia, nj, mmx, nmx, iw3
    integer(4), intent(inout) :: nact
    integer(4), intent(in)    :: msrf(0:,0:)
    integer(4), intent(in)    :: kh(0:)
    integer(4), intent(inout) :: jtask(:,:), jlows(:,:), jhghs(:,:)
    integer(4), intent(inout), optional :: jnsrf(:,:)

    integer(4) :: iweven, iw, i, j, jl, jsum, jsl, js, tsum, jsum2
    integer(4) :: jdist(nmx)

    !- even split:
    iweven = iw3/nj

    !- distribution of wet points in j-direction:
    do j=1,nmx
      iw = 0
      do i=1,mmx
        if (msrf(i,j) > 0) iw = iw + kh(msrf(i,j))
      enddo
      jdist(j) = iw
    enddo

    !- find j-slices:
    tsum  = 0
    jsum  = 0
    jsum2 = 0
    jl    = 1
    jsl   = 1
    do j=1,nmx
      jsum = jsum + jdist(j)
      if (j < nmx) jsum2 = jsum + jdist(j+1)
      ! either we have enough wetpoints, or we run out of 
      ! j-lines or of tasks
      if (jsum >= iweven .or. j == nmx .or. jsl == nj .or.                     &
          (j < nmx .and. jsum2 > iweven .and. jsum2-iweven > iweven-jsum) ) then
        jtask(jsl,ia) = jsum
        jlows(jsl,ia) = jl
        jhghs(jsl,ia) = j
        if (jsl == nj) then
          jtask(jsl,ia) = iw3 - tsum
          jhghs(jsl,ia) = nmx
          exit
        endif
        jl    = j + 1
        jsl   = jsl + 1
        tsum  = tsum + jsum
        jsum  = 0
        jsum2 = 0
      endif
    enddo
    nact = max(jsl, nact)

    if (present(jnsrf)) then
      !- find No of srf wetpoints in each slice:
      do js=1,jsl
        if (jtask(js,ia) <= 0) cycle
        iw = 0
        do j=jlows(js,ia),jhghs(js,ia)
          do i=1,mmx
            if (msrf(i,j) > 0) iw = iw + 1
          enddo
        enddo
        jnsrf(js,ia) = iw
      enddo
    endif

  end subroutine decompose_optimal_1xN


  subroutine decompose_optimal_Mx1(ni, mmx, msrf, kh, jlow, jhgh, maxtt,       &
                                   ilows, ihghs, offset_penalty, idist_out,    &
                                   isl_out)
    implicit none

    integer(4), intent(in)    :: ni, mmx, jlow, jhgh, maxtt
    integer(4), intent(in)    :: msrf(0:,0:)
    integer(4), intent(in)    :: kh(0:)
    logical,    intent(in)    :: offset_penalty
    integer(4), intent(inout) :: ilows(:), ihghs(:)
    integer(4), intent(out), optional :: idist_out(:), isl_out

    integer(4) :: iweven, iw, i, j, il, isum, isl, tsum, isum2
    integer(4) :: idist(mmx), rounding_offset, penalty_scal

    !- scaling factor controlling the penality we add when an islice is not
    !  equal to iweven (the even split)
    if (offset_penalty) then
      penalty_scal = 1
    else
      penalty_scal = 0
    endif

    !- even split:
    iweven = maxtt/ni

    !- distribution of wet points in i-direction:
    do i=1,mmx
      iw = 0
      do j=jlow,jhgh
        if (msrf(i,j) > 0) iw = iw + kh(msrf(i,j))
      enddo
      idist(i) = iw
    enddo

    !- find i-slices:
    tsum  = 0
    isum  = 0
    isum2 = 0
    il    = 1
    isl   = 1
    rounding_offset = 0
    do i=1,mmx
      isum = isum + idist(i)
      if (i < mmx) isum2 = isum + idist(i+1)
      ! out of i-lines, put remainding points on the last mpi task
      if (isl == ni) then
        ilows(isl) = il
        ihghs(isl) = mmx
        exit
      endif

      ! Do we have enough wet points to split the domain here?

      ! If offset_penalty = true, then we will try to keep track on how far from
      ! the evensplit we are. That is, if we have used to many wet points on
      ! task 1-N compared to the even split, then we will motivate task N+1 to
      ! make a smaller spit than iweven.
      if ( (isum >= iweven .or.                             &
         (i < mmx .and. isum2 > iweven .and.                                   &
         isum2-iweven-rounding_offset*penalty_scal > iweven-isum)) ) then

        ilows(isl)      = il
        ihghs(isl)      = i
        rounding_offset = rounding_offset + iweven-isum
        il    = i + 1
        isl   = isl + 1
        tsum  = tsum + isum
        isum  = 0
        isum2 = 0
      elseif (i==mmx) then
        ilows(isl) = il
        ihghs(isl) = i
      endif
    enddo

    if (present(idist_out)) idist_out(:) = idist(:)
    if (present(isl_out))   isl_out      = isl
  end subroutine decompose_optimal_Mx1

!===============================================================================
  subroutine decompose_write_MxN(nproci, nprocj, narea, msrf, kh, iw3,         &
                                 ilows, ihghs, jlows, jhghs, ntaskmax, iu06,   &
                                 tcost_max, tcost_min, hcost_max, hcost_min)

    use cmod_mem,  only : cmi2, cmi1
    use io_subs,   only : io_new_unit
    use exits,     only : exitme
    use constants, only : decomp_coeff, decomp_version

    implicit none 

    integer(4), intent(in) :: nproci, nprocj, narea, ntaskmax
    integer(4), intent(in) :: iu06
    type(cmi2), intent(in) :: msrf(:)
    type(cmi1), intent(in) :: kh(:)
    integer(4), intent(in) :: ilows(:,:), ihghs(:,:), iw3
    integer(4), intent(in) :: jlows(:,:), jhghs(:,:)
    real(8),    intent(in) :: tcost_max(:), tcost_min(:), hcost_max(:),        &
                              hcost_min(:)

    logical,    save :: FirstEntry = .true.
    integer(4), save :: luns
    integer(4)       :: lund, ios, i, j, ip, jp, ia, ntmx, nt
    integer(4)       :: il, iu, jl, ju, iw
    integer(4)       :: n_w, n_n, n_e, n_s, n_nw, n_ne, n_se, n_sw
    integer(4)       :: iwp(nproci,nprocj,narea), itt(nproci*nprocj,1:2,narea)
    integer(4)       :: ntt(nproci,nprocj,narea)
    character(1)     :: c1
    character(2)     :: c2
    character(3)     :: c3
    character(4)     :: c4
    character(256)   :: cfn, cfn2
    real(8)          :: chgh, clow, cha(narea), cla(narea), riw, hratio

    !- on first entry, open stats file:
    if (FirstEntry) then
      luns = io_new_unit()
      if     (ntaskmax < 10   ) then
        write(c1,'(i1)') ntaskmax
        cfn  = 'decompo_stats_'//c1//'.txt'
      elseif (ntaskmax < 100  ) then
        write(c2,'(i2)') ntaskmax
        cfn  = 'decompo_stats_'//c2//'.txt'
      elseif (ntaskmax < 1000 ) then
        write(c3,'(i3)') ntaskmax
        cfn  = 'decompo_stats_'//c3//'.txt'
      elseif (ntaskmax < 10000) then
        write(c4,'(i4)') ntaskmax
        cfn  = 'decompo_stats_'//c4//'.txt'
      else
        call exitme(1,'Cannot handle more than 9999 by 9999 tasks, sorry')
      endif
      open (unit=luns, file=trim(cfn), status='replace', iostat=ios)
      if (ios /= 0) call exitme(1,'Cannot create decompo_stats file')
      
      write(luns,*) 'Creating decomposition files using coeffs:'
      write(luns,'(a8,es12.5)') 'c_base =', decomp_coeff(1) 
      write(luns,'(a8,es12.5)') 'c_wet  =', decomp_coeff(2)
      write(luns,'(a8,es12.5)') 'c_halo =', decomp_coeff(3)
      FirstEntry = .false.
    endif

    !- on every entry, open decompo file:
    lund = io_new_unit()
    if     (nproci < 10   ) then
      write(c1,'(i1)') nproci
      cfn  = 'mpi_decompo_'//c1//'x'
    elseif (nproci < 100  ) then
      write(c2,'(i2)') nproci
      cfn  = 'mpi_decompo_'//c2//'x'
    elseif (nproci < 1000 ) then
      write(c3,'(i3)') nproci
      cfn  = 'mpi_decompo_'//c3//'x'
    elseif (nproci < 10000) then
      write(c4,'(i4)') nproci
      cfn  = 'mpi_decompo_'//c4//'x'
    endif
    if     (nprocj < 10   ) then
      write(c1,'(i1)') nprocj
      cfn2 = trim(cfn)//c1//'.txt'
    elseif (nprocj < 100  ) then
      write(c2,'(i2)') nprocj
      cfn2 = trim(cfn)//c2//'.txt'
    elseif (nprocj < 1000 ) then
      write(c3,'(i3)') nprocj
      cfn2 = trim(cfn)//c3//'.txt'
    elseif (nprocj < 10000) then
      write(c4,'(i4)') nprocj
      cfn2 = trim(cfn)//c4//'.txt'
    endif
    open (unit=lund, file=trim(cfn2), status='replace', iostat=ios)
    if (ios /= 0) call exitme(1,'Cannot create decompo file')

    !---------------------------------------------------------------------------
    !   run through this nproci by nprocj decompo, find pruned decompo,
    !   obtain statistics, and dump the decompo.
    itt(:,:,:) = -1
    ntt(:,:,:) = -1
    ntmx = 0
    chgh = 0.0_8
    clow = 0.0_8
    do ia=1,narea
      clow = max( clow, real(iw3,8))
    enddo
    do ia=1,narea
      nt = 0
      cha(ia) = 0.0_8
      cla(ia) = real(iw3,8)
      do jp=1,nprocj
        do ip=1,nproci
          iw = 0
          do j=jlows(jp,ia),jhghs(jp,ia)
            do i=ilows(ip,ia),ihghs(ip,ia)
              iw = iw + kh(ia)%p(msrf(ia)%p(i,j))
            enddo
          enddo
          if (iw > 0) then
            nt = nt + 1
            iwp(ip,jp,ia) = iw
            itt(nt,1,ia)  = ip
            itt(nt,2,ia)  = jp
            ntt(ip,jp,ia) = nt
            riw     = real(iw,8)
            chgh    = max( chgh, riw )
            clow    = min( clow, riw )
            cha(ia) = max( cha(ia), riw )
            cla(ia) = min( cla(ia), riw )
          else
            iwp(ip,jp,ia) = 0
          endif
        enddo
      enddo
      ntmx = max( ntmx, nt )
    enddo
    if (ntmx > 0) then
      write(luns,'(a50)') '=================================================='
      write(luns,'(a18,2i5)')   'Auto-decompo:     ',nproci, nprocj
      write(luns,'(a18,i5)')    'non-empty tasks:  ',ntmx
!fixme: shoud have a nicer format:
      write(luns,*)             'global Cmax/Cmin: ',clow, chgh, chgh/clow
      do ia=1,narea
        write(luns,*)           'area, Cmax/Cmin:  ',ia, cla(ia), cha(ia),     &
                                 cha(ia)/cla(ia)
        if (.not. decomp_version == 1) then
          if (hcost_min(ia) > zero) then
            hratio = hcost_max(ia)/hcost_min(ia)
          else
            hratio = zero
          endif
          write(luns,*)         'area, cost:       ',ia, tcost_min(ia),        &
                                 tcost_max(ia), tcost_max(ia)/tcost_min(ia)
          write(luns,*)         'area, hcost:      ',ia, hcost_min(ia),        &
                                 hcost_max(ia), hratio
        endif

        write(luns,'(a18,2i5)') 'area, empty tasks:',ia,ntmx-maxval(ntt(:,:,ia))
      enddo

      write(lund,'(I5)') ntmx
    else
      write(iu06,*) 'Cannot decompose to zero tasks'
      write(iu06,*) ' ... continue with next decompo'
      close(lund)
      return
    endif

    do nt=1,ntmx
        iw = 0
        do ia=1,narea
          if (itt(nt,1,ia) <= 0) cycle
          ip = itt(nt,1,ia)
          jp = itt(nt,2,ia)
          iw = iw + iwp(ip,jp,ia)
        enddo
        if (iw <= 0) cycle
        do ia=1,narea
          if (itt(nt,1,ia) > 0) then
            ip = itt(nt,1,ia)
            jp = itt(nt,2,ia)

            ! define task:
            il = ilows(ip,ia)
            iu = ihghs(ip,ia)
            jl = jlows(jp,ia)
            ju = jhghs(jp,ia)

            ! define neighbour tasks:
            write(lund,'(6i5)') nt, ia, il, iu, jl, ju
            if (jp > 1) then
              n_w = ntt(ip,jp-1,ia)
            else
              n_w = -1
            endif
            if (ip > 1) then
              n_n = ntt(ip-1,jp,ia)
            else
              n_n = -1
            endif
            if (jp < nprocj) then
              n_e = ntt(ip,jp+1,ia)
            else
              n_e = -1
            endif
            if (ip < nproci) then
              n_s = ntt(ip+1,jp,ia)
            else
              n_s = -1
            endif
            write(lund,'(10x,4i5)') n_w, n_n, n_e, n_s

            ! define neighbour-corner tasks:
            if (ip > 1 .and. jp > 1) then
              n_nw = ntt(ip-1,jp-1,ia)
            else
              n_nw = -1
            endif
            if (ip > 1 .and. jp < nprocj) then
              n_ne = ntt(ip-1,jp+1,ia)
            else
              n_ne = -1
            endif
            if (ip < nproci .and. jp < nprocj) then
              n_se = ntt(ip+1,jp+1,ia)
            else
              n_se = -1
            endif
            if (ip < nproci .and. jp > 1) then
              n_sw = ntt(ip+1,jp-1,ia)
            else
              n_sw = -1
            endif
            write(lund,'(10x,4i5)') n_nw, n_ne, n_se, n_sw

          else
            write(lund,'(6i5)') nt, ia, 0, 0, 0, 0
            write(lund,'(10x,4i5)')    -1,-1,-1,-1
            write(lund,'(10x,4i5)')    -1,-1,-1,-1
          endif  
        enddo
    enddo
    !---------------------------------------------------------------------------

    !- close files
    close(lund)
    if (nproci == ntaskmax .and. nprocj == ntaskmax) close(luns)

  end subroutine decompose_write_MxN


end module auto_decompo

