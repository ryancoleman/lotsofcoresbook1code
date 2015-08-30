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

program emain

  !- Modules -------------------------------------------------------------------
  use io_subs,            only : io_new_unit
  use dmi_mpi_global,     only : iu06
  use dmi_mpi,            only : dmpi_init, dmpi_decompose, dmpi_debug, dd,    &
                                 decompose_gen, mpi_size
  use exits,              only : exitme
  use dmi_omp,            only : domp_init
  use e_timer,            only : timer, timer_init, timer_destroy, timer_print
  use constants,          only : nprocj, nproci
  use params_n,           only : narea, mmx, nmx, kmx, iw2, iw3, kh,           &
                                 msrf,mcol,mainarea,nthreads,hz,hz_f,getspecs, &
                                 niter,ladv, ldiff
  use local_arrays_n,     only : iw2_l, iw3_l, ind_l, kh_l, khu_l, khv_l,      &
                                 h_old_l, h_new_l, eddyh_l,dispt_l, disps_l,   &
                                 t_l, s_l, hz_l, h_l, hu_l, hv_l,              &
                                 ui_l, vi_l, wi_l, dx_l, dy_l,msrf_l,mcol_l,   &
                                 uvdam_l, cosphi_l, alloc_local_arrays, kh_lf, &
                                 mcol_lf,ind_lf,msrf_lf
  use tflow_simd_srf_col, only : tflow_int_simd, tflow_ft_simd
  use init_n,             only : init_local_vars, Validate
  use arrays_n,           only : idx, krz, krq, kru, krv, bndz, bndq, rwqk,    &
                                 AllocArrays
  use io_control,         only : ios_init, ios_usage

  !- Directives ----------------------------------------------------------------
  implicit none
  integer(4) :: i,ia,lun,n3s,ios
  character(5), parameter :: v = '(a53)'
  character(256)     :: fnam

  !- Initialize MPI ------------------------------------------------------------
  call dmpi_init(.true.)

  !- Initialize openmp ---------------------------------------------------------
  call domp_init(nthreads)

  !- Initialize timers ---------------------------------------------------------
  call timer_init(.true.,nthreads,10)

  call timer(1,'init')

  !- Obtain the basic parameters and mcol/msrf/kh/hz 
  call getspecs ()

  call ios_usage(iu06, .false., 0, 0, .true.)

  call ios_init()

  !- Set some MPI module vars based on model specs -----------------------------
  call dmpi_init(2, 0, mainarea, 4, nthreads)

  !- Allocate memory for main arrays -------------------------------------------
  do ia=1,narea
    call AllocArrays(ia,nthreads) 
  enddo

  !- Prepare the MPI decomposition ---------------------------------------------
  if (nprocj>1) then
    ! generate and exit
    call decompose_gen(narea, mmx, nmx, kmx, iw2, iw3, msrf, kh,nproci, nprocj)
  else
    !... reading the task decomposition from a file:
    call dmpi_decompose(narea, mmx, nmx, kmx, iw2, iw3, msrf, kh)
  endif
  
  !- Print some MPI debug info -------------------------------------------------
  call dmpi_debug(narea)

  !  ... and now we are ready to allocate task-local arrays:
  call Alloc_Local_Arrays(narea,nthreads)

  do ia=1,narea
    idx(ia)%p(:,:)=-1
  enddo

  do ia=1,narea
    call init_local_vars(ia,kh(ia)%p,kh_lf(ia)%p,hz(ia)%p,                     &
            hz_f(ia)%p,mcol(ia)%p,mcol_lf(ia)%p,msrf(ia)%p,msrf_lf(ia)%p,      &
            ind_lf(ia)%p,cosphi_l(ia)%p)
  enddo

!$OMP PARALLEL DEFAULT(shared) PRIVATE(ia)
  do ia=1,narea
    call init_local_vars(ia,iw2_l(ia),kh_lf(ia)%p,mcol_lf(ia)%p,msrf_lf(ia)%p, &
                         ind_lf(ia)%p,hz_f(ia)%p,                              &
                         kh_l(ia)%p, mcol_l(ia)%p,msrf_l(ia)%p,ind_l(ia)%p,    &
                         hz_l(ia)%p, khu_l(ia)%p, khv_l(ia)%p,                 &
                         h_l(ia)%p, h_old_l(ia)%p, h_new_l(ia)%p,t_l(ia)%p,    &
                         s_l(ia)%p,  ui_l(ia)%p, vi_l(ia)%p, wi_l(ia)%p,       &
                         hu_l(ia)%p, hv_l(ia)%p,                               &
                         dispt_l(ia)%p, disps_l(ia)%p,                         &
                         eddyh_l(ia)%p,uvdam_l(ia)%p,idx(ia)%p) 
  enddo
!$OMP END PARALLEL

  if (mpi_size==1) then
    call Validate()
  endif

!$OMP PARALLEL DEFAULT(shared) PRIVATE(ia)
  do ia=1,narea
    call tflow_ft_simd(mcol_l(ia)%p,iw2_l(ia),kh_l(ia)%p,  &
                  idx(ia)%p, iw3_l(ia))
  enddo
!$OMP END PARALLEL
  call timer(2,'init')
  call timer_print('TS - Initialisation took: ','init',1)
  ia=1
  if (ladv) then
    write(iu06,*) 'Advection kernel with 2 tracers and ',                     &
                   niter,' iterations '
    call timer(1,'tflow_adv')
!$OMP PARALLEL DEFAULT(shared) PRIVATE(i)
    do i = 1,niter
      ! do the advection:
      call tflow_int_simd(kmx(ia),msrf_lf(ia)%p,mcol_l(ia)%p,ind_l(ia)%p,      &
                     iw2_l(ia), dx_l(ia), dy_l(ia), 0.0_8, 5.0_8,              &
                     kh_l(ia)%p, khu_l(ia)%p, khv_l(ia)%p,                     &
                     h_l(ia)%p, hu_l(ia)%p, hv_l(ia)%p,                        &
                     h_old_l(ia)%p, h_new_l(ia)%p,                             &
                     ui_l(ia)%p, vi_l(ia)%p, wi_l(ia)%p,                       &
                     t_l(ia)%p, s_l(ia)%p,                                     &
                     dispt_l(ia)%p, disps_l(ia)%p, eddyh_l(ia)%p,              &
                     0, krz(ia)%p, 0, krq(ia)%p, rwqk(ia)%p,                   &
                     0, kru(ia)%p, 0 , krv(ia)%p,                              &
                     bndz(ia)%p, bndq(ia)%p, cosphi_l(ia)%p,                   &
                     .true., .false., idx(ia)%p, dd(ia)%halo2, ia,             &
                     .false., uvdam_l(ia)%p)
    enddo
!$OMP END PARALLEL
    call timer(2,'tflow_adv')
    call timer_print('TS - tflow_adv took: ','tflow_adv',1)
  endif

  if (ldiff) then
    write(iu06,*) 'Diffusion kernel with 14 tracers and ',                     &
                   20*niter,' iterations '
    call timer(1,'tflow_diff')
!$OMP PARALLEL DEFAULT(shared) PRIVATE(i)
    do i = 1,20*niter
      ! do the diffusion:
      call tflow_int_simd(kmx(ia), msrf_l(ia)%p,mcol_l(ia)%p,ind_l(ia)%p,      &
                     iw2_l(ia), dx_l(ia), dy_l(ia), 0.0_8, 5.0_8,              &
                     kh_l(ia)%p, khu_l(ia)%p, khv_l(ia)%p,                     &
                     h_l(ia)%p, hu_l(ia)%p, hv_l(ia)%p,                        &
                     h_old_l(ia)%p, h_new_l(ia)%p,                             &
                     ui_l(ia)%p, vi_l(ia)%p, wi_l(ia)%p,                       &
                     t_l(ia)%p, s_l(ia)%p,                                     &
                     dispt_l(ia)%p, disps_l(ia)%p, eddyh_l(ia)%p,              &
                     0, krz(ia)%p, 0, krq(ia)%p, rwqk(ia)%p,                   &
                     0, kru(ia)%p, 0 , krv(ia)%p,                              &
                     bndz(ia)%p, bndq(ia)%p, cosphi_l(ia)%p,                   &
                     .false., .true., idx(ia)%p, dd(ia)%halo2, ia,             &
                     .false., uvdam_l(ia)%p)
    enddo
!$OMP END PARALLEL
    call timer(2,'tflow_diff')
    call timer_print('TS - tflow_diff took: ','tflow_diff',1)
  endif

  ! generate a file with output fields to ensure binary identical results
    if (mpi_size==1) then
      lun = io_new_unit()
      fnam = 'tflow_output.bin'
      open(lun, file=trim(fnam), form='unformatted', action='write',           &
         status='unknown', access='stream', iostat=ios)
      if (ios /= 0) then
        write(iu06,*) 'Error: Cannot open binary file '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      n3s=iw3_l(1)
      call timer(1,'IOtime')
      write(lun,iostat=ios) t_l(1)%p,s_l(1)%p
      call timer(2,'IOtime')
      if (ios /= 0) then
        write(iu06,*) 'Error: Writing failed for binary file '//trim(fnam)//' exit'
        write(iu06,*) 'Error: The iostat failure was:: ', ios
        call exitme(1)
      endif
      close(lun)
      call Validate()
    endif
  call timer_print('TS - IO time took: ','IOtime',1)
  call timer_destroy()
  call exitme(0, 'END TIME IN HBM' )
end program emain
