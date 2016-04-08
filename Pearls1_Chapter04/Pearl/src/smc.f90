subroutine smc()

  use advance_module
  use chemistry_module
  use derivative_stencil_module
  use initialize_module
  use layout_module
  use multifab_module
  use omp_module
  use probin_module
  use runtime_init_module
  use smcdata_module
  use threadbox_module
  use time_module
  use variables_module

  implicit none

  integer :: init_step, istep, ioproc

  real(dp_t) :: dt, courno
  real(dp_t)  , pointer     :: dx(:)

  real(dp_t) :: wt1, wt2, wt_advance, wt_fb, wt_c2p, wt_c, wt_t, wt_hd
  type(layout)   :: la
  type(multifab) :: U

  call runtime_init()

  call stencil_init()

  call chemistry_init()

  call init_variables()

  call initialize_from_scratch(la,dt,courno,dx,U)

  call build_smcdata(la)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! print processor and grid info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (parallel_IOProcessor()) then
     print *, ' '     
     print *, 'number of MPI processes = ', parallel_nprocs()
     print *, 'number of threads       = ', omp_get_max_threads()
     print *, ' '
     print *, 'number of boxes         = ', nboxes(la)
     print *, ' '
  end if
  
  init_step = 1 

  if ( parallel_IOProcessor()) then
     print*,""
     print*,"BEGIN MAIN EVOLUTION LOOP"
     print*,""
  end if

  wt1 = parallel_wtime()

  if ( (max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0) ) then

     do istep = init_step, max_step

        if (parallel_IOProcessor()) then
           print*,'Advancing time step',istep,'time = ',time
        end if

        call advance(U,dt,courno,dx,istep)

        time = time + dt

        if ( parallel_IOProcessor() ) then
           print *, 'End of step', istep,'time = ', time
           print*, ""
        end if

        if (parallel_IOProcessor() .and. verbose .ge. 2) then
           call flush(6)
        end if

        ! have we reached the stop time?
        if (stop_time >= 0.d0) then
           if (time >= stop_time) then
              goto 999
           end if
        end if

     end do

999  continue

  end if

  wt2 = parallel_wtime()

  call destroy_smcdata()
  call destroy_threadbox()

  call destroy(U)

  call destroy(la)

  call chemistry_close()

  call runtime_close()

  deallocate(dx)

  ioproc = parallel_IOProcessorNode()
  call parallel_reduce(wt_advance, wt2-wt1        , MPI_MAX, proc=ioproc)
  call parallel_reduce(wt_fb     , wt_fillboundary, MPI_MAX, proc=ioproc)
  call parallel_reduce(wt_c2p    , wt_ctoprim     , MPI_MAX, proc=ioproc)
  call parallel_reduce(wt_c      , wt_chemterm    , MPI_MAX, proc=ioproc)
  call parallel_reduce(wt_t      , wt_transprop   , MPI_MAX, proc=ioproc)
  call parallel_reduce(wt_hd     , wt_hypdiff     , MPI_MAX, proc=ioproc)

  if (parallel_IOProcessor()) then
     print*, ' '
     print*, 'SMC Advance Time = ', wt_advance
     print*, '  fill boundary Time = ', wt_fb
     print*, '  ctoprim       Time = ', wt_c2p
     print*, '  chemterm      Time = ', wt_c
     print*, '  transprop     Time = ', wt_t
     print*, '  hyper & diff  Time = ', wt_hd
     print*, ' '
  end if

end subroutine smc
