!-------------------------------------------------------------------------------
! Subroutine : FinalSave
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Generates the final data file for the simulation in the serial D3Q7/D3Q19
! Zheng-Shu-Chew multiphase LBM, which contains:
!
! - Input parameters
! - Estimated memory usage
! - Pressure difference between the inside and the outside of the drop
! - Error in the verification of Laplace's Law for the pressure
! - Mass conservation factor
! - Effective drop radius
! - Maximum velocity in the domain
! - Performance in MLUPS
! - Performance in FLOPS (approximate value)
! - Execution timings
!
! The parameters calculated are defined below.
!
! data(1) : Volume
! data(2) : NodesIn
! data(3) : NodesOut
! data(4) : Pin
! data(5) : POut
! data(6) : Memory for distribution functions
! data(7) : Memory for auxiliary arrays (velocity, gradients, etc ...)
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
! Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
!
! This file is part of MP-LABS.
!
! MP-LABS is free software: you can redistribute it and/or modify it under the
! terms of the GNU GPL version 3 or (at your option) any later version.
!
! MP-LABS is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! MP-LABS, in the file COPYING.txt. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------

 SUBROUTINE FinalSave

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Timers
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 USE MPI
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, mg, nthreads
 INTEGER :: SaveFiles, StatFiles, IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: Pin, Pout, Ro, R, Pdif, Perr, Umax, Ref, Vef, Vol
 REAL(KIND = DBL) :: distroMem, auxMem, mpiMem, totalMem, memUnit, Mb
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiSq
 REAL(KIND = DBL) :: GFLOPS, GFLOPS_CORE, GFLOPS_IO
 REAL(KIND = DBL) :: NodesIn, NodesOut, UmaxLocal, UmaxGlobal
 REAL(KIND = DBL) :: Vef_test, Umax_test, Pdif_test
 REAL(KIND = DBL), DIMENSION(1:9) :: dataLocal, dataGlobal

!$OMP PARALLEL
 nthreads = omp_get_num_threads()
!$OMP END PARALLEL

! Initialize
 Ro = bubbles(1,4)
 UmaxGlobal = 0.D0
 UmaxLocal  = 0.D0
 Umax       = 0.D0
 Vol        = 0.D0
 NodesIn    = 0.D0
 NodesOut   = 0.D0
 Pin        = 0.D0
 Pout       = 0.D0
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Calculate pressure inside and outside the bubble, maximum velocity in the
! domain and effective radius
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX
       
       m  = i + NXG*( j + NYG*k )
       mg = i + NXG*( j + NYG*k ) + now

       rhon  = g0(m)   + g1(mg)  + g2(mg)  + g3(mg)  + g4(mg)  + g5(mg)  + g6(mg)  & 
             + g7(mg)  + g8(mg)  + g9(mg)  + g10(mg) + g11(mg) + g12(mg) + g13(mg) & 
             + g14(mg) + g15(mg) + g16(mg) + g17(mg) + g18(mg)
       invRhon = 1.D0/rhon
       phin    = phi(m)
       phin2   = phin*phin

       gradPhiSq = gradPhiX(m)*gradPhiX(m) &
                 + gradPhiY(m)*gradPhiY(m) &
                 + gradPhiZ(m)*gradPhiZ(m)

! Calculate the pressure
       pressure = alpha*(3.D0*phin2*phin2 - 2.D0*phiStar2*( phin2 - phiStar2 ) ) &
                - kappa*( phin*lapPhi(m) - 0.5D0*gradPhiSq ) + invCs_sq*rhon

! Calculate the velocity
           ux = ( g1(mg) - g2(mg) + g7(mg) - g8(mg) + g9(mg) - g10(mg) + g11(mg) &
              -   g12(mg) + g13(mg) - g14(mg) )*invRhon
           uy = ( g3(mg) - g4(mg) + g7(mg) - g8(mg) - g9(mg) + g10(mg) + g15(mg) &
              -   g16(mg) + g17(mg) - g18(mg) )*invRhon
           uz = ( g5(mg) - g6(mg) + g11(mg) - g12(mg) - g13(mg) + g14(mg) &
              +   g15(mg) - g16(mg) - g17(mg) + g18(mg) )*invRhon

       R =  DSQRT( (DBLE(i+xl)-bubbles(1,1))**2 &
         +         (DBLE(j+yl)-bubbles(1,2))**2 &
         +         (DBLE(k+zl)-bubbles(1,3))**2 )

       IF ( R < (Ro - IntWidth) ) THEN
         dataLocal(2) = dataLocal(2) + 1.D0
         dataLocal(4) = dataLocal(4) + pressure
       ELSE IF ( R > (Ro + IntWidth) ) THEN
         dataLocal(3) = dataLocal(3) + 1.D0
         dataLocal(5) = dataLocal(5) + pressure
       END IF

       IF ( phi(m) >= 0.D0 ) dataLocal(1) = dataLocal(1) + 1.D0

       Umax = DSQRT( ux*ux + uy*uy + uz*uz )
       IF ( Umax > UmaxLocal ) UmaxLocal = Umax

     END DO
   END DO
 END DO

! Calculate memory use (Mb)
 Mb      = 1.D0/( 1024.D0*1024.D0 )
 memUnit   = NXG*NYG*NZG
 dataLocal(6) = 8.D0*52.D0*memUnit
 dataLocal(7) = 8.D0*memUnit
 dataLocal(8) = 8.D0*( 8.D0*( xsize + ysize + zsize) &
              + 4.D0*( xsize5 + ysize5 + zsize5 )    &
              + 16.D0*( xedge + yedge + zedge ) )

! Estimate performance (GFLOPS)
#ifdef LBMIO
 SaveFiles = MaxStep / tSave
 StatFiles = MaxStep / tStat
#else
 SaveFiles = 0
 StatFiles = 0
#endif
 GFLOPS_CORE  = ( DBLE(NXG)*DBLE(NYG)*DBLE(NZG)*414.D0                &
              + ( DBLE(xsize) + DBLE(ysize) + DBLE(zsize) )*10.D0  &
              + ( DBLE(xedge) + DBLE(yedge) + DBLE(zedge) )*20.D0 )*DBLE( MaxStep )
 GFLOPS_IO    = DBLE(NX)*DBLE(NY)*DBLE(NZ)*( DBLE(SaveFiles)*135 &
              + DBLE(StatFiles)*57 )
 dataLocal(9) = GFLOPS_CORE + GFLOPS_IO

! Gather global data
 CALL MPI_REDUCE( dataLocal, dataGlobal, 9, MPI_DOUBLE_PRECISION, MPI_SUM, &
                  master, MPI_COMM_VGRID, MPI_ERR )
 CALL MPI_REDUCE( UmaxLocal, UmaxGlobal, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                  master, MPI_COMM_VGRID, MPI_ERR )

! Save data to file from the master node only
 IF( vproc == master ) THEN
   Vol = dataGlobal(1)

! Calculate compliance with Laplace Law
   Pin  = dataGlobal(4)/dataGlobal(2)
   Pout = dataGlobal(5)/dataGlobal(3)
   Pdif = Pin - Pout
   Perr = 0.5D0*(2.D0*sigma/Ro - Pdif)*Ro/sigma

! Calculate phase conservation
   Ref = ( Vol*invPi )**(1.D0/3.D0)
   Vef = Vol*invInitVol

! Calculate Global memory use (Mb)
   Mb      = 1.D0/( 1024.D0*1024.D0 )
   distroMem = dataGlobal(6)*Mb
   auxMem    = dataGlobal(7)*Mb
   mpiMem    = dataGlobal(8)*Mb
   totalMem  = distroMem + auxMem + mpiMem

! Calculate Global performance (GFLOPS)
   GFLOPS    = dataGlobal(9)*1.0D-9/tEvol

   OPEN(UNIT = 10, FILE = "runlog.out", STATUS = "NEW", POSITION = "APPEND", &
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(A)')'*** Multiphase Zheng-Shu-Chew LBM 3D Simulation ***'
     WRITE(10,'(A)')'*** Optimized Hybrid (MPI/OMP) Implementation   ***'
     WRITE(10,*)
     WRITE(10,'(A)')'INPUT PARAMETERS'
     WRITE(10,'(A,I9)')'Number of MPI Tasks    = ',nprocs
     WRITE(10,'(A,I9)')'Number of OMP Threads  = ',nthreads
     WRITE(10,'(A,I9)')'Total Iterations       = ',MaxStep+1
     WRITE(10,'(A,I9)')'Relaxation Iterations  = ',RelaxStep
     WRITE(10,'(A,I9)')'Length in X Direction  = ',xmax
     WRITE(10,'(A,I9)')'Length in Y Direction  = ',ymax
     WRITE(10,'(A,I9)')'Length in Z Direction  = ',zmax
     WRITE(10,'(A,ES15.5)')'Interface Width        = ',IntWidth
     WRITE(10,'(A,ES15.5)')'Interface Tension      = ',sigma
     WRITE(10,'(A,ES15.5)')'Interface Mobility     = ',Gamma
     WRITE(10,'(A,ES15.5)')'RhoL    = ',rhoL
     WRITE(10,'(A,ES15.5)')'RhoH    = ',rhoH
     WRITE(10,'(A,ES15.5)')'TauRho  = ',tauRho
     WRITE(10,'(A,ES15.5)')'TauPhi  = ',tauPhi
     WRITE(10,*)
     WRITE(10,'(A)')'MEMORY USAGE (Mb)'
     WRITE(10,'(A,ES15.5)')'Distributions     = ',distroMem
     WRITE(10,'(A,ES15.5)')'Auxiliary Arrays  = ',auxMem
     WRITE(10,'(A,ES15.5)')'Total Memory Used = ',totalMem
     WRITE(10,*)
     WRITE(10,'(A)')'TIMINGS (sec)'
     WRITE(10,'(A,ES15.5)')'Setup      = ',tSetup
     WRITE(10,'(A,ES15.5)')'IO         = ',tIO
     WRITE(10,'(A,ES15.5)')'Relaxation = ',tRelax
     WRITE(10,'(A,ES15.5)')'Evolution  = ',tEvol
     WRITE(10,'(A,ES15.5)')'Total      = ',tRun
     WRITE(10,*)
     WRITE(10,'(A)')'PERFORMANCE DATA'
     WRITE(10,'(A,ES15.5)')'GFLOPS     = ',GFLOPS
     WRITE(10,'(A,ES15.5)')'MLUPS      = ',MaxStep*DBLE(xmax)*DBLE(ymax)*DBLE(zmax)*0.000001/tEvol
     WRITE(10,*)
     WRITE(10,'(A)')'OUTPUT RESULTS'
     WRITE(10,'(A,ES19.9)')'Effective Radius   = ',Ref
     WRITE(10,'(A,ES19.9)')'Phase Conservation = ',Vef
     WRITE(10,'(A,ES19.9)')'(Pin - Pout)       = ',Pdif
     WRITE(10,'(A,ES19.9)')'Laplace Error      = ',Perr
     WRITE(10,'(A,ES19.9)')'Parasitic Velocity = ',UmaxGlobal
     WRITE(10,*)
     WRITE(10,'(A)')'***       Simulation Finished Succesfully       ***'
     CLOSE(UNIT = 10)
   ELSE
     CALL MemAlloc(2)
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: unable to open output file 'runlog.out'."
   END IF

! Write basic performance data to stdout
    WRITE(*,*)
    WRITE(*,*)'=== PERFORMANCE DATA'
    WRITE(*,'(A,ES15.5)')' GFLOPS           = ',GFLOPS
    WRITE(*,'(A,ES15.5)')' MLUPS            = ',MaxStep*DBLE(xmax)*DBLE(ymax)*DBLE(zmax)*0.000001/tEvol
    WRITE(*,'(A,ES15.5)')' Time (Evolution) = ',tEvol
    WRITE(*,'(A,ES15.5)')' Time (Total)     = ',tRun
 END IF

 RETURN
 END SUBROUTINE FinalSave
