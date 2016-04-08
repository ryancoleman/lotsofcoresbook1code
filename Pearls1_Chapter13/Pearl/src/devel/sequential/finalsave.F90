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
 USE OMP_LIB
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, m, mg, nthreads
 INTEGER :: SaveFiles, StatFiles, IO_ERR
 REAL(KIND = DBL) :: Pin, Pout, Ro, R, Pdif, Perr, Umax, Ref, Vef, Vol
 REAL(KIND = DBL) :: distroMem, auxMem, totalMem, memUnit, Mb
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiSq
 REAL(KIND = DBL) :: GFLOPS, GFLOPS_CORE, GFLOPS_IO
 REAL(KIND = DBL) :: NodesIn, NodesOut, UmaxGlobal
 REAL(KIND = DBL) :: Vef_test, Umax_test, Pdif_test

!$OMP PARALLEL
 nthreads = omp_get_num_threads()
!$OMP END PARALLEL

! Initialize
 Ro = bubbles(1,4)
 UmaxGlobal = 0.D0
 Umax       = 0.D0
 Vol        = 0.D0
 NodesIn    = 0.D0
 NodesOut   = 0.D0
 Pin        = 0.D0
 Pout       = 0.D0

! Calculate pressure inside and outside the bubble, maximum velocity in the
! domain and effective radius
 DO k = 1, NZ
   DO j = 1, NY
     DO i = 1, NX

   rhon    = g(i,j,k,0,now)  + g(i,j,k,1,now)  + g(i,j,k,2,now)  + g(i,j,k,3,now)  + g(i,j,k,4,now)  &
           + g(i,j,k,5,now)  + g(i,j,k,6,now)  + g(i,j,k,7,now)  + g(i,j,k,8,now)  + g(i,j,k,9,now)  &
           + g(i,j,k,10,now) + g(i,j,k,11,now) + g(i,j,k,12,now) + g(i,j,k,13,now) + g(i,j,k,14,now) &
           + g(i,j,k,15,now) + g(i,j,k,16,now) + g(i,j,k,17,now) + g(i,j,k,18,now)
   invRhon = 1.D0/rhon
       phin    = phi(i,j,k)
       phin2   = phin*phin

       gradPhiSq = gradPhiX(i,j,k)*gradPhiX(i,j,k) &
                 + gradPhiY(i,j,k)*gradPhiY(i,j,k) &
                 + gradPhiZ(i,j,k)*gradPhiZ(i,j,k)

! Calculate the pressure
       pressure = alpha*(3.D0*phin2*phin2 - 2.D0*phiStar2*( phin2 - phiStar2 ) ) &
                - kappa*( phin*lapPhi(i,j,k) - 0.5D0*gradPhiSq ) + invCs_sq*rhon

! Calculate the velocity
   ux = ( g(i,j,k,1,now)  - g(i,j,k,2,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  + g(i,j,k,9,now)  &
      -   g(i,j,k,10,now) + g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) - g(i,j,k,14,now) &
       )*invRhon
   uy = ( g(i,j,k,3,now)  - g(i,j,k,4,now)  + g(i,j,k,7,now)  - g(i,j,k,8,now)  - g(i,j,k,9,now)  &
      +   g(i,j,k,10,now) + g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) - g(i,j,k,18,now) &
       )*invRhon
   uz = ( g(i,j,k,5,now)  - g(i,j,k,6,now)  + g(i,j,k,11,now) - g(i,j,k,12,now) - g(i,j,k,13,now) &
      +   g(i,j,k,14,now) + g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) + g(i,j,k,18,now) &
       )*invRhon

       R =  DSQRT( (DBLE(i+xl)-bubbles(1,1))**2 &
         +         (DBLE(j+yl)-bubbles(1,2))**2 &
         +         (DBLE(k+zl)-bubbles(1,3))**2 )

       IF ( R < (Ro - IntWidth) ) THEN
         NodesIn = NodesIn + 1.D0
         Pin     = Pin + pressure
       ELSE IF ( R > (Ro + IntWidth) ) THEN
         NodesOut = NodesOut + 1.D0
         Pout     = Pout + pressure
       END IF

       IF ( phi(i,j,k) >= 0.D0 ) Vol = Vol + 1.D0

       Umax = DSQRT( ux*ux + uy*uy + uz*uz )
       IF ( Umax > UmaxGlobal ) UmaxGlobal = Umax

     END DO
   END DO
 END DO

! Calculate memory use (Mb)
 Mb      = 1.D0/( 1024.D0*1024.D0 )
 memUnit   = NXG*NYG*NZG
 distroMem = 8.D0*52.D0*memUnit*Mb
 auxMem    = 8.D0*5.D0*memUnit*Mb
 totalMem  = distroMem + auxMem

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
 GFLOPS       = (GFLOPS_CORE + GFLOPS_IO)*1.0D-9/tEvol

! Calculate compliance with Laplace Law
   Pin  = Pin/NodesIn
   Pout = Pout/NodesOut
   Pdif = Pin - Pout
   Perr = 0.5D0*(2.D0*sigma/Ro - Pdif)*Ro/sigma

! Calculate phase conservation
   Ref = ( Vol*invPi*0.75D0 )**(1.D0/3.D0)
   Vef = Vol*invInitVol

   OPEN(UNIT = 10, FILE = "runlog.out", STATUS = "NEW", POSITION = "APPEND", &
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(A)')'*** Multiphase Zheng-Shu-Chew LBM 3D Simulation        ***'
     WRITE(10,'(A)')'*** SEQUENTIAL-DEVEL Implementation for Intel Xeon Phi ***'
     WRITE(10,*)
     WRITE(10,'(A)')'INPUT PARAMETERS'
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
     STOP "Error: unable to open output file 'runlog.out'."
   END IF

! Write basic performance data to stdout
    WRITE(*,*)
    WRITE(*,*)'=== PERFORMANCE DATA'
    WRITE(*,'(A,ES15.5)')' GFLOPS           = ',GFLOPS
    WRITE(*,'(A,ES15.5)')' MLUPS            = ',MaxStep*DBLE(xmax)*DBLE(ymax)*DBLE(zmax)*0.000001/tEvol
    WRITE(*,'(A,ES15.5)')' Time (Evolution) = ',tEvol
    WRITE(*,'(A,ES15.5)')' Time (Total)     = ',tRun

 RETURN
 END SUBROUTINE FinalSave
