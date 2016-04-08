!-------------------------------------------------------------------------------
! Subroutine : FinalDump
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Save relevant data at the end of the simulation run to file 'final.out'
!> @details
!! Generates the final data file for the simulation in the parallel D3Q7/D3Q19
!! Zheng-Shu-Chew multiphase LBM, which contains:
!!
!! - Input parameters
!! - Estimated memory usage
!! - Pressure difference between the inside and the outside of the drop
!! - Error in the verification of Laplace's Law for the pressure
!! - Mass conservation factor
!! - Effective drop radius
!! - Maximum velocity in the domain
!!
!! Parallel implementation using MPI. Variables with the "Local" ending refer to
!! quantities calculated locally in the current processor, vproc. Variables with
!! the ending "Global" refer to quantities calculated on the complete domain.
!! The parameters stored in the MPI exchange buffers are defined below.
!!
!! @param data(1) : Volume
!! @param data(2) : NodesIn
!! @param data(3) : NodesOut
!! @param data(4) : Pin
!! @param data(5) : POut
!!
!! @param mem(1) : Memory for distribution functions
!! @param mem(2) : Memory for auxiliary arrays (velocity, gradients, etc ...)
!! @param mem(3) : Memory for MPI buffers

!-------------------------------------------------------------------------------
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

 SUBROUTINE FinalDump

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams,   ONLY : invCs_sq, g
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, k, ie, iw, jn, js, kt, kb
 INTEGER :: IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: Pin, Pout, Ro, R, Pdif, Perr, Umax, Ref, Vef, Vol
 REAL(KIND = DBL) :: distroMem, auxMem, mpiMem, totalMem, memUnit, Mb
 REAL(KIND = DBL) :: ux, uy, uz, phin, phin2, pressure, rhon, invRhon
 REAL(KIND = DBL) :: gradPhiX, gradPhiY, gradPhiZ, gradPhiSq, lapPhi
 REAL(KIND = DBL) :: in00, in01, in02, in03, in04, in05, in06, in07, in08, in09
 REAL(KIND = DBL) :: in10, in11, in12, in13, in14, in15, in16, in17, in18
 REAL(KIND = DBL) :: UmaxLocal, UmaxGlobal
 REAL(KIND = DBL), DIMENSION(1:3) :: memLocal, memGlobal
 REAL(KIND = DBL), DIMENSION(1:5) :: dataLocal, dataGlobal


! Initialize
 Ro = bubbles(1,4)
 UmaxLocal     = 0.D0
 UmaxGlobal    = 0.D0
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Calculate pressure inside and outside the bubble, maximum velocity in the
! domain and effective radius
 DO k = zl, zu
   DO j = yl, yu
     DO i = xl, xu

       rhon    = SUM( g(i,j,k,:,now) )
       invRhon = 1.D0/rhon
       phin    = phi(i,j,k)
       phin2   = phin*phin

! Identify neighbours
       ie = ni(i,j,k,1)
       iw = ni(i,j,k,2)
       jn = ni(i,j,k,3)
       js = ni(i,j,k,4)
       kt = ni(i,j,k,5)
       kb = ni(i,j,k,6)

! Nodal values of the order parameter
       in00 = phi(i ,j ,k )
       in01 = phi(ie,j ,k )
       in02 = phi(iw,j ,k )
       in03 = phi(i ,jn,k )
       in04 = phi(i ,js,k )
       in05 = phi(i ,j ,kt)
       in06 = phi(i ,j ,kb)
       in07 = phi(ie,jn,k )
       in08 = phi(iw,js,k )
       in09 = phi(ie,js,k )
       in10 = phi(iw,jn,k )
       in11 = phi(ie,j ,kt)
       in12 = phi(iw,j ,kb)
       in13 = phi(ie,j ,kb)
       in14 = phi(iw,j ,kt)
       in15 = phi(i ,jn,kt)
       in16 = phi(i ,js,kb)
       in17 = phi(i ,jn,kb)
       in18 = phi(i ,js,kt)

! Laplacian of the order parameter
       lapPhi = ( in07 + in08 + in09 + in10 + in11 + in12 + in13 + in14 &
              + in15 + in16 + in17 + in18 + 2.0D0*( in01 + in02 + in03  &
              + in04 + in05 + in06 - 12.D0*in00 ) )*inv6

! Components of the order parameter gradient
       gradPhiX  = ( 2.0D0*( in01 - in02 ) + in07 - in08 + in09 - in10 &
                 + in11 - in12 + in13 - in14 )*inv12

       gradPhiY  = ( 2.0D0*( in03 - in04 ) + in07 - in08 + in10 - in09 &
                 + in15 - in16 + in17 - in18 )*inv12

       gradPhiZ  = ( 2.0D0*( in05 - in06 ) + in11 - in12 + in14 - in13 &
                 + in15 - in16 + in18 - in17 )*inv12

       gradPhiSq = gradPhiX*gradPhiX + gradPhiY*gradPhiY + gradPhiZ*gradPhiZ

! Calculate the pressure
       pressure = alpha*(3.D0*phin2*phin2 - 2.D0*phiStar2*( phin2 - phiStar2 ) ) &
                - kappa*( phin*lapPhi - 0.5D0*gradPhiSq ) + invCs_sq*rhon

! Calculate the velocity
       ux = ( g(i,j,k, 1,now) - g(i,j,k, 2,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) + g(i,j,k, 9,now) - g(i,j,k,10,now) &
          +   g(i,j,k,11,now) - g(i,j,k,12,now) + g(i,j,k,13,now) &
          - g(i,j,k,14,now) )*invRhon

       uy = ( g(i,j,k, 3,now) - g(i,j,k, 4,now) + g(i,j,k, 7,now) &
          -   g(i,j,k, 8,now) - g(i,j,k, 9,now) + g(i,j,k,10,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) + g(i,j,k,17,now) &
          - g(i,j,k,18,now) )*invRhon

       uz = ( g(i,j,k, 5,now) - g(i,j,k, 6,now) + g(i,j,k,11,now) &
          -   g(i,j,k,12,now) - g(i,j,k,13,now) + g(i,j,k,14,now) &
          +   g(i,j,k,15,now) - g(i,j,k,16,now) - g(i,j,k,17,now) &
          + g(i,j,k,18,now) )*invRhon

       R =  DSQRT( (DBLE(i)-bubbles(1,1))**2 + (DBLE(j)-bubbles(1,2))**2 &
         + (DBLE(k)-bubbles(1,3))**2 )

       IF ( R < (Ro - IntWidth) ) THEN
         dataLocal(2) = dataLocal(2) + 1.D0
         dataLocal(4) = dataLocal(4) + pressure
       ELSE IF ( R > (Ro + IntWidth) ) THEN
         dataLocal(3) = dataLocal(3) + 1.D0
         dataLocal(5) = dataLocal(5) + pressure
       END IF

       IF ( phi(i,j,k) >= 0.D0 ) dataLocal(1) = dataLocal(1) + 1.D0

       Umax = DSQRT( ux*ux + uy*uy + uz*uz )
       IF ( Umax > UmaxLocal ) UmaxLocal = Umax

     END DO
   END DO
 END DO

! Gather global data
 CALL MPI_ALLREDUCE(dataLocal, dataGlobal, 5, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_ALLREDUCE(UmaxLocal, UmaxGlobal, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_VGRID, MPI_ERR)
 Vol = dataGlobal(1)

! Calculate compliance with Laplace Law
 Pin  = dataGlobal(4)/dataGlobal(2)
 Pout = dataGlobal(5)/dataGlobal(3)
 Pdif = Pin - Pout
 Perr = 0.5D0*(2.D0*sigma/Ro - Pdif)*Ro/sigma

! Calculate phase conservation
 Ref = ( Vol*invPi )**(1.D0/3.D0)
 Vef = Vol*invInitVol

! Estimate memory usage (Mb)
 Mb      = 1.D0/( 1024.D0*1024.D0 )
 memUnit = ( NX + 2 )*( NY + 2)*( NZ + 2)

 memLocal(1) = 8.D0*52.D0*memUnit
 memLocal(2) = 8.D0*memUnit + 4.D0*6.D0*memUnit
 memLocal(3) = 8.D0*( 8.D0*( xsize + ysize + zsize) &
             + 4.D0*( xsize5 + ysize5 + zsize5 )    &
             + 16.D0*( xedge + yedge + zedge ) )

 CALL MPI_ALLREDUCE(memLocal, memGlobal, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_VGRID, MPI_ERR)
 distroMem = memGlobal(1)*Mb
 auxMem    = memGlobal(2)*Mb
 mpiMem    = memGlobal(3)*Mb
 totalMem  = distroMem + auxMem + mpiMem

! Save data to file from the master node only
 IF( vproc == master ) THEN
   OPEN(UNIT = 10, FILE = "final.out", STATUS = "NEW", POSITION = "APPEND", &
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(A)')'*** Multiphase Zheng-Shu-Chew LBM 3D Simulation ***'
     WRITE(10,'(A)')'*** Standard Implementation (MPI Parallel)      ***'
     WRITE(10,*)
     WRITE(10,'(A)')'INPUT PARAMETERS'
     WRITE(10,'(A,I9)')'Total Iterations       = ',MaxStep+1
     WRITE(10,'(A,I9)')'Relaxation Iterations  = ',RelaxStep
     WRITE(10,'(A,I9)')'Length in X Direction  = ',xmax
     WRITE(10,'(A,I9)')'Length in Y Direction  = ',ymax
     WRITE(10,'(A,I9)')'Length in Z Direction  = ',zmax
     WRITE(10,'(A,I9)')'Number of X Partitions = ',mpi_xdim
     WRITE(10,'(A,I9)')'Number of Y Partitions = ',mpi_ydim
     WRITE(10,'(A,I9)')'Number of Z Partitions = ',mpi_zdim
     WRITE(10,'(A,I9)')'Total Number of CPUs   = ',mpi_xdim*mpi_ydim*mpi_zdim
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
     WRITE(10,'(A,ES15.5)')'MPI buffer Arrays = ',mpiMem
     WRITE(10,'(A,ES15.5)')'Total Memory Used = ',totalMem
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
     STOP "Error: unable to open output file 'final.out'."
   END IF
 END IF

 RETURN
 END SUBROUTINE FinalDump
