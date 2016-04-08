!-------------------------------------------------------------------------------
! Subroutine : FinalDump
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Save relevant data at the end of the simulation run to file 'final.out'
!> @details
!! Generates the final data file for the simulation in the parallel dual grid
!! D2Q9 Lee-Lin multiphase LBM, which contains:
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
 USE NTypes, ONLY : DBL
 USE Domain
 USE FluidParams
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j
 INTEGER :: IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: Pin, Pout, Ro, R, Pdif, Perr, Umax, Ref, Vef, Vol
 REAL(KIND = DBL) :: distroMem, auxMem, mpiMem, totalMem, Mb
 REAL(KIND = DBL) :: memUnitF1, memUnitF2, memUnitF3, memUnitG1, memUnitG2
 REAL(KIND = DBL) :: UmaxLocal, UmaxGlobal
 REAL(KIND = DBL), DIMENSION(1:3) :: memLocal, memGlobal
 REAL(KIND = DBL), DIMENSION(1:5) :: dataLocal, dataGlobal


! Initialize
 Ro = bubbles(1,3)
 UmaxLocal     = 0.D0
 UmaxGlobal    = 0.D0
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Calculate pressure difference, maximum velocity and effective radius
 DO j = yl_f, yu_f
   DO i = xl_f, xu_f

     R =  DSQRT((DBLE(i)-bubbles(1,1))**2 + (DBLE(j)-bubbles(1,2))**2)
     IF ( R < (Ro - IntWidth) ) THEN
       dataLocal(2) = dataLocal(2) + 1.D0
       dataLocal(4) = dataLocal(4) + p_f(i,j)
     ELSE IF ( R > (Ro + IntWidth) ) THEN
       dataLocal(3) = dataLocal(3) + 1.D0
       dataLocal(5) = dataLocal(5) + p_f(i,j)
     END IF

     IF ( rho_f(i,j) < rhoStar ) dataLocal(1) = dataLocal(1) + 1.D0

     Umax = DSQRT( u_f(i,j,1)*u_f(i,j,1) + u_f(i,j,2)*u_f(i,j,2) )
     IF ( Umax > UmaxLocal ) UmaxLocal = Umax

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
 Perr = (sigma/Ro - Pdif)*Ro/sigma

! Calculate phase conservation
 Ref = DSQRT( Vol*invPi )
 Vef = Vol*invInitVol

! Estimate memory usage (Mb)
 Mb       = 1.D0/( 1024.D0*1024.D0 )
 memUnitG1 = ( NX + 2 )*( NY + 2 )
 memUnitG2 = ( NX + 4 )*( NY + 4 )
 memUnitF1 = ( NX_f + 2 )*( NY_f + 2 )
 memUnitF2 = ( NX_f + 4 )*( NY_f + 4 )
 memUnitF3 = ( NX_f + 6 )*( NY_f + 6 )

 memLocal(1) = 8.D0*18.D0*( memUnitG1 + memUnitF1 )
 memLocal(2) = 8.D0*( 7.D0*memUnitG1 + 9.D0*memUnitF1 + memUnitF2 + memUnitF3 ) &
             + 4.D0*4.D0*( memUnitG2 + memUnitF2 )
 memLocal(3) = 8.D0*( 8.D0*( xsize2 + ysize2 + xsize3 + ysize3 + xsize3_f + ysize3_f) + 144.D0 )

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
     WRITE(10,'(A)')'*** Multiphase Lee-Lim LBM 2D Simulation    ***'
     WRITE(10,'(A)')'*** Dual Grid Implementation (MPI Parallel) ***'
     WRITE(10,*)
     WRITE(10,'(A)')'INPUT PARAMETERS'
     WRITE(10,'(A,I9)')'Total Iterations      = ',MaxStep+1
     WRITE(10,'(A,I9)')'Relaxation Iterations = ',RelaxStep+1
     WRITE(10,'(A,I9)')'Length in X Direction = ',xmax
     WRITE(10,'(A,I9)')'Length in Y Direction = ',ymax
     WRITE(10,'(A,ES15.5)')'Interface Width       = ',0.5D0*IntWidth
     WRITE(10,'(A,ES15.5)')'Interfacial Tension   = ',0.5D0*sigma
     WRITE(10,'(A,ES15.5)')'TauL   = ',tauL
     WRITE(10,'(A,ES15.5)')'TauH   = ',tauH
     WRITE(10,'(A,ES15.5)')'RhoL   = ',rhoL
     WRITE(10,'(A,ES15.5)')'RhoH   = ',rhoH
     WRITE(10,*)
     WRITE(10,'(A)')'MEMORY USAGE (Mb)'
     WRITE(10,'(A,ES15.5)')'Distributions     = ',distroMem
     WRITE(10,'(A,ES15.5)')'Auxiliary Arrays  = ',auxMem
     WRITE(10,'(A,ES15.5)')'MPI Buffer Arrays = ',mpiMem
     WRITE(10,'(A,ES15.5)')'Total Memory Used = ',totalMem
     WRITE(10,*)
     WRITE(10,'(A)')'OUTPUT RESULTS (RELAXATION)'
     WRITE(10,'(A,ES19.9)')'Effective Radius   = ',0.5D0*Ref
     WRITE(10,'(A,ES19.9)')'Phase Conservation = ',Vef
     WRITE(10,'(A,ES19.9)')'(Pin - Pout)       = ',2.D0*Pdif
     WRITE(10,'(A,ES19.9)')'Laplace Error      = ',Perr
     WRITE(10,'(A,ES19.9)')'Parasitic Velocity = ',UmaxGlobal
     WRITE(10,*)
     WRITE(10,'(A)')'***    Simulation Finished Succesfully      ***'
     CLOSE(UNIT = 10)
   ELSE
     CALL MemAlloc(2)
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: unable to open output file 'final.out'."
   END IF
 END IF

 RETURN
 END SUBROUTINE FinalDump
