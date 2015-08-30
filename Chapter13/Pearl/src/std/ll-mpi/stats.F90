!-------------------------------------------------------------------------------
! Subroutine : Stats
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Calculate intermediate simulation results and save to file 'stats.out'
!> @details
!! Calculate average velocity, mass conservation factor, effective radius of the
!! drop, pressure difference between the inside and the outside of the drop and
!! the error with respect to the analytical value given by Laplace's equation,
!! and write them to file "stats.out" for the parallel D2Q9 Lee-Lin multiphase
!! LBM.
!!
!! The effective radius is calculated assuming the drop is a perfect circle with
!! area given by A = Pi*R*R.
!!
!! The pressure inside and the pressure outside of the drop are calculated as
!! the average pressures inside and outside the drop, excluding the interface
!! area.
!!
!! Parallel implementation using MPI. Variables with the "Local" ending refer to
!! quantities calculated locally in the current processor, vproc. Variables with
!! the ending "Global" refer to qunatities calculated on the complete domain.
!! The parameters stored in the MPI excahnge buffers are defined below.
!!
!! @param data(1) : Volume
!! @param data(2) : Ux
!! @param data(3) : Uy
!! @param data(4) : Pin
!! @param data(5) : NodesIn
!! @param data(6) : Pout
!! @param data(7) : NodesOut

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

 SUBROUTINE Stats

!  Common Variables
 USE NTypes,      ONLY : DBL
 USE Domain,      ONLY : invInitVol, invPi, iStep, tCall, xl, xu, yl, yu
 USE FluidParams, ONLY : bubbles, Convergence, eps, IntWidth, p, rho, rhoStar, sigma, u
 USE MPIParams,   ONLY : master, MPI_COMM_VGRID, vproc
 USE MPI
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j
 INTEGER :: IO_ERR, MPI_ERR
 REAL(KIND = DBL) :: Pdif, Perr, Pin, Pout, R, Ref, Ro, Vef, Vol, invVol, Ux, Uy
 REAL(KIND = DBL), DIMENSION(1:7) :: dataLocal, dataGlobal


! Initialization
 Ro = bubbles(1,3)
 dataLocal(:)  = 0.D0
 dataGlobal(:) = 0.D0

! Loop through all nodes inside the drop for the current processor, vproc
 DO j = yl, yu
   DO i = xl, xu

! Analyze mass conservation and average velocity
     IF ( rho(i,j) < rhoStar ) THEN
       dataLocal(1) = dataLocal(1) + 1.D0
       dataLocal(2) = dataLocal(2) + u(i,j,1)
       dataLocal(3) = dataLocal(3) + u(i,j,2)
     END IF

! Analyze pressure avoiding the interface region
     R =  DSQRT( ( DBLE(i)-bubbles(1,1) )**2 + ( DBLE(j)-bubbles(1,2) )**2 )
     IF ( R < (Ro - IntWidth) ) THEN
       dataLocal(4) = dataLocal(4) + p(i,j)
       dataLocal(5) = dataLocal(5) + 1.D0
     ELSE IF ( R > (Ro + IntWidth) ) THEN
       dataLocal(6) = dataLocal(6) + p(i,j)
       dataLocal(7) = dataLocal(7) + 1.D0
     END IF

   END DO
 END DO

! Gather global information
 CALL MPI_ALLREDUCE(dataLocal, dataGlobal, 7, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_VGRID, MPI_ERR)
 Vol    = dataGlobal(1)
 invVol = 1.D0/Vol

! Save the initial volume during the first time step
 IF( iStep == 0 ) invInitVol = invVol

! Define average velocity of the drop and effective radius and volume
 Ux  = dataGlobal(2)*invVol
 Uy  = dataGlobal(3)*invVol
 Ref = DSQRT( Vol*invPi )
 Vef = Vol*invInitVol

! Global value of the average pressure inside and outside the drop
 Pin  = dataGlobal(4)/dataGlobal(5)
 Pout = dataGlobal(6)/dataGlobal(7)
 Pdif = Pin - Pout
 Perr = (sigma/Ro - Pdif)*Ro/sigma

! Save velocity and effective radius and volume data from the master node only
 IF( vproc == master ) THEN
   OPEN(UNIT = 10, FILE = "stats.out", STATUS = "UNKNOWN", POSITION = "APPEND",&
        IOSTAT = IO_ERR)
   IF ( IO_ERR == 0 ) THEN
     WRITE(10,'(I9,6ES19.9)')iStep,Ux,Uy,Vef,Ref,Pdif,Perr
     CLOSE(UNIT = 10)
   ELSE
     CALL MemAlloc(2)
     CALL MPI_FINALIZE(MPI_ERR)
     STOP "Error: Unable to open output file 'stats.out'."
   END IF
 END IF

! Analyze convergence - Done over a certain period of time ( 10 x tStats ) to
! ensure that long range fluctuations do not give false positives
 Convergence(tCall) = Pdif
 IF ( MOD(tCall,11) == 0 ) THEN
   eps   = 0.D0
   tCall = 1
   DO i = 2, 11
     eps = eps + DABS( Convergence(i) - Convergence(i-1) )
   END DO
   eps = eps*0.1D0/Pdif
 ELSE
   tCall = tCall + 1
 END IF

 RETURN
 END SUBROUTINE Stats

