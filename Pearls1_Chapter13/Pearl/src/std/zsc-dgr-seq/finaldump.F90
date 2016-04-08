!-------------------------------------------------------------------------------
! Subroutine : FinalDump
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Save relevant data at the end of the simulation run to file 'final.out'
!> @details
!! Generates the final data file for the simulation in the dual grid D2Q5/D2Q9
!! Zheng-Shu-Chew multiphase LBM, which contains:
!!
!! - Input parameters
!! - Estimated memory usage
!! - Pressure difference between the inside and the outside of the drop
!! - Error in the verification of Laplace's Law for the pressure
!! - Mass conservation factor
!! - Effective drop radius
!! - Maximum velocity in the domain

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
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, nodesIn, nodesOut
 INTEGER :: IO_ERR
 REAL(KIND = DBL) :: Pin, Pout, Ro, R, Pdif, Perr, Uloc, Umax, Ref, Vef, Vol
 REAL(KIND = DBL) :: distroMem, auxMem, totalMem, memUnitF, memUnitG

! Initialize
 Pin  =  0.D0
 Pout =  0.D0
 Vol  =  0.D0
 Umax = -1.D0
 nodesIn  = 0
 nodesOut = 0
 Ro = bubbles(1,3)

! Pressure inside and outside the bubble, maximum velocity and effective radius
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

     R =  DSQRT( ( DBLE(i)-bubbles(1,1) )**2 + ( DBLE(j)-bubbles(1,2) )**2 )
     IF ( R < (Ro - IntWidth) ) THEN
       Pin = Pin + p_f(i,j)
       nodesIn = nodesIn + 1
     ELSE IF ( R > (Ro + IntWidth) ) THEN
       Pout = Pout + p_f(i,j)
       nodesOut = nodesOut + 1
     END IF

     Uloc = DSQRT( u_f(i,j,1)*u_f(i,j,1) + u_f(i,j,2)*u_f(i,j,2) )
     IF ( Uloc > Umax ) Umax = Uloc
     IF ( phi_f(i,j) >= 0.D0 ) Vol = Vol + 1.D0

   END DO
 END DO

! Calculate compliance with Laplace Law
 Pin  = Pin/DBLE(nodesIn)
 Pout = Pout/DBLE(nodesOut)
 Pdif = Pin - Pout
 Perr = (sigma/Ro - Pdif)*Ro/sigma

! Calculate phase conservation
 Ref = DSQRT( Vol*invPi )
 Vef = Vol*invInitVol

! Estimate memory usage (Mb)
 memUnitG = NX*NY/( 1024.D0*1024.D0 )
 memUnitF = NX_f*NY_f/( 1024.D0*1024.D0 )
 distroMem = 8.D0*( 18.D0*memUnitG + 10.D0*memUnitF )
 auxMem    = 8.D0*( 3.D0*( memUnitG + memUnitF ) + 4.D0*memUnitF ) &
           + 4.D0*( 4.D0*( memUnitG + memUnitF ) )
 totalMem  = distroMem + auxMem

! The effective radius, sigma and IntWidth are re-scaled to the momentum grid
 OPEN(UNIT = 10, FILE = "final.out", STATUS = "NEW", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(10,'(A)')'*** Multiphase Zheng-Shu-Chew LBM 2D Simulation ***'
   WRITE(10,'(A)')'*** Dual Grid Implementation (Serial)           ***'
   WRITE(10,*)
   WRITE(10,'(A)')'INPUT PARAMETERS'
   WRITE(10,'(A,I9)')'Total Iterations      = ',MaxStep+1
   WRITE(10,'(A,I9)')'Relaxation Iterations = ',RelaxStep+1
   WRITE(10,'(A,I9)')'Length in X Direction = ',xmax
   WRITE(10,'(A,I9)')'Length in Y Direction = ',ymax
   WRITE(10,'(A,ES15.5)')'Interface Width   = ',0.5D0*IntWidth
   WRITE(10,'(A,ES15.5)')'Interface Tension = ',0.5D0*sigma
   WRITE(10,'(A,ES15.5)')'Interface Mobility= ',Gamma
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
   WRITE(10,'(A)')'OUTPUT RESULTS (RELAXATION)'
   WRITE(10,'(A,ES19.9)')'Effective Radius   = ',0.5D0*Ref
   WRITE(10,'(A,ES19.9)')'Phase Conservation = ',Vef
   WRITE(10,'(A,ES19.9)')'(Pin - Pout)       = ',Pdif
   WRITE(10,'(A,ES19.9)')'Laplace Error      = ',Perr
   WRITE(10,'(A,ES19.9)')'Parasitic Velocity = ',Umax
   WRITE(10,*)
   WRITE(10,'(A)')'***       Simulation Finished Succesfully       ***'
   CLOSE(UNIT = 10)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: unable to open output file 'final.out'."
 END IF

 RETURN
 END SUBROUTINE FinalDump
