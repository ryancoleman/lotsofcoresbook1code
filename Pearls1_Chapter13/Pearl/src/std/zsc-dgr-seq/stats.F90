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
!! and write them to file "stats.out" for the dual grid D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM.
!!
!! The effective radius is calculated assuming the drop is a perfect circle with
!! area given by A = Pi*R*R.
!!
!! The pressure inside and the pressure outside of the drop are calculated as
!! the average pressures inside and outside the drop, excluding the interface
!! area

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
 USE Domain,      ONLY : invInitVol, invPi, iStep, tCall, xmax_f, xmin_f, ymax_f, ymin_f
 USE FluidParams, ONLY : bubbles, Convergence, eps, IntWidth, p_f, phi_f, sigma, u_f
 IMPLICIT NONE

!  Local Variables
 INTEGER :: i, j, nodesIn, nodesOut
 INTEGER :: IO_ERR
 REAL(KIND = DBL) :: Pin, Pout, Pdif, Perr, Ro, R, Ref, Ux, Uy, Vef, Vol

! Initialize
 Ux   = 0.D0
 Uy   = 0.D0
 Vol  = 0.D0
 Pin  = 0.D0
 Pout = 0.D0
 nodesIn  = 0
 nodesOut = 0
 Ro = bubbles(1,3)

! Loop through all nodes inside the bubble
 DO j = ymin_f, ymax_f
   DO i = xmin_f, xmax_f

     IF ( phi_f(i,j) >= 0.D0 ) THEN
       Ux  = Ux  + u_f(i,j,1)
       Uy  = Uy  + u_f(i,j,2)
       Vol = Vol + 1.D0
     END IF

     R =  DSQRT((DBLE(i)-bubbles(1,1))**2 + (DBLE(j)-bubbles(1,2))**2)
     IF (R < (Ro - IntWidth)) THEN
       Pin = Pin + p_f(i,j)
       nodesIn = nodesIn + 1
     ELSE IF (R > (Ro + IntWidth)) THEN
       Pout = Pout + p_f(i,j)
       nodesOut = nodesOut + 1
     END IF

   END DO
 END DO

! Save the initial volume during the first time step
 IF( iStep == 0 ) invInitVol = 1.D0/Vol

! Average velocity and pahse conservation (re-scale Ref to momentum grid)
 Ux  = Ux/Vol
 Uy  = Uy/Vol
 Ref = 0.5D0*DSQRT( Vol*invPi )
 Vef = Vol*invInitVol

! Compliance with Laplace's Law
 Pin  = Pin/DBLE(nodesIn)
 Pout = Pout/DBLE(nodesOut)
 Pdif = Pin - Pout
 Perr = (sigma/Ro - Pdif)*Ro/sigma

 OPEN(UNIT = 10, FILE = "stats.out", STATUS = "UNKNOWN", POSITION = "APPEND", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   WRITE(10,'(I9,6ES19.9)')istep,Ux,Uy,Vef,Ref,Pdif,Perr
   CLOSE(UNIT = 10)
 ELSE
   CALL MemAlloc(2)
   STOP "Error: Unable to open output file 'stats.out'."
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
