!-------------------------------------------------------------------------------
! Subroutine : Parameters
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Read input parameters and define simulation constants
!> @details
!! Read input parameters from files "properties.in" and "discrete.in" and define
!! constants for the simulation in the dual grid D2Q9 Lee-Lin multiphase LBM.

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

 SUBROUTINE Parameters

! Common Variables
 USE Domain
 USE FluidParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i
 INTEGER :: IO_ERR

! Read input data from file "properties.in"
 OPEN(UNIT = 10, FILE = "properties.in", STATUS = "OLD", ACTION = "READ", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   READ(10,*)
   READ(10,*) Maxstep
   READ(10,*)
   READ(10,*) tStat, tDump
   READ(10,*)
   READ(10,*) xmin, xmax, ymin, ymax
   READ(10,*)
   READ(10,*) rhoL, rhoH
   READ(10,*)
   READ(10,*) tauL, tauH
   READ(10,*)
   READ(10,*) IntWidth, sigma
   READ(10,*)
   READ(10,*) pConv
   CLOSE(UNIT = 10)
 ELSE
   STOP "Error: Unable to open input file 'properties.in'."
 END IF


!  Read bubble positions from file "discrete.in"
 OPEN(UNIT = 10, FILE = "discrete.in", STATUS = "OLD", ACTION = "READ", &
      IOSTAT = IO_ERR)
 IF ( IO_ERR == 0 ) THEN
   READ(10,*)
   READ(10,*) nBubbles
   READ(10,*)
   ALLOCATE( bubbles(1:nBubbles,1:3) )
   DO i = 1, nBubbles
     READ(10,*) bubbles(i,1), bubbles(i,2), bubbles(i,3)
   END DO
   CLOSE(UNIT=10)
 ELSE
   STOP "Error: Unable to open input file 'discrete.in'."
 END IF

! Scale parameters for order parameter mesh
 sigma    = 2.D0*sigma
 IntWidth = 2.D0*IntWidth
 DO i = 1, nBubbles
   bubbles(i,1) = 2*bubbles(i,1) - 1
   bubbles(i,2) = 2*bubbles(i,2) - 1
   bubbles(i,3) = 2*bubbles(i,3)
 END DO

! Constants needed in
 beta    = 12.D0*sigma/( IntWidth*( (rhoH - rhoL)**4 ) )
 kappa   = 1.5D0*sigma*IntWidth/( (rhoH - rhoL)**2 )
 beta4   = 4.0D0*beta
 kappaG  = 0.25D0*kappa
 kappaEf = 8.D0*inv12*kappaG
 kappa_6 = inv6*kappa
 rhostar = 0.5D0*( rhoL + rhoH )
 tauRhoStar = (tauH - tauL)/(rhoH - rhoL)

! Domain Size (order parameter mesh)
 xmin_f = xmin
 ymin_f = ymin
 xmax_f = 2*xmax - 1
 ymax_f = 2*ymax - 1

! Redefine MaxStep (iterate between 0 and MaxStep-1)
 MaxStep = MaxStep - 1


 RETURN
 END SUBROUTINE Parameters
