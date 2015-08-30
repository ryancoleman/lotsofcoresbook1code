!-------------------------------------------------------------------------------
! Subroutine : Parameters
! Revision   : 1.0 (2008-06-15)
! Author     : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Read input parameters and define simulation constants
!> @details
!! Read input parameters from files "properties.in" and "discrete.in" and define
!! constants for the simulation in the dual grid D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM.

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
 USE LBMParams
 IMPLICIT NONE

! Local Variables
 INTEGER :: i
 INTEGER :: IO_ERR


! Read parameter data
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
     READ(10,*) tauRho, tauPhi
     READ(10,*)
     READ(10,*) IntWidth, sigma, Gamma
     READ(10,*)
     READ(10,*) pConv
     CLOSE(UNIT = 10)
 ELSE
   STOP "Error: Unable to open input file 'properties.in'."
 END IF

!  Read bubble positions
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
     CLOSE(UNIT = 10)
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

! Fluid properties
 phistar    = 0.5D0*(rhoH - rhoL)
 phistar2   = phistar*phistar
 phistar4   = phistar2*phistar2
 eta        = 1.D0/(tauPhi + 0.5D0)
 eta2       = 1.D0 - eta
 invEta2    = 0.5D0/eta
 invTauPhi  = 1.D0/tauPhi
 invTauRho  = 1.D0/tauRho
 invTauRho2 = 1.D0 - 0.5D0*invTauRho

! Chemical potential parameters
 alpha   = 0.75D0*sigma/(IntWidth*phistar4)
 alpha4  = 4.D0*alpha
 kappa   = 0.50D0*alpha*(IntWidth*phistar)**2.D0
 kappaG  = 0.25D0*kappa
 kappa_6 = kappa*inv6

! Set domain limits for order parameter mesh
 xmin_f = xmin
 xmax_f = 2*xmax - 1
 ymin_f = ymin
 ymax_f = 2*ymax - 1

! Modified distribution weights for use in CollisionG
 Eg0T = Eg0C*invTauRho2
 Eg1T = Eg1C*invTauRho2
 Eg2T = Eg2C*invTauRho2

! Redefine MaxStep (iterate between 0 and MaxStep-1)
 MaxStep = MaxStep - 1

 RETURN
 END SUBROUTINE Parameters
